"""Roof coverage and repaired shared FDM scaffold, with explicit rejection logs."""
from __future__ import annotations
import math
import numpy as np
import trimesh
from scipy.spatial import cKDTree
from shapely.geometry import Point,GeometryCollection
from shapely.ops import unary_union
from shapely.strtree import STRtree
from . import support_v2_primitives as prim
from .support_custom_qc import canonical_stl
from .support_v2_mesh import support_stl
from .support_detection import bottom_points,layers,section_polygon,unsupported_layer
from .support_v2_layers import analyse,residual_support_area

def manager(mesh):
    m=trimesh.collision.CollisionManager();m.add_object('fixed',mesh);return m

def candidate_roofs(det,field,boxfield,vascular,box,branches,profile,gap=.2,topology=None):
    w=profile['line_width'];h=profile['layer_height'];vman=manager(vascular);bman=manager(box)
    cp=np.vstack([b['xyz'] for b in branches]);meta=[(b,i) for b in branches for i in range(len(b['xyz']))];tree=cKDTree(cp)
    seeds=[(np.array(p['xyz']),p,'vascular') for p in det['raw'] if p['surface_hit_valid'] and 'LAYER_UNSUPPORTED' in p['detectors_triggered']]
    # Add lower-envelope samples from every actual required polygon, including
    # locations that SWC point proximity alone can miss.
    for f,kind,obj in [(field,'vascular',vascular),(boxfield,'box',box)]:
        positions=[]
        for row,poly in zip(f['rows'],f['polygons']):
            if row['classification']!='SUPPORT_REQUIRED':continue
            pitch=6*w if kind=='vascular' else 16*w;lo=np.array(poly.bounds[:2]);hi=np.array(poly.bounds[2:]);positions.append([*row['representative_xy'],row['z_mm']+.01])
            if poly.area>pitch*pitch:
                for x in np.arange(lo[0]+pitch/2,hi[0],pitch):
                    for y in np.arange(lo[1]+pitch/2,hi[1],pitch):
                        if poly.intersects(Point(x,y)):positions.append([x,y,row['z_mm']+.01])
        if positions:
            hits=bottom_points(obj,np.asarray(positions))
            for xyz in hits:
                if np.isfinite(xyz).all() and xyz[2]>gap+3*h:seeds.append((xyz,None,kind))
    result=[];seen=set()
    for xyz,source,kind in seeds:
        key=(kind,*np.round(xyz/max(w,.2)).astype(int))
        if key in seen:continue
        seen.add(key);_,index=tree.query(xyz);b,i=meta[index];t=np.asarray(b['tangent'][i]);angle=math.atan2(t[1],t[0])
        length=prim.quantize(8*w/(1+3*abs(t[2])+b['curvature'][i]*b['radius'][i]),w,2) if kind=='vascular' else 20*w
        width=prim.quantize(min(5*w,max(3*w,1.5*b['radius'][i])),w) if kind=='vascular' else 20*w
        target=gap if kind=='vascular' else -.02;top=xyz.copy();roof,footprint=prim.roof_pad(top,length,width,3*h,target,angle if kind=='vascular' else 0.)
        valid=True;distance=None
        if kind=='vascular':
            for _ in range(31):
                distance=0. if vman.in_collision_single(roof) else float(vman.min_distance_single(roof))
                if distance>=gap-.002:break
                top[2]-=h/2;roof,footprint=prim.roof_pad(top,length,width,3*h,target,angle)
            valid=distance>=gap-.002 and not bman.in_collision_single(roof)
        elif vman.in_collision_single(roof):valid=False
        z=float(roof.bounds[1,2]);record=dict(roof_id=f'RP{len(result):04}',geometry_class='port' if kind=='vascular' and b['is_port'] else kind,
            branch_id=b['branch_id'] if kind=='vascular' else 'BOX',xyz=xyz.tolist(),local_tangent=t.tolist(),local_radius_mm=float(b['radius'][i]),
            curvature=float(b['curvature'][i]),arc_position=float(b['arc'][i]),height_mm=z,distance_to_bifurcation=source['distance_to_bifurcation'] if source else None,
            roof_direction=[float(np.cos(angle)),float(np.sin(angle))],length_mm=length,width_mm=width,thickness_mm=3*h,angle_rad=angle if kind=='vascular' else 0.,
            gap_mm=gap if kind=='vascular' else 0.,top_z_mm=z,measured_min_gap_mm=distance,valid=valid,
            rejection_reason='' if valid else 'ROOF_GAP_OR_COLLISION',required_region_ids=[],removal_accessibility='TO_BE_CHECKED_IN_CASTING_FRAME',
            _mesh=roof,_footprint=footprint)
        if kind=='vascular' and topology is not None and len(topology.get('bifurcations',[])):
            record['distance_to_bifurcation']=float(np.linalg.norm(np.array(topology['bifurcations'])-xyz,axis=1).min())
            record['bifurcation_distance_method']='Euclidean distance to transformed SWC bifurcation node'
        result.append(record)
    return result

def select_roofs(candidates,fields,profile):
    regions=[];polys=[]
    for f in fields:
        for row,poly in zip(f['rows'],f['polygons']):
            if row['classification']=='SUPPORT_REQUIRED':regions.append(row);polys.append(poly)
    if not polys:return [],dict(required_area_mm2=0,remaining_area_mm2=0)
    tree=STRtree(polys);parts=[]
    for roof in candidates:
        pieces={}
        if roof['valid']:
            fp=roof['_footprint'].buffer(profile['line_width']/2)
            for idx in tree.query(fp,predicate='intersects'):
                row=regions[idx];dz=row['z_mm']-roof['top_z_mm']
                if -.5*profile['layer_height']<=dz<=roof['gap_mm']+2*profile['layer_height']:
                    intersection=polys[idx].intersection(fp)
                    if intersection.area>1e-6:pieces[int(idx)]=intersection
        parts.append(pieces);roof['required_region_ids']=[regions[i]['region_id'] for i in pieces]
    remaining=list(polys);selected=[];available={i for i,p in enumerate(parts) if p}
    while available:
        scored=[]
        for i in available:
            gain=sum(remaining[j].intersection(poly).area for j,poly in parts[i].items())
            if gain>.001:
                r=candidates[i];damage=1 if r.get('distance_to_bifurcation') is not None and r['distance_to_bifurcation']<2 else 0
                scored.append((-gain,damage,r['curvature'],r['roof_id'],i))
        if not scored:break
        i=min(scored)[-1];selected.append(i);available.remove(i)
        for j,poly in parts[i].items():remaining[j]=remaining[j].difference(poly)
    # Delete only if the full per-layer polygon coverage remains unchanged.
    for i in list(reversed(selected)):
        redundant=True
        for j,poly in parts[i].items():
            other=unary_union([parts[k][j] for k in selected if k!=i and j in parts[k]])
            if poly.difference(other).area>1e-5:redundant=False;break
        if redundant:selected.remove(i);candidates[i]['status']='REDUNDANT_SUPPORT_PAD'
    # A short printable span between two roofs need not become a continuous
    # support wall. Remove a pad only if the *whole* layer deficit remains as
    # buildable with multi-chord bridge evidence and the conservative limit.
    previous=[]
    for field in fields:
        for row in field['rows']:
            if row['classification']=='SUPPORT_REQUIRED':
                n=row.get('layer_index',0)
                previous.append(field['layer_polygons'][n-1] if n>0 and 'layer_polygons' in field else GeometryCollection())
    expanded=STRtree([p.buffer(profile['bridge_limit_mm']+profile['line_width']) for p in polys]);eligible={j:set() for j in range(len(polys))};affects={i:[] for i in selected};footprints={i:candidates[i]['_footprint'].buffer(profile['line_width']/2) for i in selected}
    for i in selected:
        r=candidates[i]
        for j in expanded.query(footprints[i],predicate='intersects'):
            dz=regions[j]['z_mm']-r['top_z_mm']
            if -.5*profile['layer_height']<=dz<=r['gap_mm']+2*profile['layer_height']:eligible[int(j)].add(i);affects[i].append(int(j))
    def missing(j,indices):
        anchors=unary_union([footprints[k] for k in indices]);return residual_support_area(polys[j],anchors,profile,previous[j])[0].area
    baseline={j:missing(j,eligible[j]) for j in eligible};bridge_pruned=0
    for i in reversed(selected.copy()):
        trial={j:missing(j,eligible[j]-{i}) for j in affects[i]}
        if all(value<=baseline[j]+.00001 for j,value in trial.items()):
            selected.remove(i);candidates[i]['status']='REDUNDANT_WITH_TWO_ANCHORED_SHORT_BRIDGES';bridge_pruned+=1
            for j,value in trial.items():eligible[j].remove(i);baseline[j]=value
    return [candidates[i] for i in selected],dict(required_area_mm2=sum(p.area for p in polys),remaining_area_mm2=sum(baseline.values()),bridge_aware_redundant_pads_removed=bridge_pruned,selection='Deterministic full polygon supportability; bifurcation penalty; two-anchor bridge-aware redundant pad removal without loss of coverage')

def self_buildability(mesh,foundation,profile,bridge_limit=None):
    field=analyse(mesh,profile,'support',foundation,bridge_limit=bridge_limit)
    failures=[r for r in field['rows'] if r['classification']=='SUPPORT_REQUIRED' and r['area_mm2']>.005 and r['span_length_mm']>profile['line_width']]
    return dict(passed=not failures,unsupported_regions=len(failures),rows=failures,
        bridge_regions=[r for r in field['rows'] if r['classification']=='BRIDGEABLE'],bridge_calibration=profile['bridge_calibration_status'])

def prismatic_layer_check(mesh,profile):
    """Exact invariant layer sections of a constant-width sheared prism.

    Full scaffold union is still audited at every production layer. This fast
    primitive check uses actual mesh sections, not an inclination threshold.
    """
    h=profile['layer_height'];lo,hi=mesh.bounds[:,2]
    if hi-lo<2*h:return False
    levels=[lo+h*.5,lo+h*1.5,hi-h*1.5,hi-h*.5]
    polygons=[section_polygon(p) for p in mesh.section_multiplane([0,0,0],[0,0,1],levels)]
    return all(unsupported_layer(polygons[i+1],polygons[i],h,profile['line_width'])[0].area<.005 for i in [0,2])

def scaffold(roofs,vascular,box,profile,restore,inner,audit_folder=None):
    w=profile['line_width'];h=profile['layer_height'];vman=manager(vascular);bman=manager(box)
    pieces=[];records=[];nodes=[];edges=[];roots=[];cuts=[];rejected=[];attempts=[];approved=[];connected=[]
    def add(shape,kind,**meta):
        identity=f'P{len(records):04}';pieces.append(shape);records.append(dict(primitive_id=identity,kind=kind,volume_mm3=float(shape.volume),bbox_mm=shape.bounds.tolist(),**meta));return identity
    def legal(shape,clearance=.02,allow_box=False):
        if vman.in_collision_single(shape) or vman.min_distance_single(shape)<clearance-.01:return False,'VASCULAR_COLLISION_OR_GAP'
        if not allow_box and bman.in_collision_single(shape):return False,'BOX_COLLISION'
        return True,''
    def root_options(a,width):
        candidates=[]
        for shift in [(0,0),(2*w,0),(-2*w,0),(0,2*w),(0,-2*w),(6*w,0),(-6*w,0),(0,6*w),(0,-6*w)]:
            xy=a[:2]+shift;origin=np.r_[xy,a[2]-.03];hits,_,face=box.ray.intersects_location([origin],[[0,0,-1]],multiple_hits=True)
            surfaces=[(0.,'BUILD_PLATE')]
            for point,f in zip(hits,face):
                if point[2]<a[2]-.4 and box.face_normals[f,2]>.05:
                    cast=np.r_[point,1]@np.asarray(restore).T
                    domain='BOX_BOTTOM' if abs(cast[2]-inner[0,2])<.2 else 'BOX_SIDE_WALL'
                    surfaces.append((float(point[2]),domain))
            for z,domain in surfaces:
                start=np.r_[xy,z+7*h+prim.quantize(max(width-3*w,0)/2/.65,h)-.01]
                if start[2]>=a[2]-.1:continue
                trunk=prim.pillar(start,a,width);ok,why=legal(trunk,profile['support_xy_gap'])
                if not ok:continue
                foot,approved_foot,r=prim.breakaway_foot(xy,z,width,profile,domain,f'R{len(roots):04}')
                if vman.in_collision_single(foot):continue
                # Any intentional box intersection must remain inside the foot.
                # By construction, all material outside the approved pad is
                # exactly this neck/flare mesh. FCL rejects its box contact
                # conservatively without thousands of redundant Booleans.
                # The exported entire assembly still receives volume Booleans.
                if bman.in_collision_single(prim.foot_above_approved_pad(xy,z,width,profile,domain)):continue
                candidates.append((float(trunk.volume+foot.volume),domain,foot,approved_foot,r,trunk))
            if candidates:break
        return sorted(candidates,key=lambda a:(a[0],a[1]))
    def commit_root(option,load):
        _,domain,foot,approved_foot,r,trunk=option;r['root_id']=f'R{len(roots):04}';r['supported_load']=load
        roots.append(r);add(foot,'BREAKAWAY_FOOT',root_id=r['root_id'],domain=domain);add(trunk,'VERTICAL_PILLAR',root_id=r['root_id'],supported_load=load,width_mm=r['width_mm'])
        approved.append(approved_foot);cuts.append(dict(neck_id=f'N{len(cuts):04}',xyz=r['neck_xyz'],width_mm=r['neck_width_mm']+.08,height_mm=.06,role='ROOT_BREAK',root_id=r['root_id']))
    active=[]
    for roof in roofs:
        center=np.r_[np.array(roof['_footprint'].centroid.coords[0]),roof['top_z_mm']-roof['thickness_mm']]
        head_width=3*w;flare_h=prim.quantize(max(roof['length_mm']-head_width,roof['width_mm']-head_width)/2/.6,h)
        shoulder=center-[0,0,flare_h];bottom=shoulder-[0,0,3*h]
        head=prim.buttress(shoulder,center,head_width,[roof['length_mm'],roof['width_mm']],roof['angle_rad']);neck=prim.neck((shoulder+bottom)/2,head_width,3*h,roof['angle_rad'])
        ok,why=legal(trimesh.util.concatenate([head,neck]),.02)
        if not ok or bottom[2]<.5:
            rejected.append(dict(roof_id=roof['roof_id'],reason=why or 'ROOF_TOO_LOW_FOR_BREAK_NECK'));continue
        active.append(dict(roof=roof,a=bottom+[0,0,.005],head=head,neck=neck,center=(shoulder+bottom)/2))
    def commit_head(item):
        roof=item['roof'];add(roof['_mesh'],'ROOF_PAD',roof_id=roof['roof_id'],length_mm=roof['length_mm'],width_mm=roof['width_mm'],gap_mm=roof['gap_mm'])
        if roof['geometry_class']=='box':approved.append(roof['_mesh'])
        add(item['head'],'BUTTRESS',roof_id=roof['roof_id']);add(item['neck'],'BREAKAWAY_SEGMENT_NECK',roof_id=roof['roof_id'],width_mm=3*w)
        cuts.append(dict(neck_id=f'N{len(cuts):04}',xyz=item['center'].tolist(),width_mm=3*w*np.sqrt(2)+.08,height_mm=.05,role='ROOF_BREAK',roof_id=roof['roof_id']))
        connected.append(roof)
    # Prefer shared rails for nearly level stretches, then repaired Y graphs.
    active.sort(key=lambda p:(p['roof']['branch_id'],p['roof']['arc_position'],p['roof']['roof_id']))
    merges=0;rails=0;bridges=0
    while active:
        if len(active)%100==0:print('V2 GRAPH remaining',len(active),'merges',merges,'roots',len(roots),flush=True)
        first=active.pop(0);a=first['a'];peers=sorted(range(len(active)),key=lambda j:np.linalg.norm(active[j]['a']-a));joined=False
        nearby=[j for j in peers if np.linalg.norm(active[j]['a'][:2]-a[:2])<=profile['bridge_limit_mm']+8*w and abs(active[j]['a'][2]-a[2])<.15 and active[j]['roof']['branch_id']==first['roof']['branch_id']]
        if len(nearby)>=2:
            group=[first,active[nearby[0]],active[nearby[1]]];coords=np.array([q['a'] for q in group]);axis=np.argmax(np.ptp(coords[:,:2],axis=0));order=np.argsort(coords[:,axis]);left=coords[order[0]].copy();right=coords[order[-1]].copy();z=min(coords[:,2])-.5;left[2]=right[2]=z
            width=8*w;span=float(np.linalg.norm(right[:2]-left[:2])-width)
            options=[root_options(left,width),root_options(right,width)]
            rail=prim.horizontal_member(left,right,width,4*h)
            posts=[prim.pillar(np.r_[q['a'][:2],z],q['a'],width) for q in group]
            if span<=profile['bridge_limit_mm'] and all(options) and all(legal(p)[0] for p in [rail,*posts]):
                for opts in options:commit_root(opts[0],1.5)
                add(rail,'LONGITUDINAL_RAIL',span_mm=span,bridge_state='VALIDATED_BRIDGE_CANDIDATE',validation_scope='Two anchored ends and conservative software limit; physical coupon pending');rails+=1;bridges+=1
                for q,post in zip(group,posts):add(post,'VERTICAL_PILLAR',supported_load=1,width_mm=width);commit_head(q)
                for j in sorted(nearby[:2],reverse=True):active.pop(j)
                joined=True
        if joined:continue
        for j in peers[:6]:
            if j>=len(active):continue
            second=active[j];b=second['a'];distance=np.linalg.norm(a[:2]-b[:2])
            if distance>10 or np.linalg.norm(a-b)>14:continue
            base_z=min(a[2],b[2])-distance/(2*.5)-.4
            # Each hierarchy level is recorded, including every rejection.
            trials=[('LOWER_MERGE_Z',base_z-d,np.zeros(2),6*w) for d in [0,.5,1.]]
            trials += [('MOVE_MERGE_HORIZONTALLY',base_z-1,np.array(shift)*w*2,6*w) for shift in [(1,0),(-1,0),(0,1),(0,-1)]]
            trials += [('INCREASE_TRUNK_WIDTH',base_z-1,np.zeros(2),q*w) for q in [8,10]]
            trials += [('INSERT_SHORT_VERTICAL_PILLAR',base_z-2,np.zeros(2),8*w)]
            for action,z,shift,width in trials:
                merge=np.r_[(a[:2]+b[:2])/2+shift,z]
                row=dict(roof_a=first['roof']['roof_id'],roof_b=second['roof']['roof_id'],action=action,merge_xyz=merge.tolist(),width_mm=width)
                if z<1.8:attempts.append({**row,'accepted':False,'reason':'INSUFFICIENT_ROOT_HEIGHT'});continue
                arms=[prim.branch(merge,a,width),prim.branch(merge,b,width)];ok=all(legal(s,profile['support_xy_gap'])[0] for s in arms)
                if not ok:attempts.append({**row,'accepted':False,'reason':'COLLISION'});continue
                if not all(prismatic_layer_check(arm,profile) for arm in arms):attempts.append({**row,'accepted':False,'reason':'SELF_PRINTABILITY'});continue
                options=root_options(merge,prim.quantize(width*np.sqrt(2),w))
                if not options:attempts.append({**row,'accepted':False,'reason':'NO_COLLISION_FREE_ROOT'});continue
                attempts.append({**row,'accepted':True,'reason':''});commit_root(options[0],2)
                for end,shape in zip([a,b],arms):add(shape,'SELF_SUPPORTING_BRANCH',width_mm=width,supported_load=1,start=merge.tolist(),end=end.tolist());edges.append(dict(start=merge.tolist(),end=end.tolist(),kind='SELF_SUPPORTING_BRANCH',width_mm=width))
                nodes.append(dict(kind='Y_MERGE',xyz=merge.tolist(),supported_load=2));merges+=1
                commit_head(first);commit_head(second);active.pop(j);joined=True;break
            if joined:break
            # Before independent pillars, explicitly consider an anchored short
            # bridge and split shared subtrees. Their eligibility is recorded.
            attempts.append(dict(roof_a=first['roof']['roof_id'],roof_b=second['roof']['roof_id'],action='INSERT_SHORT_BRIDGE',accepted=False,reason='SPAN_OR_LEVEL_NOT_ELIGIBLE',span_mm=float(distance)))
            attempts.append(dict(roof_a=first['roof']['roof_id'],roof_b=second['roof']['roof_id'],action='SPLIT_SHARED_SUBTREES',accepted=False,reason='PAIR_HAS_NO_SMALLER_SHARED_SUBTREE'))
        if joined:continue
        width=prim.quantize(max(6*w,a[2]/18),w);options=root_options(a,width)
        if options:commit_root(options[0],1);commit_head(first)
        else:rejected.append(dict(roof_id=first['roof']['roof_id'],reason='NO_COLLISION_FREE_ROOT_AFTER_REPAIR_HIERARCHY'))
    graph=dict(nodes=nodes,edges=edges,roots=roots,segment_necks=cuts,primitives=records,merge_count=merges,rail_count=rails,bridge_count=bridges,
        roof_count=len(connected),root_count=len(roots),merge_attempts=attempts,rejected=rejected,
        topology_warning='TOPOLOGY_DEGENERATED_TO_INDEPENDENT_PILLARS' if connected and merges==0 and rails==0 and len(roots)>=len(connected)*.9 else None)
    if audit_folder:
        import pickle
        from .support_inputs import save
        save(audit_folder/'support_graph_before_export.json',graph)
        with (audit_folder/'scaffold_local_cache.pkl').open('wb') as f:pickle.dump((graph,connected,pieces,approved),f)
    model=support_stl(trimesh.boolean.union(pieces,engine='manifold')) if pieces else trimesh.Trimesh()
    graph['numerical_support_cleanup']=model.metadata.get('support_only_numerical_cleanup_tolerance_mm',0.)
    return model,graph,connected,pieces,approved
