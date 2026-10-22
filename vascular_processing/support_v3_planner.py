"""Interval-first rib planning using the unchanged V2 detection functions."""
from __future__ import annotations
import hashlib, json, pickle, time
from pathlib import Path
import numpy as np
import trimesh
from scipy.spatial import cKDTree
from shapely.geometry import Point, box as rectangle
from .support_inputs import save, table, transform, digest
from .support_v2_inputs import oriented
from .support_v2_orientation import fits, rigid
from .support_custom_detection import fused_detect
from .support_v2_layers import analyse
from .support_detection import bottom_points, contiguous_runs
from .support_custom_qc import intersection_volume, canonical_stl
from .support_v2_mesh import support_stl, exact_stl
from .support_v3_geometry import curtain, union, cut_section

POLICY = dict(architecture='FDM_SEGMENTED_RIB_SUPPORT_V3', detection='FROZEN_DETECTION_BASE',
    line_candidates=[1,2,3], gap_mm=.20, gap_sensitivity_mm=[.15,.20,.25], bridge_limit_mm=0.,
    bridge_policy='NO_UNCALIBRATED_BRIDGE', probe_diameter_mm=3., segment_search_mm=[3.,5.,8.],
    coverage_area_tolerance_mm2=1e-5, self_region_tolerance_mm2=.005,
    physical_status='NO PHYSICAL COUPON AVAILABLE', implementation_revision=1)


def fingerprint(inputs, orientation, profile):
    detector = {p.name:digest(p) for p in [Path(__file__).with_name(n) for n in
                ['support_detection.py','support_custom_detection.py','support_freefloating.py','support_v2_layers.py']]}
    value = dict(input_sha=inputs['hashes'], orientation=orientation, profile_hash=hashlib.sha256(json.dumps(profile,sort_keys=True).encode()).hexdigest(),
                 policy_hash=hashlib.sha256(json.dumps(POLICY,sort_keys=True).encode()).hexdigest(), frozen_detector=detector)
    return hashlib.sha256(json.dumps(value,sort_keys=True).encode()).hexdigest(), value


def critical_intervals(det, field, branches, topology, inner):
    """Associate full required polygon vertices with SWC arc coordinates.

    Fusion remains V2. This adds the layer-required coverage envelope, including
    short spans that the long-span detector alone does not label critical.
    """
    points=np.vstack([b['xyz'] for b in branches]); tree=cKDTree(points)
    codes=np.concatenate([np.full(len(b['xyz']),i) for i,b in enumerate(branches)])
    local=np.concatenate([np.arange(len(b['xyz'])) for b in branches])
    assigned={i:[] for i in range(len(branches))}; masks=[np.zeros(len(b['xyz']),bool) for b in branches]
    for ri,(row,p) in enumerate(zip(field['rows'],field['polygons'])):
        if row['classification']!='SUPPORT_REQUIRED':continue
        coords=np.asarray(p.exterior.coords)
        q=np.c_[coords,np.full(len(coords),row['z_mm'])]
        ids=tree.query(q)[1]
        for code in np.unique(codes[ids]):
            idx=local[ids[codes[ids]==code]];masks[code][idx.min():idx.max()+1]=True
            assigned[int(code)].append(ri)
    lookup={b['branch_id']:i for i,b in enumerate(branches)}
    for interval in det['critical_intervals']:
        i=lookup[interval['branch_id']];b=branches[i]
        masks[i]|=(b['arc']>=interval['arc_start_mm'])&(b['arc']<=interval['arc_end_mm'])
    result=[]
    for bi,b in enumerate(branches):
        # Join nearby samples to produce intervals rather than one support per
        # polygon representative point. Keep true unsupported gaps separate.
        mask=masks[bi]
        for a,z in contiguous_runs(~mask):
            if a>0 and z<len(mask)-1 and b['arc'][z]-b['arc'][a]<1.:
                mask[a:z+1]=True
        for a,z in contiguous_runs(mask):
            a=max(0,a-2);z=min(len(mask)-1,z+2)
            if a==z:continue
            xyz=b['xyz'][a:z+1];length=float(b['arc'][z]-b['arc'][a]);radius=float(np.mean(b['radius'][a:z+1]))
            bif=np.asarray(topology['bifurcations']);near_bif=bool(len(bif) and cKDTree(bif).query(xyz)[0].min()<radius*2)
            wall=min(float(np.min(xyz[:,:2]-inner[0,:2])),float(np.min(inner[1,:2]-xyz[:,:2])))
            kind='PORT_REGION' if b['is_port'] else 'BIFURCATION_REGION' if near_bif else 'WALL_ADJACENT_REGION' if wall<3 else 'LONG_TUBULAR_SPAN' if length>=3 else 'SHORT_LOCAL_OVERHANG'
            selected=[]
            for ri in assigned[bi]:
                r=field['rows'][ri]
                if cKDTree(xyz).query([*r['representative_xy'],r['z_mm']])[0]<radius*3+2:selected.append(ri)
            result.append(dict(interval_id=f'V3I{len(result):03}',branch_id=b['branch_id'],is_port=bool(b['is_port']),semantic='PORT' if b['is_port'] else 'VESSEL',
                classification=kind,arc_start_mm=float(b['arc'][a]),arc_end_mm=float(b['arc'][z]),arc_length_mm=length,
                local_tangent=np.mean(b['tangent'][a:z+1],axis=0).tolist(),local_curvature=float(np.mean(b['curvature'][a:z+1])),local_radius_mm=radius,
                height_range_mm=[float(xyz[:,2].min()),float(xyz[:,2].max())],gravity_relation_abs_tangent_z=float(np.mean(abs(b['tangent'][a:z+1,2]))),
                required_polygon_ids=[field['rows'][i]['region_id'] for i in selected],required_polygon_wkt=[field['polygons'][i].wkt for i in selected],
                bifurcation_penalty=near_bif,branch_index=bi,start_index=a,end_index=z))
    return result


def detect_candidate(folder,row,inputs,profile,cfg,mold,vascular,box,branches,topology):
    folder.mkdir(parents=True,exist_ok=True);m=np.array(row['print_transform']);key,detail=fingerprint(inputs,row,profile)
    model=mold.copy();model.apply_transform(m);v=vascular.copy();v.apply_transform(m);b=box.copy();b.apply_transform(m)
    bs,top=oriented(branches,topology,m)
    corners=np.array(np.meshgrid(*np.array(inputs['inner_bounds']).T)).reshape(3,-1).T;corners=transform(corners,m);inner=np.array([corners.min(0),corners.max(0)])
    cache=folder/'frozen_detection_cache.pkl';manifest=folder/'cache_fingerprint.json'
    if cache.exists() and manifest.exists() and json.loads(manifest.read_text())['key']==key:
        with cache.open('rb') as f:det,field,boxfield,intervals,bs=pickle.load(f)
        print('V3 DETECTOR CACHE',row['orientation_id'],flush=True)
    else:
        det=fused_detect(v,b,bs,top,inner,profile,cfg,folder)
        field=analyse(v,profile,'vascular',b,folder,bridge_limit=0.)
        boxfield=analyse(b,profile,'box',None,folder,bridge_limit=0.)
        intervals=critical_intervals(det,field,bs,top,inner)
        with cache.open('wb') as f:pickle.dump((det,field,boxfield,intervals,bs),f)
        save(manifest,dict(key=key,**detail))
    table(folder/'critical_intervals.csv',intervals)
    return model,v,b,bs,top,inner,det,field,boxfield,intervals


def planned_coverage(field,support,profile):
    """Separate design-stage audit; never substitutes for actual extrusion."""
    from .support_detection import section_polygon
    from shapely.ops import unary_union
    levels=sorted({round(r['z_mm']-d,6) for r in field['rows'] for d in [.20,.30,.40] if r['z_mm']-d>0})
    sections=support.section_multiplane([0,0,0],[0,0,1],levels)
    polys={z:section_polygon(s) for z,s in zip(levels,sections)};rows=[];area=miss=0.
    for r,p in zip(field['rows'],field['polygons']):
        if r['classification']!='SUPPORT_REQUIRED':continue
        cover=unary_union([polys.get(round(r['z_mm']-d,6),Point().buffer(0)) for d in [.20,.30,.40]])
        missing=p.difference(cover.buffer(profile['line_width']/2));u=float(missing.area);area+=p.area;miss+=u
        rows.append(dict(**r,uncovered_area_mm2=u,covered=u<=POLICY['coverage_area_tolerance_mm2']))
    return dict(percent=100*(1-miss/area) if area else 100.,uncovered_region_count=sum(not r['covered'] for r in rows),
                all_required_covered=all(r['covered'] for r in rows),rows=rows,method='DESIGN_SECTIONS_ONLY_NOT_ACTUAL_TOOLPATH')


def plan(folder,row,inputs,profile,cfg,mold,vascular,box,branches,topology,lines=2):
    start=time.monotonic()
    model,v,b,bs,top,inner,det,field,boxfield,intervals=detect_candidate(folder,row,inputs,profile,cfg,mold,vascular,box,branches,topology)
    w,h=profile['line_width'],profile['layer_height'];ribs=[];rejected=[];chooser=[]
    for interval in intervals:
        branch=bs[interval['branch_index']];a,z=interval['start_index'],interval['end_index'];xyz=branch['xyz'][a:z+1].copy()
        xy=xyz[:,:2];distance=np.r_[True,np.linalg.norm(np.diff(xy,axis=0),axis=1)>.025];xyz=xyz[distance];xy=xyz[:,:2]
        if len(xy)<2:continue
        bottom=bottom_points(v,xyz)
        good=np.isfinite(bottom).all(axis=1)
        if not good.all():bottom[~good]=xyz[~good]-[0,0,interval['local_radius_mm']]
        tops=np.floor((bottom[:,2]-profile['roof_gap_mm'])/h)*h
        roots=bottom_points(b,bottom-[0,0,profile['roof_gap_mm']+.1])
        bases=np.where(np.isfinite(roots[:,2]),roots[:,2],0.)
        # No FDM roof is allowed to root on another blood vessel. Collision is
        # checked on each interval, not silently accepted by a final union.
        seglen=3. if interval['local_curvature']>.12 or interval['bifurcation_penalty'] else 5. if interval['local_curvature']>.04 else 8.
        roofwidth=np.ceil((1.6*interval['local_radius_mm']+2*w)/w)*w
        active=tops-bases>max(2.,roofwidth)
        direction=np.gradient(xy,axis=0);normal=np.c_[-direction[:,1],direction[:,0]]/np.maximum(np.linalg.norm(direction,axis=1)[:,None],1e-9)
        # The vessel inside a wall is already supported by that wall. A rib
        # must stop before entering it, not become a welded removal fragment.
        for sign in [-1,0,1]:
            probe=np.c_[xy+sign*normal*roofwidth/2,tops-.1]
            active &= ~b.contains(probe)
        # Discontinuous root domains cannot be connected by a sloped foot.
        jumps=np.flatnonzero(abs(np.diff(bases))>h)
        for j in jumps:active[max(0,j-1):min(len(active),j+3)]=False
        for run,(aa,zz) in enumerate(contiguous_runs(active)):
            if zz<=aa:continue
            ident=interval['interval_id']+f'_{run}'
            try:
                rib=curtain(xy[aa:zz+1],tops[aa:zz+1],bases[aa:zz+1],profile,lines,roofwidth,seglen,3.)
                overlap=intersection_volume(rib['mesh'],v);lowering=0.
                while overlap>1e-5 and lowering<.4:
                    lowering=round(lowering+h,6)
                    rib=curtain(xy[aa:zz+1],tops[aa:zz+1]-lowering,bases[aa:zz+1],profile,lines,roofwidth,seglen,3.)
                    overlap=intersection_volume(rib['mesh'],v)
                offset=0.
                if overlap>1e-5:
                    # A bifurcation/lower crossing may block a vertical plane.
                    # Search neighbouring parallel curtains, preserving the
                    # anatomy and keeping a near-tangent narrow ridge.
                    for delta in [w,-w,2*w,-2*w,3*w,-3*w,4*w,-4*w]:
                        shifted=xy[aa:zz+1]+normal[aa:zz+1]*delta
                        for extra_gap in [0.,h,2*h]:
                            trial=curtain(shifted,tops[aa:zz+1]-extra_gap,bases[aa:zz+1],profile,lines,roofwidth,seglen,3.)
                            trial_overlap=intersection_volume(trial['mesh'],v)
                            if trial_overlap<=1e-5:
                                rib=trial;overlap=trial_overlap;offset=delta;lowering=extra_gap;break
                        if overlap<=1e-5:break
                if overlap>1e-5:
                    rejected.append(dict(interval_id=ident,reason='VERTICAL_RIB_INTERSECTS_LOWER_VESSEL',intersection_mm3=overlap));continue
                rib.update(interval_id=interval['interval_id'],rib_id=ident,branch_id=interval['branch_id'],is_port=interval['is_port'],collision_clearance_lowering_mm=lowering,parallel_offset_mm=offset)
                ribs.append(rib)
                chooser.append(dict(interval_id=ident,selected='RIB',alternatives=['PILLAR','RIB + PILLAR','SHORT RAIL','SHORT BRIDGE','TREE/MERGE'],
                    rationale='Continuous critical arc; vertical thin wall preferred; no free span or long inclined member',priority=['FDM_BUILDABILITY','COVERAGE','REMOVABILITY','MATERIAL']))
            except (ValueError,RuntimeError) as exc:rejected.append(dict(interval_id=ident,reason=str(exc)))
    print('V3 RIB GEOMETRY',row['orientation_id'],len(intervals),'intervals',len(ribs),'curtains',len(rejected),'rejected',flush=True)
    try:support=support_stl(union([r['mesh'] for r in ribs]))
    except ValueError:
        # Isolate the problematic interval instead of either exporting an
        # invalid STL as production or abandoning all other measured evidence.
        accepted=[];support=trimesh.Trimesh()
        for rib in ribs:
            try:
                trial=support_stl(union([support,rib['mesh']]))
            except ValueError as exc:
                rib['mesh'].export(folder/(rib['rib_id']+'_precision_failure.ply'))
                rejected.append(dict(interval_id=rib['rib_id'],reason='SUPPORT_UNION_PRECISION_FAILURE',detail=str(exc)))
                continue
            support=trial;accepted.append(rib)
        ribs=accepted
        if not ribs:raise ValueError('NO_EXPORTABLE_INTERVAL_RIBS')
    assembly_raw=union([model,support]);assembly_error=None
    try:assembly=exact_stl(assembly_raw)
    except ValueError as exc:
        # Keep the indexed Boolean as diagnostic evidence. Never smooth or
        # simplify the frozen anatomical part to repair support export seams.
        assembly=assembly_raw;assembly_error=str(exc);assembly.export(folder/'assembly_indexed.ply')
    for name,mesh in [('mold_without_support.stl',model),('vascular_with_ports.stl',v),('box.stl',b),('custom_support_v3.stl',support),('mold_with_custom_support_v3.stl',assembly)]:mesh.export(folder/name)
    exports={}
    for name,key in [('support_ribs.stl','ribs'),('support_roofs.stl','roofs'),('support_cut_necks.stl','necks')]:
        part=union([r[key] for r in ribs]);part.export(folder/name);exports[name]=dict(empty=not len(part.faces),volume_mm3=float(part.volume))
    for name,items in [('support_exterior.stl',[r['mesh'] for r in ribs if r['is_port']]),('support_pillars.stl',[])]:
        part=union(items);part.export(folder/name);exports[name]=dict(empty=not len(part.faces),volume_mm3=float(part.volume),empty_reason='No auxiliary pillars required by this topology' if not items else None)
    self_field=analyse(support,profile,'support',b,folder,bridge_limit=0.)
    self_rows=[r for r in self_field['rows'] if r['classification']=='SUPPORT_REQUIRED']
    cov=planned_coverage(field,support,profile);bcov=planned_coverage(boxfield,support,profile)
    cuts=[];designed_necks=union([c['intended'] for r in ribs for c in r['cuts']])
    for rib in ribs:
        for i,c in enumerate(rib['cuts']):
            c['cut_id']=f"{rib['rib_id']}_C{i:03}"
            audit=cut_section(support,designed_necks,c['origin'],c['normal'],c['aperture'])
            cuts.append(dict(cut_id=c['cut_id'],rib_id=rib['rib_id'],kind=c['kind'],origin=np.asarray(c['origin']).tolist(),normal=np.asarray(c['normal']).tolist(),probe_target=c['probe_target'],**audit))
    with (folder/'rib_geometry_cache.pkl').open('wb') as f:pickle.dump(ribs,f)
    table(folder/'cut_points.csv',cuts);table(folder/'topology_chooser.csv',chooser);table(folder/'rejected_intervals.csv',rejected)
    ribrows=[dict(rib_id=r['rib_id'],interval_id=r['interval_id'],branch_id=r['branch_id'],is_port=r['is_port'],segment_count=r['segment_count'],width_mm=r['width_mm'],line_count=lines,roof_width_mm=r['roof_width_mm'],window_count=len(r['windows']),windows=r['windows'],length_mm=float(np.linalg.norm(np.diff(r['path_xy'],axis=0),axis=1).sum()),bbox_mm=r['mesh'].bounds.tolist()) for r in ribs]
    table(folder/'rib_segments.csv',ribrows);save(folder/'rib_segments.json',ribrows)
    save(folder/'print_transform_v3.json',dict(transform_4x4=row['print_transform'],scale_factor=1.,rigid=rigid(row['print_transform'])))
    save(folder/'casting_restore_transform_v3.json',dict(transform_4x4=row['casting_restore_transform'],frame='FROZEN_CASTING',opening_axis=[0,0,1]))
    result=dict(candidate_id=row['orientation_id'],folder=str(folder),architecture=POLICY['architecture'],orientation=row,scale_factor=1.,
        critical_intervals=len(intervals),rib_segments=sum(r['segment_count'] for r in ribs),roof_ridges=len(ribs),pillar_count=0,root_count=sum(len(r['feet'].split()) for r in ribs),
        geometry_parameters=dict(body_lines=lines,cut_neck_lines=int(profile.get('cut_neck_lines',1))),
        window_count=sum(len(r['windows']) for r in ribs),bridge_policy=POLICY['bridge_policy'],free_bridge_count=0,long_inclined_branch_count=0,line_count=lines,
        line_candidates=[1,2,3],line_selection_reason='Two-line geometry baseline; real slicer preservation required before selection',
        support_volume_mm3=float(support.volume),component_count=len(support.split()),build_volume_fit=fits(assembly.vertices,np.eye(4),profile['build_bounds_mm']),
        vascular_intersection_mm3=intersection_volume(support,v),self_buildability=dict(passed=not self_rows,unsupported_region_count=len(self_rows),rows=self_rows,bridge_credit=False,area_tolerance_mm2=.005),
        planned_coverage={k:val for k,val in cov.items() if k!='rows'},box_coverage={k:val for k,val in bcov.items() if k!='rows'},
        cuts_valid=all(c['passed'] for c in cuts),cut_count=len(cuts),rejected_intervals=rejected,part_exports=exports,
        assembly_export_error=assembly_error,slice={},removal={},physical_gap_status='PHYSICAL_GAP_UNCALIBRATED',physical_breakaway='UNVALIDATED',planning_runtime_seconds=time.monotonic()-start)
    table(folder/'planned_coverage.csv',cov['rows']);save(folder/'candidate.json',result)
    print('V3 PLAN',row['orientation_id'],'coverage',round(cov['percent'],3),'self',len(self_rows),'cuts',result['cuts_valid'],flush=True)
    return result
