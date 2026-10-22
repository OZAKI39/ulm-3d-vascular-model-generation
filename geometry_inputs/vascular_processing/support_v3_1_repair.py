"""Local line-quantized rib edits, screened before acceptance for removal."""
from __future__ import annotations
import copy,pickle,json
from pathlib import Path
import numpy as np
import trimesh
import shapely
from shapely.geometry import box as rectangle
from .support_inputs import transform,table,save
from .support_v3_geometry import curtain,sweep,union,cut_section
from .support_v3_removal import probe_access
from .support_v2_removal import paths_for,cavity_overlap
from .support_removal import swept_path
from .support_v2_mesh import support_stl
from .support_custom_qc import intersection_volume
from .support_v3_1_residual import quantized_width,top_down_layers


def extrude(poly,lo,hi):
    parts=[]
    if hi<=lo+1e-6:return trimesh.Trimesh()
    for p in shapely.get_parts(shapely.make_valid(poly)):
        if p.geom_type!='Polygon' or p.area<1e-9:continue
        m=trimesh.creation.extrude_polygon(p,hi-lo);m.apply_translation([0,0,lo]);parts.append(m)
    return union(parts)


def candidate_preflight(trial,others,vascular,box,restore,inner,profile,other_necks=None):
    """Full changed fragments and finite 3-mm probe, before any acceptance.

    Neighbouring uncut supports remain obstacles: conservative, explicitly
    not a proof that a rejected edit has no possible removal sequence.
    """
    if not trial['mesh'].is_volume:return dict(passed=False,reason='NON_VOLUME_REPAIR')
    hit=intersection_volume(trial['mesh'],vascular)
    if hit>1e-5:return dict(passed=False,reason='FROZEN_VASCULAR_INTERSECTION',volume_mm3=hit)
    overlap=trimesh.boolean.intersection([trial['mesh'],box],engine='manifold')
    extra=trimesh.boolean.difference([overlap,trial['feet']],engine='manifold') if len(overlap.faces) else overlap
    if len(extra.faces) and abs(extra.volume)>1e-4:return dict(passed=False,reason='ILLEGAL_BOX_INTERSECTION',volume_mm3=abs(float(extra.volume)))
    full=union([trial['mesh'],others]);intended=union([c['intended'] for c in trial['cuts']]+([other_necks] if other_necks is not None else []))
    cut_checks=[cut_section(full,intended,c['origin'],c['normal'],c['aperture']) for c in trial['cuts']]
    if not all(c['passed'] for c in cut_checks):return dict(passed=False,reason='CUT_SECTION_CROSSES_UNWEAKENED_SUPPORT',cut_checks=cut_checks)
    def casting(m):
        result=m.copy();result.apply_transform(restore);return result
    v,b,whole=casting(vascular),casting(box),casting(full)
    if np.max(trial['base'])<.001 and cavity_overlap(casting(trial['mesh']).bounds,inner):
        from .support_v3_1_release import plate_release
        release=plate_release(union([vascular,box]),trial['mesh'])
        if not release['passed']:return dict(passed=False,reason='PLATE_RELEASE_FAILED',plate_release=release)
    tool=[]
    for c in trial['cuts']:
        target=transform([c['origin']],restore)[0];xy=trial['path_xy'];j=int(np.argmin(np.linalg.norm(xy-np.asarray(c['origin'])[:2],axis=1)))
        tangent=xy[min(j+1,len(xy)-1)]-xy[max(0,j-1)]
        normal=restore[:3,:3]@np.r_[-tangent[1],tangent[0],0.]
        remaining=trimesh.boolean.difference([whole,casting(c['intended'])],engine='manifold')
        check=probe_access(target,normal,trial['width_mm']/2,b,[v,b,remaining],3.)
        tool.append(dict(cut_id=c.get('cut_id'),**check))
        if not check['passed']:return dict(passed=False,reason='CUT_POINT_TOOL_ACCESS_FAILED',cut_checks=cut_checks,tool=tool)
    cut=union([c['cutter'] for c in trial['cuts']]);remaining=trimesh.boolean.difference([trial['mesh'],cut],engine='manifold') if len(cut.faces) else trial['mesh']
    parts=[casting(p) for p in remaining.split() if p.volume>1e-6];approved=casting(trial['approved_retained_stubs'])
    manager=trimesh.collision.CollisionManager();manager.add_object('vascular',v);manager.add_object('box',b)
    if len(others.faces):manager.add_object('neighbouring_supports',casting(others))
    movable=[]
    for p in parts:
        delta=trimesh.boolean.difference([p,approved],engine='manifold')
        if (not len(delta.faces) or abs(delta.volume)<2e-4) and intersection_volume(p,b)>1e-6:continue
        movable.append(p)
    movable=sorted(movable,key=lambda m:-m.bounds[1,2])
    for index,p in enumerate(movable):manager.add_object(f'own_fragment_{index}',p)
    checks=[]
    for index,p in enumerate(movable):
        manager.remove_object(f'own_fragment_{index}')
        internal=cavity_overlap(p.bounds,inner);local,paths=paths_for(p,b,not internal)
        initial=paths[0][1][0];rise=b.bounds[1,2]-p.bounds[0,2]+5
        _,_,axes=np.linalg.svd(p.vertices[:,:2]-p.vertices[:,:2].mean(0),full_matrices=False)
        for direction in (axes[1],-axes[1]):
            for distance in (.44,.88,1.5,2.5,4.):
                moved=initial.copy();moved[:2,3]+=distance*direction;end=moved.copy();end[2,3]+=rise
                paths.append((f'RIB_NORMAL_{distance}_{direction.tolist()}',[initial,moved,end]))
        found=None;attempts=[]
        for name,poses in paths:
            check=swept_path(local,[],poses,static_manager=manager);attempts.append(dict(name=name,**check))
            if check['passed']:
                end=transform(local.vertices,poses[-1]);exited=bool(end[:,2].min()>b.bounds[1,2]) if internal else bool(np.any(end.min(0)>b.bounds[1]) or np.any(end.max(0)<b.bounds[0]))
                if exited:found=dict(name=name,transforms=[x.tolist() for x in poses]);break
        checks.append(dict(fragment=index,volume_mm3=float(p.volume),passed=found is not None,path=found,attempts=attempts))
        if found is None:return dict(passed=False,reason='TRAPPED_FRAGMENT_PREDICTED',fragments=checks,tool=tool,cut_checks=cut_checks)
    return dict(passed=True,reason='LOCAL_REMOVAL_AND_TOOL_PROXY_PASSED',fragments=checks,tool=tool,cut_checks=cut_checks,
                scope='Changed rib against frozen anatomy/box and uncut neighbouring ribs; final ordered full-scaffold audit still mandatory')


def rebuilt(rib,profile,**kw):
    result=curtain(kw.get('xy',rib['path_xy']),kw.get('top',rib['top']),kw.get('base',rib['base']),profile,
                   kw.get('lines',round(rib['width_mm']/profile['line_width'])),kw.get('roof_width',rib['roof_width_mm']),
                   kw.get('segment_length',rib['segment_length_mm']),kw.get('window_width',3.))
    for k in ('rib_id','interval_id','branch_id','is_port'):result[k]=rib[k]
    for j,c in enumerate(result['cuts']):c['cut_id']=f"{result['rib_id']}_C{j:03}"
    return result


def move_cut_up(rib,index,delta,profile):
    """Physically etch a two-line notch, move its finite cutter and audit it.

    Lower material is NOT declared a retained foot. Full fragment screening
    will reject any newly trapped box-connected remainder.
    """
    result=copy.deepcopy(rib);old=result['cuts'][index];origin=np.array(old['origin'],float);origin[2]+=delta
    w,h=profile['line_width'],profile['layer_height'];xy=rib['path_xy']
    near=int(np.argmin(np.linalg.norm(xy-origin[:2],axis=1)))
    if origin[2]+2*h>=rib['top'][near]:raise ValueError('NOTCH_HAS_NO_ROOF_CLEARANCE')
    shell=sweep(xy,origin[2]-2*h,origin[2]+2*h,rib['roof_width_mm']+4*w)
    neck=sweep(xy,origin[2]-2*h-.01,origin[2]+2*h+.01,2*w)
    relief=trimesh.boolean.difference([shell,neck],engine='manifold')
    result['mesh']=trimesh.boolean.difference([rib['mesh'],relief],engine='manifold')
    result['ribs']=trimesh.boolean.difference([rib['ribs'],relief],engine='manifold')
    result['necks']=union([rib['necks'],trimesh.boolean.intersection([result['mesh'],neck],engine='manifold')])
    cutter=old['cutter'].copy();cutter.apply_translation([0,0,delta]);old.update(origin=origin,cutter=cutter,probe_target=origin.tolist())
    old['intended']=trimesh.boolean.intersection([result['mesh'],cutter],engine='manifold')
    if not len(old['intended'].faces) or old['intended'].volume<1e-7:raise ValueError('RELOCATED_CUT_MISSES_MATERIAL')
    return result


def historical_repairs(ribs,record,inputs,profile,out):
    restore=np.asarray(record['orientation']['casting_restore_transform']);folder=Path(record['folder'])
    v=trimesh.load_mesh(folder/'vascular_with_ports.stl');b=trimesh.load_mesh(folder/'box.stl');inner=np.array(inputs['inner_bounds']);w=profile['line_width']
    previous=record['removal'];actions=[];diagnoses=[]
    trapped=[x for x in previous.get('fragments',[]) if x['status']=='TRAPPED_FRAGMENT']
    # Diagnose every historical fragment by spatial overlap with restored ribs.
    for frag in trapped:
        bounds=np.array(frag['bounds_mm']);ids=[]
        for i,r in enumerate(ribs):
            vertices=transform(r['mesh'].vertices,restore)
            if np.all(vertices.max(0)>bounds[0]) and np.all(vertices.min(0)<bounds[1]):ids.append(i)
        diagnoses.append(dict(fragment_id=frag['fragment_id'],bounds_mm=frag['bounds_mm'],associated_ribs=[ribs[i]['rib_id'] for i in ids],
                              reason='All saved finite swept paths collided or lacked certified clearance; curved vessel/adjacent scaffold blocks the searched exits',
                              prior_attempt_count=len(frag.get('attempts',[]))))
        for i in ids:
            rib=ribs[i];others=union([r['mesh'] for j,r in enumerate(ribs) if j!=i])
            other_necks=union([c['intended'] for j,r in enumerate(ribs) if j!=i for c in r['cuts']])
            options=['SHORTEN_RIB','SPLIT_RIB','MOVE_CUT_WINDOW','MOVE_CUT_UPWARD','INCREASE_WINDOW','CHANGE_ROOT_DOMAIN']
            for option in options:
                action=dict(issue_id=frag['fragment_id'],rib_id=rib['rib_id'],repair=option,accepted=False)
                try:
                    if option=='SHORTEN_RIB':
                        if len(rib['path_xy'])<5:raise ValueError('RIB_TOO_SHORT')
                        trial=rebuilt(rib,profile,xy=rib['path_xy'][1:-1],top=rib['top'][1:-1],base=rib['base'][1:-1])
                        trial['mesh']=union([trial['mesh'],rib['roofs']]);trial['roofs']=rib['roofs'].copy()
                    elif option=='SPLIT_RIB':trial=rebuilt(rib,profile,segment_length=quantized_width(1.5,w))
                    elif option=='MOVE_CUT_WINDOW':trial=rebuilt(rib,profile,segment_length=quantized_width(rib['segment_length_mm']+2*w,w))
                    elif option=='MOVE_CUT_UPWARD':trial=move_cut_up(rib,0,4*w,profile)
                    elif option=='INCREASE_WINDOW':trial=rebuilt(rib,profile,window_width=quantized_width(4.,w))
                    else:
                        from .support_detection import bottom_points
                        tangent=np.gradient(rib['path_xy'],axis=0);normal=np.c_[-tangent[:,1],tangent[:,0]]/np.maximum(np.linalg.norm(tangent,axis=1)[:,None],1e-9)
                        moved=rib['path_xy']+normal*4*w;roots=bottom_points(b,np.c_[moved,rib['top']-.2]);base=np.where(np.isfinite(roots[:,2]),roots[:,2],0.)
                        if np.allclose(base,rib['base'],atol=.01):raise ValueError('NO_DIFFERENT_LEGAL_ROOT_IN_LOCAL_SEARCH')
                        trial=rebuilt(rib,profile,xy=moved,base=base);trial['mesh']=union([trial['mesh'],rib['roofs']]);trial['roofs']=rib['roofs'].copy()
                    check=candidate_preflight(trial,others,v,b,restore,inner,profile,other_necks);action['preflight']=check
                    if check['passed']:ribs[i]=trial;action['accepted']=True
                except (ValueError,RuntimeError) as exc:action['reason']=str(exc)
                actions.append(action)
                if action['accepted']:break
    # Eight is not hardcoded: enumerate the actual baseline's failed probes.
    toolpath=folder/'tool_access.json'
    if not toolpath.exists() and (folder/'nonplate_removal/tool_access.json').exists():toolpath=folder/'nonplate_removal/tool_access.json'
    tools=json.loads(toolpath.read_text()) if toolpath.exists() else []
    for bad in [x for x in tools if not x['passed']]:
        matched=False
        for i,rib in enumerate(ribs):
            for j,c in enumerate(rib['cuts']):
                if c.get('cut_id')!=bad['cut_id']:continue
                matched=True;others=union([r['mesh'] for k,r in enumerate(ribs) if k!=i])
                other_necks=union([c['intended'] for k,r in enumerate(ribs) if k!=i for c in r['cuts']])
                for steps in (1,2,4,8):
                    action=dict(issue_id=bad['cut_id'],rib_id=rib['rib_id'],repair='RELOCATE_NOTCH_TOWARD_OPEN_TOP',delta_mm=steps*w,accepted=False)
                    try:
                        trial=move_cut_up(rib,j,steps*w,profile);check=candidate_preflight(trial,others,v,b,restore,inner,profile,other_necks);action['preflight']=check
                        if check['passed']:ribs[i]=trial;action['accepted']=True
                    except (ValueError,RuntimeError) as exc:action['reason']=str(exc)
                    actions.append(action)
                    if action['accepted']:break
        if not matched:actions.append(dict(issue_id=bad['cut_id'],repair='RELOCATE_NOTCH_TOWARD_OPEN_TOP',accepted=False,reason='Cut topology changed by preceding fragment edit; final probe audit required'))
    save(out/'historical_fragment_diagnoses.json',diagnoses);save(out/'historical_removal_repairs.json',actions);table(out/'historical_removal_repairs.csv',actions)
    return ribs,actions


def local_variant(rib,grain,action,profile,actual):
    """Modify one grain neighbourhood, never dilate the whole support field."""
    w,h=profile['line_width'],profile['layer_height'];p=shapely.from_wkt(grain['footprint_wkt']);lo,hi=np.array(grain['bbox'])
    result=copy.deepcopy(rib);region=p.buffer(2*w,join_style='mitre')
    if action=='EXTEND_ROOF':
        wide=rebuilt(rib,profile,roof_width=rib['roof_width_mm']+w)
        crop=extrude(region,max(0.,lo[2]-profile['roof_gap_mm']-max(2.,rib['roof_width_mm'])),hi[2])
        addition=trimesh.boolean.intersection([wide['roofs'],crop],engine='manifold')
        result['roofs']=union([rib['roofs'],addition])
    elif action=='WIDEN_RIB_ONE_LINE':
        wide=rebuilt(rib,profile,lines=round(rib['width_mm']/w)+1,roof_width=rib['roof_width_mm']+w)
        crop=extrude(region,max(0.,lo[2]-4.),hi[2]);addition=trimesh.boolean.intersection([wide['mesh'],crop],engine='manifold')
        result['ribs']=union([rib['ribs'],addition])
    elif action=='EXTEND_RIB_SEGMENT':
        xy=rib['path_xy'].copy();tops=rib['top'].copy();bases=rib['base'].copy()
        target=np.mean(np.asarray(grain['centroid_path'])[:,:2],axis=0);at_start=np.linalg.norm(target-xy[0])<np.linalg.norm(target-xy[-1]);j=0 if at_start else -1
        direction=(xy[0]-xy[1]) if at_start else (xy[-1]-xy[-2]);direction/=np.linalg.norm(direction)
        if at_start:xy=np.vstack([xy[0]+w*direction,xy]);tops=np.r_[tops[0],tops];bases=np.r_[bases[0],bases]
        else:xy=np.vstack([xy,xy[-1]+w*direction]);tops=np.r_[tops,tops[-1]];bases=np.r_[bases,bases[-1]]
        extended=rebuilt(rib,profile,xy=xy,top=tops,base=bases)
        addition=trimesh.boolean.intersection([extended['mesh'],extrude(region,float(bases.min()),float(tops.max())+.01)],engine='manifold')
        result['ribs']=union([rib['ribs'],addition])
    else:
        # Downward closure grows the ORIGINAL rib, not a second support tree.
        footprint=p.buffer(w,join_style='mitre');start=lo[2]-h if grain['component']=='SUPPORT_SELF' else lo[2]-profile['roof_gap_mm']
        propagation=top_down_layers(footprint,start,actual,profile,float(np.min(rib['base'])))
        # A self residual at layer n needs material in n-1. Extending the
        # patch through n would create a new equally unsupported outer edge
        # (especially with first-layer elephant-foot compensation).
        base=propagation['base_z_mm'];top=max(start,hi[2]-(h if grain['component']=='SUPPORT_SELF' else profile['roof_gap_mm']))
        if action=='LOCAL_PILLAR':
            center=np.array(p.centroid.coords[0]);side=quantized_width(max(3*w,max(np.ptp(np.asarray(p.bounds).reshape(2,2),axis=0))),w)
            footprint=rectangle(*(center-side/2),*(center+side/2))
        addition=extrude(footprint,max(0.,base),top)
        result['ribs']=union([rib['ribs'],addition]);result['propagation']=propagation
    if not len(addition.faces):raise ValueError('LOCAL_REPAIR_DOES_NOT_REACH_GRAIN')
    result['mesh']=union([rib['mesh'],addition])
    if abs(result['mesh'].volume-rib['mesh'].volume)<1e-7:raise ValueError('NO_GEOMETRY_CHANGE_AT_GRAIN')
    # Additions below a designed root neck can form a wider legal foot. Never
    # exempt upper material or uncut interior bodies as retained stubs.
    rootcuts=[c for c in rib['cuts'] if c['kind']=='ROOT_NECK']
    if rootcuts:
        low=float(min(np.asarray(c['origin'])[2]-.031 for c in rootcuts))
        bbox=addition.bounds;extent=bbox[1]-bbox[0]+2.;extent[2]=max(.001,low-bbox[0,2]+1.)
        clip=trimesh.creation.box(extent,transform=trimesh.transformations.translation_matrix([*bbox.mean(0)[:2],low-extent[2]/2]))
        foot=trimesh.boolean.intersection([addition,clip],engine='manifold')
        if len(foot.faces):result['feet']=union([rib['feet'],foot]);result['approved_retained_stubs']=union([rib['approved_retained_stubs'],foot])
    result['repair_added_volume_mm3']=float(result['mesh'].volume-rib['mesh'].volume)
    return result


def repair_grains(ribs,grains,record,inputs,profile,actual,out):
    (out/'accepted_repairs').mkdir(parents=True,exist_ok=True)
    folder=Path(record['folder']);v=trimesh.load_mesh(folder/'vascular_with_ports.stl');b=trimesh.load_mesh(folder/'box.stl')
    restore=np.asarray(record['orientation']['casting_restore_transform']);inner=np.array(inputs['inner_bounds']);actions=[]
    for grain in grains:
        target=np.mean(np.asarray(grain['centroid_path']),axis=0)
        choices=sorted(range(len(ribs)),key=lambda i:(ribs[i]['branch_id']!=grain['branch'],np.min(np.linalg.norm(ribs[i]['path_xy']-target[:2],axis=1))))[:2]
        accepted=False
        for i in choices:
            rib=ribs[i];others=union([r['mesh'] for j,r in enumerate(ribs) if j!=i])
            other_necks=union([c['intended'] for j,r in enumerate(ribs) if j!=i for c in r['cuts']])
            priority=['EXTEND_ROOF','WIDEN_RIB_ONE_LINE','EXTEND_RIB_SEGMENT','LOCAL_VERTICAL_RIB','LOCAL_PILLAR'] if grain['component']!='SUPPORT_SELF' else ['TOP_DOWN_THIN_WEB','WIDEN_RIB_ONE_LINE']
            for name in priority:
                row=dict(grain_id=grain['grain_id'],component=grain['component'],rib_id=rib['rib_id'],repair=name,line_width_mm=profile['line_width'],
                         width_increment_lines=1 if name in ('EXTEND_ROOF','WIDEN_RIB_ONE_LINE','EXTEND_RIB_SEGMENT') else None,
                         local_buffer_increment_lines=1 if name in ('TOP_DOWN_THIN_WEB','LOCAL_VERTICAL_RIB','LOCAL_PILLAR') else None,accepted=False)
                try:
                    trial=local_variant(rib,grain,name,profile,actual)
                    check=candidate_preflight(trial,others,v,b,restore,inner,profile,other_necks);row['preflight']=check
                    if check['passed']:
                        ribs[i]=trial;row['accepted']=True;row['added_volume_mm3']=trial.get('repair_added_volume_mm3');accepted=True
                        trial['mesh'].export(out/'accepted_repairs'/f"{grain['grain_id']}_{rib['rib_id']}.stl")
                except (ValueError,RuntimeError) as exc:row['reason']=str(exc)
                actions.append(row)
                if accepted:break
            if accepted:break
        save(out/'repair_progress.json',actions)
        if len(actions)%20==0:print('V3.1 LOCAL REPAIRS',len(actions),'trials',sum(x['accepted'] for x in actions),'accepted',flush=True)
    save(out/'repair_preflight.json',actions);table(out/'repair_actions.csv',actions)
    return ribs,actions
