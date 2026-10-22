"""Designed fragments, casting-top opening, and a collision-checked sequence."""
import json,time
from pathlib import Path
import numpy as np
import trimesh
from .support_inputs import save,table
from .support_custom_qc import canonical_stl,intersection_volume
from .support_removal import swept_path,removal_paths
from .support_v2_primitives import loft
from .support_v2_mesh import support_stl

def fragment_support(support,necks):
    cutters=[loft(np.array(n['xyz'])-[0,0,n['height_mm']/2],np.array(n['xyz'])+[0,0,n['height_mm']/2],n['width_mm']) for n in necks]
    cut=trimesh.boolean.union(cutters,engine='manifold') if cutters else trimesh.Trimesh()
    remaining=trimesh.boolean.difference([support,cut],engine='manifold') if cutters else support
    parts=[]
    for part in remaining.split(only_watertight=False):
        if part.volume<=1e-5:continue
        try:part=support_stl(part)
        except ValueError as exc:
            if not part.is_volume:raise
            # Preserve the indexed Boolean surface for collision diagnostics.
            # STL welding can collapse distinct seam vertices. This fallback
            # is explicitly a failed export gate, never a successful STL QC.
            part.metadata['fragment_stl_precision_failure']=str(exc)
        parts.append(part)
    return parts

def break_design_audit(support,necks,allowed_necks):
    if not necks:return dict(passed=True,cut_volume_outside_designed_necks_mm3=0.)
    cutters=[loft(np.array(n['xyz'])-[0,0,n['height_mm']/2],np.array(n['xyz'])+[0,0,n['height_mm']/2],n['width_mm']) for n in necks]
    cut=trimesh.boolean.intersection([support,trimesh.boolean.union(cutters,engine='manifold')],engine='manifold')
    allowed=trimesh.boolean.union(allowed_necks,engine='manifold') if allowed_necks else trimesh.Trimesh()
    extra=trimesh.boolean.difference([cut,allowed],engine='manifold') if len(allowed.faces) else cut
    volume=abs(float(extra.volume)) if len(extra.faces) else 0.
    return dict(passed=volume<1e-5,cut_volume_outside_designed_necks_mm3=volume,physical_break_force_validated=False,
                meaning='Detect a neighbouring unweakened support filling or crossing a designed break section; a virtual cut alone is not proof of breakability')

def approved_retained_stub(part,graph,roofs,profile):
    """Recognize designed feet/pads below/above the corresponding break cut.

    Being small and touching the box is not sufficient to exempt a fragment.
    """
    lo,hi=part.bounds;tol=.0002
    for root in graph['roots']:
        if root['domain']=='BUILD_PLATE':continue
        xy=np.array(root['xyz'][:2]);half=root['foot_width_mm']/2
        if np.all(lo[:2]>=xy-half-tol) and np.all(hi[:2]<=xy+half+tol) and hi[2]<=root['neck_xyz'][2]-.029+tol and lo[2]>=root['xyz'][2]-.031:
            return 'DESIGNED_BOX_ROOT_FOOT_BELOW_BREAK'
    for roof in roofs:
        if roof['geometry_class']!='box':continue
        cut=next((n for n in graph['segment_necks'] if n.get('roof_id')==roof['roof_id']),None)
        if not cut:continue
        xy=np.array(roof['xyz'][:2]);radius=np.hypot(roof['width_mm'],roof['length_mm'])/2
        if np.all(lo[:2]>=xy-radius-tol) and np.all(hi[:2]<=xy+radius+tol) and lo[2]>=cut['xyz'][2]+.024-tol and hi[2]<=roof['top_z_mm']+tol:
            return 'DESIGNED_SACRIFICIAL_BOX_ROOF_ABOVE_BREAK'
    return None

def paths_for(part,box,external):
    local,pivot,standard=removal_paths(part,box.bounds[1,2])
    initial=np.eye(4);initial[:3,3]=pivot;rise=float(box.bounds[1,2]-part.bounds[0,2]+5)
    paths=[standard[0]]
    for distance in [2.,4.]:
        for axis,sign in [(0,-1),(0,1),(1,-1),(1,1)]:
            moved=initial.copy();moved[axis,3]+=sign*distance;raised=moved.copy();raised[2,3]+=rise
            paths.append((f'CASTING_LATERAL_{axis}_{sign*distance:+g}_THEN_TOP',[initial,moved,raised]))
    if external:
        for axis in [0,1,2]:
            for sign in [-1,1]:
                final=initial.copy();final[axis,3]+=sign*(np.ptp(box.bounds[:,axis])+np.ptp(part.bounds[:,axis])+5)
                paths.append((f'OUTSIDE_{axis}_{sign:+d}',[initial,final]))
    paths+=standard[1:]
    # Tilting a wide fragment can lower its bottom beyond the nominal rise.
    # A clear partial lift is not a completed removal.
    for name,poses in paths:
        if name.startswith('OUTSIDE_'):continue
        final=poses[-1]
        bottom=(local.vertices@final[:3,:3].T+final[:3,3])[:,2].min()
        final[2,3]+=max(0.,box.bounds[1,2]+5.-bottom)
    return local,paths

def cavity_overlap(bounds,inner):
    """Conservative domain test: straddling pieces count as internal."""
    bounds=np.asarray(bounds);inner=np.asarray(inner)
    return bool(np.all(bounds[1]>inner[0]+.001) and np.all(bounds[0]<inner[1]-.001))

def opening_exit_audit(folder,result,inner,box_bounds):
    """Check complete final exit separately from collision-free motion."""
    folder=Path(folder);box_bounds=np.asarray(box_bounds);rows=[]
    for row in result['segments']:
        if row['status']!='REMOVAL_PATH_FOUND':continue
        stem=folder/'removal_fragments'/('fragment_'+row['segment_id'][1:])
        indexed=stem.with_suffix('.indexed.ply')
        mesh=trimesh.load_mesh(indexed if indexed.exists() else stem.with_suffix('.stl'))
        poses=np.asarray(row['path']['transforms']);local=mesh.vertices-poses[0,:3,3]
        final=local@poses[-1,:3,:3].T+poses[-1,:3,3]
        internal=cavity_overlap(row['bounds_mm'],inner)
        top_exit=bool(final[:,2].min()>box_bounds[1,2])
        outside=bool(np.any(final.max(0)<box_bounds[0]) or np.any(final.min(0)>box_bounds[1]))
        rows.append(dict(segment_id=row['segment_id'],conservatively_internal=internal,passed=top_exit if internal else outside,final_min_z_mm=float(final[:,2].min()),full_top_exit=top_exit,complete_outside_exit=outside))
    return dict(passed=all(r['passed'] for r in rows),checked_fragments=len(rows),rows=rows,method='Straddling AABBs count as internal; accepted internal paths must end wholly above casting opening, independent of print +Z')

def run_removal(record,inputs):
    folder=Path(record['folder']);support=trimesh.load_mesh(folder/'custom_support.stl');graph=json.loads((folder/'support_graph_v2.json').read_text());restore=np.array(record['orientation']['casting_restore_transform']);began=time.monotonic()
    vascular=trimesh.load_mesh(inputs['core']);box=trimesh.load_mesh(inputs['box']);inner=np.array(inputs['inner_bounds']);parts=fragment_support(support,graph['segment_necks'])
    import pickle
    from .support_v2_primitives import neck
    allowed=[];cache=folder/'scaffold_local_cache.pkl'
    if cache.exists():
        with cache.open('rb') as f:cached_graph,_,primitives,_=pickle.load(f)
        allowed=[p for p,r in zip(primitives,cached_graph['primitives']) if r['kind']=='BREAKAWAY_SEGMENT_NECK']
    allowed += [neck(root['neck_xyz'],root['neck_width_mm'],root['neck_height_mm']) for root in graph['roots']]
    break_audit=break_design_audit(support,graph['segment_necks'],allowed)
    fragment_dir=folder/'removal_fragments';fragment_dir.mkdir(exist_ok=True);movable={};all_parts={};rows=[];fixed=trimesh.collision.CollisionManager();fixed.add_object('vascular_and_ports',vascular);fixed.add_object('casting_box',box)
    roofs=json.loads((folder/'roofs.json').read_text())
    for i,part in enumerate(parts):
        retained_reason=approved_retained_stub(part,graph,roofs,{})
        precision_error=part.metadata.get('fragment_stl_precision_failure');part.apply_transform(restore)
        try:part=support_stl(part)
        except ValueError as exc:
            if not part.is_volume:raise
            precision_error=str(exc)
        part.export(fragment_dir/f'fragment_{i:04}.stl')
        if precision_error:part.export(fragment_dir/f'fragment_{i:04}.indexed.ply')
        overlap=intersection_volume(part,box);external=not cavity_overlap(part.bounds,inner)
        row=dict(segment_id=f'S{i:04}',volume_mm3=float(part.volume),bounds_mm=part.bounds.tolist(),external=external,coordinate_frame='FROZEN_CASTING',status='PENDING',path=None,attempts=[],fragment_stl_precision_failure=precision_error)
        all_parts[row['segment_id']]=part
        if overlap>1e-6 and retained_reason:
            row['status']='APPROVED_FOOT_STUB_REMAINS_ON_BOX';row['retained_reason']=retained_reason;fixed.add_object(row['segment_id'],part)
        else:movable[row['segment_id']]=part;fixed.add_object(row['segment_id'],part)
        rows.append(row)
    order=[];pending=[r for r in rows if r['status']=='PENDING'];pending.sort(key=lambda r:(not r['external'],-r['bounds_mm'][1][2],r['volume_mm3'],r['segment_id']))
    # Retry blocked fragments after neighbours have actually been removed.
    for iteration in range(3):
        progress=False
        for row in pending:
            if row['status']!='PENDING':continue
            identity=row['segment_id'];part=movable[identity];fixed.remove_object(identity);local,paths=paths_for(part,box,row['external'])
            for name,poses in paths:
                audit=swept_path(local,[],poses,static_manager=fixed)
                row['attempts'].append(dict(round=iteration,path=name,**audit))
                if audit['passed']:
                    row.update(status='REMOVAL_PATH_FOUND',path=dict(name=name,transforms=[m.tolist() for m in poses],coordinate_frame='FROZEN_CASTING'),minimum_clearance_mm=audit['min_clearance_mm'],sequence_index=len(order));order.append(identity);progress=True;break
            if row['status']=='PENDING':fixed.add_object(identity,part)
        if not progress:break
    for row in pending:
        if row['status']=='PENDING':row['status']='TRAPPED_REMOVAL_FRAGMENT'
    result=dict(fragment_count=len(movable),retained_foot_stubs=sum(r['status']=='APPROVED_FOOT_STUB_REMAINS_ON_BOX' for r in rows),trapped_pieces=sum(r['status']=='TRAPPED_REMOVAL_FRAGMENT' for r in rows),sequence=order,segments=rows,
        frame='FROZEN_CASTING',opening_axis=[0,0,1],other_fragments_are_obstacles=True,break_design_audit=break_audit,precision_failed_fragments=sum(bool(r['fragment_stl_precision_failure']) for r in rows),physical_removal_validated=False,runtime_seconds=time.monotonic()-began)
    result['opening_exit_audit']=opening_exit_audit(folder,result,inner,box.bounds)
    save(folder/'support_removal_plan.json',result);save(folder/'removal_sequence.json',dict(sequence=order,frame='FROZEN_CASTING'));table(folder/'support_segments.csv',rows)
    for roof in roofs:
        point=np.r_[roof['xyz'][:2],roof['top_z_mm']-roof['thickness_mm']/2,1.]@restore.T;matched=[]
        for row in rows:
            part=all_parts[row['segment_id']]
            if np.all(point[:3]>=part.bounds[0]-.001) and np.all(point[:3]<=part.bounds[1]+.001) and part.contains([point[:3]])[0]:matched.append(row)
        roof['removal_fragment_ids']=[r['segment_id'] for r in matched]
        roof['removal_accessibility']=';'.join(sorted({r['status'] for r in matched})) if matched else 'NO_UNIQUE_FRAGMENT_MATCH'
    save(folder/'roofs.json',roofs);table(folder/'selected_roofs.csv',roofs)
    text=['# 分段拆除顺序','', '先恢复浇注姿态。先在设计窄颈处断开支撑，再按以下顺序取出。内部碎片必须经过模型开口；外部碎片也进行了血管、端口、盒体及未拆碎片的扫掠碰撞检查。物理断裂力和夹具空间尚未验证。','']
    text.extend(f'{i+1}. {identity}' for i,identity in enumerate(order));text.extend(['',f"未找到路径的碎片：{result['trapped_pieces']} 个。完整坐标和失败路径见 support_removal_plan.json。"]);(folder/'support_removal_plan.md').write_text('\n'.join(text)+'\n')
    print('V2 REMOVAL',record['candidate_id'],result['fragment_count'],'trapped',result['trapped_pieces'],flush=True);return result
