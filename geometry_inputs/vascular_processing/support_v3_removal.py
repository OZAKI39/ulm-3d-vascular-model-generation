"""Finite cuts, separate 3 mm tool-access proxy and ordered fragment removal."""
from __future__ import annotations
import pickle, time, hashlib, json
from pathlib import Path
import numpy as np
import trimesh
from .support_inputs import save, table, transform, digest
from .support_v3_geometry import union
from .support_v2_removal import paths_for, cavity_overlap
from .support_removal import swept_path
from .support_custom_qc import intersection_volume


def probe_access(target, side, halfwidth, box, obstacles, diameter=3.):
    """Probe tip reaches the exposed notch side, not the neck's solid centre.

    A 0.05 mm separation is numerical clearance, not a virtual jaw allowance.
    No support other than the intended removed neck is ignored. The probe is
    a cylinder, explicitly not a verified representation of cutters/pliers.
    """
    target=np.asarray(target,float);side=np.asarray(side,float);side/=max(np.linalg.norm(side),1e-12)
    manager=trimesh.collision.CollisionManager()
    for i,m in enumerate(obstacles):
        if len(m.faces):manager.add_object(str(i),m)
    attempts=[]
    for sign in [-1,1]:
        end=target+sign*side*(diameter/2+halfwidth+.05)
        for offset in [np.zeros(2),np.array([3.,0]),np.array([-3.,0]),np.array([0,3.]),np.array([0,-3.])]:
            start=np.r_[end[:2]+offset,box.bounds[1,2]+5.]
            direction=end-start;length=float(np.linalg.norm(direction))
            if length<=0:continue
            matrix=trimesh.geometry.align_vectors([0,0,1],direction);matrix[:3,3]=(start+end)/2
            probe=trimesh.creation.cylinder(radius=diameter/2,height=length,sections=32,transform=matrix)
            hit=manager.in_collision_single(probe)
            attempts.append(dict(start=start.tolist(),end=end.tolist(),collision=bool(hit)))
            if not hit:return dict(passed=True,status='GEOMETRIC_TOOL_ACCESS_FOUND',diameter_mm=diameter,
                target=target.tolist(),start=start.tolist(),end=end.tolist(),side_clearance_mm=.05,attempts=attempts,
                claim='Tool access is a geometric proxy only.',physical_tool_validated=False)
    return dict(passed=False,status='CUT_POINT_TOOL_ACCESS_FAILED',diameter_mm=diameter,attempts=attempts,physical_tool_validated=False)


def run(record,inputs):
    began=time.monotonic();folder=Path(record['folder']);restore=np.array(record['orientation']['casting_restore_transform'])
    fingerprint=hashlib.sha256(json.dumps(dict(support=digest(folder/'custom_support_v3.stl'),cuts=digest(folder/'rib_geometry_cache.pkl'),
        vascular=digest(inputs['core']),box=digest(inputs['box']),restore=restore.tolist(),
        code={n:digest(Path(__file__).with_name(n)) for n in ['support_v3_removal.py','support_removal.py','support_v2_removal.py']}),sort_keys=True).encode()).hexdigest()
    cached=folder/'support_removal_plan.json'
    if cached.exists():
        previous=json.loads(cached.read_text())
        if previous.get('cache_fingerprint')==fingerprint:
            print('V3 REMOVAL CACHE',record['candidate_id'],flush=True);return previous
    with (folder/'rib_geometry_cache.pkl').open('rb') as f:ribs=pickle.load(f)
    support=trimesh.load_mesh(folder/'custom_support_v3.stl');cutters=[c['cutter'] for r in ribs for c in r['cuts']]
    remaining=trimesh.boolean.difference([support,union(cutters)],engine='manifold') if cutters else support.copy()
    feet=union([r['approved_retained_stubs'] for r in ribs]);feet.apply_transform(restore)
    v=trimesh.load_mesh(inputs['core']);box=trimesh.load_mesh(inputs['box']);inner=np.array(inputs['inner_bounds'])
    fullsupport=support.copy();fullsupport.apply_transform(restore)
    tool=[]
    for rib in ribs:
        for c in rib['cuts']:
            target=transform([c['origin']],restore)[0]
            idx=int(np.argmin(np.linalg.norm(rib['path_xy']-np.asarray(c['origin'])[:2],axis=1)))
            a=max(0,idx-1);b=min(len(rib['path_xy'])-1,idx+1);t= rib['path_xy'][b]-rib['path_xy'][a]
            normal=np.r_[-t[1],t[0],0.];normal=restore[:3,:3]@normal
            intended=c['intended'].copy();intended.apply_transform(restore)
            obstacle=trimesh.boolean.difference([fullsupport,intended],engine='manifold') if len(intended.faces) else fullsupport
            access=probe_access(target,normal,rib['width_mm']/2,box,[v,box,obstacle])
            tool.append(dict(cut_id=c['cut_id'],**access))
    table(folder/'tool_access.csv',tool);save(folder/'tool_access.json',tool)
    fragment_dir=folder/'removal_fragments';fragment_dir.mkdir(exist_ok=True)
    manager=trimesh.collision.CollisionManager();manager.add_object('vascular_ports',v);manager.add_object('casting_box',box)
    parts={};rows=[]
    for i,part in enumerate(remaining.split(only_watertight=False)):
        if part.volume<1e-6:continue
        part.apply_transform(restore);identity=f'F{i:04}';part.export(fragment_dir/f'{identity}.stl')
        residual=trimesh.boolean.difference([part,feet],engine='manifold')
        retained=bool(len(residual.faces)==0 or abs(float(residual.volume))<2e-4) and intersection_volume(part,box)>1e-6
        row=dict(fragment_id=identity,bounds_mm=part.bounds.tolist(),volume_mm3=float(part.volume),internal=cavity_overlap(part.bounds,inner),
            status='APPROVED_ROOT_FOOT_RETAINED' if retained else 'PENDING',path=None,attempts=[],frame='FROZEN_CASTING',export_closed=part.is_volume)
        parts[identity]=part;rows.append(row);manager.add_object(identity,part)
    pending=[r for r in rows if r['status']=='PENDING'];pending.sort(key=lambda r:(r['internal'],-r['bounds_mm'][1][2],r['volume_mm3']))
    sequence=[]
    for iteration in range(3):
        progress=False
        for row in pending:
            if row['status']!='PENDING':continue
            identity=row['fragment_id'];part=parts[identity];manager.remove_object(identity)
            local,paths=paths_for(part,box,not row['internal'])
            initial=paths[0][1][0];rise=box.bounds[1,2]-part.bounds[0,2]+5.
            # Release normal to the actual fragment's principal curtain axis,
            # in addition to global X/Y and tilt routes retained from V2.
            _,_,axes=np.linalg.svd(part.vertices[:,:2]-part.vertices[:,:2].mean(0),full_matrices=False)
            extra=[]
            for direction in [axes[1],-axes[1]]:
                for distance in [1.5,2.5,4.]:
                    moved=initial.copy();moved[:2,3]+=direction*distance;end=moved.copy();end[2,3]+=rise
                    extra.append((f'RIB_NORMAL_{direction.tolist()}_{distance}_THEN_TOP',[initial,moved,end]))
            # A short release from the root notch can allow a downward move
            # before passing a curved vessel, without ever teleporting through
            # the root foot. Every leg has the same continuous swept audit.
            for angle in np.arange(0,360,30):
                direction=np.array([np.cos(np.radians(angle)),np.sin(np.radians(angle))])
                for release in [.8,1.2]:
                    lateral=initial.copy();lateral[:2,3]+=direction*release
                    low=lateral.copy();low[2,3]-=1.5
                    moved=low.copy();moved[:2,3]+=direction*2.5
                    end=moved.copy();end[2,3]+=rise+1.5
                    extra.append((f'RELEASE_{angle}_{release}_LOWER_THEN_TOP',[initial,lateral,low,moved,end]))
            paths=extra+paths
            for name,poses in paths:
                check=swept_path(local,[],poses,static_manager=manager)
                row['attempts'].append(dict(iteration=iteration,path=name,**check))
                if check['passed']:
                    end=transform(local.vertices,poses[-1]);exited=bool(end[:,2].min()>box.bounds[1,2]) if row['internal'] else bool(np.any(end.max(0)<box.bounds[0]) or np.any(end.min(0)>box.bounds[1]))
                    if not exited:continue
                    row.update(status='REMOVAL_PATH_FOUND',path=dict(name=name,transforms=[p.tolist() for p in poses],full_top_exit=bool(end[:,2].min()>box.bounds[1,2])),sequence_index=len(sequence))
                    sequence.append(identity);progress=True;break
            if row['status']=='PENDING':manager.add_object(identity,part)
        if not progress:break
    for row in pending:
        if row['status']=='PENDING':row['status']='TRAPPED_FRAGMENT'
    result=dict(fragment_count=len(pending),internal_fragment_count=sum(r['internal'] for r in pending),retained_root_feet=len(rows)-len(pending),
        trapped_pieces=sum(r['status']=='TRAPPED_FRAGMENT' for r in pending),tool_access_passed=all(t['passed'] for t in tool),tool_access_failed=sum(not t['passed'] for t in tool),
        sequence=sequence,fragments=rows,cut_point_count=len(tool),other_fragments_are_obstacles=True,casting_opening_axis=[0,0,1],
        all_internal_exit_top=all(r['path'] and r['path']['full_top_exit'] for r in pending if r['internal']),
        physical_breakaway='UNVALIDATED',tool_access_claim='Tool access is a geometric proxy only.',runtime_seconds=time.monotonic()-began,cache_fingerprint=fingerprint)
    save(folder/'support_removal_plan.json',result);table(folder/'fragments.csv',rows);table(folder/'removal_paths.csv',[dict(fragment_id=r['fragment_id'],**a) for r in rows for a in r['attempts']])
    print('V3 REMOVAL',record['candidate_id'],len(pending),'fragments',result['trapped_pieces'],'trapped',result['tool_access_failed'],'inaccessible cuts',flush=True)
    return result
