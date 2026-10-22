"""Individual curtain release, with actual neighbouring segments as obstacles."""
from __future__ import annotations
import pickle,time
from pathlib import Path
import numpy as np
import trimesh
from scipy.spatial.transform import Rotation
from .support_inputs import save,table,transform
from .support_rfus_geometry import disconnect_base
from .support_rfus_geometry import union,boolean
from .support_v3_removal import probe_access
from .support_removal import swept_path
from .support_custom_qc import intersection_volume


def piece_topology(segment):
    """A group of disconnected slivers is not one graspable curtain."""
    part=disconnect_base(segment)
    pieces=part.split(only_watertight=False)
    positive=[p for p in pieces if p.volume>1e-8]
    return dict(positive_components=len(positive),closed=bool(part.is_volume),
                removable_as_one_component=bool(part.is_volume and len(positive)==1),
                component_volumes_mm3=[float(p.volume) for p in positive],
                numerical_zero_volume_mm3=1e-8)


def topology_audit(record):
    """Supplement movement results with material connectivity and contact."""
    folder=Path(record['folder'])
    with (folder/'segment_geometry.pkl').open('rb') as f:segments=pickle.load(f)
    rows=[dict(segment_id=s['segment_id'],**piece_topology(s)) for s in segments]
    contacts=[]
    for i,s in enumerate(segments):
        manager=trimesh.collision.CollisionManager();manager.add_object('first',s['mesh'])
        for t in segments[i+1:]:
            if np.any(s['mesh'].bounds[1]+1e-4<t['mesh'].bounds[0]) or np.any(t['mesh'].bounds[1]+1e-4<s['mesh'].bounds[0]):continue
            d=float(manager.min_distance_single(t['mesh']))
            if d<=1e-6:contacts.append(dict(first=s['segment_id'],second=t['segment_id'],distance_mm=d))
    audit=dict(all_independent=all(r['removable_as_one_component'] for r in rows) and not contacts,
               disconnected_or_invalid_segments=sum(not r['removable_as_one_component'] for r in rows),
               touching_segment_pairs=contacts,rows=rows,
               method='Post-cut material components and FCL nominal segment separation; no physical cohesion claim')
    save(folder/'segment_topology.json',audit);table(folder/'segment_topology.csv',rows)
    record['segment_topology']=audit
    return record


def paths_for(segment, part, top):
    pivot=part.bounds.mean(0);local=part.copy();local.apply_translation(-pivot)
    initial=np.eye(4);initial[:3,3]=pivot;rise=float(top-part.bounds[0,2]+5)
    up=initial.copy();up[2,3]+=rise;paths=[('PLUS_Z',[initial,up])]
    n=np.mean(segment['lateral'],axis=0);n/=max(np.linalg.norm(n),1e-12)
    for sign in [1,-1]:
        for distance in [1.,1.5,2.,3.]:
            lateral=initial.copy();lateral[:2,3]+=sign*n*distance;end=lateral.copy();end[2,3]+=rise
            paths.append((f'LATERAL_{sign*distance:+g}_THEN_TOP',[initial,lateral,end]))
    for axis in ['x','y']:
        for angle in [10,-10,20,-20]:
            tilted=initial.copy();tilted[:3,:3]=Rotation.from_euler(axis,angle,degrees=True).as_matrix()
            # Keep the rotated base above the original cut plane; never rotate
            # through the box bottom to manufacture an apparent path.
            tiltmin=transform(local.vertices,tilted)[:,2].min();tilted[2,3]+=max(0,part.bounds[0,2]-tiltmin)
            end=tilted.copy();end[2,3]+=rise
            paths.append((f'TILT_{axis}_{angle:+g}_THEN_TOP',[initial,tilted,end]))
            for sign in [1,-1]:
                release=tilted.copy();release[:2,3]+=sign*n*2.;end=release.copy();end[2,3]+=rise
                paths.append((f'TILT_{axis}_{angle:+g}_LATERAL_{sign*2}_TOP',[initial,tilted,release,end]))
    return local,paths


def run(context,record):
    began=time.monotonic();folder=Path(record['folder']);c=context
    with (folder/'segment_geometry.pkl').open('rb') as f:segments=pickle.load(f)
    v,box=c['vascular'],c['box'];manager=trimesh.collision.CollisionManager()
    manager.add_object('vascular_and_ports',v);manager.add_object('frozen_box',box)
    for s in segments:manager.add_object(s['segment_id'],s['mesh'])
    rows=[];parts={};retained={};sequence=[];tools=[]
    for s in segments:
        part=disconnect_base(s);parts[s['segment_id']]=part
        # Leave a finite cut kerf: lower retained material stops 0.03 mm below
        # the intended plane, removable material starts 0.03 mm above it.
        upper=s.copy();upper['neck_cut_z']-=.06
        retained[s['segment_id']]=boolean([s['mesh'],disconnect_base(upper)],'difference')
        rows.append(dict(segment_id=s['segment_id'],internal=s['internal'],branch_id=s['branch_id'],arc_start_mm=s['arc_start_mm'],arc_end_mm=s['arc_end_mm'],status='PENDING',handle_accessible=False,base_accessible=False,attempts=[],removal_path=None,physical_force_validated=False))
    order=sorted(range(len(segments)),key=lambda i:(segments[i]['internal'],-segments[i]['mesh'].bounds[1,2],segments[i]['arc_length_mm']))
    for wave in range(len(segments)+1):
        progress=False
        for i in order:
            s=segments[i];r=rows[i];sid=s['segment_id']
            if r['status']!='PENDING':continue
            other=[v,box]+[segments[j]['mesh'] if rows[j]['status']!='REMOVAL_PATH_FOUND' else retained[segments[j]['segment_id']] for j in range(len(segments)) if j!=i]
            side=np.r_[np.mean(s['lateral'],axis=0),0.]
            handle=None
            if s.get('handle_target'):
                handle=probe_access(s['handle_target'],side,s['handle_width_mm']/2,box,other+[s['mesh']],c['cfg']['probe_diameter_mm'])
            basepoint=np.r_[s['path_xy'][len(s['path_xy'])//2],s['neck_cut_z']]
            base=probe_access(basepoint,side,s['neck_width_mm']/2,box,other+[s['mesh']],c['cfg']['probe_diameter_mm']) if not handle or not handle['passed'] else dict(passed=False,status='HANDLE_STRATEGY_SELECTED')
            r.update(handle_accessible=bool(handle and handle['passed']),base_accessible=base['passed'],tool_strategy='HANDLE_CUT_PEEL' if handle and handle['passed'] else 'BASE_CUT' if base['passed'] else 'NO_ACCESSIBLE_STRATEGY')
            tools.append(dict(segment_id=sid,wave=wave,handle=handle,base=base,geometric_proxy=True))
            manager.remove_object(sid)
            foot=retained[sid]
            if len(foot.faces):manager.add_object(sid+'_retained',foot)
            part=parts[sid];local,paths=paths_for(s,part,box.bounds[1,2])
            for name,poses in paths:
                check=swept_path(local,[],poses,static_manager=manager)
                r['attempts'].append(dict(wave=wave,path=name,**check))
                if check['passed'] and transform(local.vertices,poses[-1])[:,2].min()>box.bounds[1,2]:
                    r['removal_path']=dict(name=name,transforms=[p.tolist() for p in poses],full_top_exit=True)
                    if r['handle_accessible'] or r['base_accessible'] or not r['internal']:
                        r.update(status='REMOVAL_PATH_FOUND',sequence_index=len(sequence));sequence.append(sid);progress=True
                    break
            if r['status']=='PENDING':
                if len(foot.faces):manager.remove_object(sid+'_retained')
                manager.add_object(sid,s['mesh'])
        if not progress:break
    for r in rows:
        if r['status']=='PENDING':r['status']='INACCESSIBLE_SEGMENT' if r['removal_path'] else 'TRAPPED_SEGMENT'
    # Nominally separated curtains which fused in the exported union cannot
    # be treated as independent removable units. Report this as its own gate.
    support=trimesh.load_mesh(folder/'geometry/custom_full_underbody_support.stl')
    fused=[]
    for i,s in enumerate(segments):
        for t in segments[i+1:]:
            if np.any(s['mesh'].bounds[1]<t['mesh'].bounds[0]) or np.any(t['mesh'].bounds[1]<s['mesh'].bounds[0]):continue
            overlap=intersection_volume(s['mesh'],t['mesh'])
            if overlap>1e-5:fused.append(dict(first=s['segment_id'],second=t['segment_id'],intersection_mm3=overlap))
    result=dict(segment_count=len(rows),trapped_segments=sum(r['status']=='TRAPPED_SEGMENT' for r in rows),inaccessible_segments=sum(r['internal'] and not(r['handle_accessible'] or r['base_accessible']) for r in rows),accessible_handles=sum(r['handle_accessible'] for r in rows),internal_segments=sum(r['internal'] for r in rows),sequence=sequence,segments=rows,fused_segments=fused,all_segments_removable=all(r['status']=='REMOVAL_PATH_FOUND' for r in rows) and not fused,remaining_segments_are_obstacles=True,cuttable_base_neck_physically_calibrated=False,tool_proxy_diameter_mm=c['cfg']['probe_diameter_mm'],runtime_seconds=time.monotonic()-began)
    save(folder/'segment_removal.json',result);table(folder/'segment_removal.csv',rows);save(folder/'handle_accessibility.json',tools)
    print('RFUS REMOVAL',record['candidate_id'],'trapped',result['trapped_segments'],'inaccessible',result['inaccessible_segments'],flush=True)
    return result
