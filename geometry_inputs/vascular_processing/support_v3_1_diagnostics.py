"""Locate real colliders for each baseline trapped fragment and failed probe."""
import json,pickle
from pathlib import Path
import numpy as np
import trimesh
from .support_inputs import save,table,transform


def internal_root_statistics(record,inputs):
    from .support_v3_1_validation import solid_components
    from .support_custom_qc import intersection_volume
    folder=Path(record['folder']);mesh=trimesh.load_mesh(folder/'custom_support_v3.stl')
    bounds=np.array(inputs['inner_bounds']);cavity=trimesh.creation.box(bounds[1]-bounds[0],transform=trimesh.transformations.translation_matrix(bounds.mean(0)))
    cavity.apply_transform(np.array(record['orientation']['print_transform']))
    rows=[]
    for i,part in enumerate(solid_components(mesh)):
        volume=intersection_volume(part,cavity)
        if volume<=1e-6:continue
        rows.append(dict(component_id=i,root_domain='BUILD_PLATE' if part.bounds[0,2]<.001 else 'NONPLATE',inside_cavity_volume_mm3=volume))
    plate=[x for x in rows if x['root_domain']=='BUILD_PLATE'];total=sum(x['inside_cavity_volume_mm3'] for x in rows)
    result=dict(internal_components=len(rows),internal_plate_components=len(plate),internal_nonplate_components=len(rows)-len(plate),
                component_fraction=len(plate)/len(rows) if rows else None,
                cavity_volume_fraction=sum(x['inside_cavity_volume_mm3'] for x in plate)/total if total else None,
                method='Actual exported connected material components, real cavity-volume intersection, minimum print Z for bed root',rows=rows)
    record['internal_root_statistics']=result;save(folder/'internal_root_statistics.json',result);return result


def baseline_colliders(folder,inputs,out):
    record=json.loads((folder/'candidate.json').read_text());restore=np.array(record['orientation']['casting_restore_transform'])
    v=trimesh.load_mesh(inputs['core']);box=trimesh.load_mesh(inputs['box'])
    support=trimesh.load_mesh(folder/'custom_support_v3.stl');support.apply_transform(restore)
    rem=json.loads((folder/'support_removal_plan.json').read_text());fragments=[]
    for row in rem['fragments']:
        if row['status']!='TRAPPED_FRAGMENT':continue
        mesh=trimesh.load_mesh(folder/'removal_fragments'/f"{row['fragment_id']}.stl")
        manager=trimesh.collision.CollisionManager();manager.add_object('FROZEN_VASCULAR_AND_PORTS',v);manager.add_object('FROZEN_BOX',box)
        for other in rem['fragments']:
            if other['fragment_id']==row['fragment_id']:continue
            path=folder/'removal_fragments'/f"{other['fragment_id']}.stl"
            if path.exists():manager.add_object('SUPPORT_'+other['fragment_id'],trimesh.load_mesh(path))
        collision_samples=[]
        for rise in [0.,.05,.1,.2,.5,1.,2.,4.,8.,20.,40.]:
            pose=np.eye(4);pose[2,3]=rise
            hit,names=manager.in_collision_single(mesh,transform=pose,return_names=True)
            if hit:collision_samples.append(dict(upward_mm=rise,colliders=sorted(names)))
        fragments.append(dict(fragment_id=row['fragment_id'],bounds_mm=row['bounds_mm'],volume_mm3=row['volume_mm3'],
            collision_localization=collision_samples,full_swept_attempts=len(row['attempts']),
            scope='Discrete samples locate obstacles only; original bounded FCL sweeps determine path pass/fail'))
    with (folder/'rib_geometry_cache.pkl').open('rb') as f:ribs=pickle.load(f)
    cuts={c['cut_id']:c for r in ribs for c in r['cuts']};tools=json.loads((folder/'tool_access.json').read_text());toolrows=[]
    for t in tools:
        if t['passed']:continue
        cut=cuts[t['cut_id']];neck=cut['intended'].copy();neck.apply_transform(restore)
        remaining=trimesh.boolean.difference([support,neck],engine='manifold')
        manager=trimesh.collision.CollisionManager();manager.add_object('FROZEN_VASCULAR_AND_PORTS',v);manager.add_object('FROZEN_BOX',box);manager.add_object('UNCUT_SUPPORT',remaining)
        attempts=[]
        for a in t['attempts']:
            start,end=np.array(a['start']),np.array(a['end']);delta=end-start
            matrix=trimesh.geometry.align_vectors([0,0,1],delta);matrix[:3,3]=(start+end)/2
            proxy=trimesh.creation.cylinder(radius=1.5,height=np.linalg.norm(delta),sections=32,transform=matrix)
            hit,names=manager.in_collision_single(proxy,return_names=True)
            attempts.append(dict(start=a['start'],end=a['end'],colliders=sorted(names),collision=bool(hit)))
        toolrows.append(dict(cut_id=t['cut_id'],origin_casting=transform([cut['origin']],restore)[0].tolist(),attempts=attempts))
    save(out/'baseline_fragment_collision_localization.json',fragments);table(out/'baseline_fragment_collision_localization.csv',fragments)
    save(out/'baseline_cut_collision_localization.json',toolrows);table(out/'baseline_cut_collision_localization.csv',toolrows)
    return dict(trapped_fragments=len(fragments),inaccessible_cuts=len(toolrows))
