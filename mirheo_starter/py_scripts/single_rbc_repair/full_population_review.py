"""Independent CPU membership checks of complete stored native fluid states."""
import argparse
import json
from pathlib import Path
import time

import h5py
import numpy as np

from py_scripts.fluid_physics.common import now,sha256_file,write_json
from py_scripts.single_rbc_benchmark.analysis import read_csv
from py_scripts.single_rbc_benchmark.physics import read_off
from .geometry_checks import checked_mesh,membership


def inspect(run, metrics, output):
    start=time.perf_counter();run=Path(run).resolve();spec=json.loads((run/'spec.json').read_text())
    c=spec['config'];L=c['geometry']['periodic_length'];band=c['repair']['containment_band']
    _,faces=read_off(spec['mesh']);screen=read_csv(metrics);native=run/'native_relaxation'
    initial_ids=None;frames=[];sources={Path(metrics).resolve(),run/'spec.json',Path(spec['mesh'])}
    for mesh in sorted(native.glob('rbc*.h5')):
        index=int(mesh.stem[3:]);step=index*spec['sample_steps'];saved=screen[screen['step']==step]
        assert len(saved)==1 and int(saved[0]['legacy_intersection_count'])==0, 'GEOMETRY_NOT_CREDIBLE_FOR_MEMBER_COUNT'
        files=[mesh]+[native/f'{pv}{index:05d}.h5' for pv in ('outer','inner')];states=[]
        for p in files:
            sources.add(p)
            with h5py.File(p) as f:states.append({k:f[k][:] for k in ('id','position')})
        v,g=checked_mesh(states[0]['id'].ravel(),states[0]['position'].astype(float),faces,L)
        assert g['closed'] and g['minimum_face_area']>0
        current_ids=[np.sort(s['id'].ravel()) for s in states[1:]]
        if initial_ids is None:initial_ids=current_ids
        same_ids=all(np.array_equal(x,y) for x,y in zip(current_ids,initial_ids))
        p=np.concatenate([s['position'] for s in states[1:]]).astype(float)
        ids=np.concatenate([s['id'].ravel() for s in states[1:]])
        assert len(np.unique(ids))==len(ids)
        inside=np.arange(len(p))>=len(states[1]['id'])
        p[:,:2]+=L*np.rint((v.mean(0)[:2]-p[:,:2])/L)
        near_box=np.all((p>=v.min(0)-band)&(p<=v.max(0)+band),axis=1)
        selected=np.flatnonzero(near_box|inside)
        check=membership(p[selected],inside[selected],v,faces,band)
        events=[]
        for j in np.flatnonzero(check['confirmed_mismatch']):
            i=int(selected[j]);tri=int(check['nearest_triangle'][j])
            events.append(dict(particle_id=int(ids[i]),tag='inner' if inside[i] else 'outer',position=p[i].tolist(),
                distance=float(check['distance'][j]),triangle=tri,triangle_ids=faces[tri].tolist(),triangle_coordinates=v[faces[tri]].tolist()))
        frames.append(dict(step=step,elapsed_time_star=step*spec['dt'],strain=0.,total_fluid=len(p),inner_particles=int(inside.sum()),
            same_ids_in_each_species=same_ids,tested_exact_points=len(selected),outer_points_proven_outside_AABB=int((~near_box&~inside).sum()),
            confirmed_mismatch_count=len(events),uncertain_count=int(check['uncertain'].sum()),ambiguous_rays=int(check['ambiguous'].sum()),events=events))
    first=next((x for x in frames if x['confirmed_mismatch_count']),None)
    result=dict(recorded_at=now(),task=run.name,scope='All registered fluid particles at stored native beforeForces frames only, excluding unsaved endpoint and unsaved intervening steps. Every inner point and outer point within the expanded mesh AABB is tested; remaining outer points are geometrically outside that closed mesh.',
        near_surface_band=band,frames=frames,first_confirmed_membership=first,
        total_confirmed_mismatched_point_frames=sum(x['confirmed_mismatch_count'] for x in frames),
        constant_species_ids=all(x['same_ids_in_each_species'] for x in frames),
        source_sha256={str(p):sha256_file(p) for p in sorted(sources)},cpu_analysis_s=time.perf_counter()-start,qualified_speedup=None)
    write_json(output,result)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',required=True);p.add_argument('--metrics',required=True);p.add_argument('--output',required=True);a=p.parse_args()
    r=inspect(a.run,a.metrics,a.output)
    print(json.dumps({k:v for k,v in r.items() if k not in ('source_sha256','frames')},ensure_ascii=False,indent=2))
