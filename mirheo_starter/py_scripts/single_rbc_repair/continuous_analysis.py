"""Convert native, same-step HDF5 observations to the existing CPU review format.

Raw solver output is read-only. Derived CSV/NPZ files carry source hashes and
the original execution record; their creation is not a new solver repeat.
"""
import argparse
import csv
import json
from pathlib import Path
import shutil
import time

import h5py
import numpy as np

from py_scripts.fluid_physics.common import now, sha256_file, write_json
from .runtime_analysis import inspect_run


def convert(run, output):
    start=time.perf_counter();run=Path(run).resolve();output=Path(output).resolve()
    derived=output/'derived_observations';derived.mkdir(parents=True,exist_ok=False)
    spec=json.loads((run/'spec.json').read_text());c=spec['config'];pad=c['geometry']['wall_layer'];H=c['geometry']['gap']
    assert spec['continuous_observation'] and spec['steps']==0
    sources={run/'spec.json',run/'execution.json'}
    for name in ('spec.json','execution.json','completion.json'):
        if (run/name).exists():shutil.copyfile(run/name,derived/name);sources.add(run/name)
    native=run/'native_relaxation';chosen=None;populations=[]
    def read(p):
        sources.add(p)
        with h5py.File(p) as f:
            ids=f['id'][:].reshape(-1);order=np.argsort(ids)
            assert len(np.unique(ids))==len(ids)
            result={k:np.asarray(v[:])[order] for k,v in f.items()}
        result['id']=ids[order]
        return result
    names=['vertices.csv','moments.csv']
    handles=[(derived/name).open('x') for name in names]
    try:
        vertex,moment=map(csv.writer,handles)
        vertex.writerow(['step','phase','time_star','strain','vertex','x','y','z','fx','fy','fz'])
        moment.writerow(['step','phase','time_star','strain','N','inner_N','temperature','max_speed','wall_crossings'])
        for mesh_file in sorted(native.glob('rbc*.h5')):
            index=int(mesh_file.stem[3:]);step=index*spec['sample_steps']
            assert step<spec['prep_steps']
            rbc=read(mesh_file);states=[read(native/f'{name}{index:05d}.h5') for name in ('outer','inner')]
            assert np.array_equal(rbc['id'],np.arange(len(rbc['id'])))
            if chosen is None:
                chosen=[s['id'][np.arange(0,len(s['id']),max(1,len(s['id'])//96))[:96]].copy() for s in states]
            p=np.concatenate([s['position'] for s in states]).astype(float);vel=np.concatenate([s['velocity'] for s in states]).astype(float)
            ids=np.concatenate([s['id'] for s in states]);assert len(np.unique(ids))==len(ids)
            probes=[]
            for tag,(st,selection) in enumerate(zip(states,chosen)):
                take=np.searchsorted(st['id'],selection)
                assert np.all(take<len(st['id'])) and np.array_equal(st['id'][take],selection), 'PROBE_SPECIES_CHANGED'
                probes.extend([[tag,int(st['id'][i]),*(st['position'][i].astype(float)-[0,0,pad])] for i in take])
            v=rbc['position'].astype(float)-[0,0,pad];force=rbc['forces'][:,:3].astype(float)
            t=step*spec['dt']-c['protocol']['relaxation_time']
            for i,(x,f) in enumerate(zip(v,force)):vertex.writerow([step,'relaxation',t,0,i,*x,*f])
            bins=round(H);bin_id=np.clip(((p[:,2]-pad)/H*bins).astype(int),0,bins-1)
            n=np.bincount(bin_id,minlength=bins);means=np.array([np.bincount(bin_id,weights=vel[:,k],minlength=bins)/np.maximum(n,1) for k in range(3)]).T
            temperature=float(np.sum((vel-means[bin_id])**2)*c['dpd']['mass']/(3*max(1,len(vel)-np.count_nonzero(n))))
            crossings=int(np.count_nonzero((p[:,2]<pad-1e-5)|(p[:,2]>H+pad+1e-5)))
            moment.writerow([step,'relaxation',t,0,len(ids),len(states[1]['id']),temperature,float(np.linalg.norm(vel,axis=1).max()),crossings])
            populations.append(dict(step=step,outer_N=len(states[0]['id']),inner_N=len(states[1]['id'])))
            np.savez_compressed(derived/f'fluid_probes_relaxation_{step:08d}.npz',probes=np.array(probes),strain=0,phase='relaxation',global_step=step,membrane_vertex_ids=rbc['id'])
            np.savez_compressed(derived/f'membrane_state_relaxation_{step:08d}.npz',ids=rbc['id'],positions=v,velocities=rbc['velocity'].astype(float),phase='relaxation',global_step=step,strain=0)
        # The native beforeForces callback omits the final endpoint. The worker
        # captured that endpoint after its single run; append the genuine files.
        for writer,name in zip((vertex,moment),names):
            p=run/name;sources.add(p)
            with p.open() as f:
                rows=list(csv.reader(f));writer.writerows(rows[1:])
        for pattern in ('fluid_probes_*.npz','membrane_state_*.npz','*00000.log'):
            for p in run.glob(pattern):shutil.copyfile(p,derived/p.name);sources.add(p)
    finally:
        for f in handles:f.close()
    record=dict(recorded_at=now(),raw_run=str(run),derived_directory=str(derived),derived_only=True,new_solver_repeat=False,
                observation_phase='Positions/velocities at beforeForces of indexed step; native saved force channel is lagged and not used for same-step force causality.',
                constant_inner_N=len({p['inner_N'] for p in populations})==1,populations=populations,
                source_sha256={str(p):sha256_file(p) for p in sorted(sources)},cpu_conversion_s=time.perf_counter()-start)
    write_json(output/'conversion.json',record)
    return inspect_run(derived,output/'review')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',required=True);p.add_argument('--output',required=True);a=p.parse_args()
    result=convert(a.run,a.output)
    print(json.dumps({k:v for k,v in result.items() if k not in ('execution','source_sha256')},ensure_ascii=False,indent=2))
