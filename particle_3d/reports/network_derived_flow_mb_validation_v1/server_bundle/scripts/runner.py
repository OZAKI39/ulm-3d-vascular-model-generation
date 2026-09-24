"""Input/initial-state adapter for frozen P9-A.1; no dynamics implementation."""
import argparse, gzip, hashlib, json, os, resource, socket, sys, time
from pathlib import Path
from copy import deepcopy
from dataclasses import replace
from types import SimpleNamespace
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import pyvista as pv

ENV=None
INITIAL=None
IDENTITY=None
ROOT=None

def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def safe(v):
    if isinstance(v,np.ndarray):return v.tolist()
    if isinstance(v,np.generic):return v.item()
    raise TypeError(type(v).__name__)
def dump(p,v):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(v,default=safe,indent=2,allow_nan=False)+'\n')

def make_environment(root,label):
    from particle_3d.audit import read_frozen
    from particle_3d.field import FrozenFEMField
    from particle_3d.wall_geometry import WallGeometry
    from particle_3d.validation_boundary import ValidationBoundaryClassifier
    from particle_3d.hydrodynamic_resistance import viscosity_from_frozen
    fem=root/'formal_3D_flow_solver/FEM_SimVascular'
    provenance,mesh,old,boundaries=read_frozen(fem)
    contract=json.loads((root/'data/old_new_flow_contract.json').read_text())
    path=root/'inputs'/f'{label}.vtu'
    assert sha(path)==contract['sha256'][label]
    flow=pv.read(path);field=FrozenFEMField.from_grids(mesh,flow)
    wall=WallGeometry.from_frozen(fem)
    classifier=ValidationBoundaryClassifier({k:v for k,v in boundaries.items() if k.startswith('OUTLET_')})
    mu,viscosity=viscosity_from_frozen(fem)
    return SimpleNamespace(provenance=provenance,mesh=mesh,field=field,boundaries=boundaries,wall=wall,classifier=classifier,mu=mu,viscosity=viscosity)

def mb_job(spec):
    event,folder=spec
    from particle_3d.particle9a_motion import Particle9AStepper
    from particle_3d.particle9a1_audit import instrument
    from particle_3d.particle82a_integration import integrate_admitted
    from particle_3d.particle6_stepper import bind_query_dependency
    from particle_3d.particle8_replay import canonical_hash
    rows=[];created=[];start=time.time();cpu=time.process_time()
    frozen=INITIAL[str(event['particle_id'])]
    def factory(particles,*args,**kwargs):
        assert len(particles)==1
        p=particles[0]
        assert np.array_equal(p.position,frozen['position_m']) and np.array_equal(p.q,frozen['q'])
        # Explicit immutable cohort state. No step, force, tolerance or geometry changes.
        p=replace(p,velocity=np.array(frozen['velocity_m_s']),omega=np.array(frozen['omega_s_inv']))
        stepper=Particle9AStepper([p],*args,gradient_provider=lambda x:ENV.field.sample(x).velocity_gradient_s_inv,**kwargs)
        created.append(stepper)
        return instrument(stepper,enabled=True,rows=rows)
    integration=bind_query_dependency(integrate_admitted,{'environment':lambda:ENV,'SavedTrajectoryStepper':factory})
    before=canonical_hash(event);result=integration(deepcopy(event),output=Path(folder))
    assert canonical_hash(event)==before and len(created)==1
    result.update(model='P9A1',dynamics='P9A1',paired_identity=IDENTITY,diagnostics_enabled=True,
        planar_statistics=created[0].planar_statistics,
        receipt=dict(hostname=socket.gethostname(),pid=os.getpid(),started_unix_s=start,ended_unix_s=time.time(),
                     cpu_seconds=time.process_time()-cpu,peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
    stem=Path(folder)/'trajectories'/f'mb_{event["particle_id"]:06d}'
    dump(stem.with_suffix('.json'),result)
    with gzip.open(stem.with_suffix('.audit.jsonl.gz'),'wt') as f:
        for row in rows:f.write(json.dumps(row,default=safe,allow_nan=False)+'\n')
    return dict(bubble_id=event['particle_id'],completed=result['completed'],outlet=result['exit_outlet'],
                end_reason=result['end_reason'],failure_detail=result['failure_detail'],wall_seconds=time.time()-start)

def point_job(spec):
    event,folder=spec
    row=ENV.native_point.trace(event['birth_center_m'],step_m=.2e-6,error=1e-11,horizon_m=2e-3)
    path=row.pop('path');stem=Path(folder)/f'point_{event["particle_id"]:06d}';stem.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(stem.with_suffix('.npz'),path=path)
    row.update(bubble_id=event['particle_id'],initial_center_m=event['birth_center_m'],last_time_s=float(path[-1,0]),
               path_sha256=sha(stem.with_suffix('.npz')),paired_identity=IDENTITY)
    dump(stem.with_suffix('.json'),row)
    return dict(bubble_id=event['particle_id'],outlet=row['outlet'],end_reason=row['end_reason'])

def main():
    global ENV,INITIAL,IDENTITY,ROOT
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--label',choices=['OLD','NEW'],default='NEW')
    p.add_argument('--mode',choices=['mb','point'],default='mb');p.add_argument('--workers',type=int,required=True)
    p.add_argument('--ids',default='all');p.add_argument('--output',required=True)
    a=p.parse_args();ROOT=a.root.resolve();sys.path.insert(0,str(ROOT/'particle_3d/src'))
    source=json.loads((ROOT/'data/code_snapshot_manifest.json').read_text())
    for relative,digest in source.items():assert sha(ROOT/relative)==digest,relative
    ENV=make_environment(ROOT,a.label)
    INITIAL=json.loads((ROOT/'data/paired_initial_states.json').read_text())
    events=json.loads((ROOT/'data/paired_events.json').read_text())
    if a.ids!='all':
        ids=set(map(int,a.ids.split(',')));events=[e for e in events if e['particle_id'] in ids];assert len(events)==len(ids)
    assert 1<=len(events)<=30 and a.workers<=6
    out=ROOT/'outputs'/a.output;assert not out.exists(),'Never overwrite previous experiment';out.mkdir(parents=True)
    start=time.time();IDENTITY=dict(flow_sha256=sha(ROOT/'inputs'/f'{a.label}.vtu'),flow_label=a.label,
        cohort_manifest_sha256=sha(ROOT/'data/paired_cohort_manifest.csv'),events_sha256=sha(ROOT/'data/paired_events.json'),
        initial_states_sha256=sha(ROOT/'data/paired_initial_states.json'),source_snapshot_sha256=sha(ROOT/'data/code_snapshot_manifest.json'),
        runner_sha256=sha(__file__),frozen_science=True)
    config=dict(identity=IDENTITY,workers=a.workers,ids=[e['particle_id'] for e in events],mode=a.mode,started_unix_s=start,
        dt_s=.00025,horizon_s=1.5,initial_state_policy='Exact saved OLD velocities/omega and full state, then unchanged P9-A.1 integration')
    dump(out/'config.json',config)
    if a.mode=='point':
        from particle_3d.particle82_point_native import NativePointTracer
        ENV.native_point=NativePointTracer(ENV)
    results=[]
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        tasks=[pool.submit(mb_job if a.mode=='mb' else point_job,(e,str(out))) for e in events]
        for task in as_completed(tasks):
            row=task.result();results.append(row);print(json.dumps(row),flush=True)
    dump(out/'completed.json',dict(config,results=sorted(results,key=lambda r:r['bubble_id']),wall_seconds=time.time()-start,completed_unix_s=time.time()))
if __name__=='__main__':main()
