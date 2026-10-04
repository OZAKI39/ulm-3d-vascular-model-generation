"""Independent deterministic replay with unchanged numerical kernels and frozen marks.

One engine contains exactly one particle, so even CPU fallback cannot introduce
pair interactions. Numerical GPU dispatch can still be batched across engines.
"""
from pathlib import Path
import argparse,json,hashlib,sys,time,inspect,traceback
import numpy as np
import torch
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,x):
    tmp=Path(str(p)+'.pending');tmp.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n');tmp.replace(p)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--base',type=Path,required=True);ap.add_argument('--candidate',required=True);ap.add_argument('--dt',type=float,default=.0005);ap.add_argument('--validate-dispatch',action='store_true');a=ap.parse_args()
    b=a.base;s=b/'mesh_particle_validation_v2';sys.path.insert(0,str(b/'population_v1/source/particle_3d/src'))
    from particle_3d.population_v1.environment import load_environment
    from particle_3d.population_v1.gpu_engine import GPUPopulationEngine
    from particle_3d.population_v1.gpu_batch import prefetch
    from particle_3d.population_v1.output import PopulationOutput
    from particle_3d.population_v1.fixed_window import censor_remaining
    torch.set_num_threads(1)
    co=json.loads((s/'00_baseline/TRAJECTORY_COHORT.json').read_text());spec=json.loads((b/'population_v1/data/environment.json').read_text())
    if a.candidate!='M0':
        case=s/'03_cfd_runs'/a.candidate;assert json.loads((case/'postprocess/STUDY_FLOW.json').read_text())['status']=='PASS'
        spec['mesh']=str(case/'SV_MESH/mesh-complete.mesh.vtu');spec['flow']=str(case/'frozen_flow/steady_flow.vtu')
        spec['boundaries']={k:str(case/'SV_MESH/mesh-surfaces'/(k+'.vtp')) for k in spec['boundaries']}
        spec['sha256']={k:sha(p) for k,p in dict(mesh=spec['mesh'],flow=spec['flow'],**spec['boundaries']).items()}
    assert a.dt==.0005 or (s/'05_particle_replay/GATE_MESH_CONVERGENCE.json').exists(),'No timestep study before mesh gate'
    if a.dt!=.0005:
        gate=json.loads((s/'05_particle_replay/GATE_MESH_CONVERGENCE.json').read_text())
        assert gate['status']=='PASS' and gate['selected_candidate']==a.candidate,'Timestep replay must use selected refined candidate'
        assert gate['flow_sha256']==spec['sha256']['flow'] and gate['cohort_sha256']==sha(s/'00_baseline/TRAJECTORY_COHORT.json'),'Selected frozen field/cohort changed'
    folder=s/('05_particle_replay' if a.dt==.0005 else '06_dt_convergence')/(a.candidate+'_dt'+format(a.dt,'.8f')+('_dispatch_validation' if a.validate_dispatch else ''))
    import fcntl
    resource_slot=(s/'HEAVY_COMPUTE.lock').open('a')
    print(a.candidate,'WAIT_HEAVY_COMPUTE_SLOT',flush=True)
    fcntl.flock(resource_slot,fcntl.LOCK_EX)
    folder.mkdir(exist_ok=False);dump(folder/'ENVIRONMENT.json',spec)
    (folder/'executed_replay.py').write_bytes(Path(__file__).read_bytes())
    env=load_environment(spec);items=[]
    for row in co['particles']:
        event=dict(species='MB',particle_id=row['particle_id'],radius_m=row['radius_m'],q=row['q'])
        particle,status,detail=env.checker.check(event,np.array(row['position_m']),{})
        if particle is None:
            dump(folder/'GATE_ADMISSION.json',dict(status='FAIL',particle_id=row['particle_id'],reason=status,detail=detail));return
        items.append((row,particle,detail))
    dump(folder/'GATE_ADMISSION.json',dict(status='PASS',count=len(items),original_positions_unchanged=True,source_redrawn=False))
    # Retain production prefetch numerical implementation verbatim. Its only
    # altered line is the realization-ID uniqueness guard: frozen original ID 0
    # is retained on all independent singleton engines. Explicit singleton and
    # unique particle-ID assertions replace that bookkeeping guard.
    src=inspect.getsource(prefetch);guard="    if len({e.realization_id for e in live})!=len(live):raise ValueError('Distinct realization IDs required')"
    assert src.count(guard)==1
    adapted=src.replace(guard,"    assert all(e.pool.n==1 for e in live) and len({next(iter(e.births)) for e in live})==len(live)")
    ns=dict(prefetch.__globals__);exec(compile(adapted,'<study singleton prefetch>','exec'),ns);dispatch=ns['prefetch']
    (folder/'executed_prefetch.py').write_text(adapted)
    code={str(p.relative_to(b/'population_v1/source/particle_3d/src')):sha(p) for p in (b/'population_v1/source/particle_3d/src/particle_3d').rglob('*.py')}
    identity=dict(schema='FROZEN_SINGLETON_COHORT_REPLAY_V2',candidate=a.candidate,cohort_sha256=sha(s/'00_baseline/TRAJECTORY_COHORT.json'),
        environment=spec,dt_s=a.dt,observation_end_s=co['absolute_observation_end_s'],pairs=False,one_particle_per_engine=True,
        source_redrawn=False,source_contract_sha256=co['source_contract_sha256'],original_realization_id_preserved=True,
        production_prefetch_sha256=hashlib.sha256(src.encode()).hexdigest(),adapted_prefetch_sha256=hashlib.sha256(adapted.encode()).hexdigest(),
        source_files_sha256=code,replay_source_sha256=sha(__file__),GPU=torch.cuda.get_device_name(),torch_version=torch.__version__,CUDA=torch.version.cuda)
    dump(folder/'IDENTITY.json',identity)
    engines=[];outputs=[];gpu_field=None;gpu_wall=None
    for row,particle,detail in items:
        pid=row['particle_id'];out=None if a.validate_dispatch else PopulationOutput(folder/('particle_%04d'%pid),dict(identity,particle_id=pid))
        e=GPUPopulationEngine(env,dt=a.dt,pairs=False,realization_id=0,recorder=out,gpu_field=gpu_field,gpu_wall=gpu_wall,cache_path=s/'05_particle_replay'/('GPU_FIELD_'+a.candidate+'.npz'))
        gpu_field=e.field;gpu_wall=e.wall;e.time=row['birth_time_s'];e.counts['source_event_count']=1
        event=dict(row,status='ACCEPTED_BIRTH',clearance_m=detail.get('wall_gap_m',detail.get('wall_gap_lower_bound_m',0.)),admission_detail=detail)
        e.insert_accepted(particle,event)
        if out:out.source(event)
        engines.append(e);outputs.append(out)
    def target(e):
        ratio=e.time/e.dt;tick=(int(np.floor(ratio+8*abs(np.spacing(ratio))))+1)*e.dt
        assert tick>e.time
        return min(tick,co['absolute_observation_end_s'])
    if a.validate_dispatch:
        references=[]
        for row,particle,detail in items:
            e=GPUPopulationEngine(env,dt=a.dt,pairs=False,realization_id=0,gpu_field=gpu_field,gpu_wall=gpu_wall);e.time=row['birth_time_s'];e.counts['source_event_count']=1
            e.insert_accepted(particle,dict(row,clearance_m=detail.get('wall_gap_m',detail.get('wall_gap_lower_bound_m',0.))));references.append(e)
        errors=[]
        for iteration in range(8):
            dispatch(engines)
            for e,r in zip(engines,references):
                end=target(e);e.advance(end);r.advance(end);assert e.pool.n==r.pool.n
                err=float(torch.max(abs(e.pool.x[:e.pool.n]-r.pool.x[:r.pool.n])));errors.append(err);assert err<1e-12
        dump(folder/'GATE_DISPATCH.json',dict(status='PASS',singleton_count=len(items),nominal_steps_per_particle=8,max_position_difference_m=max(errors),numerical_prefetch_unchanged=True));return
    validation=s/'05_particle_replay'/('M0_dt0.00050000_dispatch_validation/GATE_DISPATCH.json');assert json.loads(validation.read_text())['status']=='PASS'
    start=time.time();last=0.;iterations=0;completed=set();problems=[]
    while True:
        active=[e for e in engines if e.pool.n and e.time<co['absolute_observation_end_s']]
        if not active:break
        dispatch(active)
        for e in active:
            try:e.advance(target(e))
            except Exception as error:
                # This is an explicit new numerical failure, never silently discard.
                pid=next(iter(e.births));row=e.births[pid];terminal=dict(particle_id=pid,source_event_id=row['source_event_id'],realization_id=0,birth_time_s=row['birth_time_s'],physical_time_s=e.time,residence_time_s=e.time-row['birth_time_s'],status='SOLVER_FAILURE',outlet=None,terminal_reason='UNHANDLED_REPLAY_INTEGRATION_EXCEPTION',error=str(error),traceback=traceback.format_exc())
                e.terminals.append(terminal);e.recorder.terminal(terminal);e.counts['solver_failed_count']+=1;e.pool.n=0;problems.append(terminal)
        iterations+=1
        for i,(e,out) in enumerate(zip(engines,outputs)):
            if e.pool.n and e.time==co['absolute_observation_end_s']:censor_remaining(e,e.time)
            if not e.pool.n and i not in completed:
                e.flush_output();out.close(dict(accounting=e.accounting(),terminal=e.terminals[0],backend_counts=dict(e.backend_counts)));completed.add(i)
        if time.time()-last>15:
            dump(folder/'PROGRESS.json',dict(elapsed_s=time.time()-start,iterations=iterations,completed=len(completed),total=len(engines),min_active_time_s=min((e.time for e in engines if e.pool.n),default=None),max_active_time_s=max((e.time for e in engines if e.pool.n),default=None)));last=time.time()
    terminals=[e.terminals[0] for e in engines];counts={k:sum(t['status']==k for t in terminals) for k in ['EXITED','SOLVER_FAILURE','SUPPORTED_STATIONARY','RIGHT_CENSORED_OBSERVATION_END']}
    assert len(completed)==46 and sum(counts.values())==46
    dump(folder/'SUMMARY.json',dict(status='COMPLETE_WITH_NUMERICAL_FAILURES' if counts['SOLVER_FAILURE'] else 'COMPLETE',counts=counts,terminals=terminals,wall_time_s=time.time()-start,cohort_count=46,identity_sha256=sha(folder/'IDENTITY.json'),unhandled_exceptions=problems))
    print('REPLAY_COMPLETE',a.candidate,counts,flush=True)
if __name__=='__main__':main()
