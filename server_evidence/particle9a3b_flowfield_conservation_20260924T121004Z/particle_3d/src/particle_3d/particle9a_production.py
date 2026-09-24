"""One fresh method-B cohort, P7 flux clock, P8.2A integration, P9-A stepper.

Smoke is the first 12 IDs of the same 500-track formal dataset. Successful
smoke trajectories may resume only under exactly identical source/event hashes.
"""
from pathlib import Path
from functools import partial
from collections import Counter
import argparse,json,multiprocessing as mp,os,resource,socket,time
import numpy as np
from .particle9a_provenance import (require_current_flow,source_identity,sha256,
    REPO,MODEL,ADMISSION_VERSION)
from .particle82a_admission import context,common_event,method_b
from .particle82a_integration import integrate_admitted
from .particle9a_motion import Particle9AStepper
from .particle6_stepper import bind_query_dependency
from .particle81_simulation import environment,dump,DT
from .particle8_replay import canonical_hash
from .injection_population import FluxClock,LinearProfile,ConstantMBConcentrationV0,C_MB

SEED=202609249
OUTPUT=REPO/'particle_3d/outputs/particle9a_2mmps'
_JOB=None


def production_provenance():
    flow=require_current_flow();c=context()
    return dict(**flow,**source_identity(),dynamics=MODEL,admission_model_version=ADMISSION_VERSION,
        sonovue_histogram_sha256=c.contract['histogram_sha256'],
        sonovue_sampler_sha256=c.contract['sampler_source_sha256'],master_seed=SEED)


def prepare_cohort(count=500,output=OUTPUT):
    output=Path(output);p=production_provenance();c=context()
    clock=FluxClock(LinearProfile([0.],[c.env.sampler.Q_m3_s]),ConstantMBConcentrationV0())
    target=output/'admission/birth_ledger.json'
    if target.exists():
        ledger=json.loads(target.read_text())
        if ledger['production_provenance']!=p or len(ledger['events'])!=count:
            raise ValueError('Fresh-cohort identity mismatch; no old birth/cache fallback')
        return ledger
    accepted=[];rows=[];pid=0
    while len(accepted)<count:
        pid+=1
        if pid>count*10:raise ValueError('Admission fraction unexpectedly low')
        e=common_event(pid,SEED,c);admit=method_b(e,c)
        birth=clock.time_at(pid)
        rows.append(dict(common_event=e,scheduled_time_s=birth,admission=admit))
        if not admit['accepted']:continue
        q=np.random.default_rng([SEED,pid,3]).normal(size=4);q/=np.linalg.norm(q)
        event=dict(particle_id=pid,species='MB',attempt_count=0,
            radius_m=admit['radius_m'],diameter_um=admit['diameter_um'],q=q.tolist(),
            position_seed=[SEED,pid,1],birth_time_s=birth,scheduled_time_s=birth,
            anchor_crossing_time_s=birth,anchor_m=e['anchor_m'],anchor_triangle=e['anchor_triangle'],
            birth_center_m=admit['birth_center_m'],admission_strategy='B',entry_transition_time_s=0.,
            first_diameter_um=e['first_diameter_um'],common_event_sha256=canonical_hash(e),
            admission_sha256=canonical_hash(admit),production_provenance=p)
        event['cache_key']=canonical_hash(dict(event=event,nominal_dt_s=DT))
        accepted.append(event)
    ledger=dict(schema='PARTICLE9A_2MMPS_FRESH_BIRTH_LEDGER_V1',production_provenance=p,
        master_seed=SEED,Q_in_m3_s=c.env.sampler.Q_m3_s,C_MB_m3=C_MB,
        Ndot_MB_s_inv=C_MB*c.env.sampler.Q_m3_s,scheduled_count=pid,admitted_count=count,
        admission_rejected_count=pid-count,events=accepted,
        scheduled_window_s=clock.time_at(pid),selection='FIRST_ACCEPTED_IDS_NO_OUTCOME_OR_OUTLET_SELECTION',
        size_distribution_role='EXISTING_METHOD_B_POSITION_CONDITIONED_SIZE_RETRIES_NOT_UNCONDITIONED_SONOVUE',
        clock='EXISTING_P7_FLUX_CLOCK_REGENERATED_FROM_NEW_INTEGRATED_Q',old_cache_used=False)
    dump(target,ledger);dump(output/'admission/all_scheduled_admission_records.json',rows)
    dump(output/'provenance/input_and_source.json',p)
    dump(output/'provenance/inlet_flux.json',c.env.audit)
    return ledger


def integrate_current(event,output=OUTPUT):
    # Check before integrate_admitted can take its resume branch.
    p=production_provenance()
    if event.get('production_provenance')!=p:raise ValueError('Event input/source provenance mismatch')
    env=environment()
    created=[]
    def factory(*args,**kwargs):
        stepper=Particle9AStepper(*args,gradient_provider=lambda x:env.field.sample(x).velocity_gradient_s_inv,**kwargs)
        created.append(stepper)
        return stepper
    # Private function dependency binding preserves the original integration
    # provider, guards, accepted substeps and outlet classification verbatim.
    integrate=bind_query_dependency(integrate_admitted,{'SavedTrajectoryStepper':factory})
    result=integrate(event,output=output)
    if created:
        result.update(production_provenance=p,dynamics=MODEL,
            timestep_role='REUSED_P81_NOMINAL_DT_WITH_UNCHANGED_P65_PHYSICAL_TIME_SUBDIVISION',
            planar_wall_statistics=dict(created[0].planar_statistics))
        dump(Path(output)/'trajectories'/f'mb_{event["particle_id"]:06d}.json',result)
    elif result.get('production_provenance')!=p or result.get('dynamics')!=MODEL:
        raise ValueError('Resume missing P9-A model provenance')
    return result


def worker(event):
    start=time.time();cpu=time.process_time();result=integrate_current(event,_JOB['output'])
    path=Path(_JOB['output'])/'trajectories'/f'mb_{event["particle_id"]:06d}.json'
    receipt=path.with_suffix('.receipt.json')
    if not receipt.exists():
        dump(receipt,dict(hostname=socket.gethostname(),pid=os.getpid(),worker=mp.current_process().name,
            remote_server_compute=True,start_unix_s=start,end_unix_s=time.time(),
            cpu_seconds=time.process_time()-cpu,peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            metadata_sha256=sha256(path),samples_sha256=result['samples_sha256'],
            cache_key=event['cache_key'],production_provenance=event['production_provenance']))
    else:
        r=json.loads(receipt.read_text())
        if r['metadata_sha256']!=sha256(path) or r['cache_key']!=event['cache_key']:
            raise ValueError('Resume receipt mismatch')
    return dict(particle_id=event['particle_id'],completed=result['completed'],end_reason=result['end_reason'],
        outlet=result['exit_outlet'],minimum_gap_m=result.get('minimum_original_wall_gap_m'),
        residence_time_s=result['residence_time_s'],wall_seconds=result['wall_seconds'],
        accepted_steps=result.get('accepted_steps',0),rejected_trials=result.get('rejected_trials',0),
        failure_detail=result['failure_detail'])


def main():
    global _JOB
    parser=argparse.ArgumentParser();parser.add_argument('--stage',choices=['prepare','smoke','production'],required=True)
    parser.add_argument('--workers',type=int,default=6);parser.add_argument('--output',type=Path,default=OUTPUT)
    args=parser.parse_args();start=time.time();ledger=prepare_cohort(output=args.output)
    if args.stage=='prepare':print(json.dumps({k:v for k,v in ledger.items() if k!='events'},indent=2));return
    if not str(REPO).startswith('/workspace/particle9a_2mmps_'):
        raise ValueError('Formal smoke/production must run in the isolated authorized server directory')
    _JOB=dict(output=str(args.output));events=ledger['events'][:12] if args.stage=='smoke' else ledger['events']
    if args.stage=='production':
        smoke=json.loads((args.output/'provenance/smoke_summary.json').read_text())
        if not smoke['smoke_gate_pass']:raise ValueError('Smoke gate failed; production blocked')
    rows=[]
    with mp.get_context('fork').Pool(args.workers) as pool:
        for row in pool.imap_unordered(worker,events,chunksize=1):
            rows.append(row)
            print(json.dumps(dict(progress=len(rows),total=len(events),**row)),flush=True)
            dump(args.output/'provenance'/f'{args.stage}_progress.json',dict(completed_tasks=len(rows),total=len(events),elapsed_s=time.time()-start))
    rows.sort(key=lambda r:r['particle_id'])
    failures=[r for r in rows if r['end_reason'] in ('INTEGRATION_SAFETY_STOP','NUMERICAL_REFINEMENT_FAILURE')]
    penetration=[r for r in rows if r['minimum_gap_m'] is None or r['minimum_gap_m']< -environment().wall.roundoff_m]
    summary=dict(stage=args.stage,hostname=socket.gethostname(),server_root=str(REPO),workers=args.workers,
        start_unix_s=start,end_unix_s=time.time(),wall_seconds=time.time()-start,
        generated_count=len(rows),completed_count=sum(r['completed'] for r in rows),
        outlet_counts=dict(Counter(r['outlet'] or 'NO_EXIT' for r in rows)),
        end_reason_counts=dict(Counter(r['end_reason'] for r in rows)),
        solver_or_subdivision_failure_count=len(failures),penetration_count=len(penetration),
        smoke_gate_pass=len(rows)==12 and not failures and not penetration and all(r['completed'] for r in rows) if args.stage=='smoke' else None,
        GPU_role='GPU_NOT_PERFORMANCE_CRITICAL_FOR_THIS_STAGE',rows=rows,
        production_provenance=ledger['production_provenance'])
    dump(args.output/'provenance'/f'{args.stage}_summary.json',summary)
    print(json.dumps({k:v for k,v in summary.items() if k not in ('rows','production_provenance')},indent=2),flush=True)


if __name__=='__main__':main()
