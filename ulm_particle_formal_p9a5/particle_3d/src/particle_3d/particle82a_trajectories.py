"""Remote P8.2A cohorts: accepted-ID order, no outlet/outcome selection."""
from pathlib import Path
from collections import Counter
import argparse, gzip, json, multiprocessing as mp, os, resource, socket, time
import numpy as np
from .particle82_provenance import atomic_json, sha256, require_remote
from .particle8_replay import canonical_hash
from .particle82a_admission import context
from .particle82a_integration import integrate_admitted

_JOB=None


def iter_rows(admission):
    for file in sorted((Path(admission)/'events').glob('events_*.json.gz')):
        with gzip.open(file,'rt') as stream:
            yield from json.load(stream)['rows']


def make_event(row, method):
    e=row['event'];a=row['methods'][method]
    if not a['accepted']:raise ValueError('Cannot integrate rejected event')
    q=np.random.default_rng([e['master_seed'],e['event_id'],3]).normal(size=4);q/=np.linalg.norm(q)
    interval=1/(8.5e12*context().env.sampler.Q_m3_s)
    transition=a.get('entry_transition_time_s',0.)
    return dict(particle_id=e['event_id'],species='MB',attempt_count=0,
        radius_m=a['radius_m'],diameter_um=a['diameter_um'],q=q.tolist(),
        position_seed=[e['master_seed'],e['event_id'],1],
        birth_time_s=e['event_id']*interval+transition,
        scheduled_time_s=e['event_id']*interval,anchor_crossing_time_s=e['event_id']*interval,
        anchor_m=e['anchor_m'],anchor_triangle=e['anchor_triangle'],birth_center_m=a['birth_center_m'],
        admission_strategy=method,entry_transition_time_s=transition,
        s_birth_m=a.get('s_birth_m',0.),point_tracer_basin=e['point_tracer_basin'],
        birth_point_basin=a['birth_point_basin'],first_diameter_um=e['first_diameter_um'],
        common_event_sha256=a['common_event_sha256'],
        entry_aperture_passable=e['aperture_passable'],
        representation_transition_not_finite_size_entry_dynamics=True)


def select_cohorts(admission, output, quota=5000):
    cohorts={k:[] for k in 'ABC'};prefix={};examined=0
    for row in iter_rows(admission):
        examined+=1
        for method in 'ABC':
            if len(cohorts[method])<quota and row['methods'][method]['accepted']:
                cohorts[method].append(make_event(row,method))
                if len(cohorts[method])==quota:prefix[method]=row['event']['event_id']
        if all(len(v)==quota for v in cohorts.values()):break
    if not all(len(v)==quota for v in cohorts.values()):
        raise ValueError('Common ledger has fewer than required admitted events: '+str({k:len(v) for k,v in cohorts.items()}))
    path=Path(output)/'TRAJECTORY_COHORTS.json'
    data=dict(selection='FIRST_5000_ACCEPTED_IDS_PER_METHOD_IN_COMMON_LEDGER_ORDER_NO_OUTLET_OR_OUTCOME_CONDITIONING',
        quota=quota,scheduled_prefix_by_method=prefix,common_prefix_examined=examined,
        conditional_sampling_warning='Accepted cohort statistics are conditional on each admission strategy, including C.',
        cohorts=cohorts)
    if path.exists() and json.loads(path.read_text())!=data:raise ValueError('Existing cohort selection differs')
    atomic_json(path,data);return cohorts


def initialize(output,provenance):
    global _JOB
    require_remote(provenance,provenance['hostname']);context()
    _JOB=dict(output=str(output),provenance=provenance)


def worker(task):
    method,event,subdir,dt=task
    out=Path(_JOB['output'])/subdir/method;start=time.time();cpu=time.process_time()
    result=integrate_admitted(event,output=out,dt=dt)
    meta=out/'trajectories'/f'mb_{event["particle_id"]:06d}.json'
    receipt_path=meta.with_suffix('.receipt.json')
    if not receipt_path.exists():
        p=_JOB['provenance']
        atomic_json(receipt_path,dict(hostname=socket.gethostname(),pid=os.getpid(),
            git_commit=p['source_git_commit'],config_hash=canonical_hash(result['integration_config']),
            frozen_input_hash=canonical_hash(p['frozen_input_sha256']),start_time=start,end_time=time.time(),
            first_id=event['particle_id'],count=1,method=method,REMOTE_SERVER_COMPUTE=True,
            worker=mp.current_process().name,metadata_sha256=sha256(meta),samples_sha256=result['samples_sha256'],
            worker_cpu_seconds=time.process_time()-cpu,
            peak_worker_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
    else:
        receipt=json.loads(receipt_path.read_text())
        if receipt['metadata_sha256']!=sha256(meta) or receipt['samples_sha256']!=result['samples_sha256']:
            raise ValueError('Trajectory resume receipt mismatch')
        if receipt['git_commit']!=_JOB['provenance']['source_git_commit'] or receipt['hostname']!=socket.gethostname():
            raise ValueError('Trajectory source/host mismatch')
    return dict(method=method,particle_id=event['particle_id'],end_reason=result['end_reason'],
        completed=result['completed'],exit_outlet=result['exit_outlet'],sample_hash=result['samples_sha256'],
        peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)


def run_tasks(tasks,workers,label):
    import psutil
    parent=psutil.Process();start=time.time();rows=[];peak_pss=0
    usage0=resource.getrusage(resource.RUSAGE_CHILDREN)
    with mp.get_context('fork').Pool(workers) as pool:
        for r in pool.imap_unordered(worker,tasks,chunksize=1):
            rows.append(r)
            if len(rows)%32==0 or len(rows)==len(tasks):
                pss=0
                for proc in [parent,*parent.children(recursive=True)]:
                    try:pss+=proc.memory_full_info().pss
                    except psutil.Error:pass
                peak_pss=max(peak_pss,pss)
                print(label,len(rows),'/',len(tasks),dict(Counter(x['end_reason'] for x in rows)),flush=True)
    elapsed=time.time()-start
    usage1=resource.getrusage(resource.RUSAGE_CHILDREN)
    cpu=usage1.ru_utime+usage1.ru_stime-usage0.ru_utime-usage0.ru_stime
    return sorted(rows,key=lambda r:(r['method'],r['particle_id'])),dict(workers=workers,
        trajectories=len(tasks),seconds=elapsed,trajectories_per_hour=len(tasks)/elapsed*3600,
        worker_cpu_seconds=cpu,mean_busy_cpu_cores=cpu/elapsed,
        worker_capacity_utilization_percent=100*cpu/(elapsed*workers),
        observed_peak_tree_pss_bytes=peak_pss,
        peak_worker_rss_kib=max(r['peak_rss_kib'] for r in rows))


def benchmark(cohorts):
    results=[];fingerprints=[]
    chosen_events=[(k,e) for k in 'ABC' for e in cohorts[k][:6]]
    for workers in [1,2,4,8,16]:
        rows,performance=run_tasks([(k,e,f'benchmark/w{workers}',.00025) for k,e in chosen_events],workers,'TRAJECTORY_SCALING')
        results.append(performance)
        fingerprints.append(canonical_hash([{k:v for k,v in x.items() if k!='peak_rss_kib'} for x in rows]))
    if len(set(fingerprints))!=1:raise ValueError('Trajectory parallel/serial parity failed')
    chosen=max(results,key=lambda r:r['trajectories_per_hour'])['workers']
    atomic_json(Path(_JOB['output'])/'TRAJECTORY_SCALING.json',dict(results=results,
        chosen_workers=chosen,scientific_parity=True,fingerprints=fingerprints,
        selection='FASTEST_MEASURED_MIXED_METHOD_FIXED_COHORT',
        shared_readonly_FEM='FORK_COPY_ON_WRITE'))
    return chosen


def baseline_parity(cohorts):
    from .particle81_simulation import integrate_one
    rows=[]
    for event in cohorts['A'][:4]:
        original=integrate_one(event,output=Path(_JOB['output'])/'original_p81_parity')
        current=integrate_admitted(event,output=Path(_JOB['output'])/'explicit_birth_parity')
        a=np.load(Path(_JOB['output'])/'original_p81_parity'/original['samples_path'])['samples']
        b=np.load(Path(_JOB['output'])/'explicit_birth_parity'/current['samples_path'])['samples']
        equal=np.array_equal(a,b)
        rows.append(dict(event_id=event['particle_id'],samples_byte_equal=equal,
            original_end_reason=original['end_reason'],current_end_reason=current['end_reason']))
        if not equal:raise ValueError('Original P8.1 integration vs explicit accepted birth mismatch')
    atomic_json(Path(_JOB['output'])/'P81_INTEGRATION_PARITY.json',dict(all_exact=True,rows=rows))


def main():
    p=argparse.ArgumentParser();p.add_argument('--admission',required=True);p.add_argument('--output',required=True)
    p.add_argument('--provenance',required=True);p.add_argument('--quota',type=int,default=5000)
    p.add_argument('--action',choices=['prepare','benchmark','run','sensitivity'],default='run')
    p.add_argument('--workers',type=int,default=8);a=p.parse_args()
    provenance=json.loads(Path(a.provenance).read_text());initialize(a.output,provenance)
    output=Path(a.output);output.mkdir(parents=True,exist_ok=True)
    cohorts=select_cohorts(a.admission,output,a.quota)
    if a.action=='prepare':baseline_parity(cohorts)
    elif a.action=='benchmark':baseline_parity(cohorts);benchmark(cohorts)
    elif a.action=='sensitivity':
        # Preselected first 12 admitted IDs per method; never selected by outcome.
        for factor in [2,4]:
            run_tasks([(k,e,f'sensitivity/dt_div_{factor}',.00025/factor) for k in 'ABC' for e in cohorts[k][:12]],a.workers,'TIMESTEP_SENSITIVITY')
    else:
        start=time.time();tasks=[(k,e,'formal',.00025) for k in 'ABC' for e in cohorts[k]]
        rows,performance=run_tasks(tasks,a.workers,'FORMAL_TRAJECTORIES')
        summary={k:dict(scheduled=len(cohorts[k]),completed=sum(r['completed'] for r in rows if r['method']==k),
            end_reasons=dict(Counter(r['end_reason'] for r in rows if r['method']==k))) for k in 'ABC'}
        atomic_json(output/'FULL_TRAJECTORY_RUN_RECEIPT.json',dict(hostname=socket.gethostname(),
            git_commit=provenance['source_git_commit'],REMOTE_SERVER_COMPUTE=True,
            start_time=start,end_time=time.time(),performance=performance,methods=summary))


if __name__=='__main__':main()
