"""Run after supervised packaging exits. No numerical work or next job."""
from common import *
import datetime
import psutil

def main():
    assert not (HERE/'TASK.json').exists(), 'Wait for supervisor exit'
    active=[]
    for proc in psutil.process_iter(['pid','cmdline']):
        if proc.pid==__import__('os').getpid():continue
        args=proc.info.get('cmdline') or []
        for arg in args:
            # Match actual script/executable arguments, not parent shell command strings.
            if arg==str(HERE/'supervise.py') or arg.startswith(str(HERE/'runtime/harness_')) or arg==str(HERE/'runtime/native_boundary'):
                active.append(dict(pid=proc.pid,program=arg))
    assert not active,active
    jobs=[]
    for p in sorted((HERE/'monitoring').glob('*/EXECUTION.json')):
        r=read(p)
        job={k:r.get(k) for k in ['status','exit_code','stop_reason','pid','elapsed_s','completed_UTC',
             'peak_tree_RSS_GiB','peak_cgroup_GiB','peak_working_set_GiB','initial_memory','memory_accounting','automatic_next_job']}
        job.update(name=p.parent.name,receipt_sha256=sha(p),command=r['task']['command'],
                   sample_count=len(r.get('samples',[])),max_seconds=r['task']['max_seconds'],
                   process_RSS_limit_GiB=r['task']['max_tree_RSS_GiB'],nonreclaimable_limit_GiB=r['task']['max_nonreclaimable_GiB'])
        log=p.parent/'job.log'
        job['log_sha256']=sha(log)
        if r.get('exit_code')!=0:job['log_tail']='\n'.join(log.read_text(errors='replace').splitlines()[-8:])
        initial=r['initial_memory'];samples=r.get('samples',[])
        job['OOM_event_increase']=max([s['oom_events']-initial['oom_events'] for s in samples]+[0])
        job['OOM_kill_increase']=max([s['oom_kill_events']-initial['oom_kill_events'] for s in samples]+[0])
        jobs.append(job)
    put('EXECUTION_RECEIPT.json',dict(status='ALL_TASKS_EXITED',UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        jobs=jobs,supervised_job_count=len(jobs),total_supervised_elapsed_s=sum(j['elapsed_s'] for j in jobs),
        peak_sampled_tree_RSS_GiB=max(j['peak_tree_RSS_GiB'] for j in jobs),
        active_jobs=active,TASK_present=False,automatic_next_job=False,
        known_failed_attempt='validate: 20 Hessian checks failed; preserved in ATTEMPT_01_HARNESS_RESULTS.json, corrected without tolerance relaxation',
        final_scientific_stop='Phase F: MIXED_CURVED_BOUNDARY_ASSEMBLY_BLOCKED',
        CFD_steps=0,linear_solves=0,GPU_compute_jobs=0,particle_runs=0,
        runtime_scope='Supervised process wall time only. Excludes reading, development, transfers, report review. RSS is 1-second sampled process tree, not a rigorous unsampled peak bound.',
        visual_review=read(HERE/'VISUAL_REVIEW.json')))
    excluded={'OUTPUT_SHA256.json','TASK.json'}
    outputs={p.name:sha(p) for p in sorted(HERE.iterdir()) if p.is_file() and p.name not in excluded}
    put('OUTPUT_SHA256.json',outputs)
    print('CLOSED',len(jobs),'jobs;',len(outputs),'artifacts;',sum(j['elapsed_s'] for j in jobs),'supervised seconds')

if __name__=='__main__':main()
