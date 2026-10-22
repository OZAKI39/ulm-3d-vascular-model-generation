"""Receipt packaging after supervised work has exited. Never submits a job."""
from common import *
import datetime,os,psutil,platform,scipy

def main():
    now=datetime.datetime.now(datetime.timezone.utc).isoformat();own=psutil.Process();skip={os.getpid(),*(p.pid for p in own.parents())};active=[]
    for p in psutil.process_iter(['pid','cmdline','status']):
        try:
            if p.pid not in skip and p.info['status']!=psutil.STATUS_ZOMBIE and any(str(HERE) in x for x in p.info['cmdline'] or []):active.append(dict(pid=p.pid,command=p.info['cmdline']))
        except psutil.Error:pass
    assert not active,active
    integrity={}
    for label,manifest,root in [('this_review',HERE/'INPUT_SHA256.json',Path('/')),('prior_cause_inputs',CAUSE/'INPUT_SHA256.json',BASE),('prior_field_outputs',PREVIOUS/'OUTPUT_SHA256.json',PREVIOUS)]:
        inputs=read(manifest);bad={k:sha(root/k) for k,v in inputs.items() if sha(root/k)!=v};assert not bad,bad
        integrity[label]=dict(count=len(inputs),unchanged=True,manifest_sha256=sha(manifest))
    put('FINAL_INPUT_INTEGRITY.json',dict(checked_UTC=now,groups=integrity))
    jobs=[]
    for path in sorted((HERE/'monitoring').glob('*/EXECUTION.json')):
        d=read(path);assert d['status']=='EXITED';assert read(path.parent/'LIVE.json')['status']=='EXITED';log=path.parent/'job.log'
        jobs.append(dict(name=path.parent.name,exit_code=d['exit_code'],stop_reason=d['stop_reason'],elapsed_s=d['elapsed_s'],completed_UTC=d['completed_UTC'],peak_tree_RSS_GiB=d['peak_tree_RSS_GiB'],peak_cgroup_GiB=d['peak_cgroup_GiB'],peak_nonreclaimable_GiB=max(v['nonreclaimable_upper_GiB'] for v in d['samples']),oom_event_increase=max(v['oom_events']-d['initial_memory']['oom_events'] for v in d['samples']),execution_sha256=sha(path),log_sha256=sha(log),log_tail=log.read_text()[-1800:],task=d['task']))
    assert read(HERE/'TEST_RESULTS.json')['status']=='PASS';total=sum(j['elapsed_s'] for j in jobs);assert total<=5400
    for item in read(HERE/'VISUAL_QC.json')['figures']:assert sha(HERE/item['file'])==item['sha256']
    put('EXECUTION_RECEIPT.json',dict(status='COMPLETED_ALL_OWN_JOBS_EXITED',completed_UTC=now,authoritative_project=str(BASE),new_CFD_count=0,new_pipe_CFD_count=0,new_BraVa_CFD_count=0,new_solver_builds=0,new_solver_patch=False,new_mesh_count=0,new_population_count=0,new_particle_integrations=0,reused_fields=NAMES,
       supervised_analysis_budget_s=5400,supervised_elapsed_s_including_failed_attempts=total,timing_scope='Supervised server jobs only, including any failures; excludes interactive source/literature inspection, coding, transfers, metadata packaging and Git.',jobs=jobs,active_project_processes=active,automatic_next_job=False,logs='Original monitoring/ retained on server and local, Git-ignored; job SHA and output tails included.',final_tests='TEST_RESULTS.json',integrity='FINAL_INPUT_INTEGRITY.json',decision='STOP_C',GPU_CFD='NOT_RUN after Stop C; retained production GPU execution identity unchanged',environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,RT_library='scikit-fem 12.0.2 isolated runtime install')))
    put('RESULT_IDENTITY.json',dict(status='COMPLETED_DIAGNOSTIC_REVIEW_STOP_C',created_UTC=now,local_project='/home/lzy/projects/computation_examples/brava_flow_roi_18mlmin',server_project=str(BASE),analysis_directory=str(HERE.relative_to(BASE)),task_start_local_branch='validation/gmsh-reference-v2',task_start_local_HEAD='4d4fc66585c471e7e3539a6da3a57c7c8fd01b5e',production_binary_sha256=read(HERE/'PRODUCTION_SOURCE_IDENTITY.json')['actual_binary_sha256'],production_source_commit=read(HERE/'PRODUCTION_SOURCE_IDENTITY.json')['source_commit'],source_status='SOURCE_CONFIRMED',term_rows=read(HERE/'STABILIZATION_TERMS.json')['rows'],term_status='OFFLINE_EXISTING_FIELD_EVALUATION',counterpart_status='OFFLINE_HYPOTHETICAL_TIME_CONSISTENT_EVALUATION',RT_REFERENCE_NOT_PRODUCTION_FIELD=True,higher_order='HIGHER_ORDER_FLUID_READY_TO_PROTOTYPE',new_CFD=0,new_population=0,new_mesh=0,new_trajectories=0,FLOW_FIELD_PHYSICALLY_VALIDATED=False,physiological_validation='NOT_RUN',time_consistent_pipe='NOT_RUN',time_consistent_BraVa='NOT_RUN',only_recommended_next_direction='Existing svMultiPhysics Taylor-Hood P2/P1 analytical-field and higher-order sampler validation prototype; NOT started',delivery_commit='git log -- this directory identifies independent commit; no self-referential hash in manifest',push_performed=False))
    ignore={'runtime','monitoring','__pycache__','.pytest_cache'};files=[p for p in HERE.rglob('*') if p.is_file() and not(set(p.relative_to(HERE).parts)&ignore) and p.name not in ['TASK.json','OUTPUT_SHA256.json']]
    for p in files:
        if p.suffix=='.json':json.loads(p.read_text(),parse_constant=lambda x:(_ for _ in ()).throw(ValueError(x)))
    put('OUTPUT_SHA256.json',{str(p.relative_to(HERE)):sha(p) for p in sorted(files)})
    print('FINALIZED',len(files),'files',len(jobs),'jobs',total,'s; active',active)
if __name__=='__main__':main()
