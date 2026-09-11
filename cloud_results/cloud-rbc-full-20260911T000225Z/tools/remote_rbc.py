"""Single formal run using the existing guarded process/seal/SSH workflow primitives."""
from common import *
import resource, shutil, sys, fcntl, threading, traceback, csv


def processes(group,include_runner=False):
    result=[]
    for p in Path('/proc').iterdir():
        if not p.name.isdigit():continue
        try:
            raw=(p/'stat').read_text();fields=raw[raw.rfind(')')+2:].split()
            if (int(fields[2])!=group and not (include_runner and int(p.name)==os.getpid())) or fields[0]=='Z':continue
            result.append(dict(pid=int(p.name),ppid=int(fields[1]),pgid=int(fields[2]),state=fields[0],
                               cpu_s=(int(fields[11])+int(fields[12]))/os.sysconf('SC_CLK_TCK'),rss_bytes=int(fields[21])*os.sysconf('SC_PAGE_SIZE'),
                               command=(p/'cmdline').read_bytes().replace(b'\0',b' ').decode(errors='replace')))
        except (FileNotFoundError,PermissionError,ProcessLookupError):pass
    return result


def claim(root):
    """An atomic persistent claim prevents any second launch, including after failure."""
    p=Path(root)/'STARTED'
    with p.open('x') as f:json.dump(dict(pid=os.getpid(),started_at=now()),f)


def main():
    cfg=read(sys.argv[1]);root=Path(cfg['output_directory']);job=root.name
    assert root.parent==Path('/workspace/bloodflow/cloud_runs') and root.resolve()==root and job==cfg['job_id']
    claim(root)
    resource.setrlimit(resource.RLIMIT_CORE,(0,0))
    lock=Path('/workspace/bloodflow/metadata/cloud_smoke.lock').open('a')
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB) # Reuse the existing solver concurrency lock.
    state=dict(job_id=job,status='STARTED',pid=os.getpid(),started_at=now(),protocol='CLOUD_FULL_RBC_V1',exit_code=None)
    write(root/'run_status.json',state);p=None;start=time.monotonic();terminal='FAILED';done=threading.Event();stop={};cache={}
    spec=read(root/'full_spec.json');work=root/'simulation';work.mkdir()
    env=dict(os.environ,**cfg['environment'],PATH='/workspace/bloodflow/.venv/bin:/usr/local/cuda/bin:/usr/bin:/bin',
             LD_LIBRARY_PATH='/usr/local/cuda/lib64',PYTHONPATH=cfg['cloud_project_root']+':'+str(root/'tools'),CLOUD_JOB_AUTHORIZED=cfg['plan_sha256'])
    env.pop('SINGLE_RBC_AUTHORIZED',None)
    monitor=None
    try:
        assert cfg['authorization']['plan_sha256']==cfg['plan_sha256'] and cfg['authorization']['attempts']==1 and cfg['authorization']['process_seconds']==1800
        assert sha(spec['native_library'])==spec['native_library_sha256'],'RUNTIME_BINARY_CHANGED'
        for name,h in cfg['tool_hashes'].items():assert sha(root/'tools'/name)==h,'WORKER_CHANGED'
        assert sha(root/'full_spec.json')==cfg['spec_sha256'],'SPEC_CHANGED'
        assert tree_identity('/workspace/bloodflow/.venv')==cfg['environment_before'],'PROJECT_VENV_CHANGED_AFTER_PREFLIGHT'
        for entry in cfg['python_inputs']:assert sha(entry['remote'])==entry['sha256'],'INPUT_CHANGED'
        assert shutil.disk_usage('/workspace').free-cfg['expected_output_bytes']>=5*1024**3,'DISK_HEADROOM'
        argv=['/usr/bin/mpirun','--allow-run-as-root','--bind-to','none','-np','2',cfg['python'],'-B',str(root/'tools/rbc_full_worker.py'),'--spec',str(root/'full_spec.json')]
        with (work/'console.log').open('wb') as log:
            launched=time.monotonic();p=subprocess.Popen(argv,cwd=work,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
            state.update(status='RUNNING',solver_pid=p.pid,solver_pgid=p.pid,argv=argv,solver_started_at=now());write(root/'run_status.json',state)
            def observe():
                from rbc_checks import scan_frames
                last_gpu=0;gpu=None
                with (root/'resources.jsonl').open('x',buffering=1) as out:
                    while not done.is_set():
                        try:
                            if time.monotonic()-last_gpu>=5:
                                q=subprocess.run(['nvidia-smi','--query-gpu=memory.used,memory.free,utilization.gpu','--format=csv,noheader,nounits'],capture_output=True,text=True,timeout=2)
                                gpu=q.stdout.strip() if q.returncode==0 else 'UNMEASURED';last_gpu=time.monotonic()
                            memory=Path('/sys/fs/cgroup/memory.current').read_text().strip();maximum=Path('/sys/fs/cgroup/memory.max').read_text().strip()
                            sample=dict(at=now(),elapsed_s=time.monotonic()-launched,processes=processes(p.pid,include_runner=True),gpu_whole_device=gpu,
                                        cgroup_memory_bytes=int(memory),disk_free_bytes=shutil.disk_usage('/workspace').free)
                            out.write(json.dumps(sample)+'\n')
                            if maximum!='max' and int(maximum)-int(memory)<512*1024**2:stop['reason']='CGROUP_MEMORY_LIMIT'
                            if (root/'STOP_REQUESTED').exists():stop['reason']='USER_STOP'
                            if sum(x.stat().st_size for x in root.rglob('*') if x.is_file())>cfg['expected_output_bytes']:stop['reason']='OUTPUT_SIZE_LIMIT'
                            scan_frames(work,spec,cache)
                            if any(v.get('phase')=='shear' for v in cache.values()) and not (root/'live_handoff.json').exists():
                                from rbc_report import handoff_audit
                                handoff=handoff_audit(work,spec)
                                if handoff['status']!='NOT_REACHED':
                                    write(root/'live_handoff.json',handoff)
                                    if handoff['status']=='FAILED':stop['reason']='HANDOFF_STATE_MISMATCH'
                            for item in cache.values():
                                if item.get('hard_failures'):stop['reason']='GEOMETRY_HARD_FAILURE';write(root/'live_failure.json',item);break
                            for f in work.glob('native_stats_*_fluid.csv'):
                                with f.open() as stream:rows=list(csv.DictReader(stream))
                                for row in rows:
                                    if None in row.values():continue
                                    if int(row['num_particles'])!=110592:stop['reason']='FLUID_POPULATION_CHANGED'
                                    import math
                                    if not all(math.isfinite(float(row[k])) for k in ['kBT','maxv']):stop['reason']='NONFINITE_FLUID_STATS'
                            progress=read(work/'progress.json') if (work/'progress.json').exists() else {}
                            lower={phase:max([x['phase_step'] for x in cache.values() if x['phase']==phase],default=0) for phase in ['relaxation','shear']}
                            write(root/'live_progress.json',dict(progress=progress,saved_progress_lower_bound=lower,observed_at=now(),stop_reason=stop.get('reason')))
                        except Exception as ex:
                            # Transient incomplete output is not a fabricated completed observation.
                            write(root/'monitor_last_error.json',dict(error=str(ex),observed_at=now()))
                        done.wait(.5)
            monitor=threading.Thread(target=observe,daemon=True);monitor.start()
            waited=wait_guarded(p,max(.01,1800-(time.monotonic()-launched)),lambda:shutil.disk_usage('/workspace').free,stop_reason=lambda:stop.get('reason'))
        state.update(exit_code=p.returncode,solver_process_wall_s=time.monotonic()-launched,stop_reason=waited['stop_reason'])
        if waited['stop_reason']:terminal='STOPPED'
        elif p.returncode:terminal='FAILED'
        elif (work/'completion.json').exists() and read(work/'completion.json').get('completed'):terminal='COMPLETED'
        else:terminal='STOPPED';state['stop_reason']=read(work/'stop_reason.json') if (work/'stop_reason.json').exists() else 'NO_COMPLETION_RECORD'
    except Exception as ex:
        state['error']=type(ex).__name__+': '+str(ex);(root/'failure_traceback.txt').write_text(traceback.format_exc())
    finally:
        if p is not None and p.poll() is None:
            os.killpg(p.pid,signal.SIGTERM)
            try:p.wait(timeout=2)
            except subprocess.TimeoutExpired:os.killpg(p.pid,signal.SIGKILL);p.wait()
        done.set()
        if monitor is not None:monitor.join(timeout=30)
        # All processes in the registered group belong to this task. Never kill by program name.
        left=processes(p.pid) if p is not None else []
        if left:
            os.killpg(p.pid,signal.SIGTERM);time.sleep(.2)
            if processes(p.pid):os.killpg(p.pid,signal.SIGKILL)
            time.sleep(.2);left=processes(p.pid)
        state.update(status=terminal,finished_at=now(),runner_elapsed_s=time.monotonic()-start,remaining_task_processes=left,
                     billing_note='Solver process wall seconds include CPU preparation/analysis inside the worker; not Vast rental billing or measured GPU-active seconds')
        if left:state['status']=terminal='FAILED';state['error']='LINGERING_TASK_PROCESSES'
        after=tree_identity('/workspace/bloodflow/.venv')
        protection=dict(project_venv_before=cfg['environment_before'],project_venv_after=after,unchanged=after==cfg['environment_before'],runtime_library_sha256=sha(spec['native_library']))
        write(root/'environment_protection.json',protection)
        if not protection['unchanged']:state['status']=terminal='FAILED';state['error']='PROJECT_VENV_CHANGED'
        write(root/'run_status.json',state);write(root/'execution.json',state)
        write(root/'budget_usage.json',dict(authorization=cfg['authorization'],attempts=1,process_wall_s=state.get('solver_process_wall_s'),process_limit_s=1800,
                                          old_budget_reused=False,automatic_retry=False))
        try:
            from rbc_report import analyze
            analyze(root,spec,state)
        except Exception as ex:
            write(root/'analysis_failure.json',dict(error=type(ex).__name__+': '+str(ex),simulation_must_not_rerun=True))
        if monitor is not None and monitor.is_alive():
            # Never claim files are closed while an observer still writes.
            raise RuntimeError('OBSERVER_NOT_CLOSED_RESULTS_NOT_SEALED')
        write(root/'disk_after.json',dict(free_bytes=shutil.disk_usage('/workspace').free,files=len(files(root)),bytes=sum(x['size'] for x in files(root)),deletion_performed=False))
        write(root/'result_manifest.json',dict(canonical_manifest='results_manifest.json',protocol='existing closed-file SHA256 manifest; this record points to the authoritative transport manifest'))
        seal(root,job,terminal)
    return 0 if terminal=='COMPLETED' else 2


if __name__=='__main__':sys.exit(main())
