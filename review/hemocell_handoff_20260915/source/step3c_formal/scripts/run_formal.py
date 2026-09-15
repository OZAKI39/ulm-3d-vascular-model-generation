#!/usr/bin/python3
"""One-shot formal execution. Durable logs; never restart a solver."""
import os,sys,time,json,csv,subprocess,threading,fcntl,shutil,datetime
from pathlib import Path
sys.dont_write_bytecode=True
from remote_common import sha,read,write,verify_bundle,verify_run_inputs,contract
from evaluator_service import serve_evaluations
R=Path(__file__).resolve().parents[1];S=R.parent/'20260914_stage4_batch_dispatch';P=R.parent/'toolchains/nvhpc/Linux_x86_64/26.5'
EXE=R/'build_gpu/vascular/vascularPoC'
def utc():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def memory():
    m={l.split(':')[0]:int(l.split()[1])*1024 for l in Path('/proc/meminfo').read_text().splitlines() if l.startswith(('MemAvailable:','MemTotal:'))}
    lim=Path('/sys/fs/cgroup/memory.max');cur=Path('/sys/fs/cgroup/memory.current')
    m['effective_available']=min(m['MemAvailable'],int(lim.read_text())-int(cur.read_text())) if lim.exists() and lim.read_text().strip()!='max' else m['MemAvailable']
    return m
def processes():
    out={}
    for d in Path('/proc').iterdir():
        if not d.name.isdigit():continue
        try:
            exe=os.readlink(d/'exe');name=(d/'comm').read_text().strip();s=(d/'stat').read_text().split()
            if exe.startswith('/workspace/hemocell_gpu_poc') or name in ['mpirun','nsys','ncu']:
                out[int(d.name)]={'exe':exe,'name':name,'ticks':int(s[13])+int(s[14]),'rss':int((d/'statm').read_text().split()[1])*os.sysconf('SC_PAGE_SIZE')}
        except OSError:pass
    return out
def gpu():
    return [float(x.strip()) for x in subprocess.check_output(['nvidia-smi','--query-gpu=utilization.gpu,memory.used,memory.total','--format=csv,noheader,nounits'],text=True).strip().split(',')]
def preflight():
    c=contract(R);verify_bundle(R);verify_run_inputs(R,R)
    m=memory();free=shutil.disk_usage(R).free;ps=processes()
    bad=[{'pid':pid,**v} for pid,v in ps.items() if v['name'].startswith(('vascular','cavity','inletCal','nvc','nsys','ncu','mpirun'))]
    apps=subprocess.check_output(['nvidia-smi','--query-compute-apps=pid,process_name,used_memory','--format=csv,noheader'],text=True).strip()
    milestones=set(c['execution']['field_snapshots'])|{5000,240000,700000}
    evaluations=set(range(5000,700001,5000))
    witnesses={s-c['formulas']['field_delta_steps'] for s in evaluations if s>=240000}
    # Keep EVERY evaluation and lag field, no lossy compression or evidence deletion.
    field_count=len(milestones|evaluations|witnesses)
    raw_fields=field_count*(8+182694*40)
    raw_sample_upper=len(milestones|evaluations)*200000
    reserve=384*2**20 # CSV, JSON, three ParaView fields, logs, final audits and spare space.
    required=raw_fields+raw_sample_upper+reserve
    result=dict(utc=utc(),status='PASS' if free>=required and m['effective_available']>=16*2**30 and not bad and not apps else 'BLOCKED_RESOURCE_PREFLIGHT',
      disk_free_bytes=free,disk_required_bytes=required,field_count=field_count,field_bytes=raw_fields,sampled_field_upper_bytes=raw_sample_upper,reserve_bytes=reserve,
      available_memory=m,other_compute_processes=bad,gpu_apps=apps,mpi_ranks=1,max_steps=700000,automatic_solver_retries=0)
    write(R/'provenance/RESOURCE_PREFLIGHT.json',result)
    if result['status']!='PASS':raise RuntimeError(json.dumps(result))
    return result
def step_from_history():
    h=R/'diagnostics/flow_history.csv'
    if not h.exists():return 0
    with h.open('rb') as f:
        f.seek(max(0,h.stat().st_size-8192));lines=f.read().splitlines()
    for l in reversed(lines):
        try:return int(l.split(b',')[0])
        except ValueError:pass
    return 0
def run():
    assert (R/'provenance/EXPLICIT_RETRY_AUTHORIZATION.json').exists(),'Explicit retry authorization required before any second solver launch'
    lock=(R/'provenance/vascular_run.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    assert not (R/'RUN_STARTED.json').exists(),'Formal run already claimed; no second solver'
    preflight()
    assert read(R/'provenance/STATIC_VALIDATION.json')['status']=='PASS'
    env=dict(os.environ,PATH='/usr/bin:/bin:/usr/sbin:/sbin',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',HWLOC_COMPONENTS='-opencl',CUDA_VISIBLE_DEVICES='0',STAGE4_LBM='batch',STAGE2_PROFILE_MODE='off',PYTHONDONTWRITEBYTECODE='1',LD_LIBRARY_PATH=str(R/'nvhpc/tbb/prefix/usr/lib/x86_64-linux-gnu')+':'+str(P/'compilers/lib'),TMPDIR=str(R/'tmp'),CUDA_CACHE_PATH=str(R/'tmp/cuda-cache'),XDG_CACHE_HOME=str(R/'tmp/cache'))
    env.pop('STAGE4_STRONG_CHECKPOINTS',None)
    cmd=['/usr/bin/mpirun','--allow-run-as-root','-np','1','--bind-to','core','--map-by','core','--report-bindings',str(EXE),str(R),str(R),'1.1197286861799598','700000']
    claim=dict(utc=utc(),command=cmd,from_zero=True,automatic_retries=0,binary_sha256=sha(EXE),freeze_sha256=sha(R/'FORMAL_STEP3C_GPU_FREEZE.json'),mpi_ranks=1)
    # Exclusive durable claim. No code path is allowed to remove/reuse it.
    with (R/'RUN_STARTED.json').open('x') as f:json.dump(claim,f,indent=2);f.flush();os.fsync(f.fileno())
    done=threading.Event();finish={};previous={};bindings={}
    with (R/'logs/solver.log').open('x') as log,(R/'gpu_performance_history.csv').open('x',newline='') as trace:
        names=['utc','elapsed_s','last_completed_step','gpu_utilization_percent','gpu_memory_MiB','CPU_percent','solver_RSS_bytes','available_RAM_bytes','disk_free_bytes']
        w=csv.DictWriter(trace,fieldnames=names);w.writeheader()
        evaluator_thread=threading.Thread(target=serve_evaluations,args=(R,done),daemon=True);evaluator_thread.start()
        start=time.monotonic();proc=subprocess.Popen(cmd,cwd=R,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
        def wait_child():
            finish['returncode']=proc.wait();finish['monotonic']=time.monotonic();done.set()
        th=threading.Thread(target=wait_child,daemon=True);th.start()
        while True:
            now=time.monotonic();ss={pid:v for pid,v in processes().items() if v['exe']==str(EXE)};cpu=0.
            for pid,v in ss.items():
                if pid in previous:cpu+=100*(v['ticks']-previous[pid][1])/os.sysconf('SC_CLK_TCK')/(now-previous[pid][0])
                previous[pid]=(now,v['ticks'])
                bindings[str(pid)]=sorted(os.sched_getaffinity(pid))
            try:g=gpu()
            except Exception as e:
                g=[None,None,None]
                with (R/'logs/resource_monitor_errors.log').open('a') as f:f.write(utc()+' '+str(e)+'\n')
            w.writerow(dict(utc=utc(),elapsed_s=now-start,last_completed_step=step_from_history(),gpu_utilization_percent=g[0],gpu_memory_MiB=g[1],CPU_percent=cpu,solver_RSS_bytes=sum(v['rss'] for v in ss.values()),available_RAM_bytes=memory()['effective_available'],disk_free_bytes=shutil.disk_usage(R).free));trace.flush()
            if done.is_set():break
            done.wait(5)
        th.join();rc=finish['returncode'];wall=finish['monotonic']-start
    statusfile=R/'diagnostics/solver_status.json';s=read(statusfile) if statusfile.exists() else {}
    if rc==0 and s.get('status') in ['AUTO_PASS_HUMAN_PENDING','NOT_CONVERGED_WITHIN_VALIDATION_HORIZON']:
        terminal='AUTO_CONVERGED' if s['auto_converged'] else 'NOT_CONVERGED_WITHIN_VALIDATION_HORIZON'
    elif (R/'diagnostics/GPU_NUMERICAL_INTEGRITY_FAILURE.json').exists():terminal='FAIL_GPU_NUMERICAL_INTEGRITY'
    elif (R/'diagnostics/solver_failure.txt').exists() and any(k in (R/'diagnostics/solver_failure.txt').read_text() for k in ['safety gate','FAIL_MPI_RUNTIME_CORRECTNESS']):terminal='FAIL_RUNTIME_SAFETY'
    else:terminal='RUNTIME_ERROR'
    write(R/'RUN_TERMINAL.json',dict(**claim,completed_utc=utc(),status=terminal,returncode=rc,actual_steps=s.get('timesteps',step_from_history()),end_to_end_seconds=wall,bindings=bindings,solver_status=s))
    print('FORMAL_RUN_TERMINAL',terminal,s.get('timesteps'),wall,flush=True)
    return 0 if rc==0 else 2
if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='--preflight':print(json.dumps(preflight(),indent=2))
    else:sys.exit(run())
