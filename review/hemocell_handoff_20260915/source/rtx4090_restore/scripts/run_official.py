#!/usr/bin/env python3
import os,sys,time,json,subprocess,hashlib,fcntl,csv,threading,sqlite3,signal,shutil
from pathlib import Path
W=Path(__file__).resolve().parents[1];D=W/'official_gpu_smoke';exe=W/'build/cavity3d'
def save(p,v):p.write_text(json.dumps(v,indent=2)+'\n')
lock=(W/'solver.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
assert not (D/'RUN_STARTED.json').exists(),'No automatic official solver retry'
assert not subprocess.check_output(['nvidia-smi','--query-compute-apps=pid','--format=csv,noheader'],text=True).strip()
nsys=shutil.which('nsys');assert nsys
cmd=['/usr/bin/mpirun','--allow-run-as-root','-np','1','--bind-to','core','--map-by','core','--report-bindings',nsys,'profile','--trace=cuda','--sample=none','--cpuctxsw=none','--export=sqlite','--output='+str(D/'official_cuda'),str(exe),'256','256','256','100','0.1','0','0','1']
env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',HWLOC_COMPONENTS='-opencl',CUDA_VISIBLE_DEVICES='0')
(D/'tmp').mkdir(exist_ok=True)
claim=dict(run_id=W.name,kind='official_cavity3d',max_steps=200,command=cmd,automatic_retries=0,binary_sha256=hashlib.sha256(exe.read_bytes()).hexdigest(),unix=time.time(),cuda_trace_purpose='Required direct GPU kernel verification only; no Stage4 profiler or optimization',normal_performance_benchmark=False)
with (D/'RUN_STARTED.json').open('x') as f:json.dump(claim,f,indent=2)
stop=threading.Event();rows=[]
def monitor():
 with (D/'GPU_RESOURCE_TRACE.csv').open('w',newline='') as f:
  w=csv.writer(f);w.writerow(['monotonic_seconds','gpu_utilization_percent','gpu_memory_used_MiB','gpu_memory_total_MiB'])
  while not stop.is_set():
   try:
    v=list(map(float,subprocess.check_output(['nvidia-smi','--query-gpu=utilization.gpu,memory.used,memory.total','--format=csv,noheader,nounits'],text=True,timeout=5).strip().split(',')))
    row=[time.monotonic(),*v];rows.append(row);w.writerow(row);f.flush()
   except Exception as e:print('RESOURCE_QUERY',repr(e),flush=True)
   stop.wait(.5)
t=threading.Thread(target=monitor,daemon=True);t.start();start=time.monotonic();reason=None
with (D/'solver.log').open('w') as log:
 p=subprocess.Popen(cmd,cwd=D,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
 save(D/'RUN_PID.json',dict(launcher_pid=p.pid,controller_pid=os.getpid(),start_monotonic=start))
 try:rc=p.wait(timeout=3600)
 except subprocess.TimeoutExpired:
  reason='WALL_TIME_STOP_3600';os.killpg(p.pid,signal.SIGTERM)
  try:rc=p.wait(timeout=30)
  except subprocess.TimeoutExpired:os.killpg(p.pid,signal.SIGKILL);rc=p.wait()
end=time.monotonic();stop.set();t.join(8)
log=(D/'solver.log').read_text();kernel=dict(status='UNVERIFIED');db=D/'official_cuda.sqlite'
if db.exists():
 con=sqlite3.connect('file:'+str(db)+'?mode=ro',uri=True)
 tables=[x[0] for x in con.execute("SELECT name FROM sqlite_master WHERE type='table'")]
 kn=[x for x in tables if x in ('CUPTI_ACTIVITY_KIND_KERNEL','CUPTI_ACTIVITY_KIND_CONCURRENT_KERNEL')]
 total=0;ns=0
 for table in kn:
  a,b=con.execute(f'SELECT count(*),sum(end-start) FROM {table}').fetchone();total+=a;ns+=b or 0
 kernel=dict(status='VERIFIED' if total>0 and ns>0 else 'UNVERIFIED',kernel_count=total,kernel_duration_sum_ns=ns,sqlite=str(db),tables=kn)
 con.close()
ok=rc==0 and 'Benchmark result:' in log and 'Energy check is OK' in log and 'Now running 100 warm-up iterations.' in log and 'Now running 100 benchmark iterations.' in log and kernel['status']=='VERIFIED'
r=dict(**claim,status='PASS' if ok else 'FAIL',returncode=rc,stop_reason=reason,wall_seconds=end-start,kernel_activity=kernel,gpu_utilization_peak=max((x[1] for x in rows),default=None),gpu_vram_peak_MiB=max((x[2] for x in rows),default=None),timing_includes_trace_export=True)
save(D/'RUN_TERMINAL.json',r);save(W/'OFFICIAL_GPU_SMOKE.json',r);print('OFFICIAL_TERMINAL',json.dumps(r),flush=True)
sys.exit(0 if ok else 2)
