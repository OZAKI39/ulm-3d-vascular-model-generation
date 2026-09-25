#!/usr/bin/env python3
import pathlib,subprocess,os,sys,csv,json,time,hashlib,shutil,datetime,fcntl,signal,threading
sys.dont_write_bytecode=True
R=pathlib.Path(__file__).resolve().parents[1];B=R/'frozen_bundle';C=R/'verification/GPU_NUMERICAL_COMPARISON_CONTRACT.json'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(2**20),b''):h.update(b)
 return h.hexdigest()
def utc():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def save(p,v):p.write_text(json.dumps(v,indent=2))
def processes():
 rows={}
 for p in pathlib.Path('/proc').iterdir():
  if not p.name.isdigit():continue
  try:
   s=(p/'stat').read_text().split();rows[int(p.name)]=dict(name=(p/'comm').read_text().strip(),exe=os.readlink(p/'exe'),ticks=int(s[13])+int(s[14]),rss=int((p/'statm').read_text().split()[1])*os.sysconf('SC_PAGE_SIZE'))
  except (OSError,ValueError):pass
 return rows
def memory():
 info={l.split(':')[0]:int(l.split()[1])*1024 for l in pathlib.Path('/proc/meminfo').read_text().splitlines() if l.startswith(('MemAvailable:','MemTotal:'))}
 cg=pathlib.Path('/sys/fs/cgroup/memory.max');used=pathlib.Path('/sys/fs/cgroup/memory.current')
 if cg.exists() and cg.read_text().strip()!='max':info['effective_available']=min(info['MemAvailable'],int(cg.read_text())-int(used.read_text()))
 else:info['effective_available']=info['MemAvailable']
 return info
def gpu():return [float(x.strip()) for x in subprocess.check_output(['nvidia-smi','--query-gpu=utilization.gpu,memory.used,memory.total','--format=csv,noheader,nounits'],text=True).strip().split(',')]
def main():
 horizon=int(sys.argv[1]);assert horizon in (200,1000,5000)
 name='stage4_'+str(horizon);D=R/name;backend='gpu';mode='off';rep=1;strong=horizon==200;ranks=1
 lock=(R/'solver.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
 assert not (D/'RUN_STARTED.json').exists(),'Run already claimed: no automatic retry'
 assert json.loads((R/'OFFICIAL_GPU_SMOKE.json').read_text())['status']=='PASS'
 if horizon>200:assert json.loads((R/'verification'/('STAGE4_'+str(200 if horizon==1000 else 1000)+'_GATE.json')).read_text())['status']=='PASS'
 assert sha(C)=='5e78c61974be7ca2a52b3add288858330804c8d74d47ae909fced717d1f02456'
 for rel,digest in json.loads((R/'provenance/WORKING_COPY_SHA256.json').read_text()).items():
  assert sha(R/rel)==digest,'Frozen working input changed: '+rel
 for p in (D/'contracts').iterdir():
  assert sha(p)==sha(R/f'references/rtx{horizon}/contracts'/p.name)
 exe=R/'build/stage4/vascularPoC'
 assert sha(exe)==json.loads((R/'provenance/STAGE4_BINARY.json').read_text())['sha256']
 forbidden=[dict(pid=pid,**v) for pid,v in processes().items() if v['exe'].startswith('/workspace/hemocell') and v['name'].startswith(('vascular','inletCal','cavity','nvc'))]
 apps=subprocess.check_output(['nvidia-smi','--query-compute-apps=pid,process_name,used_memory','--format=csv,noheader'],text=True).strip();g=gpu();m=memory();free=shutil.disk_usage(R).free
 pre=dict(utc=utc(),available_RAM=m,disk_free=free,GPU=g,GPU_apps=apps,other_solver_or_compiler=forbidden)
 save(R/'provenance'/f'{name}_preflight.json',pre)
 assert not forbidden and not apps and g[0]<5,pre
 assert m['effective_available']>16*2**30 and free>256*2**20,pre
 cmd=['/usr/bin/mpirun','--allow-run-as-root','-np','1','--bind-to','core','--map-by','core','--report-bindings',str(exe),str(D),str(B),'1.1197286861799598',str(horizon)]
 env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',HWLOC_COMPONENTS='-opencl',CUDA_VISIBLE_DEVICES='0',STAGE2_PROFILE_MODE='off',STAGE4_LBM='batch')
 env.pop('STAGE4_STRONG_CHECKPOINTS',None)
 if strong:env['STAGE4_STRONG_CHECKPOINTS']='1'
 env.update(TMPDIR=str(R/'tmp'),CUDA_CACHE_PATH=str(R/'tmp/cuda-cache'),XDG_CACHE_HOME=str(R/'tmp/cache'))
 claim=dict(utc=utc(),name=name,backend=backend,variant='batch',strong_checkpoints=strong,command=cmd,mpi_ranks=1,from_zero=True,total_steps=horizon,warmup=(100 if horizon==200 else horizon//10),timed=horizon-(100 if horizon==200 else horizon//10),automatic_retries=0,comparison_contract_sha256=sha(C),binary_sha256=sha(exe),run_id=R.name,timing_scope='monotonic immediately before Popen through process exit with required outputs',env={k:env[k] for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','HWLOC_COMPONENTS','CUDA_VISIBLE_DEVICES','STAGE2_PROFILE_MODE','STAGE4_LBM']})
 with (D/'RUN_STARTED.json').open('x') as f:json.dump(claim,f,indent=2)
 print('RUN_START',name,utc(),flush=True)
 rows=[];bindings={};lastMessage=0;lastStep=0;previous={};reason=None
 with (D/'logs/solver.log').open('w') as log,(D/'provenance/resource_trace.csv').open('w',newline='') as trace:
  fields=['utc','elapsed_s','scope','gpu_utilization_percent','gpu_memory_used_MiB','gpu_memory_total_MiB','cpu_utilization_percent','solver_pids','solver_RSS_bytes','available_RAM_bytes','last_completed_step'];w=csv.DictWriter(trace,fieldnames=fields);w.writeheader()
  start=time.monotonic();proc=subprocess.Popen(cmd,cwd=B,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
  save(D/'RUN_PID.json',dict(launcher_pid=proc.pid,controller_pid=os.getpid(),start_monotonic=start));finished={};done=threading.Event()
  def wait_child():
   finished['returncode']=proc.wait();finished['monotonic']=time.monotonic();done.set()
  waiter=threading.Thread(target=wait_child,daemon=True);waiter.start()
  while True:
   now=time.monotonic();ps=processes();ss={pid:v for pid,v in ps.items() if v['exe']==str(exe)};cp=0
   for pid,v in ss.items():
    if pid in previous:cp+=100*(v['ticks']-previous[pid][1])/os.sysconf('SC_CLK_TCK')/(now-previous[pid][0])
    previous[pid]=(now,v['ticks'])
    try:
     vals=pathlib.Path('/proc',str(pid),'environ').read_bytes().split(b'\0');rank=int(next(v.split(b'=')[1] for v in vals if v.startswith(b'OMPI_COMM_WORLD_RANK=')));aff=sorted(os.sched_getaffinity(pid));cores=sorted(set((int(pathlib.Path(f'/sys/devices/system/cpu/cpu{c}/topology/physical_package_id').read_text()),int(pathlib.Path(f'/sys/devices/system/cpu/cpu{c}/topology/core_id').read_text())) for c in aff));bindings[rank]=dict(pid=pid,affinity=aff,physical_cores=cores)
    except (OSError,StopIteration):pass
   h=D/'diagnostics/flow_history.csv'
   if h.exists():
    with h.open('rb') as f:f.seek(max(0,h.stat().st_size-8192));tail=f.read().splitlines()
    for line in reversed(tail):
     try:lastStep=int(line.split(b',')[0]);break
     except ValueError:pass
   g=gpu();row=dict(utc=utc(),elapsed_s=now-start,scope=name,gpu_utilization_percent=g[0],gpu_memory_used_MiB=g[1],gpu_memory_total_MiB=g[2],cpu_utilization_percent=cp,solver_pids=';'.join(map(str,ss)),solver_RSS_bytes=sum(v['rss'] for v in ss.values()),available_RAM_bytes=memory()['effective_available'],last_completed_step=lastStep);w.writerow(row);trace.flush();rows.append(row)
   if now-lastMessage>15:print('RUN_PROGRESS',name,'elapsed',round(now-start,1),'step',lastStep,'GPU',g[0],'RSS_GiB',round(row['solver_RSS_bytes']/2**30,2),flush=True);lastMessage=now
   if proc.poll() is not None:break
   if shutil.disk_usage(R).free<192*2**20:
    reason='DISK_STOP_LOSS_192_MIB';os.killpg(proc.pid,signal.SIGTERM);break
   if now-start>3600 or lastStep>horizon:
    reason='WALL_TIME_STOP_LOSS_3600S' if now-start>3600 else 'STEP_CAP_VIOLATION';os.killpg(proc.pid,signal.SIGTERM);break
   done.wait(1)
  try:rc=proc.wait(timeout=30)
  except subprocess.TimeoutExpired:os.killpg(proc.pid,signal.SIGKILL);rc=proc.wait()
  waiter.join(timeout=5);wall=finished.get('monotonic',time.monotonic())-start
 sfile=D/'diagnostics/solver_status.json';s=json.loads(sfile.read_text()) if sfile.exists() else {};tf=D/'diagnostics/solver_timing.json';timing=json.loads(tf.read_text()) if tf.exists() else None
 core_sets=[set(map(tuple,v['physical_cores'])) for v in bindings.values()];bindingOK=len(bindings)==ranks and all(len(c)==1 for c in core_sets) and all(not c.intersection(d) for i,c in enumerate(core_sets) for d in core_sets[i+1:])
 verified=rc==0 and s.get('timesteps')==horizon and s.get('runtime_safety')=='PASS' and bindingOK and timing is not None
 receipt=dict(**claim,completed_utc=utc(),returncode=rc,status='PASS' if verified else 'FAIL',stop_reason=reason,completed_steps=s.get('timesteps',lastStep),end_to_end_seconds=wall,solver_timing=timing,binding_status='PASS' if bindingOK else 'FAIL',bindings=bindings,peak_solver_RSS_bytes=max(x['solver_RSS_bytes'] for x in rows),gpu_utilization_peak=max(x['gpu_utilization_percent'] for x in rows),gpu_vram_peak_MiB=max(x['gpu_memory_used_MiB'] for x in rows))
 initfile=D/'diagnostics/initialization_timing.json'
 init=json.loads(initfile.read_text()) if initfile.exists() else {}
 initialization=(init['first_iteration_monotonic_ns']/1e9-start) if init else None
 receipt.update(initialization_seconds=initialization,initialization_details=init,
  launcher_to_main_seconds=(init['main_entry_monotonic_ns']/1e9-start) if init else None,
  shutdown_and_remaining_overhead_seconds=wall-initialization-timing['full_iteration_seconds'] if verified else None,
  normal_benchmark=(mode=='off'),horizon=horizon,repetition=rep)
 save(D/'RUN_TERMINAL.json',receipt)
 sums=[]
 for p in sorted(D.rglob('*')):
  if p.is_file() and p.name!='SHA256SUMS':sums.append(sha(p)+'  '+str(p.relative_to(D)))
 (D/'SHA256SUMS').write_text('\n'.join(sums)+'\n');print('RUN_TERMINAL',json.dumps({k:receipt[k] for k in ['name','status','returncode','completed_steps','initialization_seconds','end_to_end_seconds','solver_timing','binding_status','gpu_vram_peak_MiB']}),flush=True)
 return 0 if verified else 2
if __name__=='__main__':sys.exit(main())
