#!/usr/bin/env python3
import pathlib,subprocess,os,sys,csv,json,time,hashlib,shutil,datetime,fcntl,signal,threading
sys.dont_write_bytecode=True
R=pathlib.Path(__file__).resolve().parents[1];B=R/'frozen_bundle';P=pathlib.Path('/workspace/hemocell_gpu_poc/toolchains/nvhpc/Linux_x86_64/26.5');C=R/'verification/GPU_NUMERICAL_COMPARISON_CONTRACT.json'
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
def materialize(D):
 D.mkdir()
 for name in ['contracts','diagnostics/field_samples','logs','provenance']:(D/name).mkdir(parents=True,exist_ok=True)
 # Hardlinks only to task-owned frozen copies; never to protected formal inputs.
 for p in (B/'frozen_contracts').iterdir():
  if p.is_file():os.link(p,D/'contracts'/p.name)
 lines=(B/'frozen_contracts/step3_solver_parameters.txt').read_text().splitlines();old=list(lines);nums=lines[2].split();assert nums[-1]=='1000';nums[-1]='700000';lines[2]=' '.join(nums)
 lines[0]='frozen_inputs/step1_geometry_contract/geometry/cfd_surface_axis_aligned_inlet_m.stl';lines[1]='frozen_inputs/step2'
 (D/'contracts/solver_parameters.txt').write_text('\n'.join(lines)+'\n');(D/'contracts/monitor_parameters.txt').write_text('100 5000 119872 240000 180543\n')
 assert lines[2].split()[:-1]==old[2].split()[:-1] and lines[3:]==old[3:]
 save(D/'provenance/input_verification.json',dict(status='PASS',copied_frozen_contracts={p.name:sha(p) for p in (B/'frozen_contracts').iterdir() if p.is_file()},materialized_run_files={p.name:sha(p) for p in (D/'contracts').iterdir()},path_rebinding_only=True,legacy_parameter_file_max_steps=700000,external_PoC_hard_cap=5000,physical_numeric_tokens_changed=False,formal_inputs_hardlinked=False))
def main():
 name=sys.argv[1]
 choices={'base_correctness':('gpu','off',200,1,'baseline',True),'batch_correctness':('gpu','off',200,1,'batch',True),'cpu_sanity':('cpu','off',200,1,'cpu',False),'batch_profile':('gpu','trace',200,1,'batch',False)}
 for n,reps in [(200,3),(1000,2),(5000,2)]:
  for rep in range(1,reps+1):
   for variant in ['baseline','batch']:choices[f'{variant}_{n}_r{rep}']=('gpu','off',n,rep,variant,False)
 assert name in choices,'Only frozen Stage4 cases allowed'
 backend,mode,horizon,rep,variant,strong=choices[name];D=R/'runs'/name
 assert not D.exists(),'Existing run claim: no restart/reuse/overwrite'
 lock=(R/'vascular_run.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
 assert sha(C)=='5e78c61974be7ca2a52b3add288858330804c8d74d47ae909fced717d1f02456'
 if name not in ['base_correctness','batch_correctness']:
  assert json.loads((R/'verification/INITIAL_CORRECTNESS_GATE.json').read_text())['status']=='PASS'
 if name=='batch_correctness':assert json.loads((R/'runs/base_correctness/RUN_TERMINAL.json').read_text())['status']=='PASS'
 if horizon>=1000:assert json.loads((R/'verification/SHORT_GATE.json').read_text())['status']=='PASS'
 if horizon>=5000:assert json.loads((R/'verification/HORIZON_1000_GATE.json').read_text())['status']=='PASS'
 frozen=json.loads((R/'provenance/execution_hashes.json').read_text())
 for rel,digest in frozen.items():assert sha(R/rel)==digest,'Frozen execution changed '+rel
 for rel,digest in json.loads((R/'provenance/BASELINE_SOURCE_SHA256.json').read_text()).items():assert sha(R/'baseline'/rel)==digest
 exe=R/'build_gpu/vascular/vascularPoC' if backend=='gpu' else R.parent/'20260914_stage3_sustained/build_cpu/vascular/vascularPoC'
 assert sha(exe)==json.loads((R/'provenance/binary_hashes.json').read_text())[backend]
 ps0=processes();t0=time.monotonic();time.sleep(1);ps1=processes();elapsed=time.monotonic()-t0
 busy=[]
 for pid,v in ps1.items():
  v['cpu_percent']=100*(v['ticks']-ps0.get(pid,v)['ticks'])/os.sysconf('SC_CLK_TCK')/elapsed
  if pid!=os.getpid() and v['cpu_percent']>50:busy.append(dict(pid=pid,**v))
 forbidden=[dict(pid=pid,**v) for pid,v in ps1.items() if v['exe'].startswith('/workspace/hemocell') and (v['name'].startswith(('vascular','inletCal','cavity')) or v['name'] in ['nvc++','nvc'])]
 g=gpu();apps=subprocess.check_output(['nvidia-smi','--query-compute-apps=pid,process_name,used_memory','--format=csv,noheader'],text=True).strip();m=memory();free=shutil.disk_usage(R).free
 pre=dict(utc=utc(),available_RAM=m,disk_free=free,load_average=os.getloadavg(),GPU=g,GPU_apps=apps,other_busy_processes=busy,other_solver_or_compiler=forbidden)
 save(R/'provenance'/f'{name}_preflight.json',pre)
 assert not busy and not forbidden and not apps and g[0]<5,pre
 assert m['effective_available']>16*2**30 and free>140*2**20,pre
 materialize(D)
 ranks=12 if backend=='cpu' else 1
 cmd=['/usr/bin/mpirun','--allow-run-as-root','-np',str(ranks),'--bind-to','core','--map-by','core','--report-bindings',str(exe),str(D),str(B),'1.1197286861799598',str(horizon)]
 # Same launcher environment for both backends. Prevent unrelated OpenCL topology probing;
 # this does not disable CUDA execution or change CPU binding.
 env=dict(os.environ,PATH='/usr/bin:/bin:/usr/sbin:/sbin',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',HWLOC_COMPONENTS='-opencl',CUDA_VISIBLE_DEVICES='0',LD_LIBRARY_PATH=str(R/'nvhpc/tbb/prefix/usr/lib/x86_64-linux-gnu')+':'+str(P/'compilers/lib'))
 env['STAGE2_PROFILE_MODE']=mode
 env['STAGE4_LBM']=variant
 if strong:env['STAGE4_STRONG_CHECKPOINTS']='1'
 if mode=='trace':
  nsys=P/'profilers/13.2/Nsight_Systems/bin/nsys'
  pos=cmd.index(str(exe));cmd[pos:pos]=[str(nsys),'profile','--trace=cuda,nvtx,mpi,osrt','--mpi-impl=openmpi','--sample=none','--cpuctxsw=none','--cuda-um-cpu-page-faults=true','--cuda-um-gpu-page-faults=true','--capture-range=cudaProfilerApi','--capture-range-end=stop','--export=sqlite','--output='+str(R/'profiling/batch_101_200')]
 env.update(TMPDIR=str(R/'tmp'),CUDA_CACHE_PATH=str(R/'tmp/cuda-cache'),XDG_CACHE_HOME=str(R/'tmp/cache'))
 claim=dict(utc=utc(),name=name,backend=backend,variant=variant,strong_checkpoints=strong,command=cmd,mpi_ranks=ranks,from_zero=True,total_steps=horizon,warmup=(100 if horizon==200 else horizon//10),timed=horizon-(100 if horizon==200 else horizon//10),automatic_retries=0,comparison_contract_sha256=sha(C),binary_sha256=sha(exe),timing_scope='monotonic immediately before Popen through successful process exit and required outputs',env={k:env[k] for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','HWLOC_COMPONENTS','CUDA_VISIBLE_DEVICES','STAGE2_PROFILE_MODE']})
 save(D/'RUN_STARTED.json',claim);print('RUN_START',name,utc(),flush=True)
 rows=[];bindings={};lastMessage=0;lastStep=0;previous={};reason=None
 with (D/'logs/solver.log').open('w') as log,(D/'provenance/resource_trace.csv').open('w',newline='') as trace:
  fields=['utc','elapsed_s','scope','gpu_utilization_percent','gpu_memory_used_MiB','gpu_memory_total_MiB','cpu_utilization_percent','solver_pids','solver_RSS_bytes','available_RAM_bytes','last_completed_step'];w=csv.DictWriter(trace,fieldnames=fields);w.writeheader()
  start=time.monotonic();proc=subprocess.Popen(cmd,cwd=B,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
  finished={};done=threading.Event()
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
   if shutil.disk_usage(R).free<64*2**20:
    reason='DISK_STOP_LOSS_64_MIB';os.killpg(proc.pid,signal.SIGTERM);break
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
