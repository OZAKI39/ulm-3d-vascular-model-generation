"""One stage-local single-GPU run from an explicit immutable WSL plan."""
import json,os,re,shutil,signal,subprocess,sys,time
import xml.etree.ElementTree as ET
from pathlib import Path
from runner_remote import *
from flow_parser import parse_solver_log
name=sys.argv[1];assert re.fullmatch(r'[A-Z0-9_]+',name)
plan=json.loads((BASE/'configs/runplans'/(name+'.json')).read_text())
policy=json.loads((BASE/'configs/policy.json').read_text())
ref=json.loads((BASE/'configs/reference_freeze.json').read_text())
sv=load(plan['build_report']);w=sv.get('PETSc_build',ref['PETSc_build']);N=BASE.parent/'sv1_3n'
mode=plan['mode'];assert mode in ('smoke','window','early','full')
assert plan['MPI_ranks']==plan['GPUs']==plan['OMP_NUM_THREADS']==1
assert mode!='full' or name=='REAL_VASCULAR_GPU_ILU_REUSE_WINNER'
case=BASE/'outputs'/name;assert not case.exists(),'Run evidence must never be overwritten'
case.mkdir();shutil.copyfile(BASE/'outputs/vascular_template.xml',case/'solver.xml')
tree=ET.parse(case/'solver.xml');g=tree.find('GeneralSimulationParameters')
changes={}
for tag,value in {'Number_of_time_steps':plan['end_step'],'Continue_previous_simulation':'true' if plan['start_step'] else 'false','Increment_in_saving_VTK_files':10,'Increment_in_saving_restart_files':10}.items():
 e=g.find(tag);assert e is not None;changes[tag]=dict(before=e.text,after=str(value));e.text=str(value)
tree.write(case/'solver.xml',encoding='utf-8',xml_declaration=True)
initial_sha=None
if plan['start_step']:
 assert (plan['start_step'],plan['end_step']) in ((60,70),(10,20))
 cp=N/'outputs/REAL_VASCULAR_GPU/1-procs/stFile_060.bin' if plan['start_step']==60 else BASE.parent/policy['early_checkpoint_remote']
 expected_sha=policy['checkpoint_sha256'] if plan['start_step']==60 else policy['early_checkpoint_sha256']
 assert digest(cp)==expected_sha
 (case/'1-procs').mkdir();shutil.copyfile(cp,case/'1-procs/stFile_last.bin');initial_sha=digest(case/'1-procs/stFile_last.bin')
options=plan['PETSC_OPTIONS']
options+=' -log_view :petsc_profile.txt -log_view_gpu_time'
for p in (Path.home()/'.petscrc',case/'.petscrc',case/'petscrc'):assert not p.exists()
dt=float(g.findtext('Time_step_size'));eq=tree.find('Add_equation');maxnl=int(eq.findtext('Max_iterations'));nltol=float(eq.findtext('Tolerance'))
envrun=env();envrun.update(LD_LIBRARY_PATH=sv['runtime_library_path'],PETSC_OPTIONS=options)
launcher=json.loads((N/'configs/baseline_L_mpi_application_gate.json').read_text())['working_launcher']
command=[w['candidate_wrapper'],launcher,'-n','1',w['candidate_wrapper'],sv['executable'],'solver.xml']
log=L/(name+'.log');assert not log.exists()
pre=subprocess.run(['nvidia-smi','--query-gpu=name,driver_version,memory.used,utilization.gpu','--format=csv,noheader,nounits'],capture_output=True,text=True,timeout=10)
assert pre.returncode==0 and '4090' in pre.stdout
start_utc=time.time();start=time.monotonic();failure=None;offset=0;carry='';last_sample=-20;samples=[];last_slow_check=-2
write(name+'_started',dict(name=name,plan=plan,start_unix_s=start_utc,initial_checkpoint_sha256=initial_sha,config_sha256=digest(case/'solver.xml'),initial_GPU_snapshot=pre.stdout.strip()))
with log.open('x') as out:
 process=subprocess.Popen(command,cwd=case,env=envrun,stdin=subprocess.DEVNULL,stdout=out,stderr=subprocess.STDOUT,start_new_session=True)
 while process.poll() is None:
  time.sleep(.1)
  if time.monotonic()-start-last_sample>=(2 if mode=='smoke' else 15):
   s=subprocess.run(['nvidia-smi','--query-gpu=timestamp,memory.used,utilization.gpu','--format=csv,noheader,nounits'],capture_output=True,text=True,timeout=5)
   last_sample=time.monotonic()-start;samples.append(dict(elapsed_s=last_sample,exit_code=s.returncode,raw=s.stdout.strip()))
   write(name+'_GPU_samples',dict(samples=samples,scope='2 s smoke /15 s window device-wide samples, not exact process peak'))
  with log.open() as f:f.seek(offset);chunk=f.read();offset=f.tell()
  joined=carry+chunk
  if '\n' in joined:complete,carry=joined.rsplit('\n',1)
  else:complete='';carry=joined
  fatal=re.search(r'DIVERGED_\w+|PETSC ERROR|MPI_ERR_TYPE|MPI_ABORT|Resetting restart flag|PC failed due to[^\n]*|(?<![a-z])(?:nan|inf)(?![a-z])',complete,re.I)
  if plan['candidate']=='RA' and fatal and fatal.group().startswith('DIVERGED_'):fatal=None
  if plan['candidate']=='RA' and 'recovery=FRESH_RETRY_FAILED' in complete:fatal=re.search(r'FRESH_RETRY_FAILED',complete)
  rows=parse_solver_log(complete,dt)
  if plan['start_step'] and any(r['step']<=plan['start_step'] for r in rows['linear_solves']):
   fatal=re.search(r'.+','CHECKPOINT_NOT_LOADED')
  badnl=[r for r in rows['linear_solves'] if r['nonlinear_iteration']>=maxnl and min(r['nonlinear_Ri_over_R0'] or 0,r['nonlinear_Ri_over_R1'] or 0)>nltol]
  heartbeat=case/'controller_heartbeat'
  monitor_lost=mode=='full' and time.monotonic()-start>120 and (not heartbeat.exists() or time.time()-heartbeat.stat().st_mtime>120)
  slow=None
  over_budget=mode=='smoke' and time.monotonic()-start>policy['smoke_resource_budget_s']
  if failure is None and slow:
   failure=dict(slow,elapsed_s=time.monotonic()-start);pending=case/'STOP_SIM.pending';pending.write_text('0\n');pending.replace(case/'STOP_SIM')
  if failure is None and (fatal or rows['failed_linear_solves'] or rows['nonlinear_failure_messages'] or badnl or monitor_lost or over_budget):
   failure=dict(reason=fatal.group() if fatal else 'CONTROLLER_LOST' if monitor_lost else 'SMOKE_RESOURCE_BUDGET' if over_budget else 'SOLVER_FAILURE',elapsed_s=time.monotonic()-start)
   pending=case/'STOP_SIM.pending';pending.write_text('0\n');pending.replace(case/'STOP_SIM')
   if not monitor_lost:os.killpg(process.pid,signal.SIGTERM)
  if failure and failure['reason']!='EARLY_REJECT_SLOW' and time.monotonic()-start-failure['elapsed_s']>30 and not monitor_lost:os.killpg(process.pid,signal.SIGKILL)
 code=process.wait()
elapsed=time.monotonic()-start
raw=log.read_text(errors='replace');history=parse_solver_log(raw,dt)
record=dict(status='EXECUTED',name=name,mode=mode,candidate=plan['candidate'],plan=plan,MPI_ranks=1,GPUs=1,OMP_NUM_THREADS=1,exit_code=code,wall_time_s=elapsed,start_unix_s=start_utc,end_unix_s=time.time(),command=command,PETSC_OPTIONS=options,profiling=True,initial_checkpoint_sha256=initial_sha,start_step=plan['start_step'],end_step=plan['end_step'],dt_s=dt,config_sha256=digest(case/'solver.xml'),configuration_changes=changes,solver_sha256=digest(sv['executable']),PETSc_library_sha256=digest(Path(w['prefix'])/'lib/libpetsc.so'),log='logs/'+log.name,log_sha256=digest(log),history=history,monitor_stop=failure,results=[dict(path=str(p.relative_to(BASE)),sha256=digest(p),bytes=p.stat().st_size) for p in sorted(case.glob('1-procs/result_*.vtu'))],checkpoints=[dict(path=str(p.relative_to(BASE)),sha256=digest(p),bytes=p.stat().st_size) for p in sorted(case.glob('1-procs/stFile_[0-9]*.bin'))],GPU_samples=samples,wall_clock_resolution_s=.1)
write(name+'_execution',record)
print(json.dumps({k:record[k] for k in ('name','exit_code','wall_time_s','monitor_stop')})+' saved='+str(len(record['results'])),flush=True)
raise SystemExit(0 if code==0 and failure is None else 1)
