#!/usr/bin/env python3
"""Resume only never-launched Stage4 horizons after archived-provenance filename repair."""
from pathlib import Path
import os,json,time,subprocess,fcntl,traceback,sys,hashlib
W=Path(__file__).resolve().parents[1];phase='STAGE4_PRELAUNCH_RECHECK'
f=(W/'validation.lock').open('a');fcntl.flock(f,fcntl.LOCK_EX|fcntl.LOCK_NB)
assert not any(W.glob('stage4_*/RUN_STARTED.json')),'Existing Stage4 solver may not be retried'
assert json.loads((W/'OFFICIAL_GPU_SMOKE.json').read_text())['status']=='PASS'
assert json.loads((W/'BUILD_PROVENANCE.json').read_text())['status']=='PASS'
assert json.loads((W/'provenance/PROVENANCE_FILENAME_COLLISION_REPAIR.json').read_text())['status']=='PASS'
with (W/'STAGE4_RESUME_AFTER_PRELAUNCH_REPAIR_STARTED.json').open('x') as out:json.dump(dict(pid=os.getpid(),unix=time.time(),stage4_prior_solver_runs=0,official_smoke_already_passed=True,automatic_solver_retries=0),out)
os.environ.update(json.loads((W/'toolchain_env/ENV.json').read_text()));os.environ['PYTHONDONTWRITEBYTECODE']='1'
attempts=json.loads((W/'provenance/EXECUTION_COMMANDS.json').read_text())
def save(name,v):
 p=W/name;t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(v,indent=2)+'\n');t.replace(p)
def command(name,cmd):
 global phase
 phase=name;save('VALIDATION_STATE.json',dict(state='RUNNING',phase=phase,pid=os.getpid(),unix=time.time(),run_id=W.name));print('START',name,flush=True)
 start=time.monotonic()
 with (W/'logs'/f'{name}.log').open('w') as f:p=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
 attempts.append(dict(stage=name,command=cmd,returncode=p.returncode,seconds=time.monotonic()-start,log='logs/'+name+'.log',scope='Resume never-launched solver after provenance-only collision repair'))
 save('provenance/EXECUTION_COMMANDS.json',attempts)
 if p.returncode:raise RuntimeError(name+' failed '+str(p.returncode))
 print('PASS',name,flush=True)
try:
 for n in (200,1000,5000):
  command(f'stage4_{n}_actual_run',['python3','-u','-B',str(W/'scripts/run_stage4.py'),str(n)])
  command(f'stage4_{n}_compare',['python3','-u','-B',str(W/'scripts/compare_restore.py'),str(n)])
 r=dict(state='PASS',phase='ALL_REQUESTED_RUNS_TERMINAL',unix=time.time(),pid=os.getpid(),run_id=W.name,automatic_solver_retries=0,official_smoke_repeated=False)
except BaseException as e:
 traceback.print_exc();r=dict(state='FAIL',phase=phase,error=repr(e),traceback=traceback.format_exc(),unix=time.time(),pid=os.getpid(),run_id=W.name,automatic_solver_retries=0)
save('VALIDATION_STATE.json',r);save('provenance/VALIDATION_TERMINAL.json',r);print('VALIDATION_TERMINAL',json.dumps(r),flush=True)
p=subprocess.run(['/usr/bin/python3','-u','-B',str(W/'scripts/finalize_result.py')]);sys.exit(p.returncode)
