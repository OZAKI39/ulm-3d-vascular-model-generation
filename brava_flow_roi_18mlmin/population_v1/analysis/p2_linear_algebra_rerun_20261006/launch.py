from common import *
import subprocess,time,re

def run(case,mode):
 identity=read(BASE/'production_workflow/GPU_MANUFACTURED_RUN_IDENTITY.json');cmd=list(identity['command']);cmd[-2]=str(HERE/'runtime/svmultiphysics_audit_v4')
 cwd=HERE/'runtime'/case/'run';env=native_env();env.update(PETSC_OPTIONS=(cwd/'PETSC_OPTIONS.txt').read_text().strip(),BRAVA_ALGEBRA_MODE=mode,CUDA_VISIBLE_DEVICES='0')
 env.update(LD_PRELOAD=str(HERE/'runtime/cuda_admission.so'),BRAVA_CUDA_ADMISSION_LOCK=str(HERE/'runtime/cuda_admission.lock'),BRAVA_CUDA_ADMISSION_LOG=str(cwd/'cuda_admission.jsonl'))
 assert read(HERE/'CUDA_ADMISSION_GUARD.json')['status']=='PASS'
 print('RUN',case,mode,flush=True)
 proc=subprocess.Popen(cmd,cwd=cwd,env=env)
 if mode=='initial':
  if proc.wait()!=0:raise RuntimeError('Native initial-state read failed')
  return
 log=HERE/'monitoring'/read(HERE/'TASK.json')['name']/'job.log';offset=log.stat().st_size if log.exists() else 0
 while proc.poll() is None:
  time.sleep(.5)
  with log.open(errors='replace') as stream:stream.seek(offset);chunk=stream.read();offset=stream.tell()
  if re.search(r'DIVERGED_|PETSC ERROR|CUDA_ERROR|out of memory',chunk):
   managed().stop_tree(proc,{})
   raise RuntimeError('STOP_ON_FIRST_NATIVE_LINEAR_OR_BACKEND_FAILURE; see supervised log')
 if proc.returncode!=0:raise subprocess.CalledProcessError(proc.returncode,cmd)
