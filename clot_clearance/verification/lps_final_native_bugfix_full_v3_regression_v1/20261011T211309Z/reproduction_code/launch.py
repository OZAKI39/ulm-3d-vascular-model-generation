from common import *
import os,sys,subprocess,time,psutil
assert json.loads((O/'STATIC_PREFLIGHT.json').read_text())['status']=='PASS';assert sha(B/'liblammps.so')==EXPECTED
# Atomic one-run lock: never restart this full case automatically.
with (O/'FORMAL_V3_LAUNCHED.json').open('x') as f:json.dump({'epoch':time.time(),'formal_count':1,'run':C['run'],'final_library_sha256':EXPECTED},f,indent=2)
env=os.environ.copy();env['LD_LIBRARY_PATH']=str(B)+os.pathsep+env.get('LD_LIBRARY_PATH','');env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',PYTHONUNBUFFERED='1')
wrapper=O/'rank_wrapper.sh';wrapper.write_text('#!/bin/bash\n'+sys.executable+' '+str(Path(__file__).with_name('driver.py'))+'\nclot_rank_rc=$?\nprintf \'%s\\n\' "$clot_rank_rc" > '+str(P)+'/rank_exit_${OMPI_COMM_WORLD_RANK}.txt\nexit "$clot_rank_rc"\n')
cmd=['mpirun','--bind-to','core','-np','2','bash',str(wrapper)]
dump(O/'RUN_COMMAND.json',{'argv':cmd,'cwd':str(P),'environment':{k:env[k] for k in ['LD_LIBRARY_PATH','PYTHONPATH','MPI4PY_LIBMPI','MPI4PY_MPIABI','OMP_NUM_THREADS']},'code_sha256':{str(p):sha(p) for p in Path(__file__).parent.glob('*') if p.is_file()},'observer_sha256':sha(O/'observer.so')})
t=time.perf_counter();peak=0;cpus={}
with (P/'launcher.log').open('w') as log:
 proc=subprocess.Popen(cmd,cwd=P,env=env,stdout=log,stderr=subprocess.STDOUT);pp=psutil.Process(proc.pid)
 while proc.poll() is None:
  try:
   rss=0
   for pr in [pp]+pp.children(recursive=True):
    try:rss+=pr.memory_info().rss;ct=pr.cpu_times();cpus[pr.pid]=ct.user+ct.system
    except psutil.Error:pass
   peak=max(peak,rss)
  except psutil.Error:pass
  time.sleep(.2)
 rc=proc.wait()
dump(O/'FORMAL_PERFORMANCE.json',{'exit_code':rc,'wall_s':time.perf_counter()-t,'peak_process_tree_RSS_bytes_sampled':peak,'sum_process_CPU_s_sampled':sum(cpus.values()),'CPU_by_PID_s':cpus,'rank_exit_codes':{p.name:p.read_text().strip() for p in P.glob('rank_exit_*.txt')},'formal_runs':1})
print('FORMAL_FINISHED',rc,flush=True);sys.exit(rc)
