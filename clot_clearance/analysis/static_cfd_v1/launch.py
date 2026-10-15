from common import *
import os,sys,subprocess,time,psutil,signal
fr=json.loads((O/'STATIC_LOAD_DESIGN_FREEZE.json').read_text());assert sha(B/'liblammps.so')==EXPECTED
with (O/'CAMPAIGN_LAUNCHED.json').open('x') as fp:json.dump({'epoch':time.time(),'freeze_sha256':sha(O/'STATIC_LOAD_DESIGN_FREEZE.json'),'physical_load_cases':1},fp)
env=os.environ.copy();env['LD_LIBRARY_PATH']=str(B)+os.pathsep+env.get('LD_LIBRARY_PATH','');env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',PYTHONUNBUFFERED='1')
recs=[]
for name,cmd,budget in [('interface',['mpiexec','--allow-run-as-root','-n','2',sys.executable,str(Path(__file__).with_name('interface.py'))],fr['budget']['interface_wall_seconds']),('solve',[sys.executable,str(Path(__file__).with_name('solve.py'))],fr['budget']['solve_process_wall_seconds'])]:
 dump(O/(name+'_RUN_COMMAND.json'),{'argv':cmd,'cwd':str(R),'process_environment':{k:env.get(k) for k in ['LD_LIBRARY_PATH','PYTHONPATH','MPI4PY_LIBMPI','OMP_NUM_THREADS','OPENBLAS_NUM_THREADS']},'freeze_sha256':sha(O/'STATIC_LOAD_DESIGN_FREEZE.json')});t=time.perf_counter();peak=0;stop=None
 with (O/(name+'_launcher.log')).open('w') as log:
  p=subprocess.Popen(cmd,cwd=R,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True);pr=psutil.Process(p.pid)
  while p.poll() is None:
   try:peak=max(peak,sum(q.memory_info().rss for q in [pr]+pr.children(recursive=True)))
   except psutil.Error:pass
   if time.perf_counter()-t>budget or peak>fr['budget']['process_tree_max_RSS_bytes']:
    stop='FROZEN_WALL_OR_MEMORY_BUDGET';os.killpg(p.pid,signal.SIGTERM)
    try:p.wait(timeout=5)
    except subprocess.TimeoutExpired:os.killpg(p.pid,signal.SIGKILL)
    break
   time.sleep(.1)
  rc=p.wait()
 recs.append({'phase':name,'exit_code':rc,'wall_s':time.perf_counter()-t,'process_tree_peak_RSS_bytes_sampled':peak,'budget_stop':stop});dump(O/'CAMPAIGN_PERFORMANCE.json',{'records':recs})
 if rc:sys.exit(rc)
print('STATIC_CFD_CAMPAIGN_FINISHED',recs,flush=True)
