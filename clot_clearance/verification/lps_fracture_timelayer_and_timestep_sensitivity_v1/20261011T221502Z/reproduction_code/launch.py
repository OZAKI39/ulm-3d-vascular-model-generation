from common import *
import os,sys,subprocess,time,psutil
freeze=json.loads((O/'TIMESTEP_DESIGN_FREEZE.json').read_text());assert sha(B/'liblammps.so')==EXPECTED
with (O/'CAMPAIGN_LAUNCHED.json').open('x') as f:json.dump({'epoch':time.time(),'freeze_sha256':sha(O/'TIMESTEP_DESIGN_FREEZE.json'),'max_cases':2},f)
env=os.environ.copy();env['LD_LIBRARY_PATH']=str(B)+os.pathsep+env.get('LD_LIBRARY_PATH','');env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',PYTHONUNBUFFERED='1')
recs=[]
for case in freeze['formal_cases']:
 name=case['name'];out=O/name;cmd=[sys.executable,str(Path(__file__).with_name('driver.py')),name]
 dump(out/'RUN_COMMAND.json',{'argv':cmd,'cwd':str(R),'process_environment':{k:env.get(k) for k in ['LD_LIBRARY_PATH','PYTHONPATH','MPI4PY_LIBMPI','OMP_NUM_THREADS']},'freeze_sha256':sha(O/'TIMESTEP_DESIGN_FREEZE.json')})
 t=time.perf_counter();peak=0
 with (out/'launcher.log').open('w') as log:
  p=subprocess.Popen(cmd,cwd=R,env=env,stdout=log,stderr=subprocess.STDOUT);pr=psutil.Process(p.pid)
  while p.poll() is None:
   try:peak=max(peak,sum(q.memory_info().rss for q in [pr]+pr.children(recursive=True)))
   except psutil.Error:pass
   time.sleep(.05)
  rc=p.wait()
 rec={'case':name,'exit_code':rc,'wall_s':time.perf_counter()-t,'process_tree_peak_RSS_bytes_sampled':peak};recs.append(rec);dump(O/'CAMPAIGN_PERFORMANCE.json',{'attempted_cases':len(recs),'records':recs})
 if rc:sys.exit(rc)
print('TWO_CASES_FINISHED',recs,flush=True)
