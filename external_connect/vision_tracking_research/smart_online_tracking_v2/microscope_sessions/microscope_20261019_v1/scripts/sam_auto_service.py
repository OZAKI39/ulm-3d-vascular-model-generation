"""Single-instance, low-memory watcher for confirmed-label changes while annotating."""
from pathlib import Path
import sys,os,time,fcntl,subprocess,signal,argparse,traceback,json,hashlib
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from sam_auto_training import (RUNTIME,config,read_json,atomic_json,snapshot,snapshot_still_valid,
                               has_clients,active_model,attempt_blocks_retraining)

def serve(runtime):
 runtime=Path(runtime);runtime.mkdir(parents=True,exist_ok=True);lock=(runtime/'service.lock').open('a')
 try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
 except BlockingIOError:return
 child=None;folder=None;log=None;last=read_json(runtime/'last_attempt.json',{});pending=None;changed_at=0;start=0
 def status(state,**extra):
  atomic_json(runtime/'status.json',{'state':state,'service_pid':os.getpid(),'updated_at':time.time(),
    'active_model':active_model(runtime)['id'],'last_attempt':last,**extra})
 try:
  while has_clients(runtime):
   cfg=config(runtime)
   if child is not None:
    if not cfg['enabled']:(folder/'CANCEL').touch()
    if time.monotonic()-start>cfg['max_wall_s']+30:
     child.terminate()
     try:child.wait(5)
     except subprocess.TimeoutExpired:child.kill();child.wait()
    if child.poll() is None:
     progress=read_json(folder/'progress.json',{})
     status('PAUSING' if not cfg['enabled'] else 'RUNNING',job=str(folder),progress=progress);time.sleep(1);continue
    log.close();result=read_json(folder/'result.json',{'status':'FAILED','error':f'worker exited {child.returncode}'})
    job=read_json(folder/'job.json');valid=snapshot_still_valid(job['snapshot'])
    promoted=False
    if result.get('admission_pass') and valid and cfg['enabled']:
     info={'id':folder.name,'admitted':True,'classes':['SINGLE_BUBBLE'],'checkpoint':result['checkpoint'],
           'checkpoint_sha256':result['checkpoint_sha256'],'base_sha256':result['base_sha256'],
           'snapshot':job['snapshot'],'metrics':result['validation_metrics'],'updated_at':time.time()}
     atomic_json(runtime/'active_model.json',info);promoted=True
    last={'fingerprint':job['snapshot']['fingerprint'],'finished_at':time.time(),'job':str(folder),
          'process_wall_s':time.monotonic()-start,
          'status':result['status'],'reason':result.get('admission_reason',result.get('error')),
          'snapshot_still_valid':valid,'promoted':promoted}
    atomic_json(folder/'promotion.json',last);atomic_json(runtime/'last_attempt.json',last);child=None
   if not cfg['enabled']:status('PAUSED');time.sleep(2);continue
   snap=snapshot();active=active_model(runtime)
   if active.get('admitted') and not snapshot_still_valid(active['snapshot']):
    atomic_json(runtime/'active_model.json',{'id':'ORIGINAL_SAM_VIT_B','admitted':False,'reason':'LABEL_REVISION_INVALIDATED'})
   count=sum(r['meta'].get('reviewed_instance_counts',r['meta']['instance_counts'])['SINGLE_BUBBLE'] for r in snap['rows'] if r['record']['split']=='train')
   if count<cfg['min_train_objects']:status('WAITING_FOR_CONFIRMED_LABELS',confirmed_training_objects=count);time.sleep(5);continue
   if attempt_blocks_retraining(snap['fingerprint'],last):
    state='UP_TO_DATE' if last.get('promoted') and active_model(runtime).get('admitted') else 'CANDIDATE_HELD' if last.get('status')=='COMPLETED' else 'LAST_RUN_FAILED'
    status(state,confirmed_training_objects=count);time.sleep(5);continue
   if pending!=snap['fingerprint']:pending=snap['fingerprint'];changed_at=time.monotonic()
   cooldown=cfg['cooldown_s'] if last.get('status')=='COMPLETED' else 0
   wait=max(cfg['debounce_s']-(time.monotonic()-changed_at),cooldown-(time.time()-last.get('finished_at',0)))
   if wait>0:status('QUEUED',start_in_s=round(wait),confirmed_training_objects=count);time.sleep(2);continue
   # Do not reserve a new model when the GPU is already nearly full.
   try:
    free=int(subprocess.check_output(['nvidia-smi','--query-gpu=memory.free','--format=csv,noheader,nounits'],text=True,timeout=5).splitlines()[0])
   except Exception:free=0
   if free<3500:status('WAITING_FOR_GPU',free_gpu_mb=free);time.sleep(5);continue
   stamp=time.strftime('%Y%m%dT%H%M%S',time.gmtime());folder=ROOT/'results/sam_auto_training'/f'{stamp}_{snap["fingerprint"][:10]}'
   folder.mkdir(parents=True,exist_ok=False)
   sources=['scripts/sam_auto_train.py','scripts/sam_auto_service.py','src/sam_auto_training.py','src/sam_prompt_geometry.py','src/annotation_regions.py','src/annotation_store.py']
   identity={'head':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
    'source_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources}}
   atomic_json(folder/'job.json',{'runtime':str(runtime),'snapshot':snap,'active_model':active_model(runtime),'execution':identity})
   log=(folder/'worker.log').open('ab');start=time.monotonic()
   child=subprocess.Popen([sys.executable,'-u','-B',str(ROOT/'scripts/sam_auto_train.py'),str(folder/'job.json')],stdout=log,stderr=subprocess.STDOUT)
   status('RUNNING',job=str(folder));time.sleep(1)
 except Exception:
  status('SERVICE_ERROR',error=traceback.format_exc());print(traceback.format_exc(),flush=True)
 finally:
  if child is not None and child.poll() is None:
   (folder/'CANCEL').touch()
   try:child.wait(15)
   except subprocess.TimeoutExpired:child.terminate();child.wait(5)
  if log and not log.closed:log.close()
  status('STOPPED');lock.close()

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--runtime',default=str(RUNTIME));args=p.parse_args();serve(args.runtime)
