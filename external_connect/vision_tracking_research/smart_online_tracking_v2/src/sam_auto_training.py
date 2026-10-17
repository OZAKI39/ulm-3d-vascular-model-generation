"""Bounded, opt-out automatic SAM adaptation; no unreviewed labels or test data."""
from pathlib import Path
from contextlib import contextmanager
import json,os,time,hashlib,fcntl,subprocess,sys
import numpy as np
from annotation_store import load_annotation,sha
from annotation_regions import reviewed_boxes,review_identity

ROOT=Path(__file__).resolve().parents[1]
RUNTIME=ROOT/'data/auto_training'
DEFAULTS={'enabled':True,'max_steps':80,'max_wall_s':180,'cooldown_s':300,
          'debounce_s':15,'min_train_objects':8,'max_train_objects':64,
          'max_validation_objects':64,'min_validation_frames':2,'min_validation_objects':10,
          'learning_rate':1e-5,'seed':20261016,'context_width_px':128,'target_width_px':16,
    'class_name':'SINGLE_BUBBLE','protocol':'sam_decoder_v2_reviewed_regions'}

def atomic_json(path,value):
 path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
 tmp=path.with_name(path.name+f'.{os.getpid()}.tmp');tmp.write_text(json.dumps(value,ensure_ascii=False,indent=2));tmp.replace(path)

def read_json(path,default=None):
 try:return json.loads(Path(path).read_text())
 except (FileNotFoundError,ValueError):return default

def config(runtime=RUNTIME):
 # Numerical budgets are fixed here; the persisted user preference only enables/disables.
 return {**DEFAULTS,'enabled':bool(read_json(Path(runtime)/'config.json',{}).get('enabled',True))}

def set_enabled(enabled,runtime=RUNTIME):atomic_json(Path(runtime)/'config.json',{'enabled':bool(enabled)})

def snapshot(root=ROOT):
 root=Path(root);records=read_json(root/'data/candidate_manifest.json')['rows'];rows=[];excluded=[]
 for r in records:
  # Never even load test/quarantine masks for training or model admission.
  if r['split'] not in ('train','validation'):continue
  try:meta,masks=load_annotation(root/'data/annotations',r)
  except Exception as e:excluded.append({'id':r['id'],'reason':str(e)});continue
  if not meta or not reviewed_boxes(meta,masks):continue
  if 'SOFTWARE' in (meta.get('reviewer','')+' '+meta.get('notes','')).upper():continue
  if meta.get('split')!=r['split'] or meta.get('candidate_id')!=r['id']:
   excluded.append({'id':r['id'],'reason':'split/identity mismatch'});continue
  image=root/r['image']
  if sha(image)!=r['image_sha256']:
   excluded.append({'id':r['id'],'reason':'image SHA mismatch'});continue
  rows.append({'record':r,'meta':meta,'mask_path':str(root/'data/annotations'/r['id']/meta['mask_file'])})
 identity=[(r['record']['id'],r['record']['split'],r['record']['image_sha256'],r['meta']['mask_sha256'],review_identity(r['meta'])) for r in rows]
 code={str(p.relative_to(ROOT)):sha(p) for p in [ROOT/'scripts/sam_auto_train.py',ROOT/'src/sam_prompt_geometry.py',ROOT/'src/sam_auto_training.py',ROOT/'src/annotation_regions.py']}
 fingerprint=hashlib.sha256(json.dumps([DEFAULTS,identity,code],sort_keys=True).encode()).hexdigest()
 return {'fingerprint':fingerprint,'rows':rows,'excluded':excluded,'protocol':DEFAULTS,'created_at':time.time()}

def snapshot_still_valid(snap,root=ROOT):
 for row in snap['rows']:
  m=read_json(Path(root)/'data/annotations'/row['record']['id']/'latest.json',{})
  if m.get('mask_sha256')!=row['meta']['mask_sha256'] or review_identity(m)!=review_identity(row['meta']):return False
  try:
   if sha(row['mask_path'])!=row['meta']['mask_sha256']:return False
  except FileNotFoundError:return False
 return True

def build_samples(snap,split,root=ROOT):
 """Deterministic per-object crops; all geometry remains in original pixel units."""
 from scipy.ndimage import distance_transform_edt
 from sam_prompt_geometry import point_instance_requests
 groups=[];excluded=[];settings=snap['protocol']
 for row in snap['rows']:
  r=row['record'];m=row['meta']
  if r['split']!=split:continue
  if sha(row['mask_path'])!=m['mask_sha256']:raise ValueError('Frozen annotation checksum changed')
  with np.load(row['mask_path']) as z:masks={k:z[k].copy() for k in z.files}
  full_labels=masks['SINGLE_BUBBLE'];group=[];seen=set()
  # Context, targets and negative prompts are all contained in reviewed coverage.
  for region in reviewed_boxes(m,masks):
   ry0,rx0,ry1,rx1=region;offset=np.array([ry0,rx0])
   labels=full_labels[ry0:ry1,rx0:rx1]
   ignore=(masks['UNCERTAIN'][ry0:ry1,rx0:rx1]>0)|(masks['OUT_OF_FOCUS'][ry0:ry1,rx0:rx1]>0)
   ids=np.unique(labels);ids=ids[ids>0];points=[]
   for n in ids:
    distance=distance_transform_edt(labels==n);points.append([int(v) for v in np.unravel_index(np.argmax(distance),labels.shape)])
   for i,(n,point) in enumerate(zip(ids,points)):
    if int(n) in seen:continue
    try:plan=point_instance_requests(labels.shape,[point],[],None,settings['target_width_px'],settings['context_width_px'])[0]
    except ValueError as e:
     excluded.append({'id':r['id'],'instance':int(n),'reason':'invalid_local_context','detail':str(e),'review_box_yxyx':region});continue
    y0,x0,y1,x1=plan['crop_yxyx'];target=labels[y0:y1,x0:x1]==n
    reason=None
    if ignore[y0:y1,x0:x1].any():reason='ignore_region_in_context'
    if target.sum()!=np.count_nonzero(full_labels==n) or target[0].any() or target[-1].any() or target[:,0].any() or target[:,-1].any():reason='target_truncated'
    if reason:excluded.append({'id':r['id'],'instance':int(n),'reason':reason,'review_box_yxyx':region});continue
    neighbors=[(np.asarray(p)+offset).tolist() for j,p in enumerate(points) if j!=i and y0<=p[0]<y1 and x0<=p[1]<x1]
    group.append({'image_id':r['id'],'instance':int(n),'image':str(Path(root)/r['image']),
     'crop_yxyx':(np.asarray(plan['crop_yxyx'])+np.tile(offset,2)).tolist(),
     'point_yx':(np.asarray(point)+offset).tolist(),'negative_yx':neighbors,
     'box_yxyx':(np.asarray(plan['box_yxyx'])+np.tile(offset,2)).tolist(),
     'mask':target,'split':split,'annotation_revision':m['revision'],'review_box_yxyx':region})
    seen.add(int(n))
  groups.append(group)
 limit=settings['max_train_objects'] if split=='train' else settings['max_validation_objects'];samples=[]
 # Round-robin retains frame diversity rather than exhausting the first frame.
 for i in range(max((len(g) for g in groups),default=0)):
  for group in groups:
   if i<len(group):samples.append(group[i])
   if len(samples)>=limit:return samples,excluded
 return samples,excluded

def admission(metrics,settings=DEFAULTS):
 if metrics.get('validation_frames',0)<settings['min_validation_frames'] or metrics.get('validation_objects',0)<settings['min_validation_objects']:
  return False,'WAITING_FOR_INDEPENDENT_VALIDATION'
 if not metrics.get('finite',False):return False,'NONFINITE_VALIDATION'
 if metrics['candidate_mean_iou']<metrics['original_mean_iou']+.02:return False,'NO_CLEAR_GAIN_OVER_ORIGINAL'
 if metrics['candidate_mean_iou']<metrics['active_mean_iou']+.005:return False,'NO_GAIN_OVER_ACTIVE'
 if metrics['worst_frame_delta']<-.03:return False,'FRAME_REGRESSION'
 if metrics['candidate_failure_fraction']>metrics['active_failure_fraction']:return False,'MORE_FAILED_PROMPTS'
 return True,'VALIDATED_FOR_ASSISTED_SINGLE_BUBBLE_MASKS'

def active_model(runtime=RUNTIME):
 info=read_json(Path(runtime)/'active_model.json',{})
 return info if info.get('admitted') else {'id':'ORIGINAL_SAM_VIT_B','admitted':False}

def rollback(runtime=RUNTIME):
 atomic_json(Path(runtime)/'active_model.json',{'id':'ORIGINAL_SAM_VIT_B','admitted':False,'reason':'USER_ROLLBACK','updated_at':time.time()})

def attempt_blocks_retraining(fingerprint,last):
 # Pausing/cancelling must remain resumable without new annotations.
 return fingerprint==last.get('fingerprint') and last.get('status') in ('COMPLETED','FAILED')

@contextmanager
def gpu_slot(interactive=False,runtime=RUNTIME,check=None):
 """One GPU operation at a time; interactive jobs have priority over training steps."""
 runtime=Path(runtime);runtime.mkdir(parents=True,exist_ok=True);pending=runtime/'interactive';pending.mkdir(exist_ok=True)
 request=pending/str(os.getpid())
 if interactive:request.write_text(str(time.time()))
 handle=(runtime/'gpu.lock').open('a');locked=False
 try:
  while not locked:
   if check:check()
   waiting=False
   for p in pending.iterdir():
    try:os.kill(int(p.name),0);waiting=True
    except (ProcessLookupError,ValueError):p.unlink(missing_ok=True)
   if interactive or not waiting:
    try:fcntl.flock(handle,fcntl.LOCK_EX|fcntl.LOCK_NB);locked=True
    except BlockingIOError:pass
   if not locked:time.sleep(.04)
  yield
 finally:
  if locked:fcntl.flock(handle,fcntl.LOCK_UN)
  handle.close()
  if interactive:request.unlink(missing_ok=True)

def register_client(pid=None,runtime=RUNTIME):
 pid=pid or os.getpid();start=Path(f'/proc/{pid}/stat').read_text().split()[21]
 atomic_json(Path(runtime)/'clients'/f'{pid}.json',{'pid':pid,'start':start})

def unregister_client(pid=None,runtime=RUNTIME):
 (Path(runtime)/'clients'/f'{pid or os.getpid()}.json').unlink(missing_ok=True)

def has_clients(runtime=RUNTIME):
 live=False
 for p in (Path(runtime)/'clients').glob('*.json'):
  item=read_json(p,{})
  try:valid=Path(f'/proc/{item["pid"]}/stat').read_text().split()[21]==item['start']
  except (FileNotFoundError,KeyError):valid=False
  if valid:live=True
  else:p.unlink(missing_ok=True)
 return live

def ensure_service(runtime=RUNTIME):
 runtime=Path(runtime);runtime.mkdir(parents=True,exist_ok=True)
 with (runtime/'service.log').open('ab') as log:
  subprocess.Popen([sys.executable,'-u','-B',str(ROOT/'scripts/sam_auto_service.py'),'--runtime',str(runtime)],
                   stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
