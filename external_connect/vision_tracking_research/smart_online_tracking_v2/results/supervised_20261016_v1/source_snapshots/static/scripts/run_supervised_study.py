"""Bounded two-route study on frozen, reviewed ROIs; no test-set tuning."""
from pathlib import Path
import sys,json,time,signal,traceback,logging,copy,gc,os
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np,cv2,tifffile
from annotation_store import sha,CLASSES
from annotation_regions import reviewed_boxes
from vision_engine import VisionEngine,stats
from vision_evaluation import evaluate_instances,aggregate
from runtime import identity

FREEZE=ROOT/'results/expert_freezes/vision2_supervised_20261016_v1'
OUT=ROOT/'results/supervised_20261016_v1'

def dump(path,obj):path.write_text(json.dumps(obj,indent=2,allow_nan=False))

def run():
 OUT.mkdir(exist_ok=False);start=time.monotonic()
 frozen=json.loads((FREEZE/'frozen.json').read_text())
 for name,digest in frozen['snapshot_sha256'].items():
  if sha(FREEZE/name)!=digest:raise ValueError('Frozen file changed: '+name)
 frames=[];diameters=[]
 for row in frozen['rows']:
  with np.load(FREEZE/row['frozen_mask_file']) as z:masks={k:z[k].copy() for k in CLASSES}
  image=cv2.imread(str(FREEZE/row['frozen_image_file']))
  frames.append({'row':row,'image':image,'masks':masks,'boxes':reviewed_boxes(row['annotation'],masks)})
  if row['record']['split']=='train':
   for n in np.unique(masks['SINGLE_BUBBLE']):
    if n:diameters.append(2*np.sqrt(np.count_nonzero(masks['SINGLE_BUBBLE']==n)/np.pi))
 diameter=float(np.median(diameters))
 prereg={'scope':'Exploratory reviewed-ROI validation; not held-out test or cross-session generalization',
         'freeze':str(FREEZE),'freeze_sha256':sha(FREEZE/'frozen.json'),'execution':identity(ROOT),
         'primary_domain':'ROI crops; identical pixels for both routes','secondary_domain':'FULL_FRAME inference scored only in reviewed ROIs',
         'center_match_radius_px':6,'mask_iou':.5,'border_center_margin_px':6,
         'truncated_truth':'exclude and report; predictions overlapping excluded/ignore labels >=50% ignored',
         'candidate_class':'inherited area <=450 single candidate; larger regions unverified; not semantic group detection',
         'opencv_thresholds_train_only':[6,8,10],
         'cellpose':{'model':'cyto3 3.1.1.3','diameter_from_train_median_px':diameter,'min_size':5,'channels':[0,0],
                     'epochs':100,'batch_size':2,'learning_rate':.1,'SGD':True,'weight_decay':1e-4,'bsize':224,
                     'max_train_wall_s':180,'seed':20261016,'validation_passed_to_optimizer':False,
                     'finetune_trigger':'train center F1 <0.90 OR train mask F1 <0.80',
                     'selection':'final epoch only; prefer tuned if validation center F1 >= pretrained +0.02 and FP <= pretrained FP +2'},
         'no_test_data':True,'validation_source_frames':2,'group_and_track_accuracy':'NOT_EVALUATED',
         'automatic_control':False,'full_frame_accuracy':'NOT_EVALUATED; outside ROI unknown',
         'replay_frames':[360,419],'live_seconds_per_route':20,
         'references':['https://cellpose.readthedocs.io/en/v3.1.1.1/train.html','https://github.com/tryolabs/norfair']}
 dump(OUT/'preregistered.json',prereg);all_rows=[];raw=OUT/'predictions';raw.mkdir()
 def evaluate(engine,name,split,domain='ROI',persist=True):
  rows=[]
  for f in frames:
   r=f['row']['record']
   if r['split']!=split:continue
   full=f['masks']['SINGLE_BUBBLE'];counts={int(n):int(np.count_nonzero(full==n)) for n in np.unique(full) if n}
   if domain=='FULL_FRAME':
    tick=time.monotonic();objects,_,labels=engine.detect(f['image']);elapsed=1000*(time.monotonic()-tick)
   for i,box in enumerate(f['boxes']):
    y0,x0,y1,x1=box;crop=f['image'][y0:y1,x0:x1]
    if domain=='ROI':
     tick=time.monotonic();objs,_,labs=engine.detect(crop);ms=1000*(time.monotonic()-tick)
    else:
     labs=labels[y0:y1,x0:x1];objs=[dict(o,xy_local=[o['xy_local'][0]-x0,o['xy_local'][1]-y0]) for o in objects
                                 if x0<=o['xy_local'][0]<x1 and y0<=o['xy_local'][1]<y1];ms=elapsed
    metrics=evaluate_instances(full[y0:y1,x0:x1],labs,objs,
       (f['masks']['UNCERTAIN'][y0:y1,x0:x1]>0)|(f['masks']['OUT_OF_FOCUS'][y0:y1,x0:x1]>0),counts)
    key=f'{name}_{domain}_{r["id"]}_roi{i:02d}'
    row={'model':name,'split':split,'domain':domain,'image_id':r['id'],'roi_index':i,'box_yxyx':box,
         'detect_ms':ms,'metrics':metrics,'objects':objs,'artifact':key,'annotation_revision':f['row']['annotation']['revision']}
    if persist:np.savez_compressed(raw/(key+'.npz'),labels=labs,truth=full[y0:y1,x0:x1]);dump(raw/(key+'.json'),row)
    rows.append(row)
  return rows
 selections=[]
 for threshold in prereg['opencv_thresholds_train_only']:
  engine=VisionEngine('opencv',detector_config={'contrast_threshold':threshold})
  rows=evaluate(engine,f'opencv_threshold{threshold}','train');metric=aggregate(rows)
  selections.append({'threshold':threshold,'metrics':metric});all_rows+=rows
 chosen=max(selections,key=lambda r:(r['metrics']['center']['f1'] or 0,-r['metrics']['center']['fp'],-abs(r['threshold']-8)))
 opencv_cfg={'contrast_threshold':chosen['threshold']};dump(OUT/'opencv_selection.json',{'trials':selections,'chosen':opencv_cfg,'selection_split':'train'})
 engine=VisionEngine('opencv',detector_config=opencv_cfg)
 for split in ('train','validation'):
  for domain in ('ROI','FULL_FRAME'):all_rows+=evaluate(engine,'opencv',split,domain)
 cfg={'diameter':diameter,'min_size':5};engine=VisionEngine('cellpose',detector_config=cfg)
 # Warmup uses only a training image, before all timed inference.
 y0,x0,y1,x1=frames[0]['boxes'][0];engine.detect(frames[0]['image'][y0:y1,x0:x1])
 for split in ('train','validation'):
  for domain in ('ROI','FULL_FRAME'):all_rows+=evaluate(engine,'cellpose_pretrained',split,domain)
 base_train=aggregate([r for r in all_rows if r['model']=='cellpose_pretrained' and r['split']=='train' and r['domain']=='ROI'])
 dump(OUT/'baseline_rows.json',all_rows)
 training={'status':'NOT_RUN','reason':'predefined trigger not met'}
 if (base_train['center']['f1'] or 0)<.90 or (base_train['mask']['f1'] or 0)<.80:
  from cellpose import train
  import torch
  torch.manual_seed(20261016);np.random.seed(20261016)
  manifest=json.loads((FREEZE/'cellpose_train/export_manifest.json').read_text());imgs=[];labs=[]
  for item in manifest:
   if item['status']!='EXPORTED':continue
   imgs.append(cv2.imread(str(FREEZE/'cellpose_train'/(item['id']+'.png'))))
   labs.append(tifffile.imread(FREEZE/'cellpose_train'/(item['id']+'_masks.tif')))
  training={'status':'RUNNING','train_crops':len(imgs),'source_frames':2,'epochs':100,'test_data_used':False}
  dump(OUT/'training.json',training);tick=time.monotonic()
  def timeout(*_):raise TimeoutError('TRAIN_WALL_BUDGET_180S')
  previous=signal.signal(signal.SIGALRM,timeout);signal.alarm(180)
  try:
   path,loss,_=train.train_seg(engine.model.net,train_data=imgs,train_labels=labs,channels=[0,0],normalize=True,
           n_epochs=100,batch_size=2,learning_rate=.1,SGD=True,weight_decay=1e-4,bsize=224,
           min_train_masks=5,save_path=OUT,save_every=100,model_name='cyto3_bubble_reviewed_roi_100ep')
   training.update(status='COMPLETED',wall_s=time.monotonic()-tick,weights=str(path),weights_sha256=sha(path),
                   losses=loss.tolist(),cellpose_train_source_sha256=sha(Path(train.__file__)))
  except Exception:
   training.update(status='FAILED',wall_s=time.monotonic()-tick,error=traceback.format_exc())
  finally:signal.alarm(0);signal.signal(signal.SIGALRM,previous)
  dump(OUT/'training.json',training)
  if training['status']=='COMPLETED':
   del engine;gc.collect();torch.cuda.empty_cache()
   engine=VisionEngine('cellpose',detector_config={**cfg,'weights':training['weights']})
   engine.detect(frames[0]['image'][y0:y1,x0:x1])
   for split in ('train','validation'):
    for domain in ('ROI','FULL_FRAME'):all_rows+=evaluate(engine,'cellpose_tuned',split,domain)
 summaries={}
 for name in ('opencv','cellpose_pretrained','cellpose_tuned'):
  for split in ('train','validation'):
   for domain in ('ROI','FULL_FRAME'):
    rows=[r for r in all_rows if (r['model'],r['split'],r['domain'])==(name,split,domain)]
    if rows:summaries[f'{name}/{split}/{domain}']=aggregate(rows)
 base=summaries['cellpose_pretrained/validation/ROI']['center'];tuned=summaries.get('cellpose_tuned/validation/ROI',{}).get('center')
 use_tuned=bool(tuned and tuned['f1']>=base['f1']+.02 and tuned['fp']<=base['fp']+2)
 selected={**cfg,'weights':training['weights'],'weights_sha256':training['weights_sha256']} if use_tuned else cfg
 dump(OUT/'selected_opencv.json',opencv_cfg);dump(OUT/'selected_cellpose.json',selected)
 dump(OUT/'all_rows.json',all_rows)
 dump(OUT/'summary.json',{'metrics':summaries,'selected_cellpose':'cellpose_tuned' if use_tuned else 'cellpose_pretrained',
      'training_status':training['status'],'total_wall_s':time.monotonic()-start,'validation_is_not_test':True,
      'BUBBLE_DETECTION_VALIDATED':False,'SHORT_TRACKING_VALIDATED':False,'GROUP_SEGMENTATION_VALIDATED':False,
      'CONTROL_PERMISSION':False,'explanation':'Only two same-video validation frames; test, group and track labels unavailable.'})
 print(json.dumps({'selected_cellpose':'cellpose_tuned' if use_tuned else 'cellpose_pretrained','wall_s':time.monotonic()-start,
                   'validation_ROI':{k:v for k,v in summaries.items() if k.endswith('/validation/ROI')}},indent=2),flush=True)

if __name__=='__main__':
 logging.basicConfig(level=logging.INFO)
 try:run()
 except Exception:
  if OUT.exists():dump(OUT/'FAILED.json',{'error':traceback.format_exc()})
  raise
