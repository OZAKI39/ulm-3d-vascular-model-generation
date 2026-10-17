"""One bounded, reproducible SAM mask-decoder adaptation job from confirmed masks.

Uses micro-sam's trainable original SAM and its Dice + IoU-regression objective.
Frozen image embeddings match the annotation predictor's original-image transform.
No test-set access; missing independent validation never promotes the candidate.
"""
from pathlib import Path
import sys,json,time,hashlib,traceback,os,copy
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from sam_auto_training import atomic_json,read_json,build_samples,gpu_slot,admission
from annotation_store import sha

def run(job_path):
 import numpy as np,imageio.v3 as iio,torch
 from micro_sam.training.util import get_trainable_sam_model
 from segment_anything import SamPredictor
 from torch_em.loss import DiceLoss
 from sam_prompt_geometry import point_component,single_bubble_checks
 job=read_json(job_path);out=Path(job_path).parent;runtime=Path(job['runtime']);snap=job['snapshot'];cfg=snap['protocol']
 started=time.monotonic();steps=0;peak=0
 def check():
  if (out/'CANCEL').exists():raise InterruptedError('USER_PAUSE_OR_WINDOW_CLOSED')
  if time.monotonic()-started>cfg['max_wall_s']:raise TimeoutError('BOUNDED_WALL_TIME_EXCEEDED')
 def progress(stage,**extra):
  atomic_json(out/'progress.json',{'stage':stage,'step':steps,'max_steps':cfg['max_steps'],'wall_s':time.monotonic()-started,**extra})
 def serial(samples):return [{k:v for k,v in s.items() if k!='mask'} for s in samples]
 try:
  torch.set_num_threads(2);torch.manual_seed(cfg['seed']);np.random.seed(cfg['seed']);rng=np.random.default_rng(cfg['seed'])
  check()
  if not torch.cuda.is_available():raise RuntimeError('CUDA_NOT_AVAILABLE; automatic CPU training disabled')
  torch.cuda.set_per_process_memory_fraction(.40);torch.cuda.reset_peak_memory_stats()
  train,train_excluded=build_samples(snap,'train');val,val_excluded=build_samples(snap,'validation')
  if len(train)<cfg['min_train_objects']:raise ValueError('INSUFFICIENT_CONFIRMED_TRAINING_OBJECTS')
  atomic_json(out/'samples.json',{'train':serial(train),'validation':serial(val),'excluded':train_excluded+val_excluded,
    'scope':'prompt-assisted SINGLE_BUBBLE masks; not automatic detection accuracy',
    'prompts':'GT-derived interior points and neighboring object centers; fixed 16 px box; cached 128 px context',
    'test_data_used':False})
  base=ROOT/'models/sam_vit_b_01ec64.pth';base_sha=sha(base)
  progress('LOADING_MODEL',train_objects=len(train),validation_objects=len(val))
  with gpu_slot(runtime=runtime,check=check):
   model=get_trainable_sam_model('vit_b',checkpoint_path=str(base),device='cuda',freeze=['image_encoder','prompt_encoder'])
  model.eval();predictor=SamPredictor(model.sam)
  original={k:v.detach().cpu().clone() for k,v in model.sam.mask_decoder.state_dict().items()}
  active=job.get('active_model',{});active_state=original
  if active.get('admitted'):
   path=Path(active['checkpoint'])
   if sha(path)!=active['checkpoint_sha256']:raise ValueError('Active checkpoint checksum mismatch')
   checkpoint=torch.load(path,map_location='cpu',weights_only=True)
   if checkpoint['base_sha256']!=base_sha:raise ValueError('Active checkpoint base mismatch')
   active_state=checkpoint['decoder_state']
  # Warm-start from the current validated model; first run uses original SAM.
  model.sam.mask_decoder.load_state_dict(active_state,strict=True)
  cache={};images={}
  def prepared(sample,augment=False):
   check();key=(sample['image'],tuple(sample['crop_yxyx']));y0,x0,y1,x1=sample['crop_yxyx'];offset=np.array([y0,x0])
   if key not in cache:
    progress('ENCODING',encoded_crops=len(cache))
    if sample['image'] not in images:images.clear();images[sample['image']]=iio.imread(sample['image'])
    crop=images[sample['image']][y0:y1,x0:x1]
    with gpu_slot(runtime=runtime,check=check),torch.no_grad():
     predictor.set_image(crop);cache[key]=(predictor.features.detach().cpu(),predictor.input_size,predictor.original_size)
   features,input_size,size=cache[key]
   point=np.asarray(sample['point_yx'])-offset
   if augment:
    coords=np.argwhere(sample['mask']);point=coords[rng.integers(len(coords))]
   negatives=np.asarray(sample['negative_yx']).reshape(-1,2)-offset
   points=np.concatenate([point[None],negatives]);labels=np.r_[1,np.zeros(len(negatives))]
   lo=np.maximum(0,point-cfg['target_width_px']/2);hi=np.minimum(size,point+cfg['target_width_px']/2);box=np.r_[lo,hi]
   coords=predictor.transform.apply_coords(points[:,::-1].copy(),size)
   transformed_box=predictor.transform.apply_boxes(box[[1,0,3,2]][None],size)
   inp={'point_coords':torch.as_tensor(coords,dtype=torch.float32)[None],
        'point_labels':torch.as_tensor(labels,dtype=torch.int64)[None],
        'boxes':torch.as_tensor(transformed_box,dtype=torch.float32),'input_size':input_size,'original_size':size}
   return features,inp,point,negatives,box
  def evaluate(samples):
   rows=[]
   for s in samples:
    f,inp,point,negatives,box=prepared(s)
    with gpu_slot(runtime=runtime,check=check),torch.no_grad():
     pred=model([inp],f.to('cuda'),multimask_output=True)[0]
     masks=pred['masks'][0].gt(0).cpu().numpy();scores=pred['iou_predictions'][0].cpu().numpy()
    processed=[];checks=[]
    for m in masks:
     component,_=point_component(m,point);processed.append(component)
     checks.append(single_bubble_checks(component,[point],negatives,box,cfg['target_width_px']))
    order=sorted(range(3),key=lambda i:(not checks[i]['allowed'],-float(scores[i])))
    i=order[0];prediction=processed[i] if checks[i]['allowed'] else np.zeros_like(s['mask'])
    truth=s['mask'];iou=float((prediction&truth).sum()/max(1,(prediction|truth).sum()))
    rows.append({'image_id':s['image_id'],'instance':s['instance'],'iou':iou,'admitted':checks[i]['allowed']})
   return rows
  # Fitting diagnostics are explicitly training-set measurements, not validation.
  progress('TRAINING_FIT_BEFORE');fit_before=evaluate(train)
  optimizer=torch.optim.AdamW(model.sam.mask_decoder.parameters(),lr=cfg['learning_rate']);dice=DiceLoss(reduce_channel=None)
  model.sam.mask_decoder.train();losses=[]
  for step in range(cfg['max_steps']):
   check();s=train[int(rng.integers(len(train)))];f,inp,_,_,_=prepared(s,augment=True)
   with gpu_slot(runtime=runtime,check=check):
    optimizer.zero_grad(set_to_none=True);output=model([inp],f.to('cuda'),multimask_output=True)[0]
    prob=output['masks'].sigmoid();target=torch.as_tensor(s['mask'],device='cuda',dtype=torch.float32)[None,None]
    # Same best-of-three Dice + score-regression objective as micro-sam SamTrainer.
    mask_losses=torch.stack([dice(prob[:,i:i+1].swapaxes(0,1),target.swapaxes(0,1)).mean() for i in range(3)])
    with torch.no_grad():
     binary=prob>.5;truth=target.bool();true_iou=(binary&truth).sum((2,3))/(binary|truth).sum((2,3)).clamp_min(1)
    loss=mask_losses.min()+torch.nn.functional.mse_loss(output['iou_predictions'],true_iou)
    if not torch.isfinite(loss):raise FloatingPointError('Nonfinite loss')
    loss.backward();torch.nn.utils.clip_grad_norm_(model.sam.mask_decoder.parameters(),1.)
    optimizer.step();torch.cuda.synchronize();losses.append(float(loss.detach()))
   steps=step+1
   if steps%5==0:progress('TRAINING',loss=losses[-1])
  model.eval();candidate={k:v.detach().cpu().clone() for k,v in model.sam.mask_decoder.state_dict().items()}
  delta=float(sum((candidate[k]-active_state[k]).square().sum() for k in candidate).sqrt())
  if not np.isfinite(delta) or delta<=0:raise RuntimeError('Weights did not change')
  checkpoint=out/'decoder.pt';tmp=out/'decoder.tmp'
  torch.save({'schema':1,'base_sha256':base_sha,'decoder_state':candidate,'run_id':out.name},tmp);tmp.replace(checkpoint)
  progress('TRAINING_FIT_AFTER');fit_after=evaluate(train)
  metrics={'validation_frames':len({s['image_id'] for s in val}),'validation_objects':len(val),'finite':True}
  validation={}
  if metrics['validation_frames']>=cfg['min_validation_frames'] and len(val)>=cfg['min_validation_objects']:
   progress('VALIDATING');validation['candidate']=evaluate(val)
   model.sam.mask_decoder.load_state_dict(original);validation['original']=evaluate(val)
   if active.get('admitted'):
    model.sam.mask_decoder.load_state_dict(active_state);validation['active']=evaluate(val)
   else:validation['active']=validation['original']
   for key in validation:metrics[key+'_mean_iou']=float(np.mean([r['iou'] for r in validation[key]]))
   for key in ('candidate','active'):metrics[key+'_failure_fraction']=float(np.mean([r['iou']<.5 for r in validation[key]]))
   deltas=[]
   for frame in {s['image_id'] for s in val}:
    c=np.mean([r['iou'] for r in validation['candidate'] if r['image_id']==frame])
    a=np.mean([r['iou'] for r in validation['active'] if r['image_id']==frame]);deltas.append(float(c-a))
   metrics['worst_frame_delta']=min(deltas);metrics['finite']=all(np.isfinite(r['iou']) for rows in validation.values() for r in rows)
  passed,reason=admission(metrics,cfg)
  result={'status':'COMPLETED','admission_pass':passed,'admission_reason':reason,'steps':steps,
   'wall_s':time.monotonic()-started,'peak_cuda_allocated_mb':torch.cuda.max_memory_allocated()/2**20,
   'train_objects':len(train),'train_frames':len({s['image_id'] for s in train}),
   'trained_parameter_count':sum(p.numel() for p in model.parameters() if p.requires_grad),
   'frozen_modules':['image_encoder','prompt_encoder'],'decoder_delta_l2':delta,
   'checkpoint':str(checkpoint),'checkpoint_sha256':sha(checkpoint),'base_sha256':base_sha,
   'fingerprint':snap['fingerprint'],'losses':losses,'validation_metrics':metrics,'validation_rows':validation,
   'training_fit_before':fit_before,'training_fit_after':fit_after,'training_fit_is_generalization':False,
   'scientific_training_ready':False,'automatic_detection_validated':False,'test_data_used':False}
  atomic_json(out/'result.json',result);progress('COMPLETED',admission_reason=reason);print(json.dumps({k:v for k,v in result.items() if k not in ('losses','training_fit_before','training_fit_after','validation_rows')}),flush=True)
 except Exception as e:
  atomic_json(out/'result.json',{'status':'CANCELLED' if isinstance(e,InterruptedError) else 'FAILED','steps':steps,'wall_s':time.monotonic()-started,'error':str(e),'traceback':traceback.format_exc()})
  print(traceback.format_exc(),flush=True);return 1
 return 0

if __name__=='__main__':raise SystemExit(run(sys.argv[1]))
