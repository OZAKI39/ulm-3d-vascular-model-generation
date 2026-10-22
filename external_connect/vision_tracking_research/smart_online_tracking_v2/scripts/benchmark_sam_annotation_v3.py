"""Finite retrospective prompt test on eight confirmed TRAIN instances.

Prompts use the training annotations (simulated user hints), NOT held-out
recognition evaluation. No labels are changed and no model is trained.
"""
from pathlib import Path
import json,sys,time,hashlib,contextlib,datetime,subprocess
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'));sys.path.insert(0,str(ROOT/'src'))
from sam_worker import SamWorker
OUT=ROOT/'evidence/annotation_repair_v3_20261016/sam_comparison';OUT.mkdir(exist_ok=False)
p=ROOT/'evidence/annotation_repair_v3_20261016/labels_backup/video_000000/latest.json'
meta=json.loads(p.read_text());assert meta['expert_confirmed'] and meta['split']=='train'
with np.load(p.parent/meta['mask_file']) as z:labels=z['SINGLE_BUBBLE'].copy()
records=json.loads((ROOT/'data/candidate_manifest.json').read_text())['rows'];record=next(r for r in records if r['id']=='video_000000')
selected=[1,2,5,8,13,20,27,35]
plan={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'instance_ids':selected,
 'reference_revision':meta['revision'],'reference_mask_sha256':meta['mask_sha256'],
 'scope':'TRAIN_ONLY_RETROSPECTIVE_PROMPT_SIMULATION; not independent accuracy validation',
 'modes':['legacy_full_frame','local_point_default32','local_box_reference_bounds_plus3px'],
 'target_width_px':32,'context_width_px':128,'training_ready':False,
 'source_sha256':{name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in ['scripts/sam_worker.py','src/sam_prompt_geometry.py','scripts/benchmark_sam_annotation_v3.py']},
 'git_head':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()}
(OUT/'plan.json').write_text(json.dumps(plan,indent=2));worker=SamWorker();rows=[];start=time.monotonic()
for instance in selected:
 target=labels==instance;y,x=np.where(target);point=[float(y.mean()),float(x.mean())]
 box=[int(y.min()-3),int(x.min()-3),int(y.max()+4),int(x.max()+4)]
 for mode in plan['modes']:
  token=f'{instance:02d}_{mode}';req={'image':str(ROOT/record['image']),'image_id':record['id'],'token':token,
   'positive_yx':[point],'negative_yx':[],'output':str(OUT/(token+'.npz'))}
  if mode!='legacy_full_frame':req.update(mode='local_box',target_width_px=32,context_width_px=128)
  if mode.startswith('local_box'):req['box_yxyx']=box
  (OUT/(token+'_request.json')).write_text(json.dumps(req,indent=2))
  with (OUT/(token+'.log')).open('w') as log, contextlib.redirect_stdout(log):worker.run(req)
  with np.load(req['output']) as z:
   mask=z['mask'].astype(bool);details=json.loads(str(z['details']))
   iou=float((mask&target).sum()/max((mask|target).sum(),1));v=details['proposals'][0]
  rows.append({'instance':instance,'mode':mode,'area_px':int(mask.sum()),'reference_area_px':int(target.sum()),
               'iou_to_training_reference':iou,'allowed':v['allowed'],'wall_s':details['wall_s'],
               'cold_model':details['cold_model'],'cache_hit':details['embedding_cache_hit'],'reasons':v['reasons']})
  print(token,'area',int(mask.sum()),'IoU',round(iou,3),'s',round(details['wall_s'],3),flush=True)
summary={'plan':plan,'wall_s':time.monotonic()-start,'rows':rows,'summary':{}}
for mode in plan['modes']:
 r=[v for v in rows if v['mode']==mode];warm=[v['wall_s'] for v in r if not v['cold_model']]
 summary['summary'][mode]={'n':len(r),'allowed':sum(v['allowed'] for v in r),'area_median_px':float(np.median([v['area_px'] for v in r])),
   'median_iou_to_training_reference':float(np.median([v['iou_to_training_reference'] for v in r])),
   'warm_median_s':float(np.median(warm)),'warm_max_s':max(warm)}
(OUT/'summary.json').write_text(json.dumps(summary,indent=2))
# Actual mask overlays, no silhouette alteration or synthetic successful masks.
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
import imageio.v3 as iio
image=iio.imread(ROOT/record['image']);fig,axes=plt.subplots(4,4,figsize=(12,12))
for row,instance in enumerate(selected[:4]):
 target=labels==instance;y,x=np.where(target);cy,cx=np.mean(y),np.mean(x)
 for col,mode in enumerate(['reference']+plan['modes']):
  ax=axes[row,col];ax.imshow(image);mask=target
  if mode!='reference':
   with np.load(OUT/f'{instance:02d}_{mode}.npz') as z:mask=z['mask']
  if mask.any():ax.contour(mask,levels=[.5],colors=['#00a6ca'],linewidths=1)
  ax.set_xlim(cx-40,cx+40);ax.set_ylim(cy+40,cy-40);ax.set_xticks([]);ax.set_yticks([])
  ax.set_title(f'ID {instance} | '+['Human train reference','Old full-image point','Local point + box','Local explicit box'][col],fontsize=9)
fig.suptitle('Original SAM annotation assistance | Retrospective TRAIN prompts\nDisplay zoom only; not held-out recognition accuracy',fontsize=13)
fig.tight_layout();fig.savefig(OUT/'comparison.png',dpi=160);fig.savefig(OUT/'comparison.pdf');plt.close(fig)
print(json.dumps(summary['summary']),flush=True)
