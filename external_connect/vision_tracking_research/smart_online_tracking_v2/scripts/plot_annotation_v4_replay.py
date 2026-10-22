"""Visual audit of the reported image and prompts; no ground-truth claim."""
from pathlib import Path
import sys,json
import numpy as np
import imageio.v3 as iio
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from sam_batch_review import assemble_candidates
OUT=ROOT/'evidence/annotation_repair_v4_20261016'
req=json.loads((OUT/'reported_job/request.json').read_text());image=iio.imread(req['image']);pts=np.asarray(req['positive_yx'])
with np.load(OUT/'reported_job/candidate.npz') as z:old=z['masks'][0].copy()
with np.load(OUT/'component_replay/candidate.npz') as z:
 details=json.loads(str(z['details']));masks=[z[t['array_key']].copy() for t in details['targets']]
choices=[next((j for j,p in enumerate(t['proposals']) if p['allowed']),0) for t in details['targets']]
labels,items,rejected=assemble_candidates(masks,details['targets'],choices,[True]*10,image.shape[:2])
y0,x0,y1,x1=np.asarray(req['box_yxyx']).astype(int)+[-10,-10,10,10]
colors=plt.get_cmap('tab10')(np.arange(10))
fig,axes=plt.subplots(1,3,figsize=(11,5.6),layout='constrained')
for ax in axes:
 ax.imshow(image);ax.set_xlim(x0,x1);ax.set_ylim(y1,y0);ax.set_xticks([]);ax.set_yticks([])
 for i,(y,x) in enumerate(pts):
  ax.plot(x,y,'o',ms=3,color='lime',mec='black',mew=.3)
  ax.text(x+3,y-3,str(i+1),color='black',fontsize=8,bbox={'facecolor':'white','alpha':.65,'pad':.2,'edgecolor':'none'})
axes[0].set_title('Reported image and 10 points',fontsize=11)
overlay=np.zeros((*old.shape,4));overlay[old]=[.65,.28,.14,.6];axes[1].imshow(overlay)
axes[1].set_title(f'Old shared-object prompt\n{old.sum():,} pixels in one mask',fontsize=11)
overlay=np.zeros((*labels.shape,4))
for i in items:overlay[labels==i+1]=[*colors[i,:3],.55]
axes[2].imshow(overlay);axes[2].set_title('Independent point prompts\n9 candidates; point 7 rejected',fontsize=11)
fig.suptitle('Original SAM ViT-B | Candidate review, not an accuracy benchmark',fontsize=12)
fig.savefig(OUT/'reported_points_comparison.png',dpi=180);fig.savefig(OUT/'reported_points_comparison.pdf');plt.close(fig)
fig,axes=plt.subplots(2,5,figsize=(12,5),layout='constrained')
for i,ax in enumerate(axes.flat):
 y,x=pts[i];t=details['targets'][i];c=choices[i];yy,xx,_,_=t['crop_yxyx']
 ax.imshow(image);ax.contour(np.arange(masks[i].shape[2])+xx,np.arange(masks[i].shape[1])+yy,masks[i][c],levels=[.5],colors=['#e06000' if i in rejected else '#136cce'],linewidths=1)
 ax.plot(x,y,'+',color='lime',ms=7);ax.set_xlim(x-20,x+20);ax.set_ylim(y+20,y-20)
 ax.set_xticks([]);ax.set_yticks([]);ax.set_title(f'Point {i+1}: '+('REJECTED' if i in rejected else 'REVIEW'),fontsize=10)
fig.suptitle('Selected candidates in original pixels; boundaries require human review',fontsize=12)
fig.savefig(OUT/'per_point_review.png',dpi=180);plt.close(fig)
print(json.dumps({'old_mask_area':int(old.sum()),'admissible_candidates':len(items),'rejected_point_numbers':[i+1 for i in rejected]}))
