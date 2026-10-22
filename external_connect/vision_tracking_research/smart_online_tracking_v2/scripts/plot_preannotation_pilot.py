"""Reproducible engineering comparison; no ground-truth accuracy claims."""
from pathlib import Path
import json,sys
import numpy as np
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from skimage.segmentation import find_boundaries
ROOT=Path(__file__).resolve().parents[1]
out=ROOT/'results/preannotation_pilot_20261016'
ids=['video_000000','video_000054','video_000120'];models=['cyto3','microsam_lm','cellpose4']
fig,axes=plt.subplots(3,4,figsize=(15,11),layout='constrained')
for r,ident in enumerate(ids):
    q=json.loads((out/f'cyto3_{ident}/request.json').read_text());y0,x0,y1,x1=q['roi_yxyx']
    crop=np.asarray(Image.open(q['image']).convert('RGB'))[y0:y1,x0:x1]
    axes[r,0].imshow(crop);axes[r,0].set_title(ident+' | original')
    for c,model in enumerate(models,1):
        folder=out/f'{model}_{ident}';m=json.loads((folder/'proposals.json').read_text());labels=np.load(folder/'proposals.npz')['labels']
        overlay=crop.copy();overlay[find_boundaries(labels,mode='inner')]=[0,255,255]
        axes[r,c].imshow(overlay);axes[r,c].set_title(f'{model}: {len(m["objects"])} proposals\n{m["worker_wall_s"]:.2f} s, rejected {len(m["rejected"])}')
for ax in axes.ravel():ax.axis('off')
fig.suptitle('Fixed training ROIs | Unreviewed proposals, not detection accuracy\nOriginal pixel coordinates; no image resizing during inference',fontsize=14)
fig.savefig(out/'comparison.png',dpi=150);fig.savefig(out/'comparison.pdf');plt.close(fig)
print(out/'comparison.png')
