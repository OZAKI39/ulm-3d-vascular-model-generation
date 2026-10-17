from pathlib import Path
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
ROOT=Path(__file__).resolve().parents[1]
fig,axes=plt.subplots(2,3,figsize=(12,10),layout='constrained')
for ax,i in zip(axes.flat,[120,126,132,138,360,720]):
 im=Image.open(ROOT/f'evidence/frame_{i:05d}.png');ax.imshow(im);ax.set_xlim(680,815);ax.set_ylim(510,350);ax.set_xticks(range(680,816,20));ax.set_yticks(range(350,511,20));ax.grid(alpha=.3);ax.set_title(f'Frame {i} | visual annotation only');ax.tick_params(labelsize=9)
fig.savefig(ROOT/'evidence/annotation_keyframes.png',dpi=150)
