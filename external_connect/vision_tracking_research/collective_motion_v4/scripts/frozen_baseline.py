"""Execute unchanged F10U/F10W function under unchanged F10W environment."""
import sys,time,json,importlib.metadata as md
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from xcommon import *
from methods import piv
import cv2
cv2.setNumThreads(2)
images={i:(t,g) for i,t,g in frames([90,92,360,362,700,702])};rows=[];mon=ResourceMonitor()
for region,roi in {'historical':[650,300,256,384],'medium':[600,240,600,624],'full_visible_bbox':[480,0,840,1120]}.items():
 x,y,w,h=roi
 for first in [90,360,700]:
  ta,a=images[first];tb,b=images[first+2];a=a[y:y+h,x:x+w];b=b[y:y+h,x:x+w]
  for rep in range(15):
   t=time.perf_counter();r=piv(a,b,tb-ta,mode='local');state=json.dumps(clean(r['estimate']),allow_nan=False);elapsed=(time.perf_counter()-t)*1000
   if rep>=3:rows.append(dict(name='B0_FROZEN_ENV_AND_PIPELINE',region=region,frame_a=first,rep=rep-3,end_to_end_ms=elapsed,algorithm_ms=r['algorithm_ms'],estimate=r['estimate'],valid_vectors=r['valid_vector_count']))
  print(region,first,elapsed,flush=True)
dump(OUT/'results/frozen_baseline.json',dict(identity=identity(),rows=rows,resources=mon.close(),versions={x:md.version(x) for x in ['numpy','scipy','openpiv']},scope='Unchanged historical function and environment. No decoder in pair benchmark.'))
