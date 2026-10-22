"""After the CuPy ABI fix, prove selected kernels agree on the same windows."""
import sys,time,json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from xcommon import *
from backends import Backend,preprocess
from channel import Channel
import cv2
cv2.setNumThreads(2);images={i:g for i,t,g in frames([360,362])};channel=Channel();scan=channel.scan(images[360]);channel.scan(images[362]);a=preprocess(channel.cut(images[360]));b=preprocess(channel.cut(images[362]));rows=[]
for name in ['B2','B3','B4']:
 engine=Backend(name);regular=engine.core(a,b,actual_coordinates=True)
 for label,selection in [('all_safe',channel.safe),('activity',scan['selected']),('sparse_checker',channel.safe&(np.indices(channel.safe.shape).sum(0)%4==0))]:
  sparse=engine.core(a,b,selection=selection,actual_coordinates=True);good=selection&np.isfinite(regular['u']+regular['v']+sparse['u']+sparse['v']);diff=np.max(np.abs(np.c_[regular['u'][good]-sparse['u'][good],regular['v'][good]-sparse['v'][good]]))
  assert diff<1e-4,(name,label,diff)
  rows.append(dict(name=name,selection=label,n_selected=int(selection.sum()),max_component_difference=float(diff),regular_core_ms=regular['core_ms'],selected_core_ms=sparse['core_ms']))
dump(OUT/'results/sparse_equivalence.json',dict(identity=identity(),rows=rows,pass_=True,note='Correlation and displacement on identical selected windows. Spatial validation can differ because missing neighbours are correctly masked. Timings here are diagnostics, not warmed benchmark.'))
print('Sparse correlation equivalence PASS',rows)
