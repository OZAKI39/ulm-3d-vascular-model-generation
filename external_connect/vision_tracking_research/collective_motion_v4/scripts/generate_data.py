"""Reuse F10U real-patch renderer with preregistered fresh procedural seeds."""
import sys,time,json
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from xcommon import *
import cv2
import synthetic as historical
P=json.loads((ROOT/'configs/research_protocol.json').read_text())
def one(job):
 split,case,seed,source_id,gray=job;cv2.setNumThreads(1)
 old=historical.case_params
 def params(name):
  v,A=old(name)
  if name in ['static','low_texture']:v*=0;A*=0
  elif name=='translation_x':v=np.array([18.,0.]);A*=0
  elif name=='translation_y':v=np.array([0.,18.]);A*=0
  elif name=='large_translation':v=np.array([240.,-120.]);A*=0
  return v,A
 historical.case_params=params
 try:
  frames_,meta=historical.build_sequence(case,seed,(cv2.cvtColor(gray,cv2.COLOR_GRAY2BGR),source_id*.0333),length=16)
 finally:historical.case_params=old
 if case=='low_contrast':frames_=np.uint8(180+.15*(frames_.astype(float)-180))
 if case=='low_texture':frames_[:]=180
 folder=OUT/'data'/split;folder.mkdir(parents=True,exist_ok=True);p=folder/f'{seed}_{case}.npz'
 np.savez_compressed(p,frames=frames_,H=meta.pop('H'),points=meta.pop('points'),generator=meta.pop('instantaneous_generator'),times=np.arange(16)*.0333)
 return dict(split=split,case=case,seed=seed,source_frame=source_id,file=str(p.relative_to(OUT)),sha256=sha(p),meta=meta,known_global_affine=case not in ['local_opposite'],observable_truth=case!='low_texture')
def main():
 start=time.perf_counter();assert sha(VIDEO)==VIDEO_SHA
 assert not (OUT/'data/DATA_COMPLETE.json').exists(),'REFUSE_OVERWRITE'
 sources={i:g for i,t,g in frames(list(P['source_frames'].values()))};allrows={}
 for split,seeds in P['seeds'].items():
  jobs=[(split,case,seed,P['source_frames'][split],sources[P['source_frames'][split]]) for seed in seeds for case in P['cases']]
  with ProcessPoolExecutor(max_workers=4) as pool:rows=list(pool.map(one,jobs))
  dump(OUT/'data'/split/'index.json',rows);allrows[split]=dict(n=len(rows),index_sha256=sha(OUT/'data'/split/'index.json'));print(split,len(rows),flush=True)
 dump(OUT/'data/DATA_COMPLETE.json',dict(identity=identity(),splits=allrows,wall_s=time.perf_counter()-start,protocol_sha256=sha(ROOT/'configs/research_protocol.json'),scope=P['independence_scope']))
if __name__=='__main__':main()
