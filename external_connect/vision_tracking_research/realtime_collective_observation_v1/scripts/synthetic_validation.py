"""Frozen new seeds; mathematical core and independently measured real-patch images."""
import sys,os,ast,argparse,json,time
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from ycommon import ROOT,dump,identity,sha
import numpy as np,cv2,av
from aggregation import Aggregation
from channel import Channel
from backends import Backend
# Execute the exact pre-existing bounded splat function without importing training.
source=ROOT.parent/'collective_motion_v2/src/dataset.py';tree=ast.parse(source.read_text());node=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='splat');scope={'cv2':cv2,'np':np};exec(compile(ast.Module(body=[node],type_ignores=[]),str(source),'exec'),scope);splat=scope['splat']
CASES=['converge','translation','converge_translate','rotation','appearance_loss','local_converge','defocus','form_cluster','split','onset']

def trajectories(case,seed,t):
 rng=np.random.default_rng(seed);q0=rng.uniform(85,425,(64,2));center=np.array([256.,256.]);q=q0.copy();vel=np.zeros_like(q);focus=1.;onset=0.
 if case in ['converge','converge_translate','defocus','form_cluster','local_converge','onset']:
  rate=.25 if case!='form_cluster' else .5;onset=1 if case=='onset' else 0;tt=max(t-onset,0);sel=np.ones(64,bool) if case!='local_converge' else q0[:,0]<256
  if case=='local_converge':center=np.array([170.,256.])
  q[sel]=center+(q0[sel]-center)*np.exp(-rate*tt);vel[sel]=-rate*(q[sel]-center) if t>=onset else 0
 if case in ['translation','converge_translate']:q+=np.array([15,4])*t;vel+=np.array([15,4])
 if case=='rotation':
  a=.3*t;rot=np.array([[np.cos(a),-np.sin(a)],[np.sin(a),np.cos(a)]]);q=center+(q0-center)@rot.T;vel=(q-center)@np.array([[0,.3],[-.3,0]])
 if case=='appearance_loss':
  keep=np.arange(64)<64-int(15*t);q=q[keep];vel=vel[keep]
 if case=='split':vel[:,0]=np.where(q0[:,0]<256,-14,14);q+=vel*t
 if case=='defocus':focus=np.exp(-1.5*t) # manufactured apparent contraction plus visibility failure
 return q,vel,focus

def main():
 p=argparse.ArgumentParser();p.add_argument('--video',required=True);p.add_argument('--backend',default='B4');p.add_argument('--output',required=True);p.add_argument('--split',choices=['development','frozen_test'],default='development');a=p.parse_args();out=Path(a.output);out.mkdir(parents=True,exist_ok=True);cv2.setNumThreads(2)
 with av.open(a.video) as c:f=next(c.decode(video=0)).to_ndarray(format='gray')
 ch=Channel();q=ch.scan(f)['candidates'];local=ch.cut(f);bank=[]
 for x,y in np.rint(q).astype(int):
  if 12<=x<local.shape[1]-12 and 12<=y<local.shape[0]-12:
   patch=local[y-10:y+11,x-10:x+11].astype(np.float32);patch=np.maximum(np.percentile(patch,85)-patch,0)*np.outer(np.hanning(21),np.hanning(21))
   if patch.max()>20:bank.append(patch)
 bank=np.array(bank);np.savez_compressed(out/'real_patch_bank.npz',patches=bank)
 cfg=dict(image_wh=[512,512],crop_xywh=[0,0,512,512],channel_polygon=[[0,0],[511,0],[511,511],[0,511]],excluded_polygons=[],minimum_valid_window_fraction=1.,registration_min_response=.08,registration_max_translation_px=6)
 backend=Backend(a.backend);rows=[];allrecords=[];clips={};seeds=[94001] if a.split=='development' else [95001,95002]
 for seed in seeds:
  for case in CASES:
   math=Aggregation();measured=Aggregation();det=Channel(cfg);prev=None;records=[];frames=[];first_truth=None;first_measured=None;first_math=None
   for k,t in enumerate(np.arange(0,2.001,1/30)):
    pos,vel,focus=trajectories(case,seed,float(t));gold=math.update(float(t),pos,pos,vel,(512,512),100*focus,180,cells=True)
    canvas=np.full((512,512),180,np.float32)
    for i,xy in enumerate(pos):splat(canvas,bank[i%len(bank)],xy)
    if case=='defocus':canvas=cv2.GaussianBlur(canvas,(0,0),.3+3*t)
    image=np.clip(canvas,0,255).astype(np.uint8);scan=det.scan(image);qq=scan['candidates'];vv=np.empty((0,4))
    if prev is not None:
     r=backend.run(prev,image,1/30,do_fit=False,actual_coordinates=True);ok=r['valid']&scan['selected'];vv=np.c_[r['x'][ok],r['y'][ok],r['u'][ok]*30,r['v'][ok]*30]
    actual=measured.update(float(t),qq,vv[:,:2],vv[:,2:],(512,512),float(cv2.Laplacian(image,cv2.CV_32F).var()),float(image.mean()),cells=True)
    if gold['status']=='VISIBLE_CONVERGENCE' and first_math is None:first_math=float(t)
    if actual['status']=='VISIBLE_CONVERGENCE' and first_measured is None:first_measured=float(t)
    records.append(dict(time_s=float(t),gold=gold,measured=actual,candidates=qq,vectors=vv));prev=image
    if seed==seeds[0]:frames.append(image)
   expected=case in ['converge','converge_translate','local_converge','form_cluster','onset'];final=records[-1];ispositive=final['measured']['status']=='VISIBLE_CONVERGENCE';mathpos=final['gold']['status']=='VISIBLE_CONVERGENCE'
   rows.append(dict(seed=seed,case=case,expected_convergence=expected,math_status=final['gold']['status'],image_status=final['measured']['status'],math_correct=(mathpos==expected),image_correct=(ispositive==expected),image_false_positive=ispositive and not expected,image_rejected=final['measured']['status']=='NOT_OBSERVABLE',math_response_s=first_math,image_response_s=first_measured))
   allrecords.append(dict(seed=seed,case=case,records=records))
   if frames:clips[case]=np.array(frames)
   print(seed,case,rows[-1],flush=True)
 dump(out/'summary.json',dict(identity=identity(),label='SYNTHETIC_AGGREGATION_TEST',split=a.split,rows=rows,math_correct=sum(r['math_correct'] for r in rows),image_correct=sum(r['image_correct'] for r in rows),image_false_positive=sum(r['image_false_positive'] for r in rows),image_rejected=sum(r['image_rejected'] for r in rows),patch_source_sha=sha(a.video),splat_source_sha=sha(source),physical_acoustics=False))
 import gzip
 with gzip.open(out/'sequences.json.gz','wt') as f:from ycommon import clean;json.dump(clean(allrecords),f)
 np.savez_compressed(out/'clips.npz',**clips)
if __name__=='__main__':main()
