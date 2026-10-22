"""New frozen source blocks; reuse F10U generators and F10V feathered splatting."""
import sys,json,time,argparse
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from wcommon import *
from dataset import splat
from scipy.linalg import expm
from geometry import affine_residual
import cv2
MP=json.loads((ROOT/'configs/model_protocol.json').read_text())
def extract_banks():
    ids=[fid for s in MP['source_splits'].values() for fid in s['frames']];assert len(ids)==len(set(ids))
    sources=old.read_frames(ids,VIDEO);det=old.CandidateDetector();rows=[]
    for split,s in MP['source_splits'].items():
        x,y,w,h=s['roi']
        for fid in s['frames']:
            bgr,pts=sources[fid];roi=cv2.cvtColor(bgr[y:y+h,x:x+w],cv2.COLOR_BGR2GRAY);objects,_=det(roi);patches=[]
            for o in objects:
                px,py=np.round(o['xy']).astype(int)
                if 11<=px<245 and 11<=py<245 and o['classification']=='VISIBLE_SINGLE_CANDIDATE':
                    p=roi[py-10:py+11,px-10:px+11].astype(np.float32);p=np.maximum(np.percentile(p,85)-p,0)*np.outer(np.hanning(21),np.hanning(21))
                    if p.max()>10:patches.append(p)
            assert len(patches)>=2,(split,fid,len(patches))
            path=OUT/'data/banks'/split/f'{fid}.npz';path.parent.mkdir(parents=True,exist_ok=True);np.savez_compressed(path,patches=np.asarray(patches,np.float32),roi=roi)
            rows.append(dict(split=split,source_frame=fid,source_pts=pts,roi=s['roi'],patches=len(patches),file=str(path.relative_to(OUT)),sha256=sha(path)))
    dump(OUT/'data/banks/manifest.json',dict(identity=identity(),video_sha256=sha(VIDEO),banks=rows,source_frame_sets_disjoint=True,scope=MP['scope']))
def one(job):
    split,index,fid,case,rep=job;seed=MP['source_splits'][split]['seed_base']+index;rng=np.random.default_rng(seed);cv2.setNumThreads(1)
    bank=dict(np.load(OUT/'data/banks'/split/f'{fid}.npz'));patches=bank['patches'];roi=bank['roi'];size=256;center=np.array([128.,128.]);times=np.arange(16)*MP['dt_s']
    v,A=old.case_params(case);factor=[.6,1.,1.4][rep%3]
    if case.startswith('translation_'):
        v=np.zeros(2);v[int('_y_' in case)]=(24 if case.endswith('pos') else -24);A*=0
    if case in ['static','illumination_static','non_affine','low_texture','occluded']:v*=0;A*=0
    v*=factor;A*=factor;L=old.generator(v,A,center);H=np.array([expm(L*t) for t in times])
    n=120 if case=='dense' else (4 if case=='sparse' else 40);points=rng.uniform(28,228,(n,2))
    if case=='boundary':points=rng.uniform(-16,272,(n,2))
    if case=='symmetric':angles=np.arange(24)*np.pi/12;points=128+70*np.c_[np.cos(angles),np.sin(angles)];n=24
    patchids=rng.integers(len(patches),size=n)
    if case=='symmetric':patchids[:]=patchids[0]
    fixed=case in ['fixed_background','dense','dropout','boundary','non_affine']
    bg=cv2.GaussianBlur(roi.astype(np.float32),(0,0),1.4) if fixed else np.full((256,256),np.percentile(roi,75),np.float32)
    frames=[];allpoints=[];visible=[]
    for k,t in enumerate(times):
        q=(np.c_[points,np.ones(n)]@H[k].T)[:,:2]
        if case=='non_affine':q+=MP['non_affine_amplitude_px_per_s']*t*np.c_[np.sin(2*np.pi*points[:,1]/256),.6*np.sin(2*np.pi*points[:,0]/256)]
        canvas=bg.copy();vis=[]
        for j,xy in enumerate(q):
            if case=='dropout' and rng.random()<.4:continue
            splat(canvas,patches[patchids[j]],xy);vis.append(j)
        if case=='defocus':canvas=cv2.GaussianBlur(canvas,(0,0),4)
        if case=='illumination_static':canvas+=k*.8
        if case=='low_texture':canvas[:]=190
        if case=='occluded' and k>=13:canvas[:]=190
        canvas+=rng.normal(0,.4,canvas.shape);frames.append(np.uint8(np.clip(canvas,0,255)));allpoints.append(q);visible.append(vis)
    rel=H[15]@np.linalg.inv(H[13]);truth=old.parameters(rel[:2,:2],rel[:2,2],times[15]-times[13],center)
    res=affine_residual(allpoints[13],allpoints[15]-allpoints[13]);label=1 if case=='non_affine' else (2 if case in ['low_texture','occluded'] else 0)
    if label==1:assert res>1.5,(case,res)
    folder=OUT/'data'/split;folder.mkdir(parents=True,exist_ok=True);path=folder/f'{index:04d}.npz'
    np.savez_compressed(path,frames=frames,times=times,H=H,endpoint_H=rel,truth=truth,points=allpoints,label=label)
    return dict(index=index,split=split,source_frame=fid,case=case,replicate=rep,seed=seed,sha256=sha(path),file=str(path.relative_to(OUT)),class_id=label,nonaffine_support_residual=res,visible_ids_generator_only=visible,truth=truth,reference_interval=[times[13],times[15]])
def main():
    p=argparse.ArgumentParser();p.add_argument('--workers',type=int,default=8);a=p.parse_args();start=time.perf_counter()
    assert sha(VIDEO)==MP['material_video_sha256'];assert not (OUT/'data/DATASET_COMPLETE.json').exists(),'REFUSE_OVERWRITE_DATASET'
    extract_banks();manifest={}
    for split,s in MP['source_splits'].items():
        jobs=[]
        for fid in s['frames']:
            for rep in range(s['replicates']):
                for case in MP['cases']:jobs.append((split,len(jobs),fid,case,rep))
        with ProcessPoolExecutor(max_workers=a.workers) as pool:rows=list(pool.map(one,jobs,chunksize=2))
        dump(OUT/'data'/split/'index.json',rows);manifest[split]=dict(n=len(rows),index_sha256=sha(OUT/'data'/split/'index.json'));print(split,len(rows),flush=True)
    dump(OUT/'data/DATASET_COMPLETE.json',dict(identity=identity(),splits=manifest,wall_s=time.perf_counter()-start,scope=MP['scope'],final_test_policy='Generated and hashed before method selection; not loaded for tuning'))
if __name__=='__main__':main()
