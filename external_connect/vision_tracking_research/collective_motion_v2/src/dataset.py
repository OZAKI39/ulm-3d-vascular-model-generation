"""F10U real-patch generator extended to independently sampled causal triples.
Patch feathering, motion definitions and homogeneous generator reuse F10U.
Bounded splatting replaces full-canvas-per-particle warps for cache throughput.
"""
import os,time,json
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import cv2,numpy as np
from scipy.linalg import expm
from vcommon import P,ROOT,DATA,dump,sha,identity
from legacy import read_frames,case_params,generator
from geometry import transform,endpoint_six,apply_transform,affine_residual
from heatmaps import HeatmapDetector

def make_banks(video=None):
    video=video or os.environ.get('F10V_VIDEO',P['source_video'])
    if sha(video)!=P['source_video_sha256']:raise ValueError('SOURCE_VIDEO_SHA_MISMATCH')
    ids=[i for s in P['source_splits'].values() for i in s['frames']]
    sources=read_frames(ids,video); manifest=[];det=HeatmapDetector()
    for split,s in P['source_splits'].items():
        x,y,w,h=s['roi_xywh'];folder=DATA/'banks'/split;folder.mkdir(parents=True,exist_ok=True)
        for fid in s['frames']:
            image,pts=sources[fid];roi=cv2.cvtColor(image[y:y+h,x:x+w],cv2.COLOR_BGR2GRAY)
            objects,gray,heat=det(roi);patches=[];centers=[]
            for o in objects:
                px,py=np.round(o['xy']).astype(int)
                if 11<=px<w-11 and 11<=py<h-11 and o['classification']=='VISIBLE_SINGLE_CANDIDATE':
                    a=roi[py-10:py+11,px-10:px+11].astype(np.float32)
                    a=np.maximum(np.percentile(a,85)-a,0)*np.outer(np.hanning(21),np.hanning(21))
                    if a.max()>10:patches.append(a);centers.append([px,py])
            if len(patches)<2:raise ValueError(f'INSUFFICIENT_PATCH_BANK {split} {fid}: {len(patches)}')
            f=folder/f'{fid}.npz';np.savez_compressed(f,patches=np.array(patches,np.float32),roi=roi,centers=centers,heat=heat)
            manifest.append(dict(split=split,frame_id=fid,pts_s=pts,roi_xywh=s['roi_xywh'],patches=len(patches),sha256=sha(f),file=str(f.relative_to(DATA))))
    # Split equality is not claimed to be cross-recording independence.
    assert len(ids)==len(set(ids))
    for a,sa in P['source_splits'].items():
        for b,sb in P['source_splits'].items():
            if a>=b:continue
            ax,ay,aw,ah=sa['roi_xywh'];bx,by,bw,bh=sb['roi_xywh']
            assert ax+aw<=bx or bx+bw<=ax or ay+ah<=by or by+bh<=ay
            assert sa['seed_base']+sa['samples']<sb['seed_base'] or sb['seed_base']+sb['samples']<sa['seed_base']
    dump(DATA/'bank_manifest.json',dict(video_sha256=sha(video),material=manifest,split_leakage_check=True,scope=P['generalization_scope'],identity=identity()))
    return manifest

def splat(canvas,patch,xy):
    x,y=xy;h,w=canvas.shape;r=patch.shape[0]//2
    x0=int(np.floor(x))-r-1;y0=int(np.floor(y))-r-1
    x1=x0+patch.shape[1]+2;y1=y0+patch.shape[0]+2
    left,right=max(0,x0),min(w,x1);top,bottom=max(0,y0),min(h,y1)
    if left>=right or top>=bottom:return
    a=cv2.warpAffine(patch,np.array([[1,0,x-r-left],[0,1,y-r-top]],np.float32),(right-left,bottom-top),flags=cv2.INTER_LINEAR,borderMode=cv2.BORDER_CONSTANT)
    canvas[top:bottom,left:right]-=a

def draw_parameters(case,rng):
    v,A=case_params(case)
    if case in ['translation_x','translation_y']:
        v=np.zeros(2);v[int(case=='translation_y')]=rng.choice([-1,1])*rng.uniform(4,48);A*=0
    else:
        v=v*rng.uniform(.2,2.0)*rng.choice([-1,1],size=2)
        A=A*rng.uniform(.25,1.6)
    if case in ['static','illumination_static']:v*=0;A*=0
    return generator(v,A,np.array([128.,128.]))

def nonlinear_displacement(points,t,total):
    # A genuinely non-affine smooth field; opposite linear shear is not used.
    p=np.asarray(points);a=9*t/total
    return a*np.c_[np.sin(2*np.pi*p[:,1]/256),.6*np.sin(2*np.pi*p[:,0]/256)]

def sample(split,index,banks,detector):
    s=P['source_splits'][split];case=P['cases'][index%len(P['cases'])];seed=s['seed_base']+index
    rng=np.random.default_rng(seed);fid=s['frames'][(index//len(P['cases']))%len(s['frames'])];bank=banks[fid]
    roi=bank['roi'];patches=bank['patches'];dt=rng.choice(P['dt_choices_s'],2);times=np.r_[0,np.cumsum(dt)]
    L=draw_parameters(case,rng);L0=L.copy();L1=L.copy()
    if case=='start':L0*=0
    if case=='stop':L1*=0
    if case=='acceleration':L0*=.25;L1*=1.5
    H=np.array([np.eye(3),expm(L0*dt[0]),expm(L1*dt[1])@expm(L0*dt[0])])
    truth=endpoint_six(H[2],sum(dt),roi.shape)
    n=int(rng.integers(24,65))
    if case=='dense':n=120
    if case=='sparse':n=8
    points=rng.uniform(20,236,(n,2))
    if case=='boundary':points=rng.uniform(-12,268,(n,2))
    if case=='symmetric':
        angles=np.arange(24)*2*np.pi/24;points=128+70*np.c_[np.cos(angles),np.sin(angles)];n=24
    if case=='collinear':points=np.c_[np.linspace(30,226,24),np.full(24,128.)];n=24
    fixed=case in ['fixed_background','dense','sparse','dropout','boundary','local_blur','background_dominant','brightness_noise','non_affine']
    bg=cv2.GaussianBlur(roi.astype(np.float32),(0,0),1.4) if fixed else np.full((256,256),float(np.percentile(roi,75)),np.float32)
    contrast=rng.uniform(.8,1.5)
    if case=='background_dominant':contrast=.3
    if case=='low_texture':bg=np.full((256,256),190,np.float32);contrast=.001
    patchids=rng.integers(len(patches),size=n)
    if case=='symmetric':patchids[:]=patchids[0]
    frames=[];heats=[];allcenters=[];counts=[];allpoints=[]
    for k,t in enumerate(times):
        q=apply_transform(H[k],points)
        if case=='non_affine':q+=nonlinear_displacement(points,t,times[-1])
        canvas=bg.copy()
        for j,p in enumerate(q):
            if case=='dropout' and rng.random()<.4:continue
            splat(canvas,patches[patchids[j]]*contrast,p)
        if case=='defocus':canvas=cv2.GaussianBlur(canvas,(0,0),rng.uniform(2.5,4.5))
        if case=='local_blur':canvas[:,128:]=cv2.GaussianBlur(canvas,(0,0),3)[:,128:]
        if case=='brightness_noise':canvas=canvas*rng.uniform(.85,1.15)+rng.normal(0,3,canvas.shape)
        if case=='illumination_static':canvas+=k*10
        if case=='occluded' and k>0:canvas[:]=190
        canvas+=rng.normal(0,.4 if case!='low_texture' else .05,canvas.shape)
        img=np.uint8(np.clip(canvas,0,255));objects,_,heat=detector(img)
        center=np.zeros((192,2),np.float32);xy=np.asarray([o['xy'] for o in objects],float).reshape(-1,2)[:192];center[:len(xy)]=xy
        frames.append(img);heats.append(heat);allcenters.append(center);counts.append(len(xy));allpoints.append(q)
    grid=np.stack(np.meshgrid(np.linspace(16,240,16),np.linspace(16,240,16)),-1).reshape(-1,2)
    displacement=apply_transform(H[2],grid)-grid
    if case=='non_affine':displacement+=nonlinear_displacement(grid,times[-1],times[-1])
    residual=affine_residual(grid,displacement)
    label=1 if residual>P['nonaffine_min_grid_rms_px'] else 0
    if case in ['low_texture','occluded','collinear']:label=2
    # labels include exact known transforms even when not observable; evaluation masks class0.
    frames=np.array(frames);heats=np.array(heats)
    gray=np.array([cv2.resize(a,(128,128),interpolation=cv2.INTER_AREA) for a in frames])
    heat=np.array([np.uint8(np.clip(cv2.resize(a,(128,128))*255,0,255)) for a in heats])
    meta=dict(index=index,split=split,case=case,seed=seed,source_frame=fid,dt_s=dt.tolist(),times_s=times.tolist(),class_id=label,nonaffine_grid_rms_px=residual,
              label_basis='known nonlinear grid residual' if label==1 else ('explicit absent/degenerate material support' if label==2 else 'known global affine; difficulty retained'),truth=truth.tolist(),source_bank_sha256=sha(DATA/'banks'/split/f'{fid}.npz'))
    return dict(gray=gray,heat=heat,raw=frames,H=H.astype(np.float32),truth=truth.astype(np.float32),dt=dt.astype(np.float32),label=np.int64(label),centers=np.array(allcenters),counts=np.array(counts,np.int16)),meta

_BANKS=None;_DET=None;_SPLIT=None

def init_worker(split):
    global _BANKS,_DET,_SPLIT
    cv2.setNumThreads(1);_SPLIT=split;_DET=HeatmapDetector()
    _BANKS={fid:dict(np.load(DATA/'banks'/split/f'{fid}.npz')) for fid in P['source_splits'][split]['frames']}

def worker(i):return sample(_SPLIT,i,_BANKS,_DET)

def generate(split,workers=8,limit=None):
    folder=DATA/split;folder.mkdir(parents=True,exist_ok=True)
    n=limit or P['source_splits'][split]['samples'];start=time.perf_counter()
    shapes={'gray':(3,128,128),'heat':(3,128,128),'H':(3,3,3),'truth':(6,),'dt':(2,),'label':(),'centers':(3,192,2),'counts':(3,)}
    if split!='train':shapes['raw']=(3,256,256)
    dtypes={'gray':'uint8','heat':'uint8','raw':'uint8','label':'int64','counts':'int16'}
    arrays={k:np.lib.format.open_memmap(folder/f'{k}.npy',mode='w+',dtype=dtypes.get(k,'float32'),shape=(n,)+shape) for k,shape in shapes.items()}
    with (folder/'index.jsonl').open('w') as indexfile, ProcessPoolExecutor(max_workers=workers,initializer=init_worker,initargs=(split,)) as pool:
        for i,(data,meta) in enumerate(pool.map(worker,range(n),chunksize=16)):
            for k in arrays:arrays[k][i]=data[k]
            indexfile.write(json.dumps(meta)+'\n')
            if i%512==0:print(split,i,n,f'{time.perf_counter()-start:.1f}s',flush=True)
    for a in arrays.values():a.flush()
    record=dict(split=split,n=n,wall_s=time.perf_counter()-start,identity=identity(),files={f.name:sha(f) for f in folder.glob('*') if f.is_file() and f.name!='manifest.json'})
    dump(folder/'manifest.json',record);return record
