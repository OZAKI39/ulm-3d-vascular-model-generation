"""Real Chronos patches, analytic motion and independent, frozen visual test seeds."""
import cv2,numpy as np,av
from scipy.linalg import expm
from common import ROOT,VIDEO,P,dump,sha
from distribution import CandidateDetector

def read_frames(indices):
    indices=set(indices);out={}
    with av.open(str(VIDEO)) as c:
        for i,f in enumerate(c.decode(video=0)):
            if i in indices:out[i]=(f.to_ndarray(format='bgr24'),float(f.time))
            if i>=max(indices):break
    return out

def generator(v,A,center):
    L=np.zeros((3,3));L[:2,:2]=A;L[:2,2]=np.array(v)-np.array(A)@center
    return L

def case_params(name):
    v=np.array([18.,-12.]);A=np.zeros((2,2))
    if name in ['expansion','contraction']:v*=0;A=np.eye(2)*(.3 if name=='expansion' else -.3)
    elif name.startswith('rotation'):v*=0;r=.45 if name=='rotation_cw' else -.45;A=np.array([[0,-r],[r,0]])
    elif name=='stretch':v*=0;A=np.diag([.3,-.3])
    elif name=='shear':v*=0;A=np.array([[0,.5],[0,0]])
    elif name in ['mixed','fixed_background','dense','defocus','dropout','brightness_noise','boundary']:A=np.array([[.12,-.25],[.35,-.08]])
    elif name=='symmetric':v*=0;A=np.array([[0,-.45],[.45,0]])
    elif name=='local_opposite':v*=0
    elif name=='step':v*=0
    return v,A

def build_sequence(name, seed, source,length=None):
    rng=np.random.default_rng(seed);size=P['synthetic_size'];T=length or P['synthetic_frames'];dt=P['synthetic_dt_s'];center=np.array([size/2]*2)
    frame,pts=source;g=cv2.cvtColor(frame,cv2.COLOR_BGR2GRAY)
    # Disjoint background and patch sampling in development/test source frames.
    roi=g[300:300+size,650:650+size];objs,_=CandidateDetector()(roi)
    patches=[]
    for o in objs:
        x,y=np.round(o['xy']).astype(int)
        if 9<=x<size-9 and 9<=y<size-9 and o['classification']=='VISIBLE_SINGLE_CANDIDATE':
            p=roi[y-8:y+9,x-8:x+9].astype(float);p=np.maximum(np.percentile(p,85)-p,0)
            # Feather patch boundary; intensity/texture remains from the real image.
            win=np.outer(np.hanning(17),np.hanning(17));patches.append(p*win)
    if len(patches)<2:raise RuntimeError('insufficient real source patches')
    n=40
    if name=='dense':n=100
    if name=='sparse':n=4
    points=rng.uniform(32,size-32,(n,2))
    if name=='boundary':points=rng.uniform(-18,size+18,(n,2))
    if name=='symmetric':
        angles=np.arange(24)*2*np.pi/24;points=center+70*np.c_[np.cos(angles),np.sin(angles)];n=len(points)
    v,A=case_params(name);L=generator(v,A,center)
    frames=[];H=[];truthpoints=[];vis=[]
    backgrounds=[]
    fixed=name in ['fixed_background','dense','sparse','defocus','dropout','boundary','brightness_noise','local_opposite']
    bg=cv2.GaussianBlur(roi.astype(np.float32),(0,0),1.4) if fixed else np.full((size,size),190,np.float32)
    for k in range(T):
        tm=k*dt
        if name=='step':
            Hk=expm(generator([30.,0.],np.zeros((2,2)),center)*max(0,(k-7)*dt))
        else:Hk=expm(L*tm)
        q=(np.c_[points,np.ones(n)]@Hk.T)[:,:2]
        if name=='local_opposite':q[:,0]+=np.where(points[:,1]<size/2,26.,-26.)*tm
        canvas=bg.copy(); visible=[]
        for j,(x,y) in enumerate(q):
            if name=='dropout' and rng.random()<.45:continue
            if not (-8<x<size+8 and -8<y<size+8):continue
            patch=patches[0 if name=='symmetric' else j%len(patches)].astype(np.float32)
            patch=cv2.warpAffine(patch,np.array([[1,0,x-8],[0,1,y-8]],np.float32),(size,size),flags=cv2.INTER_LINEAR)
            canvas-=patch;visible.append(j)
        if name=='defocus':canvas=cv2.GaussianBlur(canvas,(0,0),5)
        if name=='brightness_noise':canvas+=12*np.sin(k*.3);canvas+=rng.normal(0,3,canvas.shape)
        else:canvas+=rng.normal(0,.5,canvas.shape)
        frames.append(np.uint8(np.clip(canvas,0,255)));H.append(Hk);truthpoints.append(q);vis.append(visible)
    return np.array(frames),{'case':name,'seed':seed,'dt':dt,'source_pts':pts,'patch_count':len(patches),'H':np.array(H),'points':np.array(truthpoints),'visible_ids_generator_only':vis,'global_affine_truth_available':name!='local_opposite','instantaneous_generator':L,'n':n}

def generate(split):
    source=read_frames([P['background_frames'][split]])[P['background_frames'][split]]
    folder=ROOT/'results/synthetic'/split;folder.mkdir(parents=True,exist_ok=True)
    manifest=[]
    for seed in P['synthetic_seeds'][split]:
        for name in P['cases']:
            frames,meta=build_sequence(name,seed,source);stem=f'{seed}_{name}'
            np.savez_compressed(folder/(stem+'.npz'),frames=frames,H=meta.pop('H'),points=meta.pop('points'),instantaneous_generator=meta.pop('instantaneous_generator'))
            dump(folder/(stem+'.json'),meta)
            manifest.append({'file':stem+'.npz','sha256':sha(folder/(stem+'.npz')),'case':name,'seed':seed})
    dump(folder/'manifest.json',manifest)
    return folder
