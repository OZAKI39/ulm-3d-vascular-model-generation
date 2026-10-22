import numpy as np,torch
from torch.utils.data import Dataset
from torch.nn import functional as F
from vcommon import DATA,SCALES
from heatmaps import centers_to_heatmap
import cv2

class CachedTriples(Dataset):
    def __init__(self,split,corruption=None):
        self.split=split;self.corruption=corruption
        self.a={k:np.load(DATA/split/f'{k}.npy',mmap_mode='r') for k in ['gray','heat','truth','dt','label','H','centers','counts']}
    def __len__(self):return len(self.a['label'])
    def __getitem__(self,i):
        gray=np.array(self.a['gray'][i],dtype=np.float32)/255
        heat=np.array(self.a['heat'][i],dtype=np.float32)/255
        if self.corruption:
            # Train RNG comes from DataLoader worker seed; eval corruption deterministic by sample.
            rng=np.random.default_rng(int(np.random.randint(0,2**31))) if self.corruption=='train' else np.random.default_rng(750000+i)
            mode=rng.choice(['drop','false','jitter','severe','none']) if self.corruption=='train' else self.corruption
            for k in range(3):
                xy=np.array(self.a['centers'][i,k,:int(self.a['counts'][i,k])])
                if mode in ['drop','severe']:xy=xy[rng.random(len(xy))>(.85 if mode=='severe' else .5)]
                if mode in ['false','severe']:xy=np.r_[xy,rng.uniform(0,256,(40 if mode=='severe' else 20,2))]
                if mode in ['jitter','severe']:xy=xy+rng.normal(0,5 if mode=='severe' else 2,xy.shape)
                if mode!='none':heat[k]=cv2.resize(centers_to_heatmap(xy,(256,256),rng.uniform(1.5,6) if mode=='severe' else 3),(128,128))
            if mode=='severe' and rng.random()<.2:heat*=0
        x=np.stack([(gray-.65)/.25,heat],1)
        return dict(images=torch.from_numpy(x),truth=torch.tensor(np.array(self.a['truth'][i])/SCALES,dtype=torch.float32),dt=torch.tensor(self.a['dt'][i].copy()),label=torch.tensor(int(self.a['label'][i])),H=torch.tensor(self.a['H'][i].copy()),index=i)

def loss_fn(out,batch,local=True):
    mask=batch['label']==0;err=out['mean']-batch['truth']
    if mask.any():
        reg=F.smooth_l1_loss(out['mean'][mask],batch['truth'][mask],beta=.1)
        log=out['log_sigma'][mask]
        nll=(.5*(err[mask].detach()/log.exp())**2+log).mean()
    else:reg=err.sum()*0;nll=reg
    ce=F.cross_entropy(out['logits'],batch['label'])
    loss=reg+.05*nll+.2*ce
    local_loss=reg*0
    if local and 'displacement' in out:
        xy=out['xy'];H=batch['H'][:,2];target=(torch.cat([xy,torch.ones_like(xy[:,:1])],1)[None]@H.transpose(1,2))[:,:,:2]-xy
        # Smooth non-affine field is known from the generator and supervised locally.
        nonaffine=batch['label']==1
        nonlinear=9*torch.stack([torch.sin(2*torch.pi*xy[:,1]/256),.6*torch.sin(2*torch.pi*xy[:,0]/256)],-1)
        target=target+nonaffine[:,None,None]*nonlinear
        visible=batch['label']!=2
        if visible.any():
            spatial=(xy[:,0]>26)&(xy[:,0]<230)&(xy[:,1]>26)&(xy[:,1]<230)
            local_loss=F.smooth_l1_loss(out['displacement'][visible][:,spatial],target[visible][:,spatial],beta=.5)
            loss=loss+.2*local_loss
    return loss,dict(regression=reg.detach(),nll=nll.detach(),classification=ce.detach(),local=local_loss.detach())
