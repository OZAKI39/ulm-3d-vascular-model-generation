"""Device-resident progress windows; CPU transfer only for a terminal candidate."""
import torch
from .safety import verify_original_support


class GPUSafetyMonitor:
    def __init__(self,env,dt,device):
        self.env=env;self.dt=dt;self.device=device;self.n=0;self.capacity=0;self.arrays={};self.last_tick=None
    def _resize(self,capacity):
        shapes=dict(calls=(),count=(),last_count=(),no_acceptance=(),unchanged=(),
            history=(64,3),times=(64,),eval_x=(32,3),velocity=(32,3),isolated=(32,),limiter=(64,),previous_nominal=(3,))
        for name,shape in shapes.items():
            dtype=torch.int64 if name in ['calls','count','last_count','no_acceptance','unchanged'] else torch.bool if name in ['isolated','limiter'] else torch.float64
            a=torch.zeros((capacity,*shape),device=self.device,dtype=dtype)
            if name in self.arrays:a[:self.n]=self.arrays[name][:self.n]
            self.arrays[name]=a
        self.capacity=capacity
    def insert(self,x,time):
        if self.n==self.capacity:self._resize(max(256,self.capacity*2))
        i=self.n;self.n+=1
        for a in self.arrays.values():a[i]=0
        self.arrays['history'][i]=x;self.arrays['times'][i]=time
        self.arrays['previous_nominal'][i]=x;self.arrays['count'][i]=1
    def remove(self,mask):
        keep=torch.nonzero(~mask).flatten();n=len(keep)
        for a in self.arrays.values():a[:n]=a[keep]
        self.n=n
    def provider(self):
        a=self.arrays;n=self.n;a['calls'][:n]+=1;calls=a['calls'][:n]
        reasons=torch.zeros(n,device=self.device,dtype=torch.int64);reasons[calls>128000]=4
        ix=torch.nonzero(calls%128==0).flatten()
        if len(ix):
            count=a['count'][ix];a['no_acceptance'][ix]=torch.where(count==a['last_count'][ix],a['no_acceptance'][ix]+1,torch.zeros_like(count))
            a['last_count'][ix]=count;reasons[ix[a['no_acceptance'][ix]>=4]]=3
            ready=ix[count>=64]
            if len(ready):
                positions=a['history'][ready];times=a['times'][ready]
                pad=0. if self.env.wall is None else self.env.wall.roundoff_m
                stagnant=(torch.linalg.vector_norm(positions.max(1).values-positions.min(1).values,dim=1)<=16*pad)&a['isolated'][ready].all(1)
                tiny=(times.max(1).values-times.min(1).values)<self.dt*.01
                responsible=a['limiter'][ready].any(1)
                reasons[ready[tiny & responsible]]=2
                reasons[ready[stagnant & (~tiny | responsible)]]=1
        return reasons
    def accepted(self,oldx,endx,velocity,time,isolated,limiter=None):
        a=self.arrays;n=self.n;ix=torch.arange(n,device=self.device);count=a['count'][:n]
        a['limiter'][ix,count%64]=False if limiter is None else limiter
        a['history'][ix,count%64]=endx;a['times'][ix,count%64]=time
        a['eval_x'][ix,(count-1)%32]=oldx;a['velocity'][ix,(count-1)%32]=velocity
        a['isolated'][ix,(count-1)%32]=isolated;a['count'][:n]+=1
    def nominal(self,time,x):
        import numpy as np
        ratio=time/self.dt;tick=round(ratio)
        reasons=torch.zeros(self.n,device=self.device,dtype=torch.int64)
        if abs(ratio-tick)>8*abs(np.spacing(ratio)) or tick==self.last_tick:return reasons
        self.last_tick=tick;a=self.arrays;n=self.n
        same=(x==a['previous_nominal'][:n]).all(1)
        a['unchanged'][:n]=torch.where(same,a['unchanged'][:n]+1,torch.zeros_like(a['unchanged'][:n]))
        a['previous_nominal'][:n]=x;reasons[(a['unchanged'][:n]>=32)&a['isolated'][:n].all(1)&(a['count'][:n]>=33)]=5
        return reasons
    def support(self,index,radius):
        a=self.arrays
        if int(a['count'][index])<33 or not bool(a['isolated'][index].all()):return False,[]
        order=(torch.arange(32,device=self.device)+int(a['count'][index])-1)%32
        return verify_original_support(self.env,a['eval_x'][index,order].cpu().numpy(),a['velocity'][index,order].cpu().numpy(),radius)
