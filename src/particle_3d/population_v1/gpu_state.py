"""Device-resident float64 SoA; only births/removals/output cross devices."""
import torch
import numpy as np


class GPUActivePool:
    integer_names=['particle_id','source_event_id','realization_id','tetra_id','status']
    scalar_names=['birth_time','age','radius','diameter','wall_gap']
    vector_names={'x':3,'v':3,'omega':3,'q':4}

    def __init__(self,capacity=256,device='cuda'):
        self.device=torch.device(device)
        if self.device.index is None:self.device=torch.device('cuda',torch.cuda.current_device())
        self.capacity=max(int(capacity),1);self.n=0;self.used_ids=set()
        for name in self.integer_names:setattr(self,name,torch.zeros(self.capacity,device=self.device,dtype=torch.int64))
        for name in self.scalar_names:setattr(self,name,torch.zeros(self.capacity,device=self.device,dtype=torch.float64))
        for name,width in self.vector_names.items():setattr(self,name,torch.zeros((self.capacity,width),device=self.device,dtype=torch.float64))

    def arrays(self):return {k:getattr(self,k) for k in self.integer_names+self.scalar_names+list(self.vector_names)}

    def insert(self,particle,event,time):
        key=(int(event.get('realization_id',0)),int(particle.particle_id))
        if key in self.used_ids:raise ValueError('STABLE_ID_NEVER_REUSED')
        if event['birth_time_s']!=time:raise ValueError('BIRTH_MUST_OCCUR_AT_EXACT_PHYSICAL_TIME')
        if self.n==self.capacity:
            self.capacity*=2
            for name,old in self.arrays().items():
                new=torch.zeros((self.capacity,*old.shape[1:]),device=self.device,dtype=old.dtype);new[:self.n]=old[:self.n];setattr(self,name,new)
        k=self.n;self.n+=1;self.used_ids.add(key)
        values=dict(particle_id=key[1],realization_id=key[0],source_event_id=event['source_event_id'],tetra_id=-1,status=0,
            birth_time=time,age=0.,radius=particle.radius,diameter=2*particle.radius,wall_gap=event.get('clearance_m',0.),
            x=particle.position,v=particle.velocity,omega=particle.omega,q=particle.q)
        for name,value in values.items():
            a=getattr(self,name)
            if isinstance(value,np.ndarray):value=value.copy()
            a[k]=torch.as_tensor(value,device=self.device,dtype=a.dtype)

    def remove_mask(self,remove):
        if remove.dtype!=torch.bool or len(remove)!=self.n:raise ValueError('Device Boolean removal mask required')
        keep=torch.nonzero(~remove).flatten();count=len(keep)
        for a in self.arrays().values():a[:count]=a[keep]
        self.n=count

    def device_snapshot(self,time):
        return dict(physical_time_s=float(time),**{k:a[:self.n].clone() for k,a in self.arrays().items()})
