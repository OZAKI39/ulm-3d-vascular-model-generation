"""Official TorchVision backbone plus project-specific temporal/geometry heads."""
import torch
from torch import nn
from torch.nn import functional as F
from torchvision.models import mobilenet_v3_small
from vcommon import SCALES

class FrameEncoder(nn.Module):
    def __init__(self,channels,local=False):
        super().__init__();base=mobilenet_v3_small(weights=None,width_mult=.5)
        # Keep official blocks, change only input stem channels for the required fusion.
        old=base.features[0][0]
        base.features[0][0]=nn.Conv2d(channels,old.out_channels,3,2,1,bias=False)
        self.features=base.features[:2] if local else base.features
        self.out_channels=8 if local else 288
    def forward(self,x):return self.features(x)

class Direct6(nn.Module):
    def __init__(self,mode='fused',frames=3):
        super().__init__();self.mode=mode;self.frames=frames;self.encoder=FrameEncoder(2 if mode=='fused' else 1)
        self.reduce=nn.Sequential(nn.Conv2d(288,24,1),nn.Hardswish(),nn.AdaptiveAvgPool2d((4,4)))
        # Spatial positions preserved; GAP alone would discard affine orientation.
        self.head=nn.Sequential(nn.Linear(24*16*(2*frames-1)+2,192),nn.Hardswish(),nn.Dropout(.1),nn.Linear(192,15))
    def forward(self,images,dt):
        x=images[:,:,0:1] if self.mode=='gray' else (images[:,:,1:2] if self.mode=='heat' else images)
        if self.frames==2:x=x[:,[0,2]]
        b,t,c,h,w=x.shape;f=self.reduce(self.encoder(x.reshape(b*t,c,h,w))).reshape(b,t,-1)
        temporal=torch.cat([f.flatten(1),(f[:,1:]-f[:,:-1]).flatten(1),dt/.0333],1)
        z=self.head(temporal)
        return dict(mean=z[:,:6],log_sigma=z[:,6:12].clamp(-4,2),logits=z[:,12:])

def torch_matrix_log(H,terms=24):
    # Near-identity transforms in normalized spatial coordinates, tested vs SciPy.
    E=H-torch.eye(3,device=H.device,dtype=H.dtype);power=E;L=E
    for k in range(2,terms+1):
        power=power@E;L=L+((-1)**(k+1)/k)*power
    return L

def six_from_normalized_H(H,dt):
    L=torch_matrix_log(H)/dt[:,None,None];A=L[:,:2,:2];v=L[:,:2,2]*128
    return torch.stack([v[:,0],v[:,1],A[:,0,0]+A[:,1,1],A[:,1,0]-A[:,0,1],(A[:,0,0]-A[:,1,1])/2,(A[:,0,1]+A[:,1,0])/2],1)

def weighted_affine(xy,displacement,weight,iterations=3):
    """Task-specific differentiable weighted/Huber affine adapter to torch.linalg."""
    X=torch.cat([torch.ones_like(xy[:,:1]),(xy-128)/128],1)
    X=X[None].expand(displacement.shape[0],-1,-1);w=weight.clamp_min(.001)
    for k in range(iterations):
        normal=X.transpose(1,2)@(w[:,:,None]*X)
        rhs=X.transpose(1,2)@(w[:,:,None]*displacement)
        B=torch.linalg.solve(normal+torch.eye(3,device=X.device)[None]*1e-5,rhs)
        residual=(X@B-displacement).norm(dim=-1)
        if k<iterations-1:w=weight/(1+(residual/.8)**2).sqrt()
    H=torch.eye(3,device=X.device)[None].repeat(len(B),1,1)
    H[:,:2,:2]+=B[:,1:].transpose(1,2)/128;H[:,:2,2]=B[:,0]/128
    return H,residual

class GeoMotion(nn.Module):
    def __init__(self):
        super().__init__();self.encoder=FrameEncoder(2,local=True)
        self.project=nn.Conv2d(8,24,1);self.log_temp=nn.Parameter(torch.tensor(2.5))
        self.weight_head=nn.Conv2d(24,1,1)
        self.head=nn.Sequential(nn.Linear(24*3+4,96),nn.Hardswish(),nn.Linear(96,9))
        self.register_buffer('scales',torch.tensor(SCALES,dtype=torch.float32));self.radius=3
    def match(self,a,b):
        n,c,h,w=a.shape;r=self.radius
        a=F.normalize(a,dim=1);b=F.normalize(b,dim=1)
        neighbors=F.unfold(b,2*r+1,padding=r).reshape(n,c,(2*r+1)**2,h*w)
        corr=(a.flatten(2)[:,:,None]*neighbors).sum(1)
        yy,xx=torch.meshgrid(torch.arange(-r,r+1,device=a.device),torch.arange(-r,r+1,device=a.device),indexing='ij')
        offset=torch.stack([xx.flatten(),yy.flatten()],-1).float()
        prob=F.softmax(corr*self.log_temp.exp().clamp(1,100),dim=1)
        disp=torch.einsum('bkn,kd->bnd',prob,offset)*8 # feature stride4 x resize inverse2
        entropy=-(prob*prob.clamp_min(1e-9).log()).sum(1)/torch.log(torch.tensor(float((2*r+1)**2),device=a.device))
        return disp,1-entropy
    def forward(self,images,dt):
        b,t,c,h,w=images.shape
        f=self.project(self.encoder(images.reshape(b*t,c,h,w))).reshape(b,t,24,32,32)
        d02,c02=self.match(f[:,0],f[:,2]);d01,c01=self.match(f[:,0],f[:,1]);d12,c12=self.match(f[:,1],f[:,2])
        yy,xx=torch.meshgrid(torch.arange(32,device=f.device),torch.arange(32,device=f.device),indexing='ij')
        # Conv stride4 receptive-center at input0; inverse resize mapping maps0 ->0.5rawpx.
        xy=torch.stack([xx.flatten()*8+.5,yy.flatten()*8+.5],-1).float()
        grid=torch.stack([xx,yy],-1).float()[None]+d01.reshape(b,32,32,2)/8
        grid=2*(grid+.5)/32-1
        chain=d01+F.grid_sample(d12.reshape(b,32,32,2).permute(0,3,1,2),grid,align_corners=False).flatten(2).transpose(1,2)
        displacement=(d02+chain)/2
        consistency=(d02-chain).norm(dim=-1)
        weight=torch.sigmoid(self.weight_head(f[:,0])).flatten(1)*(.05+c02)/(1+consistency)
        # Exclude padded match margins; detector map gives weak evidence, never GT centers.
        valid=((xy[:,0]>26)&(xy[:,0]<230)&(xy[:,1]>26)&(xy[:,1]<230)).float()
        contrast=F.adaptive_max_pool2d(images[:,0,1:2],(32,32)).flatten(1)
        weight=weight*valid[None]*(.1+contrast)
        # Float32 geometry deliberately outside AMP to keep small affine terms accurate.
        with torch.autocast(device_type=images.device.type,enabled=False):
            H,res=weighted_affine(xy,displacement.float(),weight.float())
            mean=six_from_normalized_H(H,dt.sum(1).float())/self.scales
        pooled=f.mean((-2,-1)).flatten(1)
        extra=torch.cat([dt/.0333,consistency.mean(1,keepdim=True),res.mean(1,keepdim=True)],1)
        z=self.head(torch.cat([pooled,extra],1))
        return dict(mean=mean,log_sigma=z[:,:6].clamp(-4,2),logits=z[:,6:],H_normalized=H,xy=xy,displacement=displacement,weight=weight,residual=res)

def build_model(name):
    if name.startswith('geo'):return GeoMotion()
    mode='gray' if name=='direct_gray' else ('heat' if name=='direct_heat' else 'fused')
    return Direct6(mode,2 if '2frame' in name else 3)
