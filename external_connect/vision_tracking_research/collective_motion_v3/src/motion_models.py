"""Thin adapters to official pretrained implementations and POT. No new backbone."""
import sys,time,importlib,types,contextlib,io,warnings
import numpy as np
import cv2
from wcommon import THIRD,MODELS,old,preprocess

def point_array(objects):return np.asarray([o['xy'] for o in objects],float).reshape(-1,2)
def sample_field(field,xy):
    xy=np.asarray(xy,np.float32).reshape(-1,2)
    if not len(xy):return np.empty((0,2))
    out=cv2.remap(np.asarray(field,np.float32),xy[:,0,None],xy[:,1,None],cv2.INTER_LINEAR,borderMode=cv2.BORDER_CONSTANT,borderValue=(float('nan'),float('nan')))[:,0]
    out[(xy[:,0]<0)|(xy[:,0]>field.shape[1]-1)|(xy[:,1]<0)|(xy[:,1]>field.shape[0]-1)]=np.nan
    return out
def input_image(gray,mode):
    if mode=='raw':return np.asarray(gray,np.float32)/255
    if mode=='inverted':return 1-np.asarray(gray,np.float32)/255
    if mode=='local':return np.clip(preprocess(gray,'local')/40,0,1)
    raise ValueError(mode)

class DenseOfficial:
    def __init__(self,name):
        import torch
        self.torch=torch;self.name=name;torch.set_num_threads(2)
        self.args=types.SimpleNamespace(channel_threshold=True,test=1,amp=False,iters=12,upsample='convex')
        if name=='LRAFT_ST':
            sys.path.insert(0,str(THIRD/'LRAFT-ST'));module=importlib.import_module('LRAFT-ST.RAFT_threshold')
            self.model=module.LRAFT_ST(self.args);self.path=MODELS/'lraft_student_dataset1.tar'
        elif name=='UnLiteFlowNet':
            sys.path.insert(0,str(THIRD/'UnLiteFlowNet-PIV'));module=importlib.import_module('src.model.models')
            self.model=module.Network();self.estimate=module.estimate;self.path=THIRD/'UnLiteFlowNet-PIV/models/UnsupervisedLiteFlowNet_pretrained.pt'
        else:raise ValueError(name)
        checkpoint=torch.load(self.path,map_location='cpu',weights_only=False)
        self.model.load_state_dict(checkpoint['model_state_dict'],strict=True);self.model=self.model.cuda().eval()
        self.parameters=sum(p.numel() for p in self.model.parameters())
    def flow(self,a,b,mode='inverted'):
        torch=self.torch;start=time.perf_counter();h,w=a.shape
        # Keep original pixels. Padding is cropped; no displacement scale ambiguity.
        hh=int(np.ceil(h/32)*32);ww=int(np.ceil(w/32)*32)
        pair=np.stack([input_image(a,mode),input_image(b,mode)])
        pair=np.pad(pair,((0,0),(0,hh-h),(0,ww-w)),mode='edge')
        tensor=torch.from_numpy(pair[None].copy()).float().cuda();torch.cuda.synchronize();prep=(time.perf_counter()-start)*1000
        tick=time.perf_counter()
        with torch.inference_mode():
            if self.name=='LRAFT_ST':f=self.model(tensor,args=self.args)[-1]
            else:f=self.estimate(tensor[:,0:1],tensor[:,1:2],self.model,train=False)
        torch.cuda.synchronize();model_ms=(time.perf_counter()-tick)*1000
        field=f[0,:,:h,:w].detach().float().cpu().permute(1,2,0).numpy()
        return field,dict(preprocess_ms=prep,model_ms=model_ms,parameters=self.parameters)
    def pair(self,images,objects,times,mode='inverted',fb_limit=1.5):
        start=time.perf_counter();a,b=images[-3],images[-1];dt=times[-1]-times[-3]
        field,timing=self.flow(a,b,mode);back,tb=self.flow(b,a,mode);xy=point_array(objects[-3]);d=sample_field(field,xy);rev=sample_field(back,xy+d);fb=np.linalg.norm(d+rev,axis=1)
        h,w=a.shape;good=np.isfinite(d).all(1)&np.isfinite(fb)&(fb<=fb_limit)&(xy[:,0]+d[:,0]>=0)&(xy[:,0]+d[:,0]<w)&(xy[:,1]+d[:,1]>=0)&(xy[:,1]+d[:,1]<h)
        tick=time.perf_counter();e=old.fit_affine(xy[good],d[good],dt,a.shape,self.name,bootstrap=0);fit_ms=(time.perf_counter()-tick)*1000
        return dict(estimate=e,vectors=np.c_[xy,d,good.astype(float)],field=field,fb_errors=fb,quality_score=float(1/(1+np.median(fb[good])+e.get('residual',0))) if good.any() and e.get('residual') is not None else 0.,timing=dict(**timing,backward_model_ms=tb['model_ms'],fit_ms=fit_ms,total_ms=(time.perf_counter()-start)*1000),causal=True)

class CoTrackerOfficial:
    def __init__(self):
        import torch
        self.torch=torch;torch.set_num_threads(2);sys.path.insert(0,str(THIRD/'co-tracker'))
        from cotracker.predictor import CoTrackerOnlinePredictor
        self.path=MODELS/'cotracker3_scaled_online.pth';self.model=CoTrackerOnlinePredictor(checkpoint=str(self.path),window_len=16,v2=False).cuda().eval()
        self.parameters=sum(p.numel() for p in self.model.parameters())
    def pair(self,images,objects,times,**kwargs):
        torch=self.torch;start=time.perf_counter();xy=point_array(objects[0]);h,w=images[0].shape
        if len(xy)<8:return dict(estimate=old.empty('INSUFFICIENT_QUERY_POINTS','CoTracker3_online'),vectors=[],quality_score=0.,timing=dict(total_ms=(time.perf_counter()-start)*1000),causal=True)
        # Bounded queries are chosen spatially by fixed detector order, never by future tracks.
        xy=xy[:192];video=torch.from_numpy(np.stack(images).copy())[None,:,None].repeat(1,1,3,1,1).float().cuda()
        queries=torch.from_numpy(np.c_[np.zeros(len(xy)),xy].astype(np.float32))[None].cuda();torch.cuda.synchronize();prep=(time.perf_counter()-start)*1000
        tick=time.perf_counter()
        with torch.inference_mode():
            self.model(video[:,:8],is_first_step=True,queries=queries,add_support_grid=True)
            tracks,visible=self.model(video,add_support_grid=True)
        torch.cuda.synchronize();model_ms=(time.perf_counter()-tick)*1000
        tr=tracks[0].float().cpu().numpy();vis=visible[0].cpu().numpy();x=tr[-3];d=tr[-1]-tr[-3]
        good=vis[-3]&vis[-1]&np.isfinite(d).all(1)&(x[:,0]>=0)&(x[:,0]<w)&(x[:,1]>=0)&(x[:,1]<h)
        tick=time.perf_counter();e=old.fit_affine(x[good],d[good],times[-1]-times[-3],(h,w),'CoTracker3_online',bootstrap=0);fit_ms=(time.perf_counter()-tick)*1000
        return dict(estimate=e,vectors=np.c_[x,d,good.astype(float)],tracks=tr,visibility=vis,query_points=xy,quality_score=float(good.mean()/(1+e['residual'])) if e.get('residual') is not None else 0.,timing=dict(preprocess_ms=prep,model_ms=model_ms,fit_ms=fit_ms,total_ms=(time.perf_counter()-start)*1000,parameters=self.parameters),causal=True,context_frames=len(images),native_step=self.model.step,first_output_span_s=float(times[-1]-times[0]),reference_interval=[float(times[-3]),float(times[-1])])

def optimal_transport(images,objects,times,method='unbalanced',distance_scale=8.,mass=.85):
    import ot
    start=time.perf_counter();oa,ob=objects[-3],objects[-1];x,y=point_array(oa),point_array(ob);dt=times[-1]-times[-3]
    if min(len(x),len(y))<8:return dict(estimate=old.empty('INSUFFICIENT_POINTS','POT_'+method),vectors=[],quality_score=0.,timing=dict(total_ms=(time.perf_counter()-start)*1000))
    spatial=ot.dist(x,y);maxd=300*dt+2
    aa=np.array([o['area'] for o in oa]);ab=np.array([o['area'] for o in ob]);pa=np.array([o['peak'] for o in oa]);pb=np.array([o['peak'] for o in ob])
    cost=spatial/distance_scale**2+.25*(np.log(aa[:,None]/ab[None,:])**2+((pa[:,None]-pb[None,:])/40)**2)
    cost[spatial>maxd**2]=100.;a=np.full(len(x),1/len(x));b=np.full(len(y),1/len(y));tick=time.perf_counter()
    messages=io.StringIO()
    with contextlib.redirect_stdout(messages),warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        if method=='unbalanced':plan,log=ot.unbalanced.sinkhorn_unbalanced(a,b,cost,.05,1.,numItermax=500,stopThr=1e-7,log=True)
        else:plan,log=ot.partial.entropic_partial_wasserstein(a,b,cost,.05,m=mass,numItermax=500,stopThr=1e-7,log=True)
    solver_warning=messages.getvalue()+'\n'+'\n'.join(str(w.message) for w in caught)
    if 'numerical error' in solver_warning.lower() or not np.isfinite(plan).all():
        return dict(estimate=old.empty('OT_NUMERICAL_FAILURE','POT_'+method),vectors=[],quality_score=0.,timing=dict(model_ms=(time.perf_counter()-tick)*1000,total_ms=(time.perf_counter()-start)*1000),solver_warning=solver_warning.strip(),solver_log=log,transported_mass=float(plan.sum()))
    model_ms=(time.perf_counter()-tick)*1000;rows=plan.sum(1);prob=plan/np.maximum(rows[:,None],1e-30);target=prob@y;d=target-x;peak=prob.max(1)
    ent=-(prob*np.log(np.maximum(prob,1e-30))).sum(1)/np.log(len(y))
    good=(rows>.2/len(x))&(peak>=.5)&(ent<=.6)&np.isfinite(d).all(1)&(np.linalg.norm(d,axis=1)<=maxd)
    tick=time.perf_counter();e=old.fit_affine(x[good],d[good],dt,images[-1].shape,'POT_'+method,bootstrap=0);fit_ms=(time.perf_counter()-tick)*1000
    return dict(estimate=e,vectors=np.c_[x,d,good.astype(float)],row_mass=rows,match_peak=peak,entropy=ent,transported_mass=float(plan.sum()),solver_log=log,solver_warning=solver_warning.strip(),quality_score=float(np.mean(peak[good])/(1+e['residual'])) if good.any() and e.get('residual') is not None else 0.,timing=dict(model_ms=model_ms,fit_ms=fit_ms,total_ms=(time.perf_counter()-start)*1000))
