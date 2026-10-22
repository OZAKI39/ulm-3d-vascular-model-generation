import time,json
import cv2,numpy as np,torch
from pathlib import Path
from networks import build_model
from geometry import endpoint_six,transform,original_transform
from heatmaps import HeatmapDetector
from vcommon import SCALES,NAMES,LIMITS,P,sha

def raw_geo_values(H_normalized,dt):
    T=np.array([[1/128,0,-1],[0,1/128,-1],[0,0,1.]])
    H=np.linalg.solve(T,np.asarray(H_normalized)@T)
    return endpoint_six(H,dt,(256,256))

def convert_shape(values,dt,shape):
    H=transform(values,dt,[128,128]);H=original_transform(H,shape,(256,256))
    return endpoint_six(H,dt,shape)

def shape_jacobian(shape):
    return np.array([convert_shape(np.eye(6)[i],.01,shape) for i in range(6)]).T

class Predictor:
    def __init__(self,checkpoint,device='cuda',amp=True):
        self.checkpoint=Path(checkpoint);c=torch.load(checkpoint,map_location='cpu',weights_only=False)
        self.name=c['config']['name'];self.model=build_model(self.name).to(device);self.model.load_state_dict(c['model']);self.model.eval();self.device=device;self.amp=amp and device.startswith('cuda');self.sha=sha(checkpoint)
        self.parameter_count=sum(p.numel() for p in self.model.parameters())
    @torch.inference_mode()
    def batch(self,images,dt):
        with torch.autocast(self.device.split(':')[0],dtype=torch.bfloat16,enabled=self.amp):out=self.model(images.to(self.device,non_blocking=True),dt.to(self.device,non_blocking=True))
        return {k:v.float().cpu().numpy() for k,v in out.items() if k in ['mean','log_sigma','logits','H_normalized']}
    def outputs(self,raw,dt):
        means=raw['mean']*SCALES;sigmas=np.exp(raw['log_sigma'])*SCALES
        logits=raw['logits'];prob=np.exp(logits-logits.max(1,keepdims=True));prob/=prob.sum(1,keepdims=True)
        valid=np.isfinite(means).all(1)
        if 'H_normalized' in raw:
            for i,H in enumerate(raw['H_normalized']):
                try:means[i]=raw_geo_values(H,float(np.sum(dt[i])))
                except (ValueError,np.linalg.LinAlgError):valid[i]=False;means[i]=np.nan
        # Lower means more trustworthy; not a calibrated error probability.
        scores=np.max(sigmas/LIMITS,axis=1)+(1-prob[:,0])*5
        scores[~valid]=np.inf
        return dict(mean=means,sigma=sigmas,probability=prob,score=scores,algebra_valid=valid)
    def frames(self,grays,heats,pts):
        dt=np.diff(pts);start=time.perf_counter()
        if len(dt)!=2 or np.min(dt)<=0:raise ValueError('THREE_CAUSAL_INCREASING_PTS_REQUIRED')
        gray=np.array([cv2.resize(g,(128,128),interpolation=cv2.INTER_AREA) for g in grays],np.float32)/255
        heat=np.array([cv2.resize(h,(128,128)) for h in heats],np.float32)
        x=np.stack([(gray-.65)/.25,heat],1)[None];preprocess_ms=(time.perf_counter()-start)*1000
        t=time.perf_counter();raw=self.batch(torch.from_numpy(x),torch.tensor(dt[None],dtype=torch.float32));inference_ms=(time.perf_counter()-t)*1000
        t=time.perf_counter();out={k:v[0] for k,v in self.outputs(raw,dt[None]).items()};shape=grays[0].shape
        if tuple(shape)!=(256,256) and out['algebra_valid']:
            out['mean']=convert_shape(out['mean'],sum(dt),shape);out['sigma']=np.abs(shape_jacobian(shape))@out['sigma']
        return out,dict(resize_normalize_ms=preprocess_ms,transfer_model_return_ms=inference_ms,geometry_ms=(time.perf_counter()-t)*1000)

def estimate_from_prediction(pred,calibration,name,weight_sha):
    label=int(np.argmax(pred['probability']));threshold=calibration.get('threshold');valid_algebra=bool(pred['algebra_valid'])
    accepted=valid_algebra and label==0 and threshold is not None and pred['score']<=threshold
    per_threshold=calibration.get('field_thresholds',[threshold]*6)
    per_valid=[accepted and t is not None and pred['score']<=t for t in per_threshold]
    # Full gate is evaluated separately; fields may be selectively available.
    status=P['observation_classes'][label]
    if label==0 and not any(per_valid):status='NOT_OBSERVABLE'
    values={n:float(pred['mean'][j]) if per_valid[j] else None for j,n in enumerate(NAMES)}
    reasons={n:None if per_valid[j] else ('ALGEBRA_INVALID' if not valid_algebra else (P['observation_classes'][label] if label!=0 else 'VALIDATION_RISK_GATE_REJECTED')) for j,n in enumerate(NAMES)}
    intervals={}
    for key,q in calibration.get('quantiles',{}).items():
        intervals[key]={n:[float(pred['mean'][j]-q[j]*pred['sigma'][j]),float(pred['mean'][j]+q[j]*pred['sigma'][j])] if per_valid[j] else None for j,n in enumerate(NAMES)}
    return dict(method=name,values=values,valid=dict(zip(NAMES,per_valid)),reasons=reasons,observability_status=status,full_six_accepted=accepted,raw_estimate=pred['mean'].tolist() if valid_algebra else None,
        uncertainty=dict(method='grouped_rank_scaled_residual',calibration_version='F10V-synthetic-groups-v1',real_domain_calibrated=False,intervals=intervals),sample_count=0,confidence_score=float(pred['score']) if np.isfinite(pred['score']) else None,classification_probabilities=pred['probability'].tolist(),weight_sha256=weight_sha)
