"""Causal channel observer. Missing/stale windows remain null; no flow filling."""
import time,collections,json
import cv2
import numpy as np
from sklearn.cluster import DBSCAN
from xcommon import ROOT,NAMES,fit_affine,empty,clean,sha
from channel import Channel
from backends import Backend,preprocess
from fit_adapter import fit_affine
REASONS={0:'VALID',1:'OUTSIDE_SAFE_WINDOW',2:'LOW_SUPPORT',3:'WARMUP_OR_FRAME_GAP',4:'MULTIRATE_DEFERRED',5:'PIV_QUALITY_REJECTED',6:'CHANNEL_REGISTRATION_FAILED',7:'BACKEND_FAILURE'}
def finite_list(a):
 a=np.asarray(a).ravel();return [v if ok else None for v,ok in zip(a.tolist(),np.isfinite(a).tolist())]

def groups_from_candidates(q,x,y,u,v,good,dt,crop,uncertainty=False):
 """DBSCAN spatial groups, then bounded observation cells for large components."""
 if len(q)<4:return []
 labels=DBSCAN(eps=55,min_samples=4,n_jobs=1).fit_predict(q);groups=[]
 for lab in np.unique(labels[labels>=0]):
  pts=q[labels==lab];cell=np.floor(pts/240).astype(int)
  for cx,cy in np.unique(cell,axis=0):
   p=pts[(cell[:,0]==cx)&(cell[:,1]==cy)]
   if len(p)<4:continue
   lo=np.maximum([cx*240,cy*240],p.min(0)-24);hi=np.minimum([(cx+1)*240,(cy+1)*240],p.max(0)+24)
   width,height=hi-lo;center=(lo+hi)/2
   support=good&(x>=lo[0])&(x<=hi[0])&(y>=lo[1])&(y<=hi[1])
   xy=np.c_[x[support]-lo[0],y[support]-lo[1]];d=np.c_[u[support],v[support]]
   est=fit_affine(xy,d,dt,(height,width),'LOCAL_CANDIDATE_GROUP',bootstrap=24 if uncertainty else 0) if dt else empty('WARMUP')
   est['reference_center_global_xy']=(center+np.array(crop[:2])).tolist()
   cov=np.cov(p.T,bias=True);eig,vec=np.linalg.eigh(cov);ratio=float(np.sqrt(eig[-1]/max(eig[0],1e-9)))
   groups.append(dict(observation_key=f'cell_{cx}_{cy}_component_{lab}',physical_id=None,scope='2D_SPATIAL_CANDIDATE_SET_NOT_CONTACT_OR_FUSION',reference_roi_global=[float(lo[0]+crop[0]),float(lo[1]+crop[1]),float(width),float(height)],center_global_xy=(p.mean(0)+np.array(crop[:2])).tolist(),visible_count=len(p),density_per_px2=len(p)/(width*height),covariance=cov.tolist(),Rg_px=float(np.sqrt(eig.sum())),axis_ratio=ratio,orientation_rad=float(np.arctan2(vec[1,-1],vec[0,-1])) if ratio>=1.15 else None,six=est))
 return groups

class Observer:
 def __init__(self,backend='B4',mode='local',gap=2,schedule=1,selected=True,uncertainty=False):
  self.channel=Channel();self.backend=Backend(backend);self.mode=mode;self.gap=gap;self.schedule=schedule;self.selected=selected;self.uncertainty=uncertainty
  self.history=collections.deque(maxlen=gap+1);self.fixed=None;self.background=None;self.last_pts=np.full(self.channel.x.shape,np.nan);self.prev_selected=np.zeros(self.channel.x.shape,bool);self.last_frame=None
  self.mask_sha=sha(ROOT/'configs/channel_mask.json')
 def process(self,fid,pts,gray,received=None,session='F10X_OFFLINE'):
  start=time.perf_counter();received=received or start;c=self.channel
  discontinuity=bool(self.last_frame is not None and (fid!=self.last_frame+1 or pts<=self.history[-1][1] or pts-self.history[-1][1]>.15))
  if discontinuity:self.history.clear();self.last_pts[:]=np.nan;self.prev_selected[:]=False
  self.last_frame=fid;s=c.scan(gray);a=c.cut(gray)
  if self.fixed is None:self.fixed=a.astype(np.float32);self.background=self.fixed.copy()
  if self.mode=='local':signal=s['signal'].astype(np.float32)
  else:signal=preprocess(a,self.mode,self.background if self.mode=='causal' else self.fixed)
  cv2.accumulateWeighted(a,self.background,.08)
  self.history.append((fid,pts,signal))
  selected=s['selected'];high=s['count']>=2;new=selected&~self.prev_selected;due=selected&(high|new|(fid%self.schedule==0)) if self.selected else c.safe.copy();self.prev_selected=selected.copy()
  shape=c.x.shape;u=np.full(shape,np.nan);v=u.copy();snr=u.copy();good=np.zeros(shape,bool);fresh=np.zeros(shape,bool);reason=np.full(shape,2,np.uint8);reason[~c.safe]=1;reason[selected&~due]=4;dt=None;timings={};failure=None
  if not s['registration']['valid']:reason[c.safe]=6
  elif len(self.history)<=self.gap:reason[c.safe]=3
  elif due.any():
   prev=self.history[0];dt=pts-prev[1]
   try:
    # Shared signal is passed as raw float, avoiding repeated Gaussian filtering.
    selection=due if self.selected else None
    r=self.backend.run(prev[2],signal,dt,mode='raw',do_fit=False,selection=selection,actual_coordinates=True)
    u,v,snr=r['u'],r['v'],r['snr'];snr=np.full(shape,np.nan) if snr is None else snr
    good=r['valid']&due&c.safe&selected;fresh=due.copy();reason[due&selected]=5;reason[good]=0;self.last_pts[due]=pts
    timings={k:r[k] for k in ['core_ms','validation_ms','preprocess_ms']}
   except Exception as e:failure=repr(e);reason[due]=7
  t=time.perf_counter();groups=groups_from_candidates(s['candidates'],c.x,c.y,u,v,good,dt,c.crop,self.uncertainty);group_ms=(time.perf_counter()-t)*1000
  age=(pts-self.last_pts)*1000
  windows=dict(x=(c.x+c.crop[0]).ravel().tolist(),y=(c.y+c.crop[1]).ravel().tolist(),vx=finite_list(np.where(good,u/dt if dt else np.nan,np.nan)),vy=finite_list(np.where(good,v/dt if dt else np.nan,np.nan)),valid=good.ravel().tolist(),fresh=fresh.ravel().tolist(),last_update_pts=finite_list(self.last_pts),age_ms=finite_list(age),stale=(~fresh).ravel().tolist(),reason_code=reason.ravel().tolist(),snr=finite_list(snr),candidate_count=s['count'].ravel().tolist(),active=selected.ravel().tolist(),safe=c.safe.ravel().tolist())
  now=time.perf_counter();state=dict(schema_version='M4-channel-3.0',source='CHRONOS_FILE_CAUSAL',session_id=session,frame_id=fid,source_PTS=pts,host_receive_timestamp=received,processing_finished_timestamp=now,ROI=c.crop,calibration_status='PIXEL_ONLY_NO_PHYSICAL_SCALE',age_ms=(now-received)*1000,stale_or_discontinuous=discontinuity,channel_mask_sha256=self.mask_sha,channel_mask_valid=s['registration']['valid'],backend=self.backend.name,device='CUDA' if self.backend.name in ['B3','B4'] else 'CPU',reason_dictionary=REASONS,window_shape=list(shape),windows=windows,groups=groups,registration=s['registration'],candidates_global_xy=(s['candidates']+np.array(c.crop[:2])).tolist(),candidate_areas=s['areas'].tolist(),excluded_candidate_count=len(s['excluded_candidates']),time_window=dict(start_pts=self.history[0][1],end_pts=pts,dt_s=dt,midpoint_age_ms=dt*500 if dt else None),mode=self.mode,schedule=self.schedule,selected=self.selected,scan_ms=s['scan_ms'],group_ms=group_ms,component_times=timings,backend_failure=failure,valid_window_fraction_channel=float(good.sum()/max(1,c.safe.sum())),fresh_valid_group_count=sum(all(g['six']['valid'].values()) for g in groups),units=['px/s','px/s','1/s','1/s','1/s','1/s'],observability_scope='2D_VISIBLE_BUBBLE_PATTERN_APPARENT_MOTION',scientific_accuracy_verified=False,hardware_control_executed=False,training_ready=False)
  state['groups']=clean(groups);state['algorithm_ms']=(time.perf_counter()-start)*1000
  return state
