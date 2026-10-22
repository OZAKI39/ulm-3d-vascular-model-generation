import time,json
from collections import deque
import cv2,numpy as np
from ycommon import ROOT
from channel import Channel
from backends import Backend
from aggregation import Aggregation
from microscope_detector import load_profile,profiled_channel
class Observer:
 def __init__(self,cfg):
  cv2.setNumThreads(2);self.cfg=cfg;self.backend=Backend(cfg['backend']);self.channel=None;self.geometry=None;self.session=None;self.history=deque(maxlen=4);self.agg=Aggregation();self.background=None
  self.detector_profile=load_profile(cfg['detector_profile']) if cfg.get('detector_profile') else None
 def reset(self):
  self.history.clear();self.agg=Aggregation();self.background=None
  if hasattr(self,'live_ref'):del self.live_ref
  if self.channel:self.channel.prev=None;self.channel.ref=None
 def geometry_for(self,p):
  path=self.cfg['replay_geometry'] if p.source_mode=='REALTIME_REPLAY' else self.cfg.get('live_geometry')
  if not path:return None
  g=json.loads(open(path).read())
  if tuple(g['image_wh'])!=p.original_resolution or p.image.shape[:2]!=(p.original_resolution[1],p.original_resolution[0]):return None
  if p.source_mode=='REALTIME_REPLAY' and self.cfg.get('source_sha256')!=g.get('source_sha256'):return None
  if p.source_mode=='LIVE_CAMERA' and (not g.get('human_verified') or g.get('source_mode')!='LIVE_CAMERA'):return None
  return g
 def analyze(self,p,reset=False):
  t=time.perf_counter();gray=p.image if p.image.ndim==2 else cv2.cvtColor(p.image,cv2.COLOR_BGR2GRAY)
  g=self.geometry_for(p)
  if g is None:return dict(valid=False,reason='MASK_NOT_VERIFIED',candidates=[],vectors=[],aggregation=None,timings={'analysis_ms':(time.perf_counter()-t)*1000},geometry=None)
  if self.detector_profile and list(p.original_resolution)!=self.detector_profile['camera_image_wh']:
   return dict(valid=False,reason='DETECTOR_PROFILE_RESOLUTION_MISMATCH',candidates=[],vectors=[],aggregation=None,timings={'analysis_ms':(time.perf_counter()-t)*1000},geometry=None)
  if g!=self.geometry:self.channel=profiled_channel(g,self.detector_profile) if self.detector_profile else Channel(g);self.geometry=g;reset=True
  if reset:self.reset()
  ch=self.channel;sc=ch.scan(gray)
  # Generic LIVE registration uses explicit user-selected fixed structure, never
  # claims the historical x=450:600 strip is valid in a new optical setup.
  registration=sc['registration'];require_reference=True
  if p.source_mode=='LIVE_CAMERA':
   box=g.get('registration_roi')
   if not box:
    registration=dict(valid=False,method='LIVE_REFERENCE_NOT_VERIFIED')
    # Interactive masks may explicitly request unregistered image-coordinate
    # observation. No wall-reference check is forged or labelled as passed.
    require_reference=g.get('motion_reference')!='IMAGE_SPACE_UNREGISTERED'
   else:
    x,y,w,h=box;ref=gray[y:y+h,x:x+w].astype(np.float32)
    if not hasattr(self,'live_ref') or reset:self.live_ref=ref.copy()
    shift,response=cv2.phaseCorrelate(self.live_ref,ref);registration=dict(valid=bool(response>=.08 and np.linalg.norm(shift)<=6),shift_xy_px=list(shift),response=response)
  local=ch.cut(gray);q=sc['candidates'];self.background=local.astype(np.float32) if self.background is None else .92*self.background+.08*local
  timestamp=p.source_pts if p.source_pts is not None else p.arrival_monotonic # only local receive intervals for LIVE; explicit approximation
  self.history.append((timestamp,local,self.background.copy()));vectors=np.empty((0,4));piv_ms=0;computed=0;validmask=np.zeros_like(ch.safe);dt=None
  eligible=[item for item in list(self.history)[:-1] if 0<timestamp-item[0]<=self.cfg['max_pts_gap_s']]
  if eligible:
   pt,prev,bg=eligible[0];dt=timestamp-pt
   if 0<dt<=self.cfg['max_pts_gap_s']:
    z=self.backend.run(prev,local,dt,mode='causal',background_a=bg,background_b=self.background,do_fit=False,actual_coordinates=True)
    validmask=z['valid']&ch.safe&sc['selected'];computed=int(ch.safe.sum());piv_ms=z['total_ms']
    if require_reference and not registration['valid']:validmask[:]=False
    vectors=np.c_[ch.x[validmask],ch.y[validmask],z['u'][validmask]/dt,z['v'][validmask]/dt]
  grad=cv2.Laplacian(local,cv2.CV_32F);focus=float(np.var(grad[ch.local.astype(bool)]));illum=float(np.mean(local[ch.local.astype(bool)]))
  near=(q[:,0]<20)|(q[:,1]<20)|(q[:,0]>local.shape[1]-20)|(q[:,1]>local.shape[0]-20) if len(q) else np.zeros(0,bool)
  ta=time.perf_counter();agg=self.agg.update(timestamp,q,vectors[:,:2],vectors[:,2:],local.shape,focus,illum,float(near.mean()) if len(q) else 0.,support_mask=ch.local)
  ams=(time.perf_counter()-ta)*1000;quality=np.clip(sc['areas']/30,0,1)
  return dict(valid=bool(len(vectors)>=8 and len(q)>=8),reason='OK' if len(vectors)>=8 and len(q)>=8 else 'LOW_PIV_SUPPORT_OR_REGISTRATION',motion_reference='IMAGE_SPACE_UNREGISTERED' if not require_reference else 'FIXED_WALL_CHECK',motion_reference_verified=bool(require_reference and registration['valid']),candidates=q,candidate_areas=sc['areas'],candidate_quality=quality,boundary_candidates=int(near.sum()),low_quality_candidates=int((quality<.5).sum()),vectors=vectors,computed_windows=computed,safe_windows=int(ch.safe.sum()),valid_windows=int(validmask.sum()),active_windows=int(sc['selected'].sum()),window_grid_xy=np.c_[ch.x.ravel(),ch.y.ravel()],window_valid=validmask.ravel(),registration=registration,focus=focus,illumination=illum,geometry=g,detector_profile_id=self.detector_profile['id'] if self.detector_profile else 'LEGACY_OPENCV',aggregation=agg,dt_s=dt,motion_time_basis='SOURCE_PTS' if p.source_pts is not None else 'LOCAL_RECEIVE_INTERVAL_APPROXIMATION',timings={'scan_ms':sc['scan_ms'],'piv_ms':piv_ms,'aggregation_ms':ams,'analysis_ms':(time.perf_counter()-t)*1000})
