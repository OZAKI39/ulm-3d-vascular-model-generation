"""Editable geometry, whole-window validity and causal whole-channel scanning."""
import json,time
import cv2
import numpy as np
from xcommon import ROOT
from openpiv import pyprocess
class Channel:
 def __init__(self,config=None,ws=48,ov=24):
  self.config=config or json.loads((ROOT/'configs/channel_mask.json').read_text());c=self.config;w,h=c['image_wh'];self.ws=ws;self.ov=ov
  self.channel=np.zeros((h,w),np.uint8);cv2.fillPoly(self.channel,[np.asarray(c['channel_polygon'],np.int32)],1)
  self.excluded=np.zeros_like(self.channel)
  for p in c['excluded_polygons']:cv2.fillPoly(self.excluded,[np.asarray(p,np.int32)],1)
  self.valid=self.channel*(1-self.excluded);self.crop=c['crop_xywh'];x,y,w,h=self.crop;self.local=self.valid[y:y+h,x:x+w]
  self.x,self.y=pyprocess.get_coordinates(self.local.shape,ws,ov,center_on_field=False)
  self.x0=np.rint(self.x-ws/2).astype(int);self.y0=np.rint(self.y-ws/2).astype(int)
  self.coverage=self.window_mean(self.local.astype(np.float32));self.safe=self.coverage>=c['minimum_valid_window_fraction']
  self.ref=None;self.prev=None;self.registration=None
 def cut(self,a):
  x,y,w,h=self.crop;return a[y:y+h,x:x+w]
 def window_mean(self,a):
  ii=cv2.integral(np.asarray(a,np.float32),sdepth=cv2.CV_64F);x,y=self.x0,self.y0;w=self.ws
  return (ii[y+w,x+w]-ii[y,x+w]-ii[y+w,x]+ii[y,x])/(w*w)
 def scan(self,gray):
  start=time.perf_counter();a=self.cut(gray);signal=cv2.subtract(cv2.GaussianBlur(a,(0,0),4),a)
  # Frozen OpenCV connected-component criteria from F10U; no material IDs.
  binary=np.uint8(signal>=10);n,labels,stats,centers=cv2.connectedComponentsWithStats(binary,8)
  s=stats[1:];q=centers[1:];keep=(s[:,4]>=8)&(s[:,4]<=450)&(s[:,2]/np.maximum(s[:,3],1)<=2.8)&(s[:,3]/np.maximum(s[:,2],1)<=2.8)&(s[:,0]>2)&(s[:,1]>2)&(s[:,0]+s[:,2]<a.shape[1]-2)&(s[:,1]+s[:,3]<a.shape[0]-2)
  allq=q[keep];allareas=s[keep,4];inside=self.local[np.rint(allq[:,1]).astype(int),np.rint(allq[:,0]).astype(int)].astype(bool)
  q=allq[inside];areas=allareas[inside];pointmap=np.zeros_like(a,np.float32)
  if len(q):np.add.at(pointmap,(np.rint(q[:,1]).astype(int),np.rint(q[:,0]).astype(int)),1)
  counts=self.window_mean(pointmap)*self.ws**2
  texture=np.sqrt(np.maximum(self.window_mean(signal.astype(np.float32)**2)-self.window_mean(signal)**2,0))
  activity=self.window_mean(cv2.absdiff(a,self.prev)) if self.prev is not None else np.zeros_like(counts)
  # A candidate alone triggers checking; motion is not required for activation.
  selected=self.safe&((counts>=1)|((texture>=2)&(activity>=.4)))
  # Registration uses a fixed structural crop, never moving bubble candidates.
  reference=gray[:,450:600].astype(np.float32);reference=cv2.resize(reference,(75,280));reference-=reference.mean()
  if self.ref is None:self.ref=reference.copy();shift=(0.,0.);response=1.
  else:shift,response=cv2.phaseCorrelate(self.ref,reference)
  dx,dy=shift[0]*2,shift[1]*4;registered=bool(response>=self.config['registration_min_response'] and np.hypot(dx,dy)<=self.config['registration_max_translation_px'])
  self.prev=a.copy();self.registration={'shift_xy_px':[dx,dy],'response':response,'valid':registered,'method':'OpenCV phaseCorrelate first-frame fixed wall strip','background_reference_velocity':None}
  return dict(signal=signal,candidates=q,areas=areas,excluded_candidates=allq[~inside],count=counts,texture=texture,activity=activity,selected=selected,registration=self.registration,scan_ms=(time.perf_counter()-start)*1000)
