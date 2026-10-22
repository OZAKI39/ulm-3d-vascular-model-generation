"""A single frame contract. Time fields never imply an exposure timestamp."""
from dataclasses import dataclass,field
import time,uuid,threading,json,struct
import numpy as np
@dataclass
class FramePacket:
 source_mode:str
 source_name:str
 session_id:str
 frame_id:int
 source_pts:float|None
 image:np.ndarray
 original_resolution:tuple
 pixel_format:str='GRAY8'
 capture_timestamp:float|None=None
 capture_timestamp_status:str='EXPOSURE_TIMESTAMP_UNKNOWN'
 host_decode_timestamp:dict=field(default_factory=dict)
 host_receive_timestamp:dict=field(default_factory=dict)
 arrival_monotonic:float=field(default_factory=time.monotonic)
 roi:list|None=None
 source_status:str='FRAME'
 connection_state:str='STREAMING'
 timestamp_quality:str='SOURCE_PTS_KNOWN_EXPOSURE_UNKNOWN'
 age_upper_bound_at_arrival_ms:float=0.
 decode_ms:float=0.
 def validate(self):
  if self.source_mode not in ['REALTIME_REPLAY','LIVE_CAMERA']:raise ValueError('SOURCE_MODE')
  if not self.session_id or not isinstance(self.frame_id,int) or self.frame_id<0:raise ValueError('FRAME_ID')
  a=self.image
  if not isinstance(a,np.ndarray) or a.dtype!=np.uint8 or a.size==0:raise ValueError('EMPTY_OR_DTYPE')
  if a.ndim not in [2,3] or (a.ndim==3 and a.shape[2]!=3):raise ValueError('SHAPE')
  if self.pixel_format not in ['GRAY8','BGR8'] or (a.ndim==2)!=(self.pixel_format=='GRAY8'):raise ValueError('PIXEL_FORMAT')
  w,h=self.original_resolution
  if w<1 or h<1 or w*h>16_000_000:raise ValueError('RESOLUTION')
  roi=self.roi or [0,0,w,h];x,y,rw,rh=roi
  if min(x,y)<0 or min(rw,rh)<1 or x+rw>w or y+rh>h or a.shape[:2]!=(rh,rw):raise ValueError('ROI_SHAPE')
  if self.source_pts is not None and (not np.isfinite(self.source_pts) or self.source_pts<0):raise ValueError('PTS')
  if not np.isfinite(self.arrival_monotonic) or self.age_upper_bound_at_arrival_ms<0:raise ValueError('CLOCK')
  return self
 def age_ms(self):return self.age_upper_bound_at_arrival_ms+1000*(time.monotonic()-self.arrival_monotonic)

class Continuity:
 def __init__(self,max_gap=.25):self.last=None;self.max_gap=max_gap
 def check(self,p):
  p.validate();reset=False;reason='OK'
  if self.last:
   s,i,t,shape=self.last
   if s!=p.session_id:reset=True;reason='SESSION_RESET'
   elif p.frame_id<=i:raise ValueError('DUPLICATE_OR_REVERSED_FRAME_ID')
   elif shape!=p.image.shape:raise ValueError('RESOLUTION_CHANGED_WITHOUT_SESSION')
   elif t is not None and p.source_pts is not None and (p.source_pts<=t or p.source_pts-t>self.max_gap):reset=True;reason='PTS_DISCONTINUITY'
  else:reset=True;reason='NEW_SESSION'
  self.last=(p.session_id,p.frame_id,p.source_pts,p.image.shape)
  return reset,reason

class LatestSlot:
 def __init__(self):self.value=None;self.lock=threading.Lock();self.event=threading.Event();self.overwritten=0
 def put(self,x):
  with self.lock:
   if self.value is not None:self.overwritten+=1
   self.value=x;self.event.set()
 def take(self):
  with self.lock:x=self.value;self.value=None;self.event.clear();return x
 def peek(self):
  with self.lock:return self.value

def read_cv2frame(stream):
 # Same framing limits as original WindowsLatestBridge, for conformance/fault tests.
 from live_adapter import read_exact
 p=read_exact(stream,12)
 if p[:8]!=b'CV2FRAME':raise ValueError('MAGIC')
 n=struct.unpack('<I',p[8:])[0]
 if not 0<n<=65536:raise ValueError('HEADER_LENGTH')
 h=json.loads(read_exact(stream,n));size=h['payload_bytes']
 if not isinstance(size,int) or not 0<=size<=16_000_000:raise ValueError('PAYLOAD_LENGTH')
 b=read_exact(stream,size)
 if h['status']=='FRAME':
  shape=h['shape']
  if len(shape)!=3 or shape[2]!=3 or min(shape)<=0 or int(np.prod(shape))!=size:raise ValueError('PAYLOAD_SHAPE')
  return h,np.frombuffer(b,np.uint8).reshape(shape)
 return h,None
