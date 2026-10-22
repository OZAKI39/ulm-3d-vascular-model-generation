"""Bounded library evaluation. No device control, RTSP reimplementation or learned model.
Library algorithms: OpenCV detection/shape/flow and Norfair association/Kalman.
Replay subclass changes only the producer; snapshot/start/stop are original Chronos.
"""
from pathlib import Path
import sys,json,time,threading,argparse,math,collections,subprocess,hashlib
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT.parent/'CHRONOS_control'))
from chronos_common import LatestFrameStream,load_config,utc_now
import av,cv2,numpy as np,psutil,zmq
from norfair import Detection,Tracker
cv2.setNumThreads(1)
CFG=json.loads((ROOT/'protocol.json').read_text())
VIDEO=json.loads((ROOT/'evidence/input_identity.json').read_text())['path']
DT=.0333

def dump(path,obj):Path(path).write_text(json.dumps(obj,indent=2,allow_nan=False))
def diststats(v):
 a=np.asarray(v,float)
 return {k:float(f(a)) for k,f in [('mean',np.mean),('p50',lambda x:np.percentile(x,50)),('p95',lambda x:np.percentile(x,95)),('max',np.max)]} if len(a) else None

class ReplayStream(LatestFrameStream):
 def __init__(self,start,end,fault=False):
  super().__init__(load_config()); self.begin=start; self.end=end; self.done=threading.Event(); self.records=[]; self.fault=fault; self.error=None
 def _run(self):
  try:
   with av.open(VIDEO) as c:
    s=c.streams.video[0]; s.codec_context.thread_count=1
    # Decode preroll strictly before the test; do not process it or expose future images.
    c.seek(int(self.begin*DT/s.time_base),stream=s,backward=True)
    it=c.decode(s); first=None
    for f in it:
     if int(round(float(f.time)/DT))>=self.begin:first=f;break
    origin=time.monotonic(); offset=0
    for index in range(self.begin,self.end):
     if self.fault and index==self.begin+35:
      with self.lock:self.latest=None;self.state='reconnecting'
      self.stop_event.wait(.35);offset+=.35
     due=origin+(index-self.begin)*DT+offset
     if self.stop_event.wait(max(0,due-time.monotonic())):break
     t=time.monotonic(); f=first if index==self.begin else next(it); bgr=f.to_ndarray(format='bgr24'); received=time.monotonic()
     fid=int(round(float(f.time)/DT)); pts=float(f.time)
     if self.fault and self.begin+65<=index<self.begin+72:pts-=.3
     rec={'frame_id':fid,'pts_s':pts,'due_mono':due,'decode_start_mono':t,'receive_mono':received,'decode_ms':(received-t)*1000,'decode_lateness_ms':(received-due)*1000}
     self.records.append(rec)
     with self.lock:
      self.frames+=1;self.latest=(fid,received,utc_now(),bgr); self.state='streaming';self.first_received=self.first_received or received; self.last_received=received
    self.state='eof'
  except Exception as e:self.error=repr(e)
  finally:self.done.set()

class Detector:
 def __init__(self,kind):
  self.kind=kind
  p=cv2.SimpleBlobDetector_Params();p.minThreshold=6;p.maxThreshold=30;p.thresholdStep=6;p.minRepeatability=2
  p.filterByColor=True;p.blobColor=255;p.filterByArea=True;p.minArea=8;p.maxArea=450;p.filterByCircularity=False;p.filterByConvexity=False;p.filterByInertia=False;p.minDistBetweenBlobs=5
  self.blob=cv2.SimpleBlobDetector_create(p)
 def __call__(self,bgr):
  x,y,w,h=CFG['roi_xywh']; g=cv2.cvtColor(bgr[y:y+h,x:x+w],cv2.COLOR_BGR2GRAY)
  signal=cv2.subtract(cv2.GaussianBlur(g,(0,0),4),g)
  mask=np.full(g.shape,255,np.uint8)
  for x1,y1,x2,y2 in CFG['fixed_structure_exclusion_xyxy']:
   mask[max(0,y1-y):min(h,y2-y),max(0,x1-x):min(w,x2-x)]=0
  signal=cv2.bitwise_and(signal,mask); binary=np.uint8(signal>=CFG['contrast_threshold'])*255
  cc,labs,stats,centers=cv2.connectedComponentsWithStats(binary,8)
  allowed=[]
  for k in range(1,cc):
   xx,yy,ww,hh,area=stats[k]
   if not(8<=area<=450 and max(ww/hh,hh/ww)<=2.8 and min(xx,yy)>2 and xx+ww<w-2 and yy+hh<h-2):continue
   peak=float(signal[yy:yy+hh,xx:xx+ww].max()); objects={'xy':[float(centers[k,0]+x),float(centers[k,1]+y)],'area':int(area),'bbox':[int(xx+x),int(yy+y),int(ww),int(hh)],'peak':peak,'label':int(k),'classification':'VISIBLE_SINGLE_CANDIDATE' if area<80 and max(ww,hh)<22 else 'BLUR_OR_OVERLAP_UNKNOWN','bubble_identity_verified':False}
   allowed.append(objects)
  if self.kind=='blob':
   keys=self.blob.detect(signal,mask); selected=[]
   for kp in keys:
    pos=np.array(kp.pt)+[x,y]
    if allowed:
     j=int(np.argmin([np.linalg.norm(pos-a['xy']) for a in allowed])); d=allowed[j]
     if np.linalg.norm(pos-d['xy'])<6 and d['label'] not in [a['label'] for a in selected]:selected.append({**d,'xy':pos.tolist()})
   allowed=selected
  return allowed,g,labs,mask

def shape_of(obj,labs):
 x,y,w,h=obj['bbox'];rx,ry,_,_=CFG['roi_xywh']; sub=np.uint8(labs[y-ry:y-ry+h,x-rx:x-rx+w]==obj['label'])
 m=cv2.moments(sub,True); n=m['m00']
 if n<3:return None
 tensor=np.array([[m['mu20'],m['mu11']],[m['mu11'],m['mu02']]])/n
 eig,vec=np.linalg.eigh(tensor);axis=vec[:,-1]
 contours,_=cv2.findContours(sub,cv2.RETR_EXTERNAL,cv2.CHAIN_APPROX_SIMPLE)
 contour=max(contours,key=cv2.contourArea).reshape(-1,2)+[x,y] if contours else np.empty((0,2))
 return {'area_px2':n,'centroid_xy_px':obj['xy'],'tensor_px2':tensor.tolist(),'aspect_ratio':float(np.sqrt(eig[1]/max(eig[0],.1))),'major_axis_deg':float(np.degrees(np.arctan2(axis[1],axis[0]))%180),'contour_xy_px':contour.tolist(),'scope':'THRESHOLDED_VISUAL_REGION_ONLY','group_shape_valid':False,'group_shape':None}

class MockController:
 """Local ZeroMQ consumer: never imports or calls Tektronix modules."""
 def __init__(self):
  self.context=zmq.Context();self.ready=threading.Event();self.rows=[];self.stop=threading.Event();self.thread=threading.Thread(target=self.loop,daemon=True);self.thread.start();self.ready.wait(3)
  self.tx=self.context.socket(zmq.PUSH);self.tx.setsockopt(zmq.SNDHWM,1);self.tx.setsockopt(zmq.LINGER,0);self.tx.connect(self.endpoint)
 def loop(self):
  rx=self.context.socket(zmq.PULL);rx.setsockopt(zmq.RCVHWM,1);rx.setsockopt(zmq.CONFLATE,1);port=rx.bind_to_random_port('tcp://127.0.0.1');self.endpoint=f'tcp://127.0.0.1:{port}';self.ready.set();last=-1
  while not self.stop.is_set():
   if not rx.poll(20):continue
   obj=rx.recv_json();now=time.monotonic(); age=(now-obj['receive_monotonic_s'])*1000
   fresh=0<=age<=CFG['stale_limit_ms'] and obj['frame_id']>last and not obj['timestamp_discontinuity']
   valid=bool(fresh and any(t['visual_valid'] for t in obj['targets']))
   self.rows.append({'frame_id':obj['frame_id'],'received_mono':now,'age_ms':age,'transfer_ms':(now-obj['send_monotonic_s'])*1000,'fresh':fresh,'valid_visual':valid,'bubble_control_admitted':False,'status':'OBSERVED' if valid else 'UNKNOWN','reason':'UNVALIDATED_BUBBLE_IDENTITY_NO_HARDWARE_PERMISSION'})
   last=max(last,obj['frame_id'])
  rx.close(0)
 def send(self,obj):
  obj['send_monotonic_s']=time.monotonic()
  try:self.tx.send_json(obj,flags=zmq.NOBLOCK);return True
  except zmq.Again:return False
 def close(self):
  time.sleep(.03);self.stop.set();self.thread.join(2);self.tx.close(0);self.context.term()

def run(method,start,end,tag,delay=0,fault=False):
 out=ROOT/'results'/tag;out.mkdir(exist_ok=False)
 detector=Detector('blob' if method=='blob_norfair' else 'cc')
 tracker=Tracker(distance_function='euclidean',distance_threshold=8,initialization_delay=2,hit_counter_max=5,past_detections_length=3)
 mock=MockController(); stream=ReplayStream(start,end,fault).start(); process=psutil.Process(); cpu0=process.cpu_times();wall0=time.monotonic()
 log=[]; messages=[]; seen={}; prevgray=None;previd=None;prevpts=None;prevactive=set();lost_events=[];births=[];flow_outputs=[];resets=[]
 while True:
  packet=stream.snapshot()
  if packet is None or packet[0]==previd:
   if stream.done.is_set():break
   time.sleep(.0005);continue
  fid,receive,utc,bgr=packet; t0=time.monotonic();rec=next(r for r in reversed(stream.records) if r['frame_id']==fid);pts=rec['pts_s']; gap=1 if previd is None else fid-previd
  discontinuity=prevpts is not None and (pts<=prevpts or abs((pts-prevpts)-gap*DT)>.01)
  if discontinuity or gap>5:
   resets.append({'frame_id':fid,'reason':'timestamp_discontinuity' if discontinuity else 'large_frame_gap'})
   tracker=Tracker(distance_function='euclidean',distance_threshold=8,initialization_delay=2,hit_counter_max=5);seen={};prevgray=None;prevactive=set()
  else:
   for _ in range(max(0,gap-1)):tracker.update(detections=[])
  tg=time.monotonic();det,g,labs,mask=detector(bgr);td=time.monotonic()
  detections=[Detection(np.array([o['xy']]),scores=np.array([1.]),data={**o,'frame_id':fid}) for o in det]
  objects=tracker.update(detections=detections);tt=time.monotonic();targets=[];active=set()
  xy=np.array([o['xy'] for o in det]);neardist={}
  if len(xy)>1:
   from scipy.spatial.distance import cdist
   dm=cdist(xy,xy);np.fill_diagonal(dm,np.inf);neardist={d['label']:float(dm[i].min()) for i,d in enumerate(det)}
  for o in objects:
   obs=o.last_detection.data;isobserved=obs['frame_id']==fid;tid=int(o.id);active.add(tid);isnew=tid not in seen
   if isnew:seen[tid]=collections.deque(maxlen=8);births.append({'frame_id':fid,'id':tid})
   hist=seen[tid];velocity=None;uncertainty=None;pred=None;cross='UNKNOWN'
   if isobserved:
    pos=np.asarray(obs['xy']);
    if hist and pts>hist[-1][0]:velocity=((pos-hist[-1][1])/(pts-hist[-1][0])).tolist()
    hist.append((pts,pos.copy()));sh=shape_of(obs,labs)
    motion=float(np.linalg.norm(hist[-1][1]-hist[0][1])) if len(hist)>1 else 0
    ambiguity=neardist.get(obs['label'],999)<9
    valid=len(hist)>=3 and motion>=2 and obs['peak']>=10 and obs['classification']=='VISIBLE_SINGLE_CANDIDATE' and not ambiguity and not discontinuity
    if len(hist)>=4:
     times=np.array([h[0] for h in hist]); positions=np.array([h[1] for h in hist]);coef=np.polyfit(times-times[-1],positions,1);fit=(times-times[-1])[:,None]*coef[0]+coef[1];sigma=float(max(1,np.sqrt(np.mean((fit-positions)**2))));uncertainty={'position_sigma_px':sigma,'calibrated':False,'method':'causal linear residual floor 1 px'}
     vy=coef[0,1];line=CFG['target_line_y_px'];tau=(line-pos[1])/vy if abs(vy)>1 else -1
     if valid and 0<tau<.25 and vy*CFG['target_direction_y']>0:pred={'status':'PREDICTED','xy_px':(pos+coef[0]*tau).tolist(),'seconds_to_line':float(tau),'time_sigma_s':float(sigma/max(abs(vy),1))};cross='PREDICTED'
    if valid and len(hist)>1 and hist[-2][1][1]<CFG['target_line_y_px']<=pos[1]:cross='OBSERVED'
    elif valid and cross!='PREDICTED':cross='NOT_CROSSED'
    classification='FIXED_OR_SLOW_VISUAL_FEATURE' if len(hist)>=6 and motion<2 else obs['classification']
   else:pos=o.estimate[0];sh=None;valid=False;classification='UNKNOWN';cross='LOST'
   targets.append({'id':tid,'status':'OBSERVED' if isobserved else 'LOST','new':isnew,'lost':not isobserved,'xy_px':pos.tolist() if isobserved else None,'predicted_xy_px':pos.tolist() if not isobserved else None,'velocity_px_per_replay_s':velocity,'visual_valid':bool(valid),'confidence_semantics':'uncalibrated deterministic quality gate','classification':classification,'bubble_identity_verified':False,'shape':sh,'prediction':pred,'uncertainty':uncertainty,'target_crossing':cross})
  for tid in prevactive-active:lost_events.append({'frame_id':fid,'id':tid,'status':'LOST'})
  ts=time.monotonic();flow=None
  if prevgray is not None and method.endswith('_lk'):
   points=cv2.goodFeaturesToTrack(prevgray,maxCorners=150,qualityLevel=.025,minDistance=8,mask=mask)
   if points is not None:
    dest,ok,err=cv2.calcOpticalFlowPyrLK(prevgray,g,points,None,winSize=(21,21),maxLevel=3)
    back,bok,_=cv2.calcOpticalFlowPyrLK(g,prevgray,dest,None,winSize=(21,21),maxLevel=3)
    keep=(ok[:,0]>0)&(bok[:,0]>0)&(np.linalg.norm(back[:,0]-points[:,0],axis=1)<1)&(err[:,0]<20)
    vectors=dest[:,0][keep]-points[:,0][keep];flow={'frame_id':fid,'n':int(keep.sum()),'xy_local':points[:,0][keep].tolist(),'vectors_px':vectors.tolist(),'median_px_per_replay_s':(np.median(vectors,axis=0)/(pts-prevpts)).tolist() if len(vectors) else None}
  if prevgray is not None and method.endswith('_farneback'):
   f=cv2.calcOpticalFlowFarneback(prevgray,g,None,.5,3,15,3,5,1.2,0)
   yy,xx=np.mgrid[8:g.shape[0]:24,8:g.shape[1]:24];keep=mask[yy,xx]>0;vectors=f[yy,xx][keep];flow={'frame_id':fid,'n':int(keep.sum()),'xy_local':np.stack([xx[keep],yy[keep]],axis=1).tolist(),'vectors_px':vectors.tolist(),'median_px_per_replay_s':(np.median(vectors,axis=0)/(pts-prevpts)).tolist()}
  tf=time.monotonic()
  if delay:time.sleep(delay)
  msg={'schema_version':1,'test_kind':'CAUSAL_REPLAY_TEST','frame_id':fid,'image_source':VIDEO,'capture_timestamp':None,'source_pts_s':pts,'time_unit':'ENCODED_REPLAY_SECONDS','receive_monotonic_s':receive,'receive_utc':utc,'process_done_monotonic_s':time.monotonic(),'timestamp_discontinuity':discontinuity,'stale':(time.monotonic()-receive)*1000>CFG['stale_limit_ms'],'targets':targets,'group_shape':None,'group_shape_status':'UNKNOWN_UNVALIDATED_SEGMENTATION','control_permission':False}
  sent=mock.send(msg);te=time.monotonic()
  log.append({'frame_id':fid,'pts_s':pts,'decode_ms':rec['decode_ms'],'frame_age_at_start_ms':(t0-receive)*1000,'detect_ms':(td-tg)*1000,'track_ms':((tt-td)+(tg-t0))*1000,'shape_quality_ms':(ts-tt)*1000,'flow_ms':(tf-ts)*1000,'serialize_send_ms':(te-tf-delay)*1000,'processing_ms':(te-t0)*1000,'receive_to_send_ms':(te-receive)*1000,'n_detections':len(det),'n_tracks':len(targets),'valid_visual_targets':sum(t['visual_valid'] for t in targets),'latest_slot_occupancy':int(stream.snapshot() is not None),'rss_mb':process.memory_info().rss/1e6,'sent':sent})
  messages.append(msg)
  if flow is not None:flow_outputs.append(flow)
  prevgray=g;previd=fid;prevpts=pts;prevactive=active
  if stream.done.is_set() and stream.snapshot()[0]==fid:break
 stream.stop();mock.close();wall=time.monotonic()-wall0;cpu1=process.cpu_times()
 decoded=len(stream.records);span=(end-start)*DT;valid=sum(v['valid_visual'] for v in mock.rows)
 summary={'execution_git_head':subprocess.run(['git','-C',str(ROOT),'rev-parse','HEAD'],capture_output=True,text=True).stdout.strip(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'roi_xywh':CFG['roi_xywh'],'tag':tag,'method':method,'test_kind':'CAUSAL_REPLAY_TEST','frames_requested':[start,end],'source_duration_s':span,'wall_s':wall,'decoded_frames':decoded,'processed_frames':len(log),'source_clock_fps':1/DT,'decode_rate_wall_hz':decoded/wall,'processed_rate_wall_hz':len(log)/wall,'mock_received':len(mock.rows),'mock_valid_visual_states':valid,'mock_valid_visual_rate_wall_hz':valid/wall,'bubble_control_admitted':0,'dropped_latest_frames':decoded-len(log),'drop_fraction':(decoded-len(log))/max(decoded,1),'cpu_one_core_pct':100*((cpu1.user-cpu0.user)+(cpu1.system-cpu0.system))/wall,'max_rss_mb':max([r['rss_mb'] for r in log],default=0),'gpu_used':False,'max_python_frame_slots':1,'decoder_error':stream.error,'timestamp_resets':resets,'lost_track_events':len(lost_events),'birth_events':len(births),'stage_ms':{k:diststats([r[k] for r in log]) for k in ['decode_ms','frame_age_at_start_ms','detect_ms','track_ms','shape_quality_ms','flow_ms','serialize_send_ms','processing_ms','receive_to_send_ms']},'transfer_ms':diststats([v['transfer_ms'] for v in mock.rows]),'receive_to_mock_ms':diststats([v['age_ms'] for v in mock.rows]),'CAMERA_TO_CONTROLLER_END_TO_END_LATENCY':'NOT_MEASURED','semantics':'Visual candidate quality gates are not calibrated bubble identity or biological accuracy.'}
 dump(out/'summary.json',summary);dump(out/'frames.json',log);dump(out/'decoder.json',stream.records);dump(out/'mock_received.json',mock.rows);dump(out/'events.json',{'births':births,'lost':lost_events});dump(out/'flow.json',flow_outputs)
 with (out/'states.jsonl').open('w') as f:
  for m in messages:f.write(json.dumps(m,allow_nan=False)+'\n')
 print(json.dumps({'tag':tag,'processed':len(log),'drops':summary['dropped_latest_frames'],'p95_ms':(summary['stage_ms']['processing_ms'] or {}).get('p95'),'valid_hz':summary['mock_valid_visual_rate_wall_hz'],'error':stream.error}),flush=True)
 return summary

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--suite',choices=['smoke','compare','full','fault']);p.add_argument('--target-roi',action='store_true');a=p.parse_args()
 if a.target_roi:CFG['roi_xywh']=[650,350,180,350]
 suffix='_target_roi' if a.target_roi else ''
 if a.suite=='smoke':run('cc_norfair',0,60,'smoke_cc')
 if a.suite=='compare':
  for method in ['cc_norfair','blob_norfair','cc_norfair_lk','cc_norfair_farneback']:
   for start in [120,360,720]:run(method,start,start+90,f'{method}_{start}{suffix}')
 if a.suite=='full':run('cc_norfair',0,937,'full_cc'+suffix)
 if a.suite=='fault':run('cc_norfair',120,240,'fault_cc'+suffix,delay=.12,fault=True)
