"""Bounded acquisition/analysis/recording, lightweight TTL publisher and UI snapshot."""
import threading,time,queue,json,gzip,uuid,traceback,subprocess
from collections import Counter,deque
from pathlib import Path
import numpy as np,psutil
from ycommon import clean,dump,stats,identity
from protocol import Continuity,LatestSlot
from sources import RecordedVideoFrameSource
from observer import Observer
from experiment_video import ExperimentVideoRecorder,DEFAULT_LAYERS
from piv_display import PIVDisplayFilter

def observation_record(result,diagnostics=False):
 """Keep scientific observations and provenance, omit default profiling detail."""
 if diagnostics:return result
 excluded={'timings','window_grid_xy','window_valid','computed_windows','safe_windows',
           'active_windows','candidate_areas','candidate_quality','piv_display_vectors'}
 return {key:value for key,value in result.items() if key not in excluded}

class Recorder:
 def __init__(self,path,capacity=32):
  self.q=queue.Queue(capacity);self.stop=threading.Event();self.dropped=0;self.written=0;self.times=deque(maxlen=6000);self.error=None;self.path=Path(path);self.thread=threading.Thread(target=self.run,name='state-recorder',daemon=True);self.thread.start()
 def submit(self,x):
  try:self.q.put_nowait(x)
  except queue.Full:self.dropped+=1
 def run(self):
  try:
   with gzip.open(self.path,'wt') as f:
    while not self.stop.is_set() or not self.q.empty():
     try:x=self.q.get(timeout=.05)
     except queue.Empty:continue
     t=time.perf_counter();f.write(json.dumps(clean(x),separators=(',',':'),allow_nan=False)+'\n');self.times.append((time.perf_counter()-t)*1000);self.written+=1
  except Exception as e:self.error=repr(e)
 def close(self):
  self.stop.set();self.thread.join(10)
  if self.thread.is_alive():raise RuntimeError('RECORDER_STOP_TIMEOUT')
class Pipeline:
 def __init__(self,cfg,output,observer_factory=Observer):
  self.cfg=cfg.copy();self.output=Path(output);self.output.mkdir(parents=True,exist_ok=True);self.factory=observer_factory;self.input=LatestSlot();self.raw=None;self.analyzed_packet=None;self.state=None;self.lock=threading.Lock();self.status='STOPPED';self.events=[];self.counter=Counter();self.stop_event=threading.Event();self.source=None;self.worker=None;self.publisher=None;self.start_time=time.monotonic();self.timings={};self.samples=deque(maxlen=3600);self.record=None;self.error=None;self.stages=deque(maxlen=300);self.ui_times=deque(maxlen=6000);self.last_ui_key=None;self.transitions=[];self.metric_history=deque(maxlen=600);self.published=None;self.geometry_revision=0
  self.video_record=None;self.video_geometry=None;self.video_layers=dict(DEFAULT_LAYERS,**cfg.get('video_layers',{}));self._stopped=False
 def set_video_layers(self,layers):
  with self.lock:self.video_layers={name:bool(layers.get(name,default)) for name,default in DEFAULT_LAYERS.items()}
 def set_live_geometry(self,path):
  """Change the analysis ROI without restarting camera reception.

  A generation check discards in-flight results from the previous mask.
  Historical files remain available even when the active mask is cleared.
  """
  if self.cfg['source_mode']!='LIVE_CAMERA':raise ValueError('LIVE_ONLY')
  geometry=json.loads(Path(path).read_text()) if path else None
  with self.lock:
   self.cfg['live_geometry']=str(path) if path else None;self.geometry_revision+=1
   self.video_geometry=geometry
   self.state=None;self.analyzed_packet=None;self.metric_history.clear()
   self.published={'valid':False,'reason':'MASK_CHANGED','velocity':None,'aggregation':None,'training_ready':False,'scientific_control_admitted':False}
   event=dict(elapsed=time.monotonic()-self.start_time,status='MASK_APPLIED' if path else 'MASK_CLEARED',geometry_revision=self.geometry_revision,path=self.cfg['live_geometry'])
   self.transitions.append(event)
  dump(self.output/'geometry_changes.json',[e for e in self.transitions if 'geometry_revision' in e])
 def source_status(self,s):
  if s!=self.status:self.transitions.append({'elapsed':time.monotonic()-self.start_time,'status':s});self.status=s
 def receive(self,p):
  try:p.validate()
  except Exception as e:self.counter['invalid_packets']+=1;self.source_status('INVALID_PACKET:'+str(e));return
  self.counter['bridge_received']+=1
  with self.lock:
   self.raw=p;revision=self.geometry_revision;geometry=self.video_geometry;layers=self.video_layers
  if self.video_record:self.video_record.submit(p,revision,geometry,layers)
  self.input.put(p)
 def start(self,source_factory=None):
  if self.worker and self.worker.is_alive():raise RuntimeError('ALREADY_RUNNING')
  self.stop_event.clear();self.start_time=time.monotonic();self.input=LatestSlot();self.raw=None;self.analyzed_packet=None;self.state=None;self.counter=Counter();self.timings={};self.samples=deque(maxlen=3600);self.error=None;self.metric_history.clear();self._stopped=False
  geometry_path=self.cfg.get('live_geometry') if self.cfg['source_mode']=='LIVE_CAMERA' else self.cfg.get('replay_geometry')
  self.video_geometry=json.loads(Path(geometry_path).read_text()) if geometry_path else None
  self.record=Recorder(self.output/'states.jsonl.gz',self.cfg.get('record_queue',32))
  if self.cfg.get('record_experiment_video',True):self.video_record=ExperimentVideoRecorder(self.output,self.cfg,self.start_time)
  self.worker=threading.Thread(target=self.work,name='observation-worker',daemon=True);self.worker.start()
  self.publisher=threading.Thread(target=self.publish,name='ttl-state-publisher',daemon=True);self.publisher.start()
  if self.cfg.get('diagnostics',False):
   self.resource_thread=threading.Thread(target=self.monitor,name='resource-monitor',daemon=True);self.resource_thread.start()
  if source_factory:self.source=source_factory(self.receive,self.source_status)
  elif self.cfg['source_mode']=='REALTIME_REPLAY':self.source=RecordedVideoFrameSource(self.cfg['video'],self.receive,self.source_status,self.cfg['speed'],self.cfg['loop'])
  else:
   from live_adapter import LiveCameraFrameSource
   self.source=LiveCameraFrameSource(self.cfg,self.receive,self.source_status,self.output)
  self.source.start();return self
 def work(self):
  continuity=Continuity(self.cfg['max_pts_gap_s']);last_geometry_revision=-1
  try:
   observer=self.factory(self.cfg);piv_display=PIVDisplayFilter()
   while not self.stop_event.is_set() or self.input.peek() is not None:
    self.input.event.wait(.03);p=self.input.take()
    if p is None:continue
    try:reset,reason=continuity.check(p)
    except ValueError as e:self.counter['continuity_rejections']+=1;self.source_status('INVALID_PACKET:'+str(e));continue
    if p.age_ms()>self.cfg['ttl_ms']:self.counter['analysis_expired']+=1;continue
    with self.lock:geometry_revision=self.geometry_revision
    reset=reset or geometry_revision!=last_geometry_revision
    last_geometry_revision=geometry_revision
    t=time.perf_counter()
    if self.cfg.get('analysis_delay_ms'):time.sleep(self.cfg['analysis_delay_ms']/1000)
    r=observer.analyze(p,reset);r.update(schema='F10Y_OBSERVATION_1.0',source_mode=p.source_mode,session_id=p.session_id,frame_id=p.frame_id,source_pts=p.source_pts,arrival_monotonic=p.arrival_monotonic,age_upper_bound_at_arrival_ms=p.age_upper_bound_at_arrival_ms,timestamp_quality=p.timestamp_quality,capture_timestamp=p.capture_timestamp,host_decode_timestamp=p.host_decode_timestamp,host_receive_timestamp=p.host_receive_timestamp,reset_reason=reason,completed_monotonic=time.monotonic(),training_ready=False,scientific_control_admitted=False)
    r['timings'].update(decode_ms=p.decode_ms,complete_processing_ms=p.decode_ms+r['timings']['analysis_ms'],queue_age_ms=(t-p.arrival_monotonic)*1000,end_to_end_software_ms=p.age_ms())
    if reset:piv_display.reset()
    r['geometry_revision']=geometry_revision
    r['piv_display_vectors']=piv_display.update(r,self.cfg['ttl_ms'])
    with self.lock:
     if geometry_revision!=self.geometry_revision:
      self.counter['discarded_old_geometry']+=1;continue
     r['geometry_revision']=geometry_revision;r['geometry_path']=self.cfg.get('live_geometry') if p.source_mode=='LIVE_CAMERA' else self.cfg.get('replay_geometry')
     self.counter['processed']+=1
     if r['valid'] and p.age_ms()<=self.cfg['ttl_ms']:self.counter['fresh_valid']+=1
     if r.get('aggregation'):
      a=r['aggregation']['regions'][0];self.metric_history.append([p.session_id,p.frame_id,p.source_pts,a['count'],a['density_per_10k_px'],a['nn_px'],a['rg_px'],a['radial_px_s'] if r['valid'] else None,self.counter['fresh_valid']/max(time.monotonic()-self.start_time,.001)])
     self.state=r;self.analyzed_packet=p
    if self.video_record:self.video_record.add_analysis(r)
    if self.cfg.get('diagnostics',False):
     for k,v in r['timings'].items():self.timings.setdefault(k,deque(maxlen=6000)).append(v)
    self.record.submit(observation_record(r,self.cfg.get('diagnostics',False)))
  except Exception as e:self.error=traceback.format_exc();self.source_status('ANALYSIS_ERROR:'+repr(e))
 def snapshot(self):
  with self.lock:s=self.state;p=self.raw;analyzed_packet=self.analyzed_packet;geometry_revision=self.geometry_revision
  now=time.monotonic();age=(now-s['arrival_monotonic'])*1000+s['age_upper_bound_at_arrival_ms'] if s else None
  fresh=bool(not self.error and s and age<=self.cfg['ttl_ms'] and self.status=='STREAMING' and p and p.session_id==s['session_id'])
  valid=bool(fresh and s['valid'])
  reason='ANALYSIS_ERROR' if self.error else 'FRESH' if valid else self.status if self.status!='STREAMING' else 'STALE' if age and age>self.cfg['ttl_ms'] else s.get('reason','WARMUP') if s else 'WAITING'
  return dict(raw=p,state=s,analyzed_packet=analyzed_packet,status=reason,source_status=self.status,fresh=fresh,valid=valid,age_ms=age,ttl_ms=self.cfg['ttl_ms'],geometry_revision=geometry_revision,counters=dict(self.counter),analysis_skipped=self.input.overwritten,recorder_dropped=self.record.dropped if self.record else 0,video_recording=self.video_record.status() if self.video_record else None)
 def publish(self):
  import zmq
  context=zmq.Context();socket=context.socket(zmq.PUB);socket.setsockopt(zmq.SNDHWM,1);socket.setsockopt(zmq.LINGER,0)
  port=socket.bind_to_random_port('tcp://127.0.0.1');self.publish_endpoint=f'tcp://127.0.0.1:{port}';dump(self.output/'publisher.json',{'endpoint':self.publish_endpoint,'ttl_ms':self.cfg['ttl_ms'],'local_only':True})
  try:
   while not self.stop_event.wait(.05):
    t=time.perf_counter();snap=self.snapshot();s=snap['state'];valid=snap['valid'];agg=s.get('aggregation') if s else None
    message={'schema':'F10Y_M4_1.0','session_id':s['session_id'] if s else None,'frame_id':s['frame_id'] if s else None,'source_mode':self.cfg['source_mode'],'valid':valid,'detection_valid':snap['fresh'],'aggregation_valid':bool(valid and agg and agg['status']!='NOT_OBSERVABLE'),'state_age_ms':snap['age_ms'],'ttl_ms':self.cfg['ttl_ms'],'published_monotonic':time.monotonic(),'clock_scope':'LOCAL_HOST_ONLY','reason':snap['status'],'velocity':agg['regions'][0]['mean_velocity'] if valid and agg else None,'aggregation':agg['status'] if valid and agg else None,'six_parameters':None,'material_id':None,'training_ready':False}
    message.update(motion_reference=s.get('motion_reference') if s else None,motion_reference_verified=bool(s and s.get('motion_reference_verified')),scientific_control_admitted=False,geometry_revision=snap['geometry_revision'])
    with self.lock:
     if snap['geometry_revision']!=self.geometry_revision:continue
     self.published=message
     try:socket.send_json(clean(message),flags=zmq.DONTWAIT)
     except zmq.Again:self.counter['publisher_dropped']+=1
    if self.cfg.get('diagnostics',False):self.counter['state_published']+=1;self.timings.setdefault('publish_ms',deque(maxlen=6000)).append((time.perf_counter()-t)*1000)
  finally:
   socket.send_json({'schema':'F10Y_M4_1.0','valid':False,'velocity':None,'aggregation':None,'reason':'STOPPED'},flags=zmq.DONTWAIT);socket.close();context.term()
 def monitor(self):
  proc=psutil.Process();proc.cpu_percent()
  while not self.stop_event.wait(1):
   now=time.monotonic();row={'elapsed_s':now-self.start_time,'rss_mib':proc.memory_info().rss/2**20,'cpu_percent':proc.cpu_percent(),'processed':self.counter['processed'],'age_ms':self.snapshot()['age_ms']}
   try:
    gpu=subprocess.check_output(['nvidia-smi','--query-gpu=utilization.gpu,memory.used','--format=csv,noheader,nounits'],text=True,timeout=.8).strip().splitlines()[0].split(',');row.update(gpu_percent=float(gpu[0]),gpu_mib=float(gpu[1]))
   except Exception:row.update(gpu_percent=None,gpu_mib=None)
   self.samples.append(row)
 def event(self,state,source='SIMULATED_EVENT',channel='ALL',note='',conditions=None):
  if state not in ['ON','OFF'] or source not in ['MANUAL_ANNOTATION','DEVICE_TELEMETRY','SIMULATED_EVENT']:raise ValueError('EVENT_SCHEMA')
  p=self.raw;e=dict(event_id=uuid.uuid4().hex,channel=channel,state=state,source=source,command_time=None,actual_execution_time=None,camera_time=None,source_pts=p.source_pts if p else None,session_id=p.session_id if p else None,host_event_monotonic=time.monotonic(),time_sync_quality='HOST_ONLY_EXPOSURE_UNSYNCHRONIZED',conditions=conditions or {},note=note,hardware_control=False)
  self.events.append(e);dump(self.output/'pzt_events.json',self.events);return e
 def stop(self):
  if self._stopped:return self.summary
  shutdown_errors=[]
  def finish(name,fn):
   try:fn()
   except Exception:shutdown_errors.append(name+': '+traceback.format_exc())
  if self.source:finish('source stop',self.source.stop)
  self.stop_event.set();self.input.event.set()
  if self.worker:finish('analysis stop',lambda:self.worker.join(15))
  if self.publisher:finish('publisher stop',lambda:self.publisher.join(3))
  if hasattr(self,'resource_thread'):finish('monitor stop',lambda:self.resource_thread.join(2))
  if self.record:finish('state log finalize',self.record.close)
  if self.video_record:finish('video finalize',self.video_record.close)
  self.source_status('STOPPED');self.published={'valid':False,'reason':'STOPPED','velocity':None,'aggregation':None,'six_parameters':None}
  if self.worker and self.worker.is_alive():shutdown_errors.append('WORKER_STOP_TIMEOUT')
  if shutdown_errors:self.error=(self.error or '')+'\n'.join(shutdown_errors)
  self.summary=self.report();dump(self.output/'summary.json',self.summary);self._stopped=True;return self.summary
 def report(self):
  span=time.monotonic()-self.start_time
  result=dict(identity=identity(),config=self.cfg,wall_s=span,status=self.status,error=self.error,counters=dict(self.counter),source_decoded=getattr(self.source,'decoded',None),source_sent=getattr(self.source,'sent',0),source_sessions=getattr(self.source,'sessions',0),bridge_overwritten=getattr(self.source,'overwritten',0),analysis_skipped=self.input.overwritten,recorder_dropped=self.record.dropped if self.record else 0,recorder_error=self.record.error if self.record else None,recorder_written=self.record.written if self.record else 0,fresh_valid_hz=self.counter['fresh_valid']/span if span else 0,analysis_hz=self.counter['processed']/span if span else 0,transitions=self.transitions,events=self.events,video_recording=self.video_record.report() if self.video_record else None,threads_closed=not(self.worker and self.worker.is_alive()))
  if self.cfg.get('diagnostics',False):result.update(record_ms=stats(self.record.times if self.record else []),stages={k:stats(v) for k,v in self.timings.items()},ui_refresh_hz=self.counter['ui_ticks']/span if span else 0,ui_render_ms=stats(self.ui_times),resources=list(self.samples))
  return result
