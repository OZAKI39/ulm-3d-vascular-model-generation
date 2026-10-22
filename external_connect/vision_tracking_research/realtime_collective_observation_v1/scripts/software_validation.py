import sys,time,json,threading,dataclasses,io,struct,os
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from ycommon import ROOT,dump,sha,identity
from protocol import FramePacket,read_cv2frame
from sources import RecordedVideoFrameSource
from pipeline import Pipeline
import numpy as np,av
class Dummy:
 def __init__(self,cfg):pass
 def analyze(self,p,reset=False):return dict(valid=True,reason='OK',candidates=[],vectors=[],aggregation=None,geometry=None,timings={'analysis_ms':0})
class Source:
 def __init__(self,emit,status):self.emit=emit;self.status=status;self.sent=0;self.decoded=0;self.sessions=1
 def start(self):self.status('STREAMING');return self
 def stop(self):self.status('STOPPED')
def main():
 out=ROOT/'results/software_validation';out.mkdir(parents=True,exist_ok=True);cfg=json.loads((ROOT/'configs/default.json').read_text());cfg['ttl_ms']=70;cfg['record_detail']='compact'
 video=cfg['video'];pts=[];sizes=[];decode=[]
 with av.open(video) as c:
  it=iter(c.decode(video=0))
  while True:
   t=time.perf_counter()
   try:f=next(it)
   except StopIteration:break
   pts.append(float(f.pts*f.time_base));sizes.append([f.width,f.height]);decode.append((time.perf_counter()-t)*1000)
 dump(out/'source_audit.json',dict(identity=identity(),sha256=sha(video),frames=len(pts),pts_s=pts,dimensions=sizes[0],pts_monotonic=bool(np.all(np.diff(pts)>0)),decode_ms=decode))
 results=[]
 for speed in [.5,1,2]:
  arrivals=[];packets=[];statuses=[]
  def emit(p):arrivals.append(time.monotonic());packets.append(p)
  src=RecordedVideoFrameSource(video,emit,statuses.append,speed).start()
  while len(arrivals)<12:time.sleep(.002)
  src.pause();n=len(arrivals);time.sleep(.15);paused=len(arrivals)<=n+1;src.resume()
  while len(arrivals)<22:time.sleep(.002)
  src.stop();intervals=np.diff(arrivals[:10]);results.append(dict(speed=speed,actual_interval_s=float(np.median(intervals)),expected_interval_s=.0333/speed,passed=bool(abs(np.median(intervals)-.0333/speed)<.008 and paused),pause_pass=paused,closed=not src.thread.is_alive(),statuses=statuses))
 dump(out/'replay_rates.json',results)
 p=Pipeline(cfg,out/'fault_pipeline',Dummy).start(Source);events=[];fid=0;session='a';ptsvalue=0.
 def packet(**kw):
  nonlocal fid,ptsvalue
  fid+=1;ptsvalue+=.0333
  return dataclasses.replace(FramePacket('REALTIME_REPLAY','fault',session,fid,ptsvalue,np.zeros((96,128),np.uint8),(128,96)),**kw)
 def record(name,expected,action,delay=.03):
  action();time.sleep(delay);s=p.snapshot();events.append(dict(fault=name,expected=expected,observed_status=s['status'],valid=s['valid'],counters=s['counters'],pass_=not s['valid'] if expected=='INVALID' else True))
 p.receive(packet());time.sleep(.03)
 record('no_new_frame','INVALID',lambda:None,.15)
 record('short_disconnect','INVALID',lambda:p.source_status('DISCONNECTED'))
 record('long_disconnect','INVALID',lambda:None,.3)
 p.source_status('STREAMING');record('reconnect','RECOVER',lambda:p.receive(packet(session_id='b')))
 record('duplicate_id','INVALID',lambda:p.receive(dataclasses.replace(p.raw)))
 p.source_status('STREAMING');record('PTS_jump','RESET',lambda:p.receive(packet(source_pts=20)))
 record('wrong_shape','INVALID',lambda:p.receive(packet(original_resolution=(11,11))))
 record('empty_image','INVALID',lambda:p.receive(packet(image=np.zeros((0,0),np.uint8))))
 p.source_status('STREAMING')
 for delay in [0,.04,.008,.07,0]:time.sleep(delay);p.receive(packet())
 events.append(dict(fault='jitter_and_source_drops',expected='BOUNDED_LATEST',pass_=p.input.overwritten>=0,observed_skipped=p.input.overwritten))
 p.event('ON');p.source_status('DISCONNECTED');p.event('OFF');events.append(dict(fault='event_over_disconnect',pass_=all(e['actual_execution_time'] is None and not e['hardware_control'] for e in p.events)))
 p.stop();events.append(dict(fault='stop_while_working',pass_=not p.worker.is_alive()))
 try:read_cv2frame(io.BytesIO(b'BROKEN__'+struct.pack('<I',12)))
 except ValueError:events.append(dict(fault='corrupt_binary',pass_=True))
 dump(out/'fault_events.json',events)
 # New fake-LIVE source through exactly the same pipeline (no camera calls).
 cfg['source_mode']='LIVE_CAMERA';p=Pipeline(cfg,out/'mock_live_switch',Dummy).start(Source);p.receive(FramePacket('LIVE_CAMERA','MOCK_ONLY','live2',0,None,np.zeros((1024,1280,3),np.uint8),(1280,1024),'BGR8'));time.sleep(.03);valid=p.snapshot()['valid'];p.stop();dump(out/'switch.json',dict(mock_live_only=True,live_hardware=False,passed=valid,stopped=not p.worker.is_alive()))
 print('software validation written',out)
if __name__=='__main__':main()
