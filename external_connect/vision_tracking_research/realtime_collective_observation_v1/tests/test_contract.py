import sys,time,io,json,struct,dataclasses,threading
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np,pytest
from protocol import FramePacket,Continuity,LatestSlot,read_cv2frame
from live_adapter import from_bridge,ConfiguredWindowsBridge,WindowsLatestBridge
from aggregation import Aggregation
from pipeline import Pipeline

def packet(fid=0,pts=0,session='a',shape=(96,128)):
 return FramePacket('REALTIME_REPLAY','test',session,fid,pts,np.zeros(shape,np.uint8),(shape[1],shape[0]),'BGR8' if len(shape)==3 else 'GRAY8')
@pytest.mark.parametrize('shape',[(96,128),(1024,1280),(1120,1476),(96,128,3)])
def test_frame_shapes(shape):assert packet(shape=shape).validate()
@pytest.mark.parametrize('change',[{'frame_id':-1},{'source_pts':-1},{'source_pts':float('nan')},{'image':np.empty((0,0),np.uint8)},{'original_resolution':(12,12)},{'pixel_format':'BGR8'},{'source_mode':'FAKE'},{'image':np.zeros((96,128,4),np.uint8)}])
def test_bad_packet(change):
 with pytest.raises(ValueError):dataclasses.replace(packet(),**change).validate()
def test_continuity():
 c=Continuity();assert c.check(packet())[0]
 with pytest.raises(ValueError):c.check(packet())
 assert c.check(packet(1,.03))[0] is False
 assert c.check(packet(2,2))[1]=='PTS_DISCONTINUITY'
 assert c.check(packet(0,0,'b'))[1]=='SESSION_RESET'
 with pytest.raises(ValueError):c.check(packet(1,.1,'b',(128,96)))
def binary(header,payload):
 h=json.dumps(header).encode();return b'CV2FRAME'+struct.pack('<I',len(h))+h+payload
class Chunked(io.BytesIO):
 def read(self,n=-1):return super().read(min(n,7))
def test_binary_chunking_and_corruption():
 a=np.zeros((10,10,3),np.uint8);h=dict(status='FRAME',shape=list(a.shape),payload_bytes=a.nbytes)
 assert read_cv2frame(Chunked(binary(h,a.tobytes())))[1].shape==a.shape
 for b in [b'BADMAGIC'+struct.pack('<I',0),b'CV2FRAME'+struct.pack('<I',999999),binary(h,b'a'),binary(dict(h,payload_bytes=2),b'ab')]:
  with pytest.raises((ValueError,EOFError)):read_cv2frame(Chunked(b))
def test_original_bridge_loop_reused():assert ConfiguredWindowsBridge._loop is WindowsLatestBridge._loop

def test_clock_no_cross_host_subtraction():
 h=dict(session='s',seq=1,original_resolution_wh=[10,10],sender_decode_mono=900000,arrival_mono=time.monotonic(),roi_xywh=[0,0,10,10],age_upper_bound_at_arrival_ms=12)
 p=from_bridge(h,np.zeros((10,10,3),np.uint8));assert 12<=p.age_ms()<30;assert p.capture_timestamp is None

def test_latest_slot():
 s=LatestSlot()
 for i in range(100):s.put(i)
 assert s.take()==99 and s.overwritten==99 and s.take() is None
class DummyObserver:
 def __init__(self,cfg):pass
 def analyze(self,p,reset=False):return dict(valid=True,reason='OK',candidates=[],vectors=[],aggregation=None,timings={'analysis_ms':1},geometry=None)
class MockSource:
 def __init__(self,emit,status):self.emit=emit;self.status=status;self.decoded=0;self.sent=0;self.sessions=1
 def start(self):self.status('STREAMING');self.emit(packet());self.sent+=1;return self
 def stop(self):self.status('STOPPED')
def config():return dict(source_mode='REALTIME_REPLAY',max_pts_gap_s=.25,ttl_ms=50,record_queue=8,analysis_delay_ms=0,record_experiment_video=False)
def test_pipeline_fault_states(tmp_path):
 p=Pipeline(config(),tmp_path,DummyObserver).start(MockSource);time.sleep(.025);assert p.snapshot()['valid']
 time.sleep(.09);assert not p.snapshot()['valid'];assert p.published['velocity'] is None
 p.receive(packet(1,.1));time.sleep(.02);p.source_status('DISCONNECTED');assert not p.snapshot()['valid']
 p.source_status('STREAMING');p.receive(packet(0,0,'new'));time.sleep(.02);assert p.snapshot()['state']['reset_reason']=='SESSION_RESET'
 p.event('ON');assert p.events[0]['source']=='SIMULATED_EVENT' and not p.events[0]['hardware_control']
 p.stop();assert not p.published['valid'];assert not p.worker.is_alive() and not p.publisher.is_alive()
def test_slow_bounded(tmp_path):
 c=config();c['analysis_delay_ms']=80;p=Pipeline(c,tmp_path,DummyObserver).start(MockSource)
 for i in range(1,100):p.receive(packet(i,i*.001))
 time.sleep(.15);p.stop();assert p.input.overwritten>90 and not p.snapshot()['valid']
def test_source_switch(tmp_path):
 p=Pipeline(config(),tmp_path/'r',DummyObserver).start(MockSource);time.sleep(.02);p.stop()
 c=config();c['source_mode']='LIVE_CAMERA';p=Pipeline(c,tmp_path/'l',DummyObserver).start(MockSource);time.sleep(.02);p.stop();assert p.cfg['source_mode']=='LIVE_CAMERA'
def test_motion_cases():
 rng=np.random.default_rng(87);base=rng.uniform(120,380,(50,2));center=base.mean(0)
 for mode in ['converge','translation','rotation','disappear','defocus']:
  a=Aggregation()
  for t in np.arange(0,1.51,.05):
   q=base.copy();vel=np.zeros_like(q);focus=100
   if mode in ['converge','defocus']:q=center+(base-center)*np.exp(-.2*t);vel=-.2*(q-center)
   if mode=='translation':q+=np.array([20,5])*t;vel[:]=[20,5]
   if mode=='rotation':theta=.3*t;rot=np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]]);q=center+(base-center)@rot.T;vel=(q-center)@np.array([[0,.3],[-.3,0]])
   if mode=='disappear':q=q[:max(8,50-int(t*20))];vel=vel[:len(q)]
   if mode=='defocus':focus=100*np.exp(-t)
   r=a.update(t,q,q,vel,(512,512),focus,100,cells=False)
  status=r['status']
  assert (status=='VISIBLE_CONVERGENCE')==(mode=='converge'),(mode,status)

def test_zmq_publisher_invalidates(tmp_path):
 import zmq
 c=config();c['ttl_ms']=250;p=Pipeline(c,tmp_path,DummyObserver).start(MockSource)
 while not hasattr(p,'publish_endpoint'):time.sleep(.005)
 ctx=zmq.Context();sock=ctx.socket(zmq.SUB);sock.setsockopt(zmq.SUBSCRIBE,b'');sock.setsockopt(zmq.RCVTIMEO,2000);sock.connect(p.publish_endpoint)
 time.sleep(.1);p.receive(packet(1,.1));msg=sock.recv_json();assert msg['schema']=='F10Y_M4_1.0'
 p.source_status('DISCONNECTED');time.sleep(.1)
 while sock.poll(0):msg=sock.recv_json()
 assert msg['valid'] is False and msg['velocity'] is None
 p.stop();sock.close();ctx.term()

def test_unverified_live_and_wrong_recording_geometry():
 import json
 from observer import Observer
 from ycommon import ROOT
 cfg=json.loads((ROOT/'configs/default.json').read_text());cfg['backend']='B1';o=Observer(cfg)
 p=FramePacket('LIVE_CAMERA','mock','a',0,None,np.zeros((1024,1280,3),np.uint8),(1280,1024),'BGR8');assert o.geometry_for(p) is None
 p=FramePacket('REALTIME_REPLAY','wrong source','a',0,0,np.zeros((1120,1476),np.uint8),(1476,1120));assert o.geometry_for(p) is None

def test_event_windows_no_acoustic_claim():
 sys.path.append(str(Path(__file__).resolve().parents[1]/'scripts'))
 from event_window_analysis import compare
 event=dict(source='SIMULATED_EVENT',source_pts=1.,session_id='s')
 states=[dict(session_id='s',source_pts=i*.1,valid=True,aggregation={'regions':[dict(count=10,nn_px=4,rg_px=8,radial_px_s=-1,density_per_10k_px=20)]}) for i in range(20)]
 r=compare([event],states)[0];assert r['before']['valid_samples']==10 and r['after']['valid_samples']==10;assert r['acoustic_attribution']=='NOT_VERIFIED'
 assert compare([dict(event,source_pts=None)],states)[0]['status']=='NOT_OBSERVABLE'
