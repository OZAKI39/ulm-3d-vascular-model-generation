"""Configurable launch + finite session supervision; original decoder/protocol loop reused."""
import sys,os,json,time,threading,subprocess,uuid
from collections import deque
from pathlib import Path
from ycommon import ROOT,dump
sys.path.insert(0,str(ROOT.parent/'smart_online_tracking_v2/src'))
from live_bridge import WindowsLatestBridge,read_exact,winpath
from protocol import FramePacket
class ConfiguredWindowsBridge(WindowsLatestBridge):
 def __init__(self,output,project,python,seconds=120):
  super().__init__(output,seconds);self.project=Path(project);self.python=python;self.records=deque(maxlen=20000)
 def start(self):
  if not 1<=self.seconds<=180:raise ValueError('FINITE_SESSION_REQUIRED')
  cfg=json.loads((self.project/'chronos_config.json').read_text())
  sys.path.insert(0,str(ROOT.parents[1]/'CHRONOS_control'))
  from chronos_common import video_options
  env=os.environ.copy();env['OPENCV_FFMPEG_CAPTURE_OPTIONS']=video_options(cfg);env['PYTHONDONTWRITEBYTECODE']='1'
  tokens=env.get('WSLENV','').split(':')
  for name in ['OPENCV_FFMPEG_CAPTURE_OPTIONS','PYTHONDONTWRITEBYTECODE']:
   if name not in tokens:tokens.append(name)
  env['WSLENV']=':'.join(tokens)
  server=ROOT.parent/'smart_online_tracking_v2/scripts/windows_frame_server.py'
  args=[self.python,'-B',winpath(server),'--source',winpath(self.project),'--summary',winpath(self.output/'windows_receiver.json'),'--seconds',str(self.seconds)]
  self.log=(self.output/'windows_stderr.log').open('wb');self.process=subprocess.Popen(args,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=self.log,env=env,bufsize=0)
  dump(self.output/'launch.json',{'args':args,'reused_loop':'WindowsLatestBridge._loop','original_receiver':'LatestFrameStream','finite_seconds':self.seconds})
  self.thread=threading.Thread(target=self._loop,daemon=True);self.thread.start();return self
 def stop(self):
  self.stop_event.set()
  if self.process:
   try:self.process.stdin.write(b'QUIT\n');self.process.stdin.flush();self.process.wait(timeout=9)
   except (OSError,subprocess.TimeoutExpired):
    self.process.terminate()
    try:self.process.wait(timeout=3)
    except subprocess.TimeoutExpired:self.process.kill();self.process.wait(timeout=2)
   self.process.stdin.close();self.process.stdout.close()
  if self.thread:self.thread.join(2)
  if self.log:self.log.close()
  dump(self.output/'bridge_status.json',{'error':self.error,'exit_code':self.process.returncode if self.process else None,'messages':list(self.records),'clock':'sender age + GET RTT + local age'})
  if self.thread and self.thread.is_alive():raise RuntimeError('BRIDGE_THREAD_STOP_TIMEOUT')

def from_bridge(h,a):
 return FramePacket('LIVE_CAMERA','Chronos original LatestFrameStream',h['session'],int(h['seq']),None,a,tuple(h['original_resolution_wh']),'BGR8',host_decode_timestamp={'clock':'windows_monotonic','value':h['sender_decode_mono'],'utc':h.get('decode_utc'),'meaning':'decode, not exposure'},host_receive_timestamp={'clock':'wsl_monotonic','value':h['arrival_mono']},arrival_monotonic=h['arrival_mono'],roi=h['roi_xywh'],timestamp_quality='EXPOSURE_UNKNOWN_AGE_UPPER_BOUND',age_upper_bound_at_arrival_ms=h['age_upper_bound_at_arrival_ms']).validate()
class LiveCameraFrameSource:
 def __init__(self,config,emit,status,output,bridge_factory=ConfiguredWindowsBridge):
  self.cfg=config;self.emit=emit;self.status=status;self.output=Path(output);self.factory=bridge_factory;self.stop_event=threading.Event();self.thread=None;self.bridge=None;self.decoded=None;self.sent=0;self.sessions=0;self.received=0;self.overwritten=0
 def start(self):self.thread=threading.Thread(target=self.run,name='live-session-supervisor',daemon=True);self.thread.start();return self
 def run(self):
  try:
   while not self.stop_event.is_set():
    self.sessions+=1;self.status('CONNECTING');last=None;started=time.monotonic()
    self.bridge=self.factory(self.output/f'bridge_{uuid.uuid4().hex}',self.cfg['windows_project'],self.cfg['windows_python'],self.cfg['bridge_session_s'])
    try:
     self.bridge.start()
     while not self.stop_event.wait(.002):
      if self.bridge.done.is_set() or time.monotonic()-started>=self.cfg['bridge_session_s']-1:break
      snap=self.bridge.snapshot()
      if snap is None:self.status('DISCONNECTED');continue
      h,a=snap;key=(h['session'],h['seq'])
      if key==last:continue
      if last and key[0]==last[0]:self.overwritten+=max(0,key[1]-last[1]-1)
      last=key;self.received+=1;self.decoded=h.get('received_frames');self.status('STREAMING');self.emit(from_bridge(h,a));self.sent+=1
    finally:self.bridge.stop()
    self.status('RECONNECTING');self.stop_event.wait(.5)
  except Exception as e:self.status('ERROR:'+repr(e)) # No replay fallback.
  finally:
   if self.stop_event.is_set():self.status('STOPPED')
 def pause(self):raise ValueError('LIVE_HAS_STOP_NOT_PAUSE')
 def resume(self):pass
 def stop(self):
  self.stop_event.set()
  if self.thread:self.thread.join(16)
  if self.thread and self.thread.is_alive():raise RuntimeError('LIVE_STOP_TIMEOUT')
