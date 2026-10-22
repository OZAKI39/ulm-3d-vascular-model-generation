"""Exercise inherited binary loop with real subprocess pipes; never import/call camera."""
import sys,time,subprocess,threading,json,os
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from ycommon import ROOT,dump
from live_adapter import WindowsLatestBridge,ConfiguredWindowsBridge,LiveCameraFrameSource
class MockBridge(ConfiguredWindowsBridge):
 def start(self):
  self.log=(self.output/'mock_stderr.log').open('wb');self.process=subprocess.Popen([sys.executable,__file__,'--server'],stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=self.log,bufsize=0);self.thread=threading.Thread(target=self._loop,daemon=True);self.thread.start();return self
if '--server' in sys.argv:
 import struct,numpy as np,uuid
 seq=0;session=uuid.uuid4().hex
 for command in sys.stdin.buffer:
  if command.strip()==b'QUIT':break
  seq+=1
  if seq>7:break # finite child EOF tests automatic supervisor rotation
  a=np.zeros((96,128,3),np.uint8);h=dict(status='FRAME',session=session,seq=seq,shape=list(a.shape),original_resolution_wh=[128,96],roi_xywh=[0,0,128,96],payload_bytes=a.nbytes,sender_frame_age_ms=2.,sender_decode_mono=123456.,receiver_state='streaming',received_frames=seq)
  raw=json.dumps(h).encode();sys.stdout.buffer.write(b'CV2FRAME'+struct.pack('<I',len(raw))+raw+a.tobytes());sys.stdout.buffer.flush();time.sleep(.01)
else:
 out=ROOT/'results/bridge_mock';out.mkdir(parents=True,exist_ok=True);packets=[];statuses=[];cfg=dict(windows_project='MOCK',windows_python=sys.executable,bridge_session_s=3)
 src=LiveCameraFrameSource(cfg,packets.append,statuses.append,out,bridge_factory=MockBridge).start();time.sleep(1.9);src.stop()
 dump(out/'summary.json',dict(real_os_pipes=True,inherited_binary_loop=True,mock_only=True,hardware_attempts=0,received=len(packets),sessions=len(set(p.session_id for p in packets)),max_age_upper_bound_ms=max(p.age_upper_bound_at_arrival_ms for p in packets),exposure_unknown=all(p.capture_timestamp is None for p in packets),closed=not src.thread.is_alive(),statuses=statuses,pass_=len(packets)>0 and len(set(p.session_id for p in packets))>=2 and not src.thread.is_alive()))
