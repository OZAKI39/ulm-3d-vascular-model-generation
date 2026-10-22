"""Hardware-free checks using real image bytes and original receiver control flow."""
import sys,time,json,hashlib,threading
from pathlib import Path
from unittest.mock import patch
ROOT=Path(__file__).resolve().parents[1];sys.dont_write_bytecode=True
from replay_benchmark import MockController,LatestFrameStream,load_config,dump
import cv2
results=[]
def check(name,condition,evidence):
 results.append({'name':name,'passed':bool(condition),'evidence':evidence});assert condition,name
m=MockController()
for fid,age,disc,valid in [(1,0,False,True),(2,.3,False,True),(2,0,False,True),(3,0,True,True),(4,0,False,False)]:
 obj={'frame_id':fid,'receive_monotonic_s':time.monotonic()-age,'timestamp_discontinuity':disc,'targets':[{'visual_valid':valid}]};sent=m.send(obj)
 deadline=time.monotonic()+.3
 while len(m.rows)<len(results)+1 and time.monotonic()<deadline:time.sleep(.002)
 # Evaluation performed after all messages below.
 time.sleep(.03)
m.close();observed=[r['valid_visual'] for r in m.rows];check('mock_rejects_stale_duplicate_discontinuous_unreliable',observed==[True,False,False,False,False],m.rows)
cfg=load_config();cfg.update(rtsp_control_relay=False,reconnect_seconds=.05,open_timeout_ms=100,read_timeout_ms=100)
s=LatestFrameStream(cfg);im=cv2.imread(str(ROOT/'evidence/frame_00120.png'));calls=[];failure=threading.Event();release=threading.Event()
class FakeCapture:
 def __init__(self,*args):self.number=len(calls);calls.append(self);self.n=0
 def isOpened(self):return True
 def get(self,key):return 30 if key==cv2.CAP_PROP_FPS else 1
 def read(self):
  self.n+=1;time.sleep(.01)
  if self.number==0 and self.n>2:failure.set();return False,None
  if self.number>0:release.wait(.2)
  return True,im.copy()
 def release(self):pass
with patch.object(cv2,'VideoCapture',FakeCapture):
 s.start();deadline=time.monotonic()+2
 while not failure.is_set() and time.monotonic()<deadline:time.sleep(.002)
 deadline=time.monotonic()+1
 while s.failures<1 and time.monotonic()<deadline:time.sleep(.002)
 check('original_receiver_clears_slot_after_read_failure',s.snapshot() is None,{'failures':s.failures,'state':s.state})
 release.set();deadline=time.monotonic()+2
 while s.frames<4 and time.monotonic()<deadline:time.sleep(.002)
 check('original_receiver_reconnects',s.connections>=2 and s.frames>=4,s.stats())
 s.stop()
check('original_receiver_stops_cleanly',not s.thread.is_alive(),{'state':s.state})
checks=[]
for row in json.loads((ROOT/'evidence/protected_sources_before.json').read_text()):
 p=ROOT.parent/row['path'];sha=hashlib.sha256(p.read_bytes()).hexdigest();checks.append({'path':row['path'],'unchanged':sha==row['sha256'],'sha256':sha})
check('protected_source_and_config_unchanged',all(r['unchanged'] for r in checks),{'files':len(checks)})
dump(ROOT/'evidence/protected_sources_after.json',checks);dump(ROOT/'evidence/interface_checks.json',results);print([(r['name'],r['passed']) for r in results])
