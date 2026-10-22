"""NO CAMERA ACCESS: synthetic BGR frames for actual Windows/WSL IPC validation."""
import os,sys,time,struct,json,threading,uuid
import numpy as np
if os.name=='nt':
 import msvcrt
 msvcrt.setmode(sys.stdin.fileno(),os.O_BINARY);msvcrt.setmode(sys.stdout.fileno(),os.O_BINARY)
session=uuid.uuid4().hex;start=time.monotonic();seq=0
watch=threading.Timer(15,lambda:os._exit(124));watch.daemon=True;watch.start()
try:
 for cmd in sys.stdin.buffer:
  if cmd.strip()==b'QUIT':break
  seq+=1;now=time.monotonic();a=np.full((1024,1280,3),seq%255,np.uint8)
  h=dict(status='FRAME',source='LIVE_CAMERA',session=session,seq=seq,shape=list(a.shape),original_resolution_wh=[1280,1024],roi_xywh=[0,0,1280,1024],payload_bytes=a.nbytes,sender_frame_age_ms=(time.monotonic()-now)*1000,sender_decode_mono=now,receiver_state='streaming',received_frames=seq,synthetic=True)
  raw=json.dumps(h).encode();sys.stdout.buffer.write(b'CV2FRAME'+struct.pack('<I',len(raw))+raw+a.tobytes());sys.stdout.buffer.flush()
  if time.monotonic()-start>10:break
finally:watch.cancel()
