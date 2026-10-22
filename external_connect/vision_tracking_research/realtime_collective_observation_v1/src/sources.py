import time,threading,uuid
import av
from protocol import FramePacket
class RecordedVideoFrameSource:
 def __init__(self,path,emit,status,speed=1,loop=False):
  if speed not in [.5,1,2]:raise ValueError('SPEED')
  self.path=path;self.emit=emit;self.status=status;self.speed=speed;self.loop=loop;self.stop_event=threading.Event();self.paused=threading.Event();self.decoded=0;self.sent=0;self.sessions=0;self.thread=None
 def start(self):self.thread=threading.Thread(target=self.run,name='frame-acquisition',daemon=True);self.thread.start();return self
 def run(self):
  try:
   while not self.stop_event.is_set():
    session=uuid.uuid4().hex;self.sessions+=1;self.status('CONNECTING');base=time.monotonic();first=None
    with av.open(self.path) as c:
     c.streams.video[0].codec_context.thread_count=2
     it=iter(c.decode(video=0));i=0
     while not self.stop_event.is_set():
      t=time.perf_counter()
      try:f=next(it)
      except StopIteration:break
      image=f.to_ndarray(format='gray');dms=(time.perf_counter()-t)*1000;self.decoded+=1;decoded_at=time.monotonic()
      pts=float(f.pts*f.time_base);first=pts if first is None else first;due=base+(pts-first)/self.speed
      while not self.stop_event.is_set():
       if self.paused.is_set():
        start=time.monotonic();self.status('PAUSED')
        while self.paused.is_set() and not self.stop_event.wait(.02):pass
        delta=time.monotonic()-start;base+=delta;due+=delta
       left=due-time.monotonic()
       if left<=0:break
       self.stop_event.wait(min(.01,left))
      if self.stop_event.is_set():break
      now=time.monotonic();self.status('STREAMING')
      p=FramePacket('REALTIME_REPLAY',str(self.path),session,i,pts,image,(image.shape[1],image.shape[0]),host_decode_timestamp={'clock':'local_monotonic','value':decoded_at,'meaning':'replay decode completed before scheduled release'},host_receive_timestamp={'clock':'local_monotonic','value':now},decode_ms=dms)
      self.emit(p);self.sent+=1;i+=1
    if not self.loop:self.status('EOF');break
    self.status('LOOP_RESET')
  except Exception as e:self.status('ERROR:'+repr(e))
  finally:
   if self.stop_event.is_set():self.status('STOPPED')
 def pause(self):self.paused.set();self.status('PAUSED')
 def resume(self):self.paused.clear()
 def stop(self):
  self.stop_event.set();self.paused.clear()
  if self.thread:self.thread.join(5)
  if self.thread and self.thread.is_alive():raise RuntimeError('SOURCE_STOP_TIMEOUT')
