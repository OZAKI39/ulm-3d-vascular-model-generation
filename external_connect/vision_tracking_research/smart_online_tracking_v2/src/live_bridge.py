"""Bounded OS-pipe bridge: one GET in flight and one latest frame in WSL.

Freshness uses a conservative age upper bound, avoiding subtraction of Windows
and Linux monotonic clocks: sender frame age + Linux GET roundtrip + local age.
"""
from pathlib import Path
import json,os,struct,subprocess,threading,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1]

def winpath(p):return subprocess.check_output(['wslpath','-w',str(p)],text=True).strip()

def read_exact(stream,n):
    chunks=[];left=n
    while left:
        block=stream.read(left)
        if not block:raise EOFError(f'IPC closed with {left} bytes missing')
        chunks.append(block);left-=len(block)
    return b''.join(chunks)

class WindowsLatestBridge:
    def __init__(self,output,seconds=90,roi=None):
        self.output=Path(output).resolve();self.output.mkdir(parents=True,exist_ok=True)
        self.seconds=seconds;self.roi=roi;self.latest=None;self.lock=threading.Lock()
        self.stop_event=threading.Event();self.done=threading.Event();self.records=[];self.error=None
        self.process=None;self.thread=None;self.log=None
    def start(self):
        cfg=json.loads(Path('/mnt/f/CHRONOS_control/chronos_config.json').read_text())
        import sys
        sys.path.insert(0,str(ROOT.parent.parent/'CHRONOS_control'))
        from chronos_common import video_options
        env=os.environ.copy();env['OPENCV_FFMPEG_CAPTURE_OPTIONS']=video_options(cfg)
        env['PYTHONDONTWRITEBYTECODE']='1'
        # WSL does not automatically forward arbitrary environment variables to
        # Windows. WSLENV explicitly passes these two values at process creation.
        tokens=[s for s in env.get('WSLENV','').split(':') if s]
        for name in ['OPENCV_FFMPEG_CAPTURE_OPTIONS','PYTHONDONTWRITEBYTECODE']:
            if name not in [v.split('/')[0] for v in tokens]:tokens.append(name)
        env['WSLENV']=':'.join(tokens)
        args=['/mnt/f/CHRONOS_control/.venv/Scripts/python.exe','-B',winpath(ROOT/'scripts/windows_frame_server.py'),
              '--summary',winpath(self.output/'windows_receiver.json'),'--seconds',str(self.seconds)]
        if self.roi is not None:args+=['--roi',*map(str,self.roi)]
        self.log=(self.output/'windows_stderr.log').open('wb')
        self.process=subprocess.Popen(args,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=self.log,env=env,bufsize=0)
        (self.output/'bridge_launch.json').write_text(json.dumps({'args':args,'capture_options':env['OPENCV_FFMPEG_CAPTURE_OPTIONS'],
           'protocol':'one request/response in flight; raw BGR; no resizing; single latest slot'},indent=2))
        self.thread=threading.Thread(target=self._loop,daemon=True);self.thread.start();return self
    def _loop(self):
        try:
            while not self.stop_event.is_set():
                request=time.monotonic();self.process.stdin.write(b'GET\n');self.process.stdin.flush()
                prefix=read_exact(self.process.stdout,12)
                if prefix[:8]!=b'CV2FRAME':raise ValueError('Invalid binary IPC magic')
                n=struct.unpack('<I',prefix[8:])[0]
                if not 0<n<=65536:raise ValueError('Invalid header size')
                header=json.loads(read_exact(self.process.stdout,n))
                size=header['payload_bytes']
                if not 0<=size<=16_000_000:raise ValueError('Invalid image payload size')
                payload=read_exact(self.process.stdout,size) if size else None
                arrived=time.monotonic();header['arrival_mono']=arrived;header['get_roundtrip_ms']=(arrived-request)*1000
                if header['status']=='FRAME':
                    shape=header['shape']
                    if len(shape)!=3 or shape[2]!=3 or int(np.prod(shape))!=size:raise ValueError('Bad frame shape')
                    header['age_upper_bound_at_arrival_ms']=header['sender_frame_age_ms']+header['get_roundtrip_ms']
                    frame=np.frombuffer(payload,np.uint8).reshape(shape)
                    with self.lock:self.latest=(header,frame)
                elif header['receiver_state']!='streaming':
                    with self.lock:self.latest=None
                self.records.append(header)
                if header['status']!='FRAME':time.sleep(.003)
        except (OSError,EOFError,ValueError) as exc:
            if not self.stop_event.is_set():self.error=repr(exc)
        finally:
            with self.lock:self.latest=None
            self.done.set()
    def snapshot(self):
        with self.lock:return self.latest
    def stop(self):
        self.stop_event.set()
        if self.thread:self.thread.join(8)
        if self.process:
            try:
                self.process.stdin.write(b'QUIT\n');self.process.stdin.flush()
                self.process.wait(timeout=10)
            except (OSError,subprocess.TimeoutExpired):
                self.process.terminate();self.process.wait(timeout=5)
            self.process.stdin.close();self.process.stdout.close()
        if self.log:self.log.close()
        (self.output/'bridge_frames.jsonl').write_text(''.join(json.dumps(r)+'\n' for r in self.records))
        (self.output/'bridge_status.json').write_text(json.dumps({'error':self.error,'exit_code':self.process.returncode,
             'frames':sum(r['status']=='FRAME' for r in self.records),'single_slot':True,
             'clock_method':'sender age + GET roundtrip conservative upper bound, no cross-clock subtraction'},indent=2))
