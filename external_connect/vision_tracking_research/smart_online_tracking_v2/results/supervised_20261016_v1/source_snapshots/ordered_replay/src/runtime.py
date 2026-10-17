"""Bounded recording and mock-only status transport; never imports device control."""
import time,json,threading,queue,subprocess,hashlib
from pathlib import Path
import cv2,zmq,psutil
from vision_engine import stats

def identity(root):
    return {'git_head':subprocess.check_output(['git','-C',str(root),'rev-parse','HEAD'],text=True).strip(),
            'git_status':subprocess.check_output(['git','-C',str(root),'status','--short'],text=True),
            'source_sha256':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest()
                             for d in ['src','scripts'] for p in sorted((root/d).glob('*.py'))},
            'protocol_sha256':hashlib.sha256((root/'protocol.json').read_bytes()).hexdigest()}

def enforce_freshness(msg,now=None,limit_ms=100):
    now=time.monotonic() if now is None else now
    age=msg.get('age_upper_bound_at_arrival_ms',float('inf'))+1000*(now-msg.get('arrival_mono',now))
    fresh=0<=age<=limit_ms and msg.get('image_available',True)
    msg['age_upper_bound_ms']=age if age!=float('inf') else None;msg['fresh']=fresh
    if not fresh:
        # Explicitly remove every location; diagnostic raw output is stored separately.
        for key in ['candidates','targets','unreviewed_region_shapes','apparent_incorporation_events']:msg[key]=[]
        msg.update(flow=None,group_shape=None,group_neighborhood=None,visual_candidate_available=False,
                   status='NO_VALID_TARGET',reason='STALE_OR_NO_FRAME')
    msg.update(biological_valid=False,control_permission=False,group_shape=None)
    return fresh

def invalid_state(session='NO_SESSION',frame_id=-1,reason='NO_IMAGE'):
    return {'session':session,'frame_id':frame_id,'arrival_mono':time.monotonic(),
            'age_upper_bound_at_arrival_ms':1000,'image_available':False,'status':'NO_VALID_TARGET',
            'reason':reason,'candidates':[],'targets':[],'unreviewed_region_shapes':[],
            'group_shape':None,'group_neighborhood':None,'flow':None,
            'visual_candidate_available':False,'biological_valid':False,'control_permission':False}

class MockConsumer:
    """One-message ZeroMQ queues. No connection to physical controllers."""
    def __init__(self):
        self.context=zmq.Context();self.rows=[];self.ready=threading.Event();self.stop=threading.Event()
        self.thread=threading.Thread(target=self.loop,daemon=True);self.thread.start();self.ready.wait(3)
        self.tx=self.context.socket(zmq.PUSH);self.tx.setsockopt(zmq.SNDHWM,1);self.tx.setsockopt(zmq.LINGER,0)
        self.tx.connect(self.endpoint);self.send_drops=0
    def loop(self):
        rx=self.context.socket(zmq.PULL);rx.setsockopt(zmq.RCVHWM,1);rx.setsockopt(zmq.CONFLATE,1)
        port=rx.bind_to_random_port('tcp://127.0.0.1');self.endpoint=f'tcp://127.0.0.1:{port}';self.ready.set()
        previous={}
        while not self.stop.is_set():
            if not rx.poll(20):continue
            start=time.monotonic();msg=rx.recv_json();now=time.monotonic()
            fresh=enforce_freshness(msg,now);session=msg.get('session','NONE');seq=msg.get('frame_id',-1)
            new=seq>previous.get(session,-1)
            if seq>=0:previous[session]=max(seq,previous.get(session,-1))
            candidate=bool(fresh and new and msg.get('visual_candidate_available',False))
            self.rows.append({'frame_id':seq,'session':session,'fresh':fresh,'new_frame':new,
                 'candidate_valid':candidate,'bubble_valid':False,'control_permission':False,
                 'status':'NO_VALID_TARGET','reason':'EXPERT_GROUND_TRUTH_PENDING' if candidate else msg.get('reason'),
                 'age_upper_bound_ms':msg['age_upper_bound_ms'],'consume_ms':(time.monotonic()-start)*1000,
                 'transport_ms':(now-msg['send_mono'])*1000,'received_mono':now})
        rx.close(0)
    def send(self,msg):
        msg['send_mono']=time.monotonic();start=time.monotonic()
        try:self.tx.send_json(msg,flags=zmq.NOBLOCK)
        except zmq.Again:self.send_drops+=1
        return (time.monotonic()-start)*1000
    def close(self):
        time.sleep(.04);self.stop.set();self.thread.join(2);self.tx.close(0);self.context.term()

class Recorder:
    def __init__(self,out):
        self.out=Path(out);self.out.mkdir(parents=True,exist_ok=True);self.q=queue.Queue(maxsize=2)
        self.rows=[];self.drops=0;self.error=None
        self.thread=threading.Thread(target=self.loop,daemon=True);self.thread.start()
    def put(self,image,record):
        try:self.q.put_nowait((image,record))
        except queue.Full:self.drops+=1
    def loop(self):
        try:
            while True:
                item=self.q.get()
                if item is None:break
                image,record=item;name=f"frame_{record['frame_id']:07d}.png"
                if not cv2.imwrite(str(self.out/name),image,[cv2.IMWRITE_PNG_COMPRESSION,1]):raise OSError(name)
                self.rows.append({**record,'image':name})
        except Exception as exc:self.error=repr(exc)
    def close(self):
        # Finite tests only; drain two queued images, never block camera thread.
        self.q.put(None);self.thread.join(10)
        (self.out/'frames.jsonl').write_text(''.join(json.dumps(r)+'\n' for r in self.rows))

class ResourceMonitor:
    def __init__(self):
        self.stop_event=threading.Event();self.rows=[];self.thread=threading.Thread(target=self.loop,daemon=True)
        self.thread.start()
    def loop(self):
        p=psutil.Process();p.cpu_percent()
        while not self.stop_event.wait(.5):
            r={'mono':time.monotonic(),'cpu_pct_one_core_basis':p.cpu_percent(),'rss_mb':p.memory_info().rss/1e6}
            try:
                raw=subprocess.check_output(['nvidia-smi','--query-gpu=utilization.gpu,memory.used','--format=csv,noheader,nounits'],text=True,timeout=2)
                u,m=raw.strip().splitlines()[0].split(',');r.update(gpu_device_util_pct=float(u),gpu_device_memory_mb=float(m))
            except Exception as e:r['gpu_query_error']=repr(e)
            self.rows.append(r)
    def close(self):self.stop_event.set();self.thread.join(3)
