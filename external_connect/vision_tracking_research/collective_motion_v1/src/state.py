"""M4 observation-only schema and bounded local consumer; no hardware imports."""
import time,threading,collections
import numpy as np
import zmq
from common import NAMES,P,clean
from affine import empty

def make_state(source,session,frame_id,pts,received,roi,dist,estimate,image_quality,window):
    now=time.monotonic()
    return clean({'schema_version':'M4-collective-1.0','source':source,'session_id':session,'frame_id':frame_id,
      'available_capture_timestamp':None,'capture_timestamp_status':'NOT_AVAILABLE',
      'host_receive_timestamp':received,'processing_finished_timestamp':now,'host_clock':'monotonic_s_same_process',
      'source_PTS':pts,'ROI':roi,'calibration_status':'PIXEL_ONLY_NO_PHYSICAL_SCALE','image_quality':image_quality,
      'distribution':dist,'motion':estimate,'time_window':window,'age_ms':(now-received)*1000,'stale_or_discontinuous':False,
      'six_parameter_units':dict(zip(NAMES,['px/s','px/s','1/s','1/s','1/s','1/s'])),
      'camera_to_controller_end_to_end_latency':'NOT_MEASURED','action_event':None,'hardware_control_executed':False,'training_ready':False,
      'scientific_accuracy_verified':False,'observability_scope':'Algebraic/quality checks; material-motion accuracy is not guaranteed by valid=true'})

def invalidate(state,reason):
    import copy
    state=copy.deepcopy(state);state['stale_or_discontinuous']=True
    state['motion']=empty(reason,state['motion']['method'])
    # Stale visible distribution must not masquerade as current observation either.
    state['distribution']['field_validity']={k:False for k in state['distribution'].get('field_validity',{})}
    for key in ['visible_count','group_center_xy','covariance_2x2','gyration_radius','distribution_axis_ratio','distribution_orientation','density_summary']:state['distribution'][key]=None
    return state

class ConsumerGate:
    def __init__(self):self.last=None;self.pts=None;self.session=None;self.latest=None
    def accept(self,obj,now=None):
        now=time.monotonic() if now is None else now;reason=None
        if obj['session_id']!=self.session:self.last=None;self.pts=None;self.session=obj['session_id']
        if self.last is not None and obj['frame_id']<=self.last:reason='DUPLICATE_OR_OLD_FRAME'
        if self.pts is not None and obj['source_PTS'] is not None and obj['source_PTS']<=self.pts:reason='PTS_DISCONTINUITY'
        age=(now-obj['host_receive_timestamp'])*1000
        if age<0 or age>P['replay_stale_limit_ms']:reason='STALE_OBSERVATION'
        if obj['stale_or_discontinuous']:reason='UPSTREAM_DISCONTINUITY'
        self.last=max(obj['frame_id'],self.last if self.last is not None else obj['frame_id'])
        if obj['source_PTS'] is not None:self.pts=max(obj['source_PTS'],self.pts if self.pts is not None else obj['source_PTS'])
        result=invalidate(obj,reason) if reason else obj
        self.latest=result
        return result,{'frame_id':obj['frame_id'],'age_ms':age,'fresh':reason is None,'reason':reason,'valid':result['motion']['valid'],'transfer_ms':max(0,(now-obj.get('send_monotonic',now))*1000)}
    def current(self,now=None):
        now=time.monotonic() if now is None else now
        if self.latest is None:return None
        return invalidate(self.latest,'NO_NEW_FRAME_EXPIRED') if (now-self.latest['host_receive_timestamp'])*1000>P['replay_stale_limit_ms'] else self.latest

def tensor_input(state):
    """NumPy adapter. Missing = NaN plus explicit mask; torch.as_tensor optional."""
    vals=state['motion']['values'];valid=state['motion']['valid']
    mask=np.array([valid[k] and vals[k] is not None for k in NAMES],bool)
    return np.array([vals[k] if m else np.nan for k,m in zip(NAMES,mask)],np.float32),mask

class MockConsumer:
    def __init__(self):
        self.ctx=zmq.Context();self.stop_event=threading.Event();self.ready=threading.Event();self.rows=[];self.states=[];self.watchdog_events=[]
        self.thread=threading.Thread(target=self._run,daemon=True);self.thread.start();self.ready.wait(2)
        self.tx=self.ctx.socket(zmq.PUSH);self.tx.setsockopt(zmq.SNDHWM,1);self.tx.setsockopt(zmq.LINGER,0);self.tx.connect(self.endpoint)
    def _run(self):
        rx=self.ctx.socket(zmq.PULL);rx.setsockopt(zmq.CONFLATE,1);rx.setsockopt(zmq.RCVHWM,1);p=rx.bind_to_random_port('tcp://127.0.0.1');self.endpoint=f'tcp://127.0.0.1:{p}';self.ready.set();gate=ConsumerGate()
        while not self.stop_event.is_set():
            if rx.poll(10):
                s,row=gate.accept(rx.recv_json());self.rows.append(row);self.states.append(s)
            else:
                current=gate.current()
                if current and current['motion']['reasons'].get('translation_vx')=='NO_NEW_FRAME_EXPIRED' and (not self.watchdog_events or self.watchdog_events[-1]['frame_id']!=current['frame_id']):
                    self.watchdog_events.append({'frame_id':current['frame_id'],'time':time.monotonic(),'reason':'NO_NEW_FRAME_EXPIRED'})
        rx.close(0)
    def send(self,state):
        state['send_monotonic']=time.monotonic()
        try:self.tx.send_json(clean(state),flags=zmq.NOBLOCK);return True
        except zmq.Again:return False
    def close(self):
        self.stop_event.wait(.03);self.stop_event.set();self.thread.join(2);self.tx.close(0);self.ctx.term()

class CausalWindow:
    def __init__(self,size):self.size=size;self.frames=collections.deque(maxlen=size);self.last_reset=None
    def add(self,fid,pts,gray,objects):
        self.last_reset=None
        if self.frames and (fid!=self.frames[-1][0]+1 or pts<=self.frames[-1][1] or pts-self.frames[-1][1]>.15):
            self.frames.clear();self.last_reset='FRAME_GAP_OR_PTS_DISCONTINUITY'
        self.frames.append((fid,pts,gray.copy(),objects))
        return len(self.frames)==self.size
