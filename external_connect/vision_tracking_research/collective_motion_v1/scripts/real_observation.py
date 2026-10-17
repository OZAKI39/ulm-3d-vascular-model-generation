"""Actual-PTS recording analysis, followed by latest-slot causal replay or live reception."""
import sys,time,json,threading,subprocess,argparse
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from common import *
import av,cv2,psutil
from distribution import CandidateDetector,distribution
from methods import piv,cpd
from affine import empty
from state import make_state,MockConsumer,CausalWindow

def estimate_pair(a,b,override=None):
    mode=override or json.loads((ROOT/'results/selection.json').read_text())['mode']
    return piv(a[2],b[2],b[1]-a[1],mode='local' if mode=='multipass' else mode,objects_a=a[3],objects_b=b[3],multipass=mode=='multipass',bootstrap=0)

def decode_roi():
    lo,hi=P['real_frame_range'];x,y,w,h=P['real_roi_xywh']
    with av.open(str(VIDEO)) as c:
        for fid,f in enumerate(c.decode(video=0)):
            if lo<=fid<hi:yield fid,float(f.time),f.to_ndarray(format='bgr24')[y:y+h,x:x+w].copy()
            if fid>=hi:break

def offline():
    start=time.perf_counter();det=CandidateDetector();history=[];rows=[];frames=[];times=[]
    for fid,pts,im in decode_roi():
        tick=time.perf_counter();objects,gray=det(im);dist,density=distribution(objects,gray.shape);history.append((fid,pts,gray,objects));methods={}
        for n in P['windows']:
            if len(history)>=n:
                r=estimate_pair(history[-n],history[-1]);methods['PIV_'+str(n)]=r
                if (fid-P['real_frame_range'][0])%5==0:
                    methods['CPD_'+str(n)]=cpd([o['xy'] for o in history[-n][3]],[o['xy'] for o in objects],pts-history[-n][1],gray.shape)
        rows.append({'frame_id':fid,'pts':pts,'objects':objects,'distribution':dist,'methods':methods,'processing_ms':(time.perf_counter()-tick)*1000,'focus_laplacian_var':float(cv2.Laplacian(gray,cv2.CV_32F).var())})
        frames.append(im);times.append(pts)
    np.savez_compressed(ROOT/'results/real_frames.npz',images=frames,pts=times)
    dump(ROOT/'results/real_observation.json',{'identity':identity(),'source':'RECORDED_CHRONOS_OFFLINE_CAUSAL_ENDPOINTS','source_video':str(VIDEO),'source_sha256':sha(VIDEO),'roi':P['real_roi_xywh'],'rows':rows,'elapsed_s':time.perf_counter()-start,'truth_six_parameters':'NOT_AVAILABLE'})
    print('offline',len(rows),round(time.perf_counter()-start,2),flush=True)

class ReplaySlot:
    def __init__(self):self.lock=threading.Lock();self.latest=None;self.done=threading.Event();self.records=[];self.error=None
    def start(self):self.thread=threading.Thread(target=self.run,daemon=True);self.thread.start();return self
    def run(self):
        try:
            origin=None;base=None
            for fid,pts,im in decode_roi():
                if origin is None:origin=time.monotonic();base=pts
                delay=origin+pts-base-time.monotonic()
                if delay>0:self.done.wait(delay)
                received=time.monotonic()
                with self.lock:self.latest=(fid,pts,received,im)
                self.records.append({'frame_id':fid,'pts':pts,'receive':received})
        except Exception as e:self.error=repr(e)
        finally:self.done.set()
    def snapshot(self):
        with self.lock:return self.latest

def replay(override=None):
    det=CandidateDetector();buffer=CausalWindow(3);mock=MockConsumer();stream=ReplaySlot().start();last=None;records=[];start=time.monotonic();proc=psutil.Process();cpu0=proc.cpu_times();resources=[]
    while True:
        q=stream.snapshot()
        if q is None or q[0]==last:
            if stream.done.is_set():break
            time.sleep(.0005);continue
        fid,pts,received,im=q;last=fid;tick=time.perf_counter();objects,gray=det(im);dist,_=distribution(objects,gray.shape)
        ready=buffer.add(fid,pts,gray,objects);r=estimate_pair(buffer.frames[0],buffer.frames[-1],override) if ready else {'estimate':empty('CAUSAL_WARMUP_OR_FRAME_GAP'),'algorithm_ms':0}
        state=make_state('CAUSAL_REPLAY_TEST','f10u_replay',fid,pts,received,P['real_roi_xywh'],dist,r['estimate'],{'focus_laplacian_var':float(cv2.Laplacian(gray,cv2.CV_32F).var())},{'frames':3,'span_s':pts-buffer.frames[0][1],'causal':True})
        if buffer.last_reset:state['stale_or_discontinuous']=True
        sent=mock.send(state);records.append({'frame_id':fid,'processing_ms':(time.perf_counter()-tick)*1000,'piv_ms':r['algorithm_ms'],'fit_ms':r['estimate'].get('fit_ms'),'state':state,'sent':sent})
        resources.append({'t':time.monotonic()-start,'rss_mb':proc.memory_info().rss/2**20})
    mock.close();elapsed=time.monotonic()-start;cpu1=proc.cpu_times()
    dump(ROOT/('results/causal_replay_fast.json' if override else 'results/causal_replay.json'),{'identity':identity(),'method_override':override,'scope':'secondary speed probe without new heldout accuracy claim' if override else 'frozen development winner','source':'CAUSAL_REPLAY_TEST','LIVE_CAMERA_PASS':False,'records':records,'consumer':mock.rows,'producer':stream.records,'error':stream.error,'elapsed_s':elapsed,
         'input_count':len(stream.records),'processed_count':len(records),'coverage_fraction':len(records)/max(1,len(stream.records)),
         'software_processing_ms':stats([r['processing_ms'] for r in records]),'piv_ms':stats([r['piv_ms'] for r in records if r['piv_ms']]),'fit_ms':stats([r['fit_ms'] for r in records if r['fit_ms'] is not None]),'transfer_ms':stats([r['transfer_ms'] for r in mock.rows]),
         'output_rate_hz':len(mock.rows)/elapsed,'valid_update_rate_hz':{k:sum(r['valid'][k] for r in mock.rows)/elapsed for k in NAMES},'unobservable_fraction':{k:1-np.mean([r['valid'][k] for r in mock.rows]) for k in NAMES},'queue_capacity':1,'queue_backlog':'bounded latest-slot overwrite; skipped frame IDs counted',
         'cpu_core_equivalents':((cpu1.user+cpu1.system)-(cpu0.user+cpu0.system))/elapsed,'resources':resources,'gpu_used_for_replay':False,'CAMERA_TO_CONTROLLER_END_TO_END_LATENCY':'NOT_MEASURED'})
    print('replay',len(records),'/',len(stream.records),'p95',stats([r['processing_ms'] for r in records]),flush=True)

def camera():
    sys.path.insert(0,str(ROOT.parent/'smart_online_tracking_v2/src'))
    from live_bridge import WindowsLatestBridge
    out=ROOT/'evidence/camera_readonly';b=WindowsLatestBridge(out,seconds=P['camera_readonly_seconds'],roi=[496,568,256,384]);records=[];last=None;buf=CausalWindow(3);det=CandidateDetector();mock=MockConsumer();start=time.monotonic();error=None
    try:
        b.start()
        while time.monotonic()-start<P['camera_readonly_seconds']+2 and not b.done.is_set():
            q=b.snapshot()
            if q is None or q[0]['seq']==last:time.sleep(.002);continue
            header,im=q;last=header['seq'];received=header['arrival_mono']-header['age_upper_bound_at_arrival_ms']/1000
            pts=header['sender_decode_mono'];objects,gray=det(im);dist,_=distribution(objects,gray.shape)
            ready=buf.add(last,pts,gray,objects);r=estimate_pair(buf.frames[0],buf.frames[-1]) if ready else {'estimate':empty('WARMUP')}
            s=make_state('LIVE_CHRONOS_RECEIVER_DECODE_CLOCK',header['session'],last,None,received,[496,568,256,384],dist,r['estimate'],{}, {'frames':3,'span_s':pts-buf.frames[0][1],'timestamp_semantics':'sender decode clock; exposure timestamps unavailable'})
            mock.send(s);records.append(s)
    except Exception as exc:error=repr(exc)
    finally:
        if b.process is not None:b.stop()
        mock.close()
    dump(ROOT/'results/camera_status.json',{'identity':identity(),'read_only':True,'frames':len(records),'states':records,'consumer':mock.rows,'LIVE_CAMERA_PASS':len(records)>0,'error':error or b.error,'elapsed_s':time.monotonic()-start,'network_settings_changed':False,'hardware_actions':False})
    print('camera frames',len(records),'error',error or b.error,flush=True)

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('mode',choices=['offline','camera','replay','replay_fast']);mode=a.parse_args().mode
    if mode=='replay_fast':replay('candidate')
    else:globals()[mode]()
