"""Full-file actual-PTS paced producer, latest-slot one and real local M4 consumer."""
import sys,time,json,threading,gzip,argparse
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from xcommon import *
from observer import Observer
from channel_state import validate_channel,to_v1,expire
from state import MockConsumer
import av,cv2
class Producer:
 def __init__(self):
  self.lock=threading.Condition();self.latest=None;self.done=False;self.rows=[];self.error=None;self.thread=threading.Thread(target=self.run,daemon=True)
 def run(self):
  try:
   with av.open(str(VIDEO)) as c:
    c.streams.video[0].codec_context.thread_count=1;it=c.decode(video=0);origin=time.monotonic();fid=0
    while True:
     t=time.monotonic()
     try:f=next(it)
     except StopIteration:break
     gray=f.to_ndarray(format='gray');decode_ms=(time.monotonic()-t)*1000;pts=float(f.pts*f.time_base);due=origin+pts
     time.sleep(max(0,due-time.monotonic()));received=time.monotonic()
     row=dict(frame_id=fid,pts_s=pts,decode_ms=decode_ms,due_mono=due,received_mono=received)
     with self.lock:self.latest=(fid,pts,gray,row);self.rows.append(row);self.lock.notify()
     fid+=1
  except Exception as e:self.error=repr(e)
  finally:
   with self.lock:self.done=True;self.lock.notify_all()
def main():
 p=argparse.ArgumentParser();p.add_argument('--strategy',default='regular');a=p.parse_args();cv2.setNumThreads(2)
 assert sha(VIDEO)==VIDEO_SHA;freeze=json.loads((SERVER/'runs/research_v1/results/CONFIGURATION_FREEZE.json').read_text());folder=OUT/'results'/('replay_'+a.strategy);folder.mkdir(exist_ok=True);assert not (folder/'summary.json').exists()
 # CUDA/import warmup from the first available image only; never future data.
 _,_,first=next(frames([0]));warm=Observer('B4',freeze['preprocess'],freeze['gap'],1,False);t=time.monotonic()
 for i in range(5):warm.process(i,i*.0333,first)
 warm_s=time.monotonic()-t
 obs=Observer('B4',freeze['preprocess'],freeze['gap'],int(a.strategy[-1]) if a.strategy!='regular' else 1,a.strategy!='regular');consumer=MockConsumer();producer=Producer();mon=ResourceMonitor();producer.thread.start();rows=[];states=[];last=-1
 while True:
  with producer.lock:
   producer.lock.wait_for(lambda:(producer.latest is not None and producer.latest[0]!=last) or producer.done,timeout=.1)
   item=producer.latest
   if item is None or item[0]==last:
    if producer.done:break
    continue
  fid,pts,g,pr=item;start=time.monotonic();state=obs.process(fid,pts,g,received=pr['received_mono'],session='F10X_CAUSAL_'+a.strategy);state=expire(state,time.monotonic());validate_channel(state)
  group=next((v['observation_key'] for v in state['groups'] if all(v['six']['valid'].values())),None)
  sent=consumer.send(to_v1(state,group));wire=json.dumps(state,separators=(',',':'),allow_nan=False);end=time.monotonic()
  full_ms=(end-start)*1000+pr['decode_ms'];fresh=(end-pr['received_mono'])*1000<=100 and not state['stale_or_discontinuous']
  rows.append(dict(frame_id=fid,pts_s=pts,full_ms=full_ms,service_ms=(end-start)*1000,decode_ms=pr['decode_ms'],age_ms=(end-pr['received_mono'])*1000,scan_ms=state['scan_ms'],algorithm_ms=state['algorithm_ms'],group_ms=state['group_ms'],component_times=state['component_times'],fresh=fresh,valid_groups=state['fresh_valid_group_count'] if fresh else 0,raw_valid_groups=state['fresh_valid_group_count'],valid_windows=sum(state['windows']['valid']) if fresh else 0,computed_windows=sum(state['windows']['fresh']),safe_windows=sum(state['windows']['safe']),selected_group_for_legacy_adapter=group,sent=sent,frame_gap=fid-last-1))
  states.append(wire);last=fid
 producer.thread.join(2);consumer.close();resources=mon.close()
 assert len(producer.rows)==937 and producer.error is None
 elapsed=producer.rows[-1]['pts_s']+.0333;processed={r['frame_id'] for r in rows};dropped=sorted(set(range(937))-processed)
 with gzip.open(folder/'records.jsonl.gz','wt',compresslevel=1) as f:f.write('\n'.join(states)+'\n')
 dump(folder/'summary.json',dict(identity=identity(),strategy=a.strategy,warmup_s=warm_s,rows=rows,producer_rows=producer.rows,producer_error=producer.error,input_frames=937,processed_frames=len(rows),dropped_frames=dropped,max_pending_queue=1,full_pipeline_ms=stats([r['full_ms'] for r in rows]),ready_pipeline_ms=stats([r['full_ms'] for r in rows if r['computed_windows']>0]),source_duration_s=elapsed,full_channel_scan_hz=len(rows)/elapsed,fresh_valid_group_motion_frame_hz=sum(r['valid_groups']>0 for r in rows)/elapsed,fresh_valid_group_states_hz=sum(r['valid_groups'] for r in rows)/elapsed,local_piv_update_frame_hz=sum(r['computed_windows']>0 for r in rows)/elapsed,consumer_rows=consumer.rows,consumer_watchdog=consumer.watchdog_events,resources=resources,record_sha256=sha(folder/'records.jsonl.gz'),live_camera=False,camera_end_to_end_latency='NOT_MEASURED',scope='Actual PTS release; only past/current frames; processing plus per-frame decoder cost; no camera. Group-state Hz counts local spatial observations, not unique physical groups.'))
 print(a.strategy,'processed',len(rows),'drops',len(dropped),'p95',stats([r['full_ms'] for r in rows]),flush=True)
if __name__=='__main__':main()
