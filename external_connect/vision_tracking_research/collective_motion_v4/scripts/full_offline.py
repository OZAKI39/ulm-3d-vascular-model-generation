"""Every source frame, actual PTS, no interpolation, four scheduling strategies."""
import sys,json,time,gzip,argparse,traceback
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from xcommon import *
from observer import Observer
from channel_state import validate_channel
import av,cv2
def main():
 p=argparse.ArgumentParser();p.add_argument('--strategies',default='regular,sparse1,sparse3,sparse6');a=p.parse_args();cv2.setNumThreads(2)
 freeze=json.loads((SERVER/'runs/research_v1/results/CONFIGURATION_FREEZE.json').read_text());assert sha(VIDEO)==VIDEO_SHA
 for strategy in a.strategies.split(','):
  folder=OUT/'results'/('offline_'+strategy);folder.mkdir(exist_ok=True);path=folder/'records.jsonl.gz';assert not path.exists(),'REFUSE_OVERWRITE'
  obs=Observer('B4',freeze['preprocess'],freeze['gap'],int(strategy[-1]) if strategy!='regular' else 1,strategy!='regular');mon=ResourceMonitor();rows=[];start=time.perf_counter()
  with gzip.open(path,'wt',compresslevel=1) as log,av.open(str(VIDEO)) as c:
   c.streams.video[0].codec_context.thread_count=1;it=c.decode(video=0);i=0
   while True:
    t=time.perf_counter()
    try:f=next(it)
    except StopIteration:break
    gray=f.to_ndarray(format='gray');decode_ms=(time.perf_counter()-t)*1000;pts=float(f.pts*f.time_base)
    state=obs.process(i,pts,gray,received=t);validate_channel(state);ts=time.perf_counter();wire=json.dumps(state,separators=(',',':'),allow_nan=False);encode_ms=(time.perf_counter()-ts)*1000;full_ms=(time.perf_counter()-t)*1000
    log.write(wire+'\n');rows.append(dict(frame_id=i,pts_s=pts,decode_ms=decode_ms,algorithm_ms=state['algorithm_ms'],full_ms=full_ms,serialization_ms=encode_ms,scan_ms=state['scan_ms'],group_ms=state['group_ms'],valid_windows=sum(state['windows']['valid']),computed_windows=sum(state['windows']['fresh']),active_windows=sum(state['windows']['active']),safe_windows=sum(state['windows']['safe']),candidates=len(state['candidates_global_xy']),n_groups=len(state['groups']),valid_groups=state['fresh_valid_group_count'],mask_valid=state['channel_mask_valid'],backend_failure=state['backend_failure']))
    if i%100==0:print(strategy,i,'full_ms',round(full_ms,2),'valid_windows',rows[-1]['valid_windows'],'valid_groups',rows[-1]['valid_groups'],flush=True)
    i+=1
  assert [r['frame_id'] for r in rows]==list(range(937))
  dump(folder/'summary.json',dict(identity=identity(),strategy=strategy,freeze=freeze,rows=rows,frame_count=len(rows),coverage=1.,file_sha256=sha(path),wall_s=time.perf_counter()-start,resources=mon.close(),full_pipeline_ms=stats([r['full_ms'] for r in rows]),scope='Full offline analysis including decode, validation and state serialization; compressed research logging timed separately in wall_s. Not paced live throughput.'))
if __name__=='__main__':main()
