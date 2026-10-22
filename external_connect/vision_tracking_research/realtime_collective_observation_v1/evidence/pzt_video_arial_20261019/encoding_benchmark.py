from pathlib import Path
import sys,time,json,tempfile,statistics
import av,numpy as np
root=Path('/home/lzy/projects/external_connect/vision_tracking_research/realtime_collective_observation_v1');sys.path.insert(0,str(root/'src'))
from experiment_video import ExperimentVideoRecorder,DEFAULT_LAYERS
from protocol import FramePacket
from ycommon import dump
out=root/'evidence/pzt_video_arial_20261019';live=json.loads((root/'evidence/pzt_video_20261019/live_check/validation.json').read_text());session=Path(live['session'])
with av.open(str(session/'experiment_raw.mp4')) as c:images=[f.to_ndarray(format='bgr24') for f in c.decode(video=0)]
rows=[json.loads(l) for l in (session/'video_frames.jsonl').read_text().splitlines()];samples={s['sample_id']:s for s in map(json.loads,(session/'pzt_readbacks.jsonl').read_text().splitlines())}
with tempfile.TemporaryDirectory(prefix='pzt_video_encoding_benchmark_') as d:
 start=time.monotonic();r=ExperimentVideoRecorder(d,dict(video_annotation_wait_ms=0,video_queue=48),start);times=[]
 for i in range(90):
  due=start+i/30;left=due-time.monotonic()
  if left>0:time.sleep(left)
  at=time.monotonic();j=i%len(images);a=images[j]
  p=FramePacket('LIVE_CAMERA','ENCODING BENCHMARK: repeated captured frames, not a new experiment','encoding_benchmark',i,None,a,(a.shape[1],a.shape[0]),'BGR8',arrival_monotonic=at)
  row=rows[j];sample=dict(samples[row['pzt_sample_id']]);sample.pop('sample_id')
  r.submit(p,0,None,DEFAULT_LAYERS,dict(frame_unix_s=row['pzt_frame_unix_s'],sample=sample));times.append((time.monotonic()-at)*1000)
 report=r.close();wall=time.monotonic()-start
 counts={}
 for name in ('experiment_raw.mp4','experiment_overlay.mp4'):
  with av.open(str(Path(d)/name)) as c:counts[name]=sum(1 for _ in c.decode(video=0))
 result=dict(scope='ENCODER_ONLY_REPEATED_LIVE_FRAMES_NOT_NEW_EXPERIMENT',requested_hz=30,submitted=90,written=report['written_frames'],dropped=report['dropped_frames'],wall_s=wall,submit_ms_p50=statistics.median(times),submit_ms_max=max(times),decoded_counts=counts,error=report['error'],source_wh=[1280,1024],overlay_wh=report['parts'][0]['overlay_encoded_wh'],analysis_included=False)
 dump(out/'encoding_benchmark.json',result);print(json.dumps(result))
 assert report['video_complete'] and all(v==90 for v in counts.values())
