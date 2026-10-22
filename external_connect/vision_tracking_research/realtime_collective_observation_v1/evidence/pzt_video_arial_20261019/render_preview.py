from pathlib import Path
import sys,json,hashlib
from fractions import Fraction
import av,cv2
root=Path('/home/lzy/projects/external_connect/vision_tracking_research/realtime_collective_observation_v1');sys.path.insert(0,str(root/'src'))
from experiment_video import render_overlay
from pzt_video import PZTVideoRenderer
from video_text import font_record
from protocol import FramePacket
from ycommon import dump
out=root/'evidence/pzt_video_arial_20261019';src=Path(json.loads((root/'evidence/pzt_video_20261019/live_check/validation.json').read_text())['session'])
rows=[json.loads(s) for s in (src/'video_frames.jsonl').read_text().splitlines()];samples={s['sample_id']:s for s in map(json.loads,(src/'pzt_readbacks.jsonl').read_text().splitlines())}
renderer=PZTVideoRenderer();target=out/'experiment_overlay_arial.mp4';container=av.open(str(target),'w');stream=None;count=0
with av.open(str(src/'experiment_raw.mp4')) as raw:
 for i,frame in enumerate(raw.decode(video=0)):
  row=rows[i];a=frame.to_ndarray(format='bgr24');p=FramePacket(row['source_mode'],row['source_name'],row['source_session_id'],row['source_frame_id'],row['source_pts_s'],a,tuple(row['original_resolution']),'BGR8',arrival_monotonic=row['arrival_monotonic'])
  assert row['analysis_frame_id'] is None,'This preview script expects the original no-mask test.'
  image=render_overlay(p,None,None,row['layers'],row['annotation_match'],row['annotation_lag_ms'],row['receive_utc'],row['elapsed_since_first_frame_s'])
  image=renderer.compose(image,dict(frame_unix_s=row['pzt_frame_unix_s'],sample=samples[row['pzt_sample_id']]))
  h,w=image.shape[:2]
  if h%2 or w%2:image=cv2.copyMakeBorder(image,0,h%2,0,w%2,cv2.BORDER_REPLICATE)
  if stream is None:
   stream=container.add_stream('libx264',rate=30);stream.width=image.shape[1];stream.height=image.shape[0];stream.pix_fmt='yuv420p';stream.time_base=Fraction(1,1000000);stream.codec_context.time_base=stream.time_base;stream.options={'crf':'18','preset':'ultrafast','tune':'zerolatency','bf':'0','threads':'2'}
  f=av.VideoFrame.from_ndarray(image,format='bgr24');f.pts=row['video_pts_us'];f.time_base=Fraction(1,1000000)
  for packet in stream.encode(f):container.mux(packet)
  if i==3:cv2.imwrite(str(out/'arial_preview.png'),image)
  count+=1
for packet in stream.encode():container.mux(packet)
container.close()
with av.open(str(target)) as c:pts=[float(f.pts*f.time_base) for f in c.decode(video=0)]
assert len(pts)==len(rows) and all(abs(t-r['video_pts_us']/1e6)<2e-6 for t,r in zip(pts,rows))
dump(out/'preview.json',dict(source=str(src),scope='Rerender of previously recorded real camera/AFG data; no new device commands',frames=count,original_pts_preserved=True,font=font_record(),video=str(target),video_sha256=hashlib.sha256(target.read_bytes()).hexdigest()))
print(json.dumps(dict(frames=count,font=font_record(),video=str(target))))
