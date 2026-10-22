import sys,time,json,gzip,argparse
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from ycommon import ROOT,dump,clean,stats,identity,sha
from observer import Observer
from protocol import FramePacket
from pipeline import Pipeline
import av,numpy as np

def main():
 p=argparse.ArgumentParser();p.add_argument('--video',required=True);p.add_argument('--output',required=True);p.add_argument('--backend',default='B4');a=p.parse_args();out=Path(a.output);out.mkdir(parents=True,exist_ok=True)
 cfg=json.loads((ROOT/'configs/default.json').read_text());cfg.update(video=a.video,backend=a.backend,replay_geometry=str(ROOT/'configs/replay_geometry.json'),source_sha256=sha(a.video))
 dump(out/'protocol.json',dict(identity=identity(),config=cfg,synthetic_frozen_seeds=[95001,95002],thresholds_frozen=True))
 observer=Observer(cfg);start=time.perf_counter();rows=[];sourcepts=[]
 with gzip.open(out/'offline.jsonl.gz','wt') as logfile,av.open(a.video) as video:
  video.streams.video[0].codec_context.thread_count=2
  for fid,f in enumerate(video.decode(video=0)):
   pts=float(f.pts*f.time_base);im=f.to_ndarray(format='gray');p=FramePacket('REALTIME_REPLAY',a.video,'offline_all_937',fid,pts,im,(im.shape[1],im.shape[0]));r=observer.analyze(p,reset=fid==0);r.update(frame_id=fid,source_pts=pts,session_id='offline_all_937');logfile.write(json.dumps(clean(r),separators=(',',':'))+'\n');rows.append(r['timings']);sourcepts.append(pts)
 dump(out/'offline_summary.json',dict(identity=identity(),source_sha256=sha(a.video),frames=len(rows),source_pts_s=sourcepts,wall_s=time.perf_counter()-start,stages={k:stats([r[k] for r in rows]) for k in rows[0]},all_frames_analyzed=True))
 for speed,name in [(1.,'replay_1x'),(2.,'stress_2x'),(.5,'debug_half')]:
  cfg['speed']=speed;cfg['loop']=False;p=Pipeline(cfg,out/name).start();start=time.monotonic()
  try:
   while time.monotonic()-start<90:
    time.sleep(.05)
    if p.status=='EOF' and p.input.peek() is None:break
    if p.error:raise RuntimeError(p.error)
  finally:r=p.stop()
  print(name,r['counters'],r['fresh_valid_hz'],flush=True)
if __name__=='__main__':main()
