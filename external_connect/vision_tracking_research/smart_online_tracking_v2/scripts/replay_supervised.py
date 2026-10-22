"""Strictly ordered two-route replay; offline processing is never called live FPS."""
from pathlib import Path
import sys,json,time,copy
import numpy as np,cv2,av
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from vision_engine import VisionEngine,stats
from runtime import identity
from annotation_store import sha
from prepare_candidates import VIDEO

def main():
 study=ROOT/'results/supervised_20261016_v1';out=study/'replay';out.mkdir(exist_ok=False);images=out/'images';images.mkdir()
 configs={r:json.loads((study/f'selected_{r}.json').read_text()) for r in ('opencv','cellpose')}
 for c in configs.values():
  if c.get('weights_sha256') and sha(c['weights'])!=c['weights_sha256']:raise ValueError('Checkpoint changed')
 # Same original-pixel crop for both, fixed before replay inference.
 roi=[650,180,550,560];x,y,w,h=roi
 (out/'preregistered.json').write_text(json.dumps({'roi_xywh':roi,'frames':[360,419],
     'selection':'Fixed view of vessel used identically for both routes; no future frame input',
     'configs':configs,'execution':identity(ROOT)},indent=2))
 engines={r:VisionEngine(r,origin=(x,y),detector_config=c) for r,c in configs.items()}
 freeze=ROOT/'results/expert_freezes/vision2_supervised_20261016_v1'
 frozen=json.loads((freeze/'frozen.json').read_text());train=next(r for r in frozen['rows'] if r['record']['split']=='train')
 warm=cv2.imread(str(freeze/train['frozen_image_file']))[y:y+h,x:x+w]
 for e in engines.values():e.detect(warm);e.reset('past_training_image_warmup')
 rows=[];started=time.monotonic()
 with av.open(str(VIDEO)) as container:
  stream=container.streams.video[0];stream.codec_context.thread_count=1
  for index,frame in enumerate(container.decode(stream)):
   if index<360:continue
   if index>419:break
   image=frame.to_ndarray(format='bgr24')[y:y+h,x:x+w];pts=float(frame.time)
   filename=f'{index:06d}.png';cv2.imwrite(str(images/filename),image)
   for route,engine in engines.items():
    header={'seq':index,'source_time_s':pts,'session':'SUPERVISED_CAUSAL_REPLAY','source':'CAUSAL_REPLAY',
            'source_time_semantics':'recorded video PTS; offline processing','arrival_mono':time.monotonic(),
            'age_upper_bound_at_arrival_ms':0,'roi_xywh':roi}
    msg,timing,_=engine.process(image,header)
    rows.append({'frame':index,'image':filename,'source_pts_s':pts,'route':route,'message':msg,'timing':timing})
 (out/'records.jsonl').write_text(''.join(json.dumps(r)+'\n' for r in rows))
 summary={'source':'CAUSAL_REPLAY','roi_xywh':roi,'source_frame_range':[360,419],'frames_per_route':len(rows)//2,
          'wall_s':time.monotonic()-started,'video_path':str(VIDEO),'video_sha256':sha(VIDEO),'execution':identity(ROOT),
          'configurations':configs,'future_frames_used':False,'interpolation':False,'accuracy':'NOT_EVALUATED_NO_HUMAN_TRACK_IDS',
          'routes':{r:{'timing_ms':{k:stats([x['timing'][k] for x in rows if x['route']==r]) for k in rows[0]['timing']},
                        'software_ids':len({t['id'] for x in rows if x['route']==r for t in x['message']['targets']}),
                        'lost_observations':sum(t['status']=='LOST' for x in rows if x['route']==r for t in x['message']['targets'])}
                    for r in engines}}
 (out/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps({'frames_per_route':len(rows)//2,'wall_s':summary['wall_s']}))

if __name__=='__main__':main()
