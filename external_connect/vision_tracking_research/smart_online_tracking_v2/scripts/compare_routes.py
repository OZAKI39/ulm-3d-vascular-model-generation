"""Paired images and strictly causal sequential replay; NOT a live fps claim."""
from pathlib import Path
import sys,time,json,hashlib
import cv2,numpy as np,av
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from vision_engine import VisionEngine,stats
from runtime import identity,ResourceMonitor
from prepare_candidates import VIDEO

def main():
    p=json.loads((ROOT/'protocol.json').read_text());budget=p['inference_budget']
    static=set(budget['static_frame_indices']);lo,hi=budget['causal_comparison_frames_inclusive']
    roi=budget['saved_video_roi_xywh'];out=ROOT/'results/paired_comparison';out.mkdir(exist_ok=False)
    (out/'execution_identity.json').write_text(json.dumps(identity(ROOT),indent=2))
    images=out/'images';images.mkdir();rows=[];monitor=ResourceMonitor();start=time.monotonic()
    models={name:VisionEngine(name,origin=roi[:2]) for name in ['opencv','cellpose']}
    warm=cv2.imread(str(ROOT/'data/candidates/video_000000.png'));warmup=[]
    for route,engine in models.items():
        for scope in ['FULL_FRAME','ROI']:
            x,y,w,h=roi;frame=warm if scope=='FULL_FRAME' else warm[y:y+h,x:x+w]
            tick=time.monotonic();engine.detect(frame);warmup.append({'route':route,'scope':scope,'ms':1000*(time.monotonic()-tick)})
    # Each decoder iteration provides only the present frame. Neither engine sees future frames.
    with av.open(str(VIDEO)) as c:
        stream=c.streams.video[0];stream.codec_context.thread_count=1
        for i,f in enumerate(c.decode(stream)):
            if i not in static and not lo<=i<=hi:continue
            full=f.to_ndarray(format='bgr24');x,y,w,h=roi;crop=full[y:y+h,x:x+w]
            if i in static:
                for scope,frame in [('FULL_FRAME',full),('ROI',crop)]:
                    for route,engine in models.items():
                        tick=time.monotonic();objs,gray,labels=engine.detect(frame);elapsed=(time.monotonic()-tick)*1000
                        name=f'static_{i:06d}_{scope}_{route}'
                        np.save(out/f'{name}_labels.npy',labels)
                        (out/f'{name}.json').write_text(json.dumps({'objects':objs,'scope':scope,'roi_xywh':[0,0,full.shape[1],full.shape[0]] if scope=='FULL_FRAME' else roi,
                            'source_pts':f.pts,'source_time_base':str(f.time_base),'source_time_s':float(f.time),'status':'UNREVIEWED'},indent=2))
                        rows.append({'kind':'STATIC_PAIRED','frame':i,'scope':scope,'route':route,'count':len(objs),'detect_ms':elapsed})
                cv2.imwrite(str(images/f'static_{i:06d}.png'),full,[cv2.IMWRITE_PNG_COMPRESSION,1])
            if lo<=i<=hi:
                frame_rows=[]
                for route,engine in models.items():
                    header={'seq':i,'session':'saved_video_test_block','source_time_s':float(f.time),'source':'CAUSAL_REPLAY_OFFLINE',
                            'source_time_semantics':'actual encoded video PTS seconds; not wall-clock realtime',
                            'arrival_mono':time.monotonic(),'age_upper_bound_at_arrival_ms':0,'roi_xywh':roi}
                    msg,timing,labels=engine.process(crop,header)
                    row={'kind':'CAUSAL_REPLAY_OFFLINE','frame':i,'scope':'ROI','route':route,'count':len(msg['candidates']),
                         'source_pts':f.pts,'source_time_base':str(f.time_base),'source_time_s':float(f.time),
                         'message':msg,**timing}
                    rows.append(row);frame_rows.append(row)
                cv2.imwrite(str(images/f'replay_{i:06d}.png'),crop,[cv2.IMWRITE_PNG_COMPRESSION,1])
            if i>=max(hi,max(static)):break
    monitor.close();wall=time.monotonic()-start
    (out/'records.jsonl').write_text(''.join(json.dumps(r)+'\n' for r in rows))
    (out/'resources.json').write_text(json.dumps(monitor.rows,indent=2))
    summary={'source_file':str(VIDEO),'source_sha256':hashlib.sha256(VIDEO.read_bytes()).hexdigest(),'roi_xywh':roi,
         'same_source_frames':True,'same_pixels':True,'no_resizing':True,'no_future_frame_in_tracking':True,
         'replay_mode':'CAUSAL_REPLAY_OFFLINE; sequential computation, not paced realtime',
         'comparison_source_frames':[lo,hi],'wall_s':wall,'warmup':warmup,'expert_accuracy':'NOT_EVALUATED',
         'scope_results':{f'{kind}/{scope}/{route}':{'n':len(part),'candidates':stats([r['count'] for r in part]),
                                                'detect_ms':stats([r['detect_ms'] for r in part]),
                                                'engine_ms':stats([r['engine_ms'] for r in part if 'engine_ms' in r])}
              for kind in ['STATIC_PAIRED','CAUSAL_REPLAY_OFFLINE'] for scope in ['FULL_FRAME','ROI'] for route in models
              if (part:=[r for r in rows if r['kind']==kind and r['scope']==scope and r['route']==route])}}
    (out/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary['scope_results'],indent=2))
if __name__=='__main__':main()
