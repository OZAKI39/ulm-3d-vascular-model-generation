"""Same-server causal latest-frame replay; no live camera/device connection."""
import sys,os,time,threading,json,argparse,collections
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np,cv2,av,torch,jsonschema
from vcommon import OUT,ROOT,P,NAMES,dump,stats,ResourceMonitor,identity,sha
from heatmaps import centers_to_heatmap
from inference import Predictor,estimate_from_prediction
from legacy import CandidateDetector,distribution,piv,empty,CausalWindow,MockConsumer
from state_v2 import new_state,to_v1,VersionedGate

class PacedRecording:
    def __init__(self,video):
        self.video=video;self.lock=threading.Lock();self.latest=None;self.done=threading.Event();self.produced=[];self.error=None;self.start=time.monotonic();self.thread=threading.Thread(target=self.run,daemon=True)
    def run(self):
        begin,end=P['real_replay']['frame_range'];origin=None;pts0=None
        try:
            with av.open(self.video) as c:
                stream=c.streams.video[0];stream.codec_context.thread_count=1
                for fid,frame in enumerate(c.decode(stream)):
                    if fid<begin:continue
                    if fid>=end:break
                    decode_start=time.monotonic();bgr=frame.to_ndarray(format='bgr24');pts=float(frame.time)
                    if origin is None:origin=time.monotonic();pts0=pts
                    due=origin+pts-pts0
                    time.sleep(max(0,due-time.monotonic()));received=time.monotonic()
                    row=dict(frame_id=fid,pts=pts,due=due,received=received,convert_ms=(received-decode_start)*1000,decode_time_scope='conversion plus scheduled wait; full decoder startup recorded separately')
                    with self.lock:self.latest=(fid,pts,received,bgr);self.produced.append(row)
        except Exception as exc:self.error=repr(exc)
        finally:self.done.set()
    def take(self,last):
        with self.lock:return self.latest if self.latest is not None and self.latest[0]!=last else None

@torch.inference_mode()
def pure_model_profile(model):
    x=torch.randn(1,3,2,128,128,device='cuda');dt=torch.full((1,2),.0333,device='cuda');times=[]
    for i in range(70):
        torch.cuda.synchronize();t=time.perf_counter()
        with torch.autocast('cuda',dtype=torch.bfloat16):out=model.model(x,dt)
        torch.cuda.synchronize()
        if i>=10:times.append((time.perf_counter()-t)*1000)
    return stats(times)

def run(method,video,frozen):
    cv2.setNumThreads(1);torch.set_num_threads(4);start=time.perf_counter();model=None;cal=None
    if method!='PIV_SINGLE':
        key=frozen['champions_by_validation_only'][method];cal=frozen['models'][key];model=Predictor(OUT/cal['checkpoint'])
    cold=time.perf_counter()-start;pure=pure_model_profile(model) if model else None
    detector=CandidateDetector();window=CausalWindow(3);heats={};stream=PacedRecording(video);gate=VersionedGate();consumer=MockConsumer();rows=[];last=None
    schema=jsonschema.Draft202012Validator(json.loads((ROOT/'configs/collective_state_v2.schema.json').read_text()))
    monitor=ResourceMonitor();torch.cuda.reset_peak_memory_stats();stream.thread.start();x,y,w,h=P['real_replay']['roi_xywh'];start=time.monotonic()
    try:
        while True:
            item=stream.take(last)
            if item is None:
                if stream.done.is_set():break
                time.sleep(.0005);continue
            fid,pts,received,bgr=item;last=fid;t0=time.perf_counter();crop=bgr[y:y+h,x:x+w]
            t=time.perf_counter();objects,gray=detector(crop);detect_ms=(time.perf_counter()-t)*1000
            t=time.perf_counter();dist,density=distribution(objects,gray.shape);distribution_ms=(time.perf_counter()-t)*1000
            ready=window.add(fid,pts,gray,objects);heat_ms=0
            if model:
                t=time.perf_counter();heats[fid]=centers_to_heatmap([o['xy'] for o in objects],gray.shape);heat_ms=(time.perf_counter()-t)*1000
                heats={k:v for k,v in heats.items() if k>=fid-2}
            estimate=empty('CAUSAL_WINDOW_NOT_READY',method);prediction=None;timing={};vectors=[]
            if ready:
                frames=list(window.frames);dt=frames[-1][1]-frames[0][1]
                if model:
                    prediction,timing=model.frames([z[2] for z in frames],[heats[z[0]] for z in frames],[z[1] for z in frames]);estimate=estimate_from_prediction(prediction,cal,method,model.sha)
                else:
                    t=time.perf_counter();r=piv(frames[0][2],frames[-1][2],dt,mode='candidate',objects_a=frames[0][3],objects_b=frames[-1][3]);timing={'piv_ms':(time.perf_counter()-t)*1000};estimate=r['estimate'];vectors=r['vectors'];estimate['observability_status']='VALID_GLOBAL_MOTION' if all(estimate['valid'].values()) else ('NON_AFFINE_MOTION' if 'NON_AFFINE_OR_INCONSISTENT_VECTORS' in estimate['reasons'].values() else 'NOT_OBSERVABLE')
            ids=[r[0] for r in window.frames];ts=[r[1] for r in window.frames]
            state=new_state('REAL_RECORDING_CAUSAL_REPLAY',method+'_20261017',fid,pts,received,[x,y,w,h],dist,estimate,{},dict(frame_ids=ids,pts_s=ts,history_span_s=ts[-1]-ts[0],estimate_reference='interval-equivalent generator; timestamp midpoint',reference_pts_s=(ts[0]+ts[-1])/2),model_metadata=dict(name=method,weight_sha256=model.sha if model else 'NO_LEARNED_WEIGHTS',calibration_version='F10V-synthetic-groups-v1' if model else 'F10U-frozen-quality-gate'))
            schema.validate(state);checked,age=gate.accept(state);consumer.send(to_v1(checked));pipeline_ms=(time.perf_counter()-t0)*1000
            rows.append(dict(frame_id=fid,pts_s=pts,objects=objects,state=checked,fresh=age['fresh'],six_valid=all(checked['motion']['valid'].values()),pipeline_ms=pipeline_ms,detect_ms=detect_ms,heatmap_ms=heat_ms,distribution_ms=distribution_ms,stages=timing,age_at_finish_ms=(time.monotonic()-received)*1000,interval_reference_age_ms=(time.monotonic()-received)*1000+(ts[-1]-ts[0])*500,vectors=vectors))
        stream.thread.join(2)
        # Explicit stale watchdog verification after the source ends.
        time.sleep(.13);expired=gate.current();stale_ok=expired is not None and expired['stale_or_discontinuous']
    finally:consumer.close()
    duration=stream.produced[-1]['pts']-stream.produced[0]['pts'] if stream.produced else 0
    resource=monitor.close();fresh=sum(r['fresh'] and r['six_valid'] for r in rows);readyrows=[r for r in rows if len(r['state']['time_window']['frame_ids'])==3]
    report=dict(method=method,identity=identity(),source_sha256=sha(video),source='REAL_RECORDING_CAUSAL_REPLAY',real_accuracy='NOT_VERIFIED',camera_end_to_end_latency='NOT_MEASURED',cold_model_load_s=cold,pure_model_ms=pure,recording_duration_s=duration,cold_decode_to_first_publish_s=stream.produced[0]['received']-stream.start if stream.produced else None,
        n_produced=len(stream.produced),n_consumed=len(rows),n_dropped=len(stream.produced)-len(rows),drop_fraction=(len(stream.produced)-len(rows))/max(1,len(stream.produced)),queue_capacity=1,n_ready=len(readyrows),n_fresh_valid=fresh,fresh_six_hz=fresh/duration if duration else 0,processed_hz=len(rows)/duration if duration else 0,expired_fraction=sum(not r['fresh'] for r in rows)/max(1,len(rows)),no_new_frame_TTL_pass=stale_ok,
        full_pipeline_ms=stats([r['pipeline_ms'] for r in readyrows]),all_frames_pipeline_ms=stats([r['pipeline_ms'] for r in rows]),detect_ms=stats([r['detect_ms'] for r in rows]),heatmap_ms=stats([r['heatmap_ms'] for r in rows]),age_ms=stats([r['age_at_finish_ms'] for r in rows]),interval_reference_age_ms=stats([r['interval_reference_age_ms'] for r in readyrows]),peak_cuda_allocated_mb=torch.cuda.max_memory_allocated()/2**20,parameter_count=model.parameter_count if model else 0,resource=resource,consumer_rows=consumer.rows,consumer_watchdog_events=consumer.watchdog_events,source_error=stream.error,producer=stream.produced,rows=rows)
    dump(OUT/'results'/f'real_replay_{method}.json',report)
    print(method,report['full_pipeline_ms'],report['fresh_six_hz'],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--video',default=os.environ.get('F10V_VIDEO',P['source_video']));a=p.parse_args();frozen=json.loads((OUT/'results/model_freeze.json').read_text())
    for method in ['PIV_SINGLE','direct_gray','direct_fused','geo_fused']:run(method,a.video,frozen)
