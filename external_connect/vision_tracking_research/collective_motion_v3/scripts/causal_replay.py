"""Full 937-frame real-PTS paced replay, latest-slot=1, no hardware."""
import sys,time,threading,json,argparse
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from wcommon import *
from benchmark import *
import cv2,av,jsonschema
class PacedRecording:
    def __init__(self):
        self.lock=threading.Lock();self.latest=None;self.done=threading.Event();self.produced=[];self.error=None;self.thread=threading.Thread(target=self.run,daemon=True)
    def run(self):
        origin=None;pts0=None
        try:
            with av.open(str(VIDEO)) as c:
                stream=c.streams.video[0];stream.codec_context.thread_count=1
                for fid,frame in enumerate(c.decode(stream)):
                    tick=time.perf_counter();bgr=frame.to_ndarray(format='bgr24');convert_ms=(time.perf_counter()-tick)*1000;pts=float(frame.time)
                    if origin is None:origin=time.monotonic();pts0=pts
                    due=origin+pts-pts0;time.sleep(max(0,due-time.monotonic()));received=time.monotonic()
                    with self.lock:
                        self.latest=(fid,pts,received,bgr);self.produced.append(dict(frame_id=fid,pts=pts,due=due,received=received,convert_ms=convert_ms,scheduler_lateness_ms=(received-due)*1000))
        except Exception as exc:self.error=repr(exc)
        finally:self.done.set()
    def take(self,last):
        with self.lock:return self.latest if self.latest is not None and self.latest[0]!=last else None
def run(method):
    cv2.setNumThreads(1);start=time.perf_counter();model=construct(method);cold=time.perf_counter()-start;weight_sha=sha(model.path) if model else 'NO_LEARNED_WEIGHTS'
    frozen=json.loads((OUT/'results/model_freeze.json').read_text()) if method!='PIV' else None
    config=frozen['selected'][method] if frozen else {};threshold=frozen['gates'][method]['threshold'] if frozen else None
    nwin=16 if method=='CoTracker3' else 3;window=old.CausalWindow(nwin);det=old.CandidateDetector();gate=VersionedGate();consumer=old.MockConsumer();stream=PacedRecording();rows=[];last=None;monitor=ResourceMonitor()
    schema=jsonschema.Draft202012Validator(json.loads((ROOT.parent/'collective_motion_v2/configs/collective_state_v2.schema.json').read_text()))
    x,y,w,h=P['roi_xywh'];stream.thread.start();start=time.monotonic()
    try:
        while True:
            item=stream.take(last)
            if item is None:
                if stream.done.is_set():break
                time.sleep(.0005);continue
            fid,pts,received,bgr=item;last=fid;t0=time.perf_counter();tick=time.perf_counter();objects,gray=det(bgr[y:y+h,x:x+w]);detect_ms=(time.perf_counter()-tick)*1000
            tick=time.perf_counter();dist,_=old.distribution(objects,gray.shape);distribution_ms=(time.perf_counter()-tick)*1000
            ready=window.add(fid,pts,gray,objects);frames=list(window.frames);ts=[z[1] for z in frames];ids=[z[0] for z in frames]
            r=dict(estimate=old.empty(window.last_reset or 'WARMUP',method),quality_score=0.,timing={},vectors=[])
            if ready:r=predict(method,model,[z[2] for z in frames],[z[3] for z in frames],ts,config)
            base=r['estimate'];estimate=apply_gate(base,r['quality_score'],threshold) if frozen else base
            # PIV remains exactly its historical gate; extra AI validation gate is separately visible.
            if all(estimate['valid'].values()):estimate['observability_status']='VALID_GLOBAL_MOTION'
            elif 'NON_AFFINE_OR_INCONSISTENT_VECTORS' in estimate['reasons'].values():estimate['observability_status']='NON_AFFINE_MOTION'
            reference=(ts[-3]+ts[-1])/2 if ready else ts[-1]
            state=new_state('REAL_RECORDING_CAUSAL_REPLAY','F10W_'+method,fid,pts,received,[x,y,w,h],dist,estimate,{},dict(frame_ids=ids,pts_s=ts,history_span_s=ts[-1]-ts[0],reference_pts_s=reference,estimate_reference='last three frames, endpoint logm; fixed ROI center'),model_metadata=dict(name=method,weight_sha256=weight_sha,calibration_version='F10U frozen gate' if not frozen else 'F10W validation empirical gate'))
            schema.validate(state);checked,age=gate.accept(state);sent=consumer.send(to_v1(checked));pipeline_ms=(time.perf_counter()-t0)*1000
            rows.append(clean(dict(frame_id=fid,pts_s=pts,ready=ready,state=checked,fresh=age['fresh'],six_valid=all(checked['motion']['valid'].values()),base_six_valid=all(base['valid'].values()),pipeline_ms=pipeline_ms,detect_ms=detect_ms,distribution_ms=distribution_ms,stages=r['timing'],quality_score=r['quality_score'],age_ms=(time.monotonic()-received)*1000,reference_age_ms=(time.monotonic()-received)*1000+1000*(pts-reference),sent=sent)))
        stream.thread.join(2);time.sleep(.13);stale=gate.current();stale_ok=bool(stale and stale['stale_or_discontinuous'])
    finally:consumer.close()
    duration=stream.produced[-1]['pts']-stream.produced[0]['pts'];ready=[r for r in rows if r['ready']];nvalid=sum(r['fresh'] and r['six_valid'] for r in rows);processed={r['frame_id'] for r in rows};dropped=[r['frame_id'] for r in stream.produced if r['frame_id'] not in processed]
    report=dict(identity=identity(),method=method,config=config,source_sha256=sha(VIDEO),source_type='PACED_RECORDING_NOT_LIVE',n_input=len(stream.produced),n_processed=len(rows),n_dropped=len(dropped),dropped_frame_ids=dropped,queue_capacity=1,queue_backlog_max=1,n_ready=len(ready),n_fresh_valid=nvalid,fresh_six_hz=nvalid/duration,base_valid_count=sum(r['base_six_valid'] for r in rows),expired_fraction=sum(not r['fresh'] for r in rows)/len(rows),full_pipeline_ms=stats([r['pipeline_ms'] for r in ready]),all_frames_pipeline_ms=stats([r['pipeline_ms'] for r in rows]),reference_age_ms=stats([r['reference_age_ms'] for r in ready]),age_ms=stats([r['age_ms'] for r in rows]),model_stage_ms=stats([r['stages'].get('model_ms',float('nan')) for r in ready]),cold_model_load_s=cold,history_frames=nwin,first_ready_pts_s=ready[0]['pts_s'] if ready else None,no_new_frame_TTL_pass=stale_ok,wall_s=time.monotonic()-start,resource=monitor.close(),producer=stream.produced,consumer=consumer.rows,source_error=stream.error,rows=rows,real_six_parameter_accuracy='NOT_VERIFIED',camera_to_controller_end_to_end_latency='NOT_MEASURED')
    dump(OUT/'results'/f'causal_replay_{method}.json',report);print(method,'input',len(stream.produced),'processed',len(rows),'freshHz',nvalid/duration,report['full_pipeline_ms'],flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--methods',nargs='+',default=['PIV','LRAFT_ST','CoTracker3','POT_unbalanced']);a=p.parse_args()
    for method in a.methods:run(method)
