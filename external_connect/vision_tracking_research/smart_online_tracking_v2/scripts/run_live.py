"""Finite true LIVE_CAMERA benchmark through the unchanged Windows receiver."""
from pathlib import Path
import argparse,sys,time,json,copy,traceback
import numpy as np,cv2,psutil
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from vision_engine import VisionEngine,stats
from live_bridge import WindowsLatestBridge
from runtime import MockConsumer,Recorder,ResourceMonitor,enforce_freshness,identity,invalid_state
from annotation_store import sha

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--route',choices=['opencv','cellpose'],required=True)
    ap.add_argument('--tag',required=True);ap.add_argument('--seconds',type=float,default=20)
    ap.add_argument('--detector-config',type=Path);a=ap.parse_args()
    protocol=json.loads((ROOT/'protocol.json').read_text());roi=protocol['inference_budget']['live_roi_xywh']
    out=ROOT/'results'/a.tag;out.mkdir(exist_ok=False)
    (out/'execution_identity.json').write_text(json.dumps(identity(ROOT),indent=2))
    config=json.loads(a.detector_config.read_text()) if a.detector_config else {}
    if config.get('weights_sha256') and sha(config['weights'])!=config['weights_sha256']:raise ValueError('Model checksum mismatch')
    monitor=ResourceMonitor();engine=VisionEngine(a.route,origin=roi[:2],detector_config=config)
    # Report startup/JIT separately. Only past, saved calibration images used here.
    warm=cv2.imread(str(ROOT/'results/bridge_smoke/first_frame.png'));x,y,w,h=roi;warm=warm[y:y+h,x:x+w]
    warms=[]
    for _ in range(3):
        t=time.monotonic();engine.detect(warm);warms.append((time.monotonic()-t)*1000)
    engine.reset('warmup_complete')
    bridge=WindowsLatestBridge(out/'bridge',seconds=min(180,int(a.seconds+35)),roi=roi).start()
    mock=MockConsumer();recorder=Recorder(out/'recorded');rows=[];last=None;first=None;start=time.monotonic()
    deadline=start+25+a.seconds;heart=time.monotonic();stale_start=0;error=None;source0=None
    try:
        while time.monotonic()<deadline and not bridge.done.is_set():
            item=bridge.snapshot();now=time.monotonic()
            if first is not None and now-first>=a.seconds:break
            if item is None or (item[0]['session'],item[0]['seq'])==last:
                if now-heart>.1:
                    empty={'session':last[0] if last else 'NO_SESSION','frame_id':last[1] if last else -1,
                           'arrival_mono':now,'age_upper_bound_at_arrival_ms':None,'image_available':False,
                           'age_upper_bound_ms':None,'status':'NO_VALID_TARGET','reason':'NO_NEW_IMAGE',
                           'candidates':[],'targets':[],'group_shape':None,'visual_candidate_available':False,
                           'control_permission':False,'biological_valid':False}
                    empty['age_upper_bound_at_arrival_ms']=1000
                    mock.send(empty);heart=now
                    if item is None:engine.reset('no_frame')
                time.sleep(.001);continue
            header,frame=item;last=(header['session'],header['seq']);heart=now
            first=first or now;source0=source0 or header['sender_decode_mono']
            header={**header,'source':'LIVE_CAMERA','source_time_s':header['sender_decode_mono'],
                    'source_time_semantics':'Windows host decode completion seconds; not sensor exposure PTS'}
            if header['age_upper_bound_at_arrival_ms']+1000*(now-header['arrival_mono'])>100:
                stale_start+=1;engine.reset('stale_frame_before_inference')
                mock.send(invalid_state(header['session'],header['seq'],'STALE_FRAME_BEFORE_INFERENCE'));continue
            tic=time.monotonic();msg,timing,labels=engine.process(frame,header)
            raw=copy.deepcopy(msg);fresh=enforce_freshness(msg)
            if not fresh:engine.reset('stale_result')
            timing['serialize_send_ms']=mock.send(msg)
            timing['software_ms']=(time.monotonic()-tic)*1000
            record={'frame_id':header['seq'],'session':header['session'],'source_time_s':header['source_time_s'],
                    'source_relative_s':header['source_time_s']-source0,'source':'LIVE_CAMERA','roi_xywh':roi,
                    'header':header,'message':raw,'output_fresh':fresh,'timing':timing}
            recorder.put(frame,record)
            rows.append({'frame_id':header['seq'],'source_time_s':header['source_time_s'],
                         'candidate_count':len(raw['candidates']),'observed_targets':sum(t['xy_px'] is not None for t in raw['targets']),
                         'lost_targets':sum(t['xy_px'] is None for t in raw['targets']),
                         'fresh_at_send':fresh,'age_upper_bound_ms':msg['age_upper_bound_ms'],**timing})
    except Exception as exc:
        error=repr(exc);(out/'failure.txt').write_text(traceback.format_exc())
    finally:
        measured=time.monotonic()-(first or start)
        mock.send(invalid_state(last[0] if last else 'NO_SESSION',last[1] if last else -1,'STREAM_ENDED_OR_RUN_STOPPED'))
        bridge.stop();mock.close();recorder.close();monitor.close()
    headers=[r for r in bridge.records if r['status']=='FRAME']
    selected=[r for r in headers if rows and rows[0]['frame_id']<=r['seq']<=rows[-1]['frame_id']]
    seqspan=rows[-1]['frame_id']-rows[0]['frame_id']+1 if rows else 0
    times=[r['sender_decode_mono'] for r in selected]
    receiver=json.loads((out/'bridge/windows_receiver.json').read_text()) if (out/'bridge/windows_receiver.json').exists() else {}
    if not rows:error=error or 'NO_LIVE_FRAMES: '+str(bridge.error)
    summary={'route':a.route,'source':'LIVE_CAMERA','detector_config':config,'wall_measured_s':measured,'startup_wait_s':None if first is None else first-start,
             'roi_xywh':roi,'original_resolution_wh':[1280,1024],'resized':False,'processed_frames':len(rows),
             'processed_hz':len(rows)/measured,'input_seq_span':seqspan,
             'input_host_decode_hz':(selected[-1]['seq']-selected[0]['seq'])/(times[-1]-times[0]) if len(times)>1 else None,
             'bridge_frames_in_processing_span':len(selected),'bridge_or_receiver_skipped':max(0,seqspan-len(selected)),
             'vision_latest_overwrites_in_processing_span':max(0,len(selected)-len(rows)),
             'unprocessed_source_frames_in_span':max(0,seqspan-len(rows)),
             'fresh_candidate_states':sum(r['candidate_valid'] for r in mock.rows),
             'fresh_candidate_hz':sum(r['candidate_valid'] for r in mock.rows)/measured,
             'validated_bubble_hz':0,'validated_group_hz':0,'stale_start_rejected':stale_start,
             'stale_result_rejected':sum(not r['fresh_at_send'] for r in rows),'mock_stale_rejected':sum(not r['fresh'] for r in mock.rows),
             'mock_no_valid_target':len(mock.rows),'mock_received':len(mock.rows),'mock_send_drops':mock.send_drops,
             'recorder_drops':recorder.drops,'recorder_error':recorder.error,'error':error,
             'stage_ms':{k:stats([r[k] for r in rows]) for k in ['detect_ms','track_ms','morphology_ms','flow_ms','engine_ms','serialize_send_ms','software_ms']},
             'get_roundtrip_ms':stats([r['get_roundtrip_ms'] for r in selected]),
             'frame_age_at_bridge_ms':stats([r['age_upper_bound_at_arrival_ms'] for r in selected]),
             'mock_transport_ms':stats([r['transport_ms'] for r in mock.rows]),
             'mock_consume_ms':stats([r['consume_ms'] for r in mock.rows]),
             'mock_age_upper_bound_ms':stats([r['age_upper_bound_ms'] for r in mock.rows if r['frame_id']>=0]),
             'warmup_detect_ms':warms,'reset_events':engine.reset_events,
             'resources':{k:stats([r[k] for r in monitor.rows if k in r]) for k in ['cpu_pct_one_core_basis','rss_mb','gpu_device_util_pct','gpu_device_memory_mb']},
             'resource_scope':'WSL vision process CPU/RAM; device-wide GPU includes annotation GUI and Windows workloads',
             'gpu_torch_peak_allocated_mb':engine.torch.cuda.max_memory_allocated()/1e6 if engine.torch else 0,
             'windows_receiver':receiver.get('receiver'),'windows_cpu_seconds':receiver.get('cpu_seconds'),
             'windows_wall_s':receiver.get('wall_s'),
             'windows_pipe_write_ms':stats([r['send_write_ms'] for r in receiver.get('messages',[]) if r['status']=='FRAME']),
             'EXPOSURE_TIMESTAMP_UNKNOWN':True,'CAMERA_TO_CONTROLLER_END_TO_END_LATENCY':'NOT_MEASURED',
             'HARDWARE_CONTROL_DISABLED':True,'biological_accuracy':'EXPERT_GROUND_TRUTH_PENDING'}
    for name,data in [('summary.json',summary),('resources.json',monitor.rows),('mock.json',mock.rows)]:
        (out/name).write_text(json.dumps(data,indent=2))
    (out/'frames.jsonl').write_text(''.join(json.dumps(r)+'\n' for r in rows))
    print(json.dumps({k:summary[k] for k in ['route','processed_frames','processed_hz','fresh_candidate_hz','stale_result_rejected','error']}))
    if error:raise RuntimeError(error)
if __name__=='__main__':main()
