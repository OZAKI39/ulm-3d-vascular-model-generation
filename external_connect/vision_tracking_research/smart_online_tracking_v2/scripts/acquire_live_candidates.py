"""Append 12 real, original-resolution live images; never overwrite existing labels."""
from pathlib import Path
import sys,time,json,hashlib,datetime
import cv2
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from live_bridge import WindowsLatestBridge
from runtime import identity
out=ROOT/'results/live_candidate_capture';out.mkdir(exist_ok=False)
(out/'execution_identity.json').write_text(json.dumps(identity(ROOT),indent=2))
bridge=WindowsLatestBridge(out,seconds=35).start();rows=[];next_sample=None;deadline=time.monotonic()+30
try:
    while len(rows)<12 and time.monotonic()<deadline and not bridge.done.is_set():
        item=bridge.snapshot()
        if item is None:time.sleep(.01);continue
        h,frame=item;now=time.monotonic()
        if next_sample is not None and now<next_sample:time.sleep(.01);continue
        name=f'live_20261016_{len(rows):03d}';p=ROOT/'data/candidates'/f'{name}.png'
        if p.exists():raise FileExistsError(p)
        cv2.imwrite(str(p),frame,[cv2.IMWRITE_PNG_COMPRESSION,3])
        digest=hashlib.sha256(p.read_bytes()).hexdigest();w,hh=frame.shape[1],frame.shape[0]
        rows.append({'id':name,'image':str(p.relative_to(ROOT)),'image_sha256':digest,
             'source_kind':'LIVE_CAMERA','source_file':None,'source_sha256':digest,
             'source_sha_semantics':'saved raw-frame PNG; live RTSP has no source file',
             'source_frame_index':h['seq'],'source_pts':None,'source_time_base':None,'source_pts_s':None,
             'host_decode_monotonic_s':h['sender_decode_mono'],'host_decode_utc':h['decode_utc'],
             'timestamp_status':'EXPOSURE_TIMESTAMP_UNKNOWN; sensor PTS not exposed by original receiver',
             'resolution_wh':[w,hh],'saved_roi_xywh':[0,0,w,hh],'resized':False,
             'split':'QUARANTINE_UNASSIGNED','time_block':h['session'],'selection_kind':'live_one_second_spacing',
             'sequence_id':None,'review_status':'UNREVIEWED','camera_configuration':h['camera_config'],
             'capture_fps':None,'objective_magnification':None,'adapter':None,'focal_plane':None,
             'scene_class':'EXPERT_SCENE_CONFIRMATION_PENDING'})
        next_sample=now+1
finally:bridge.stop()
if len(rows)!=12:raise RuntimeError(f'Only {len(rows)} images; original manifest unchanged')
p=ROOT/'data/candidate_manifest.json';manifest=json.loads(p.read_text())
assert not ({r['id'] for r in rows}&{r['id'] for r in manifest['rows']})
manifest['rows']+=rows;manifest['live_append_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
p.with_suffix('.tmp').write_text(json.dumps(manifest,indent=2));p.with_suffix('.tmp').replace(p)
(out/'capture_manifest.json').write_text(json.dumps(rows,indent=2));print('Appended',len(rows),'total',len(manifest['rows']))
