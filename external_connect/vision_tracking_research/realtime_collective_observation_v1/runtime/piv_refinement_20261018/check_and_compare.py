from pathlib import Path
import av,cv2,gzip,json,sys,importlib.util
import numpy as np
from fractions import Fraction
root=Path(__file__).resolve().parents[2];sys.path.insert(0,str(root/'src'))
from piv_display import PIVDisplayFilter
from experiment_video import render_overlay
from overlay_style import DEFAULT_LAYERS
from protocol import FramePacket
qa=Path(__file__).parent;run=qa/'ui_replay_final';session=next(run.glob('session_*'))
with gzip.open(session/'states.jsonl.gz','rt') as f:states=[json.loads(x) for x in f]
rows=[json.loads(x) for x in (session/'video_frames.jsonl').read_text().splitlines()]
manifest=json.loads((session/'video_manifest.json').read_text())
assert manifest['dropped_frames']==0 and manifest['written_frames']==manifest['submitted_frames']==len(rows)==303
assert all(not r['layers']['Regions'] and not r['layers']['Shape'] for r in rows)
assert all('piv_display_vectors' not in s for s in states)
f=PIVDisplayFilter();lookup={};prev=None;jump_raw=[];jump_smooth=[]
for s in states:
    if s.get('reset_reason') not in [None,'OK']:f.reset();prev=None
    v=f.update(s);coarse=PIVDisplayFilter().update(s)
    current={tuple(p[:2]):(p[2:],q[2:]) for p,q in zip(coarse,v)}
    if prev is not None:
        for key in current.keys()&prev.keys():
            a,b=current[key];c,d=prev[key];jump_raw.append(np.linalg.norm(a-c));jump_smooth.append(np.linalg.norm(b-d))
    prev=current
    lookup[(s['session_id'],s['frame_id'])]=dict(s,piv_display_vectors=v)
def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod
old=load('piv_old_renderer',qa/'before/experiment_video.py')
old_piv=load('piv_old_geometry',qa/'before/piv_display.py')
old_style=load('piv_old_opacity',qa/'before/overlay_style.py')
for key in ['piv_arrows','arrow_shape','draw_cv_arrows','DISPLAY_SPEC','REFERENCE_SPEED','GAIN_SECONDS','MAX_LENGTH_PX']:setattr(old,key,getattr(old_piv,key))
old.overlay_opacity=old_style.overlay_opacity
out=av.open(str(qa/'piv_style_comparison.mp4'),'w');stream=out.add_stream('libx264',rate=30)
stream.width=1584;stream.height=640;stream.pix_fmt='yuv420p';stream.time_base=Fraction(1,1000000);stream.codec_context.time_base=stream.time_base;stream.options={'crf':'20','preset':'ultrafast','threads':'2','bf':'0'}
checks={};contact=[];first_pts=None
with av.open(str(session/'experiment_raw.mp4')) as video:
    pts=[]
    for i,frame in enumerate(video.decode(video=0)):
        pts.append(float(frame.pts*frame.time_base))
        if not 65<=i<=215:continue
        r=rows[i];a=frame.to_ndarray(format='bgr24');s=lookup.get((r['analysis_session_id'],r['analysis_frame_id']))
        p=FramePacket('REALTIME_REPLAY','recorded camera data',r['source_session_id'],r['source_frame_id'],r['source_pts_s'],a,tuple(r['original_resolution']),pixel_format='BGR8',arrival_monotonic=r['arrival_monotonic'])
        geometry=s.get('geometry') if s else None
        args=(p,s,geometry,DEFAULT_LAYERS,r['annotation_match'],r['annotation_lag_ms'],r['receive_utc'],r['elapsed_since_first_frame_s'])
        left=old.render_overlay(*args);right=render_overlay(*args)
        canvas=np.full((640,1584,3),22,np.uint8)
        canvas[40:640,:792]=cv2.resize(left,(792,600),interpolation=cv2.INTER_AREA)
        canvas[40:640,792:]=cv2.resize(right,(792,600),interpolation=cv2.INTER_AREA)
        cv2.putText(canvas,'BEFORE | thick arrows, age fade',(16,27),cv2.FONT_HERSHEY_SIMPLEX,.65,(235,235,235),1,cv2.LINE_AA)
        cv2.putText(canvas,'AFTER | thin arrows, smoothed display',(808,27),cv2.FONT_HERSHEY_SIMPLEX,.65,(235,235,235),1,cv2.LINE_AA)
        if i in [100,140,180]:contact.append(canvas)
        if first_pts is None:first_pts=r['video_pts_us']
        result=av.VideoFrame.from_ndarray(canvas,format='bgr24');result.pts=r['video_pts_us']-first_pts;result.time_base=stream.time_base
        for packet in stream.encode(result):out.mux(packet)
    checks['experiment_raw.mp4']=dict(frames=len(pts),timestamp_error_s=float(np.max(np.abs(np.asarray(pts)-[r['video_pts_us']/1e6 for r in rows]))))
for packet in stream.encode():out.mux(packet)
out.close()
with av.open(str(session/'experiment_overlay.mp4')) as video:
    pts=[float(f.pts*f.time_base) for f in video.decode(video=0)]
    checks['experiment_overlay.mp4']=dict(frames=len(pts),timestamp_error_s=float(np.max(np.abs(np.asarray(pts)-[r['video_pts_us']/1e6 for r in rows]))))
assert all(v['frames']==303 and v['timestamp_error_s']<1e-6 for v in checks.values())
with av.open(str(qa/'piv_style_comparison.mp4')) as video:assert sum(1 for _ in video.decode(video=0))==151
cv2.imwrite(str(qa/'comparison_contact_sheet.png'),np.concatenate(contact,axis=0))
samples=json.loads((run/'qa_ui_samples.json').read_text());alphas=[s['alpha'] for s in samples if s['valid']]
result=dict(source='Existing Chronos recording, actual Qt replay; not a live camera validation',decoded=checks,
    recorded_frames=303,dropped_frames=0,valid_ui_samples=len(alphas),dimmed_valid_samples=sum(x<1 for x in alphas),
    consecutive_grid_velocity_change=dict(samples=len(jump_raw),raw_spatial_mean_median_px_s=float(np.median(jump_raw)),smoothed_median_px_s=float(np.median(jump_smooth)),scope='Display stability only; not a scientific accuracy improvement'),
    comparison=str(qa/'piv_style_comparison.mp4'),initial_attempt={'recorded_frames':297,'dropped_frames':5,'all_drops_source_frames':[38,46,47,49,50],'mitigation':'Bounded recording queue 24 -> 48 frames for startup burst; no TTL changes'})
(qa/'validation.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
