import json,sys,time
from pathlib import Path
root=Path(__file__).resolve().parents[2];sys.path.insert(0,str(root/'src'))
from PySide6 import QtWidgets as W,QtCore as C
from ui import Window
from ycommon import sha,stats
from overlay_style import overlay_opacity
out=Path(__file__).parent/sys.argv[1];out.mkdir(parents=True,exist_ok=False)
cfg=json.loads((root/'configs/default.json').read_text());cfg.update(source_mode='REALTIME_REPLAY',backend='B4',source_sha256=sha(cfg['video']))
app=W.QApplication([]);window=Window(cfg,out,auto=True,seconds=10);window.canvas.diagnostics=True;window.show()
records=[];saved=False
probe=C.QTimer()
def sample():
    global saved
    s=window.canvas.snap
    if not s or not s['state']:return
    p=s['display_packet'];state=s['state']
    records.append(dict(t=time.monotonic(),raw_frame=p.frame_id,analysis_frame=state['frame_id'],fresh=s['fresh'],valid=s['valid'],age_ms=s['age_ms'],alpha=overlay_opacity(s['age_ms'],s['ttl_ms']) if s['fresh'] else 0))
    if not saved and time.monotonic()-window.started>5 and s['valid']:
        window.grab().save(str(out/'ui_stable_default_layers.png'));saved=True
probe.timeout.connect(sample);probe.start(50)
app.exec();probe.stop()
(out/'qa_ui_samples.json').write_text(json.dumps(records,indent=2))
paint=stats(window.canvas.paint_times);(out/'qa_paint_ms.json').write_text(json.dumps(paint,indent=2))
r=window.pipe.summary;v=r['video_recording']
result=dict(wall_s=r['wall_s'],counters=r['counters'],recorded_frames=v['written_frames'],dropped_frames=v['dropped_frames'],matches=v['match_counts'],screenshot_saved=saved,paint_ms=paint)
(out/'qa_summary.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
