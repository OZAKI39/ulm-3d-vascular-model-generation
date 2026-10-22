import os
os.environ['QT_QPA_PLATFORM']='offscreen'
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='2'
os.environ['RAYON_NUM_THREADS']='4'
from pathlib import Path
import sys,json,time
root=Path('/home/lzy/projects/external_connect/vision_tracking_research/realtime_collective_observation_v1');sys.path.insert(0,str(root/'src'))
os.environ['F10X_SERVER']=str(root/'vendor')
from PySide6.QtWidgets import QApplication
from PySide6.QtTest import QTest
from ui import Window
from ycommon import dump
out=root/'evidence/pzt_video_20261019/live_check'
cfg=json.loads((root/'configs/default.json').read_text());cfg.update(source_mode='LIVE_CAMERA',backend='B4',bridge_session_s=45)
app=QApplication([]);w=Window(cfg,out);w.show()
report={}
try:
 deadline=time.monotonic()+6
 while w.pzt_panel.data.get('controller') in (None,'CLOSED') and time.monotonic()<deadline:QTest.qWait(50)
 w.pzt_panel.read_button.click();deadline=time.monotonic()+21
 while w.pzt_panel.pending_request and time.monotonic()<deadline:QTest.qWait(50)
 report['pzt_read_ack']=w.pzt_panel.data.get('refresh');report['pzt_before']=w.pzt_panel.data
 w.buttons['Start'].click();started=time.monotonic();first=None
 while time.monotonic()-started<25:
  QTest.qWait(50)
  if w.pipe.video_record.written and first is None:first=time.monotonic()
  if first and time.monotonic()-first>8:break
  if w.pipe.error:break
 report['elapsed_s']=time.monotonic()-started;report['frames_before_stop']=w.pipe.video_record.written
 report['last_source_status']=w.pipe.status;report['session']=str(w.pipe.output)
 w.grab().save(str(out/'live_ui.png'));w.buttons['Stop'].click()
 report['video']=w.pipe.video_record.report();report['pzt_after']=w.pzt_panel.data
finally:
 w.close();app.processEvents();dump(out/'validation.json',report)
print(json.dumps({k:report[k] for k in ['frames_before_stop','last_source_status','session','elapsed_s']}))
