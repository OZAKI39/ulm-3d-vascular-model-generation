import argparse,json,time,sys,os
from pathlib import Path
from ycommon import ROOT,dump,sha
from pipeline import Pipeline

def main(mode='REALTIME_REPLAY',headless=False):
 p=argparse.ArgumentParser();p.add_argument('--config',default=str(ROOT/'configs/default.json'));p.add_argument('--video');p.add_argument('--backend',choices=['B1','B2','B4']);p.add_argument('--speed',type=float);p.add_argument('--loop',action='store_true');p.add_argument('--seconds',type=float);p.add_argument('--output');p.add_argument('--autostart',action='store_true');p.add_argument('--capture-ui');p.add_argument('--simulate-events',action='store_true');p.add_argument('--live-geometry');p.add_argument('--windows-project');p.add_argument('--windows-python');a=p.parse_args()
 cfg=json.loads(Path(a.config).read_text());cfg['source_mode']=mode;cfg['replay_geometry']=str(ROOT/'configs/replay_geometry.json')
 for k in ['video','backend','speed','live_geometry','windows_project','windows_python']:
  if getattr(a,k) is not None:cfg[k]=getattr(a,k)
 cfg['loop']=a.loop or cfg['loop']
 if mode=='REALTIME_REPLAY':cfg['source_sha256']=sha(cfg['video'])
 out=Path(a.output or ROOT/f'runtime/{time.time_ns()}')
 if headless:
  pipe=Pipeline(cfg,out).start();start=time.monotonic()
  try:
   while not a.seconds or time.monotonic()-start<a.seconds:
    time.sleep(.05)
    if pipe.error or (pipe.status=='EOF' and pipe.input.peek() is None):break
  finally:
   r=pipe.stop();print(json.dumps({k:r[k] for k in ['wall_s','counters','analysis_skipped','fresh_valid_hz','error']}))
  if r['error']:raise RuntimeError(r['error'])
 else:
  from PySide6.QtWidgets import QApplication
  from ui import Window
  app=QApplication(sys.argv);win=Window(cfg,out,a.autostart,a.seconds,a.capture_ui);win.show()
  if a.simulate_events:
   from PySide6.QtCore import QTimer
   QTimer.singleShot(5100,lambda:win.mark_event('ON'));QTimer.singleShot(18100,lambda:win.mark_event('OFF'))
  return app.exec()
