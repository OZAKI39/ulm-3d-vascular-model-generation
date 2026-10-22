"""Real Qt controls, isolated labels; process launching mocked to avoid duplicate jobs."""
from pathlib import Path
import sys,argparse,json,traceback
from unittest.mock import patch
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from annotate import Workbench,QApplication,QTimer,QFont,napari,auto_sam
OUT=ROOT/'evidence/auto_sam_training_20261016';OUT.mkdir(exist_ok=True)
app=QApplication.instance() or QApplication([]);app.setFont(QFont('Microsoft YaHei',10))
report={'software_fixture':True,'checks':{},'PASS':False}
patches=[patch.object(auto_sam,'register_client'),patch.object(auto_sam,'unregister_client'),patch.object(auto_sam,'ensure_service'),
         patch.object(auto_sam,'set_enabled'),patch.object(auto_sam,'config',return_value={'enabled':True}),patch.object(auto_sam,'rollback')]
mocks=[p.start() for p in patches]
w=Workbench(argparse.Namespace(index=0,smoke=False,manual_only=False,no_auto_train=False,annotation_dir=str(OUT/'ui_sandbox')))
def run():
 try:
  report['checks']['default_checkbox_enabled']=w.auto_training.value and w.auto_training.enabled
  report['checks']['starts_background_service']=mocks[0].call_count==1 and mocks[2].call_count==1
  w.auto_training.value=False;report['checks']['pause_toggle']=mocks[3].call_args.args==(False,)
  w.auto_training.value=True;report['checks']['resume_toggle']=mocks[3].call_args.args==(True,)
  state={'state':'CANDIDATE_HELD','last_attempt':{'reason':'WAITING_FOR_INDEPENDENT_VALIDATION'}}
  with patch.object(auto_sam,'read_json',return_value=state):w.poll_training()
  report['checks']['validation_wait_explained']='2 张验证帧' in w.training_status.value
  w.layers['SINGLE_BUBBLE'].data[10,10]=1;w.restore_original_sam()
  report['checks']['rollback_keeps_labels']=bool(mocks[5].call_count==1 and w.layers['SINGLE_BUBBLE'].data[10,10]==1)
  w.notes.value='SOFTWARE_UI_FIXTURE; not human annotations'
  w.viewer.screenshot(path=str(OUT/'auto_training_window.png'),canvas_only=False)
  w.viewer.window._qt_window.close();report['checks']['close_unregisters_client']=mocks[1].call_count==1
  report['checks']['close_stops_poll_timers']=not w.training_timer.isActive() and not w.watch.isActive()
  report['PASS']=all(report['checks'].values())
 except Exception:report['error']=traceback.format_exc();w.closing=True;w.viewer.close()
 finally:
  for p in patches:p.stop()
  (OUT/'ui_checks.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True);app.quit()
QTimer.singleShot(1500,run);napari.run()
if not report['PASS']:raise SystemExit(1)
