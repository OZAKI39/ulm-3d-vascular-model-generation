import os,sys,json,time
os.environ.setdefault('QT_QPA_PLATFORM','offscreen')
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from PySide6.QtWidgets import QApplication
from PySide6.QtCore import QEvent
from PySide6.QtTest import QTest
from ycommon import ROOT,sha
from ui import Window

def test_ui_buttons_and_unexpected_close(tmp_path):
 app=QApplication.instance() or QApplication([]);cfg=json.loads((ROOT/'configs/default.json').read_text());cfg['source_sha256']=sha(cfg['video']);w=Window(cfg,tmp_path);w.show();app.processEvents()
 assert w.event(QEvent(QEvent.Type.User)) in [True,False] # must not override Qt event with PZT handler
 w.buttons['Start'].click();QTest.qWait(700);w.buttons['Pause'].click();QTest.qWait(150);assert not w.pipe.snapshot()['valid'];w.buttons['Resume'].click();QTest.qWait(200)
 w.mark_event('ON');assert w.pipe.events[-1]['hardware_control'] is False
 w.checks['PIV'].setChecked(False);w.buttons['Stop'].click();assert not w.pipe.worker.is_alive()
 w.buttons['Start'].click();QTest.qWait(300);w.close();assert not w.pipe.worker.is_alive();assert (tmp_path/'ui_summary.json').exists()
