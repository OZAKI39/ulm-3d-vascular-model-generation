import os
os.environ['QT_QPA_PLATFORM']='offscreen'
import sys,time,json,threading
from pathlib import Path
import numpy as np,cv2
from PySide6 import QtWidgets as W,QtCore as C
from PySide6.QtTest import QTest
ROOT=Path(__file__).resolve().parents[2]
os.environ['F10X_SERVER']=str(ROOT/'vendor')
sys.path.insert(0,str(ROOT/'src'))
from ui import Window
from pipeline import Pipeline
from protocol import FramePacket

class SyntheticSource:
    def __init__(self,emit,status):self.emit=emit;self.status=status;self.done=threading.Event();self.starts=0;self.sent=0
    def start(self):
        self.starts+=1;self.thread=threading.Thread(target=self.run);self.thread.start();return self
    def run(self):
        self.status('STREAMING');rng=np.random.default_rng(30)
        centers=np.array([(x,y) for y in range(80,430,27) for x in range(50,610,27)])+rng.integers(-4,5,(273,2))
        while not self.done.is_set():
            a=np.full((480,640),180,np.uint8)
            for x,y in centers:cv2.circle(a,(int(x+self.sent%12),int(y)),3,35,-1)
            cv2.putText(a,'SYNTHETIC SOFTWARE TEST - NO CAMERA',(15,28),cv2.FONT_HERSHEY_SIMPLEX,.55,30,1)
            p=FramePacket('LIVE_CAMERA','SYNTHETIC SOFTWARE TEST - NO CAMERA','synthetic-check',self.sent,None,a,(640,480))
            self.emit(p);self.sent+=1;self.done.wait(.10)
    def stop(self):self.done.set();self.thread.join(2);self.status('STOPPED')

out=Path(__file__).resolve().parent/f'attempt_{time.time_ns()}';out.mkdir()
cfg=json.loads((ROOT/'configs/default.json').read_text());cfg.update(source_mode='LIVE_CAMERA',backend='B4',live_geometry=None)
app=W.QApplication([]);w=Window(cfg,out/'records');w.setWindowTitle('SOFTWARE UI CHECK — synthetic input, no camera')
w.pipe=Pipeline(cfg,out/'records'/'session_synthetic').start(SyntheticSource);w.show()
def until(test,seconds=15):
    start=time.monotonic()
    while not test() and time.monotonic()-start<seconds:app.processEvents();time.sleep(.02)
    assert test(),(w.pipe.error,w.pipe.snapshot()['status'])
def draw(points):
    w.draw_mask();w.tick();app.processEvents();x,y,sw,sh,iw,ih=w.canvas.image_rect
    for px,py in points:QTest.mouseClick(w.canvas,C.Qt.LeftButton,pos=C.QPoint(round(x+px*sw/iw),round(y+py*sh/ih)))
    QTest.keyClick(w.canvas,C.Qt.Key_Return)
try:
    until(lambda:w.pipe.raw is not None)
    draw([[20,45],[620,45],[620,460],[20,460]])
    until(lambda:w.pipe.snapshot()['valid'])
    s=w.pipe.snapshot();first=dict(candidates=len(s['state']['candidates']),vectors=s['state']['valid_windows'],reference=s['state']['motion_reference'],geometry_revision=s['state']['geometry_revision'])
    w.tick();app.processEvents();w.grab().save(str(out/'mask_applied.png'))
    draw([[320,45],[620,45],[620,460],[320,460]])
    until(lambda:w.pipe.snapshot()['valid'] and w.pipe.state['geometry_revision']==2)
    s=w.pipe.snapshot();second=dict(candidates=len(s['state']['candidates']),vectors=s['state']['valid_windows'],geometry_revision=s['state']['geometry_revision'])
    assert second['candidates']<first['candidates']
    w.clear_mask();assert not w.pipe.snapshot()['valid'] and w.pipe.state is None
    w.tick();app.processEvents();w.grab().save(str(out/'mask_cleared.png'))
    assert w.pipe.source.starts==1
    report=dict(source='SYNTHETIC SOFTWARE TEST - NO CAMERA',backend='B4',first_mask=first,replacement_mask=second,clear_immediate_valid=False,camera_source_start_calls=w.pipe.source.starts,passed=True)
    (out/'software_check.json').write_text(json.dumps(report,indent=2));print(json.dumps(report))
finally:w.close()
