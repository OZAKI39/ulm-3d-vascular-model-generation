"""Interactive mask updates must affect a running worker, not just UI config."""
import os
os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
import json
import sys
import threading
import time
from pathlib import Path
from types import SimpleNamespace

import cv2
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from live_geometry import new_geometry, validate_geometry
from ycommon import ROOT, dump
from protocol import FramePacket
from pipeline import Pipeline


def cfg():
    result = json.loads((ROOT/'configs/default.json').read_text())
    result.update(source_mode='LIVE_CAMERA', backend='B1', ttl_ms=3000,record_experiment_video=False)
    return result


def image_packet(frame_id=0):
    a = np.full((480,640), 180, np.uint8)
    for y in range(60,430,36):
        for x in range(60,600,36):
            cv2.circle(a,(x+frame_id,y),3,30,-1)
    return FramePacket('LIVE_CAMERA','synthetic software test','test',frame_id,None,
        a,(640,480),arrival_monotonic=time.monotonic()+frame_id*.05)


def polygon():
    return [[20,20],[620,20],[620,460],[20,460]]


@pytest.mark.parametrize('points', [[], [[1,1],[2,2],[3,3]],
    [[-1,20],[200,20],[200,200]], [[20,20],[400,400],[20,400],[400,20]],
    [[20,20],[640,20],[500,200]]])
def test_bad_masks_are_rejected(points):
    with pytest.raises(ValueError):
        new_geometry(points,[640,480])


def test_old_mask_without_reference_loads_as_explicit_image_space():
    g = new_geometry(polygon(),[640,480])
    del g['motion_reference']
    assert validate_geometry(g,[640,480])['motion_reference']=='IMAGE_SPACE_UNREGISTERED'
    with pytest.raises(ValueError):validate_geometry(g,[1280,1024])


def test_image_motion_does_not_claim_verified_reference(tmp_path):
    from observer import Observer
    config = cfg();path=tmp_path/'geometry.json'
    g=new_geometry(polygon(),[640,480]);dump(path,g);config['live_geometry']=str(path)
    obs=Observer(config)
    class Backend:
        def run(self,*args,**kwargs):
            shape=obs.channel.safe.shape
            return dict(valid=np.ones(shape,bool),u=np.ones(shape),v=np.zeros(shape),total_ms=1.)
    obs.backend=Backend()
    first=obs.analyze(image_packet(),True)
    assert not first['valid']  # A mask cannot manufacture a second image.
    r=obs.analyze(image_packet(1))
    assert r['valid'] and r['valid_windows']>=8
    assert not r['registration']['valid'] and not r['motion_reference_verified']
    assert r['motion_reference']=='IMAGE_SPACE_UNREGISTERED'
    # Old external callers keep the strict reference requirement unless they
    # explicitly select image-space observation through the new mask workflow.
    del g['motion_reference'];dump(path,g)
    obs.analyze(image_packet(2));r=obs.analyze(image_packet(3))
    assert not r['valid'] and r['valid_windows']==0
    assert r['registration']['method']=='LIVE_REFERENCE_NOT_VERIFIED'


def test_failed_configured_wall_reference_still_rejects_motion(tmp_path):
    from observer import Observer
    config=cfg();g=new_geometry(polygon(),[640,480])
    g['registration_roi']=[0,0,20,20] # Flat patch: no trustworthy registration.
    g=validate_geometry(g,[640,480]);path=tmp_path/'g.json';dump(path,g);config['live_geometry']=str(path)
    obs=Observer(config)
    class Backend:
        def run(self,*args,**kwargs):
            s=obs.channel.safe.shape
            return dict(valid=np.ones(s,bool),u=np.ones(s),v=np.zeros(s),total_ms=1.)
    obs.backend=Backend();obs.analyze(image_packet(),True);r=obs.analyze(image_packet(1))
    assert not r['valid'] and r['valid_windows']==0


def test_clear_discards_inflight_result_without_restarting_source(tmp_path):
    entered=threading.Event();release=threading.Event()
    class Observer:
        def __init__(self,config):self.config=config
        def analyze(self,p,reset=False):
            if p.frame_id==0:entered.set();assert release.wait(4)
            return dict(valid=bool(self.config.get('live_geometry')),reason='MASK_NOT_VERIFIED',
                candidates=[],vectors=[],aggregation=None,geometry=None,timings={'analysis_ms':1.})
    class Source:
        def __init__(self,emit,status):self.emit=emit;self.status=status;self.starts=0
        def start(self):self.starts+=1;self.status('STREAMING');self.emit(image_packet());return self
        def stop(self):self.status('STOPPED')
    p=Pipeline(cfg(),tmp_path,Observer).start(Source)
    try:
        assert entered.wait(3)
        dump(tmp_path/'first.json',new_geometry(polygon(),[640,480]))
        p.set_live_geometry(tmp_path/'first.json')
        p.set_live_geometry(None)
        assert not p.snapshot()['valid'] and p.state is None
        release.set()
        deadline=time.monotonic()+3
        while not p.counter['discarded_old_geometry'] and time.monotonic()<deadline:time.sleep(.01)
        assert p.counter['discarded_old_geometry']==1 and p.state is None
        assert p.source.starts==1 and p.worker.is_alive()
        p.receive(image_packet(1))
        deadline=time.monotonic()+3
        while p.state is None and time.monotonic()<deadline:time.sleep(.01)
        assert p.state['geometry_revision']==2 and not p.state['valid']
        assert len(json.loads((tmp_path/'geometry_changes.json').read_text()))==2
    finally:release.set();p.stop()


@pytest.mark.parametrize('finish', ['enter','double_click','first_point'])
def test_draw_applies_clear_removes_and_history_survives(tmp_path,finish):
    from PySide6 import QtCore as C
    from PySide6.QtWidgets import QApplication
    from PySide6.QtTest import QTest
    from ui import Window
    app=QApplication.instance() or QApplication([])
    w=Window(cfg(),tmp_path);w.pipe=Pipeline(cfg(),tmp_path/'session_test')
    w.pipe.raw=image_packet();w.pipe.status='STREAMING';w.show();w.tick();app.processEvents()
    try:
        w.geometry_buttons['Draw LIVE mask'].click();app.processEvents()
        x,y,sw,sh,iw,ih=w.canvas.image_rect
        def qpoint(point):return C.QPoint(round(x+point[0]*sw/iw),round(y+point[1]*sh/ih))
        for point in polygon()[:3]:QTest.mouseClick(w.canvas,C.Qt.LeftButton,pos=qpoint(point))
        if finish=='double_click':QTest.mouseDClick(w.canvas,C.Qt.LeftButton,pos=qpoint(polygon()[3]))
        else:
            QTest.mouseClick(w.canvas,C.Qt.LeftButton,pos=qpoint(polygon()[3]))
            if finish=='enter':QTest.keyClick(w.canvas,C.Qt.Key_Return)
            else:QTest.mouseClick(w.canvas,C.Qt.LeftButton,pos=qpoint(polygon()[0]))
        app.processEvents()
        path=Path(w.cfg['live_geometry']);assert path.exists()
        assert w.pipe.cfg['live_geometry']==str(path) and w.pipe.geometry_revision==1
        assert not w.canvas.editing and not w.canvas.points
        assert w.canvas.active_geometry['registration_roi'] is None
        # Escape cancels a replacement; Backspace undoes a point.
        w.draw_mask();QTest.mouseClick(w.canvas,C.Qt.LeftButton,pos=qpoint([100,100]))
        QTest.keyClick(w.canvas,C.Qt.Key_Backspace);assert not w.canvas.points
        QTest.keyClick(w.canvas,C.Qt.Key_Escape);assert w.cfg['live_geometry']==str(path)
        w.geometry_buttons['Clear current mask'].click();app.processEvents()
        assert w.cfg['live_geometry'] is None and w.pipe.cfg['live_geometry'] is None
        assert w.canvas.active_geometry is None and not w.pipe.snapshot()['valid']
        assert path.exists()  # Clear does not delete historical evidence.
    finally:w.close()


def test_wsl_opener_passes_literal_unc_path(tmp_path,monkeypatch):
    import desktop
    target=tmp_path/'records with spaces & punctuation';calls=[]
    monkeypatch.setattr(desktop,'is_wsl',lambda:True)
    monkeypatch.setattr(desktop.shutil,'which',lambda name:'/bin/'+name)
    unc=r'\\wsl.localhost\Ubuntu\home\lzy\records with spaces & punctuation'
    def run(args,**kwargs):
        assert args==['/bin/wslpath','-w',str(target.resolve())]
        return SimpleNamespace(stdout=unc+'\n')
    monkeypatch.setattr(desktop.subprocess,'run',run)
    monkeypatch.setattr(desktop.subprocess,'Popen',lambda args,**kwargs:calls.append((args,kwargs)))
    ok,_=desktop.open_directory(target)
    assert ok and target.is_dir() and calls[0][0]==['/bin/explorer.exe',unc]
    assert calls[0][1].get('shell',False) is False


def test_failed_open_returns_actionable_path(tmp_path,monkeypatch):
    import desktop
    from PySide6.QtGui import QDesktopServices
    monkeypatch.setattr(desktop,'is_wsl',lambda:False)
    monkeypatch.setattr(QDesktopServices,'openUrl',lambda url:False)
    ok,message=desktop.open_directory(tmp_path)
    assert not ok and str(tmp_path.resolve()) in message
