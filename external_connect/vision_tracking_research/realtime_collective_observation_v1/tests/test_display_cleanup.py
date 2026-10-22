import os
os.environ.setdefault('QT_QPA_PLATFORM','offscreen')
import gzip
import json
from pathlib import Path
import re
import sys
import time

import numpy as np
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from piv_display import arrow_shape,piv_arrows,GAIN_SECONDS,MAX_LENGTH_PX
from pipeline import Pipeline,observation_record
from protocol import FramePacket
from ycommon import ROOT


def test_arrow_geometry_keeps_measured_direction_and_velocity():
    v=np.array([[24.,24.,25.,0.],[72.,24.,0.,-25.],[120.,24.,300.,400.]])
    original=v.copy();arrows=piv_arrows(v)
    assert len(arrows)==3 and np.array_equal(v,original)
    for (start,end,head),row in zip(arrows,v):
        delta=end-start;velocity=row[2:]
        assert np.dot(delta,velocity)>0
        assert np.allclose(delta/np.linalg.norm(delta),velocity/np.linalg.norm(velocity),atol=1e-14)
        assert np.linalg.norm(delta)==pytest.approx(min(np.linalg.norm(velocity)*GAIN_SECONDS,MAX_LENGTH_PX))
        a,b=head[1]-end,head[2]-end
        assert np.allclose(head[0],end) and abs(a[0]*b[1]-a[1]*b[0])>0
    assert arrow_shape([0,0],[0,0]) is None
    assert arrow_shape([0,0],[1e-12,0]) is None


def test_spatial_sampling_selects_real_vectors_without_averaging():
    v=np.array([[0.,0.,99.,99.],[24.,24.,12.,-3.],[25.,25.,-30.,0.]])
    arrows=piv_arrows(v)
    assert len(arrows)==1
    assert np.allclose(arrows[0][0],[24,24])
    assert np.allclose(arrows[0][1]-arrows[0][0],[12*GAIN_SECONDS,-3*GAIN_SECONDS])


def test_default_records_keep_science_and_remove_profiling(tmp_path):
    cfg=dict(source_mode='LIVE_CAMERA',backend='B1',ttl_ms=3000,max_pts_gap_s=.25,record_experiment_video=False)
    class Observer:
        def __init__(self,cfg):pass
        def analyze(self,p,reset=False):
            return dict(valid=False,reason='MASK_NOT_VERIFIED',candidates=np.array([[20.,20.]]),
                vectors=np.array([[20.,20.,3.,2.]]),aggregation=None,geometry=None,
                timings={'analysis_ms':1.,'scan_ms':.1},computed_windows=100,window_grid_xy=np.ones((3,2)))
    class Source:
        def __init__(self,emit,status):self.emit=emit;self.status=status
        def start(self):
            self.status('STREAMING');self.emit(FramePacket('LIVE_CAMERA','test','s',0,None,np.zeros((96,128),np.uint8),(128,96)));return self
        def stop(self):self.status('STOPPED')
    p=Pipeline(cfg,tmp_path,Observer).start(Source)
    try:
        deadline=time.monotonic()+3
        while p.state is None and time.monotonic()<deadline:time.sleep(.01)
        assert p.state is not None
        assert not hasattr(p,'resource_thread') and not p.samples and not p.timings
    finally:summary=p.stop()
    assert not {'resources','stages','record_ms','ui_render_ms','ui_refresh_hz'}&summary.keys()
    with gzip.open(tmp_path/'states.jsonl.gz','rt') as f:row=json.loads(next(f))
    assert not {'timings','window_grid_xy','computed_windows'}&row.keys()
    assert row['vectors']==[[20.,20.,3.,2.]] and row['candidates']==[[20.,20.]]
    assert 'arrival_monotonic' in row and 'geometry_revision' in row and 'valid' in row
    assert observation_record({'timings':{'analysis_ms':1}},diagnostics=True)['timings']['analysis_ms']==1


def test_english_interface_and_relevant_summary(tmp_path):
    from PySide6 import QtWidgets as W
    from ui import Window,observation_text
    app=W.QApplication.instance() or W.QApplication([])
    cfg=json.loads((ROOT/'configs/default.json').read_text());cfg['source_mode']='LIVE_CAMERA'
    window=Window(cfg,tmp_path);window.show();app.processEvents()
    try:
        assert window.speed.isHidden() and window.loop.isHidden()
        assert window.buttons['Pause'].isHidden() and window.buttons['Resume'].isHidden()
        assert 'Validity' not in window.checks
        assert window.checks['Smooth PIV'].isChecked()
        assert window.event_source.currentText()=='MANUAL_ANNOTATION'
        texts=[item.text() for item in window.findChildren(W.QLabel)]+[item.text() for item in window.findChildren(W.QAbstractButton)]
        assert not any(re.search(r'[\u4e00-\u9fff]',text) for text in texts)
        summary=observation_text(dict(status='FRESH',valid=True,fresh=True,state=dict(
            aggregation={'regions':[dict(count=12,mean_velocity=[12.345,-3.456])],'status':'VISIBLE_CONVERGENCE'})))
        assert '12' in summary and 'x=12.3, y=-3.5 px/s' in summary and 'Visible convergence' in summary
        assert not any(key in summary for key in ['CPU','GPU','RSS','Backend','PTS','Hz','frame:'])
        window.mode.setCurrentText('REALTIME_REPLAY');app.processEvents()
        assert not window.speed.isHidden() and not window.buttons['Pause'].isHidden()
    finally:window.close()
    persisted=json.loads((tmp_path/'ui_summary.json').read_text())
    assert 'canvas_paint_ms' not in persisted and 'ui_state_update_ms' not in persisted
