"""Render-level regressions: no brightness blinking, no stale-data persistence."""
import os
os.environ.setdefault('QT_QPA_PLATFORM','offscreen')
from dataclasses import replace
import copy
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from experiment_video import render_overlay
from overlay_style import DEFAULT_LAYERS,overlay_opacity
from protocol import FramePacket
from ycommon import ROOT


def sample():
    packet=FramePacket('REALTIME_REPLAY','synthetic render regression','session',1,.1,
        np.full((480,640),100,np.uint8),(640,480),arrival_monotonic=10.)
    geometry=dict(image_wh=[640,480],channel_polygon=[[10,10],[630,10],[630,470],[10,470]],
        excluded_polygons=[],crop_xywh=[0,0,640,480])
    state=dict(session_id='session',frame_id=1,arrival_monotonic=10.,geometry=geometry,valid=True,
        candidates=[[300.,240.]],vectors=[[320.,240.,25.,0.]],
        aggregation=dict(clusters=[dict(extent=[250,200,400,300],center=[320,250],rg_px=30)]))
    return packet,state


def render(p,s,lag,layers=None):
    return render_overlay(p,s,s['geometry'],DEFAULT_LAYERS if layers is None else layers,
        'EXACT_FRAME' if lag==0 else 'RECENT_PAST_FRAME',lag,'2026-10-18T00:00:00Z',1.)


def test_recorded_overlay_does_not_blink_between_exact_and_recent():
    p,s=sample();original=copy.deepcopy(s)
    exact=render(p,s,0)
    for lag in [33.,66.,100.,160.,190.,230.,249.]:
        recent=render(replace(p,frame_id=2),s,lag)
        assert np.array_equal(exact[100:400],recent[100:400])
    assert s==original and np.all(p.image==100) # Display cannot mutate observations.


def test_display_keeps_contrast_until_original_expiry_and_removes_invalid_motion():
    p,s=sample()
    strengths=[]
    for age in [170,190,210,230,249,250,300]:
        a=render(p,s,age)[225:255,290:355].astype(float)
        strengths.append(np.abs(a-100).sum())
    assert len(set(strengths[:-2]))==1 and strengths[0]>0
    assert strengths[-2:]==[0,0]
    assert overlay_opacity(float('nan'))==0 and overlay_opacity(-1)==0
    s['valid']=False
    invalid=render(p,s,0)
    assert np.all(invalid[230:250,315:350]==100) # No arrow while motion invalid.
    assert np.any(invalid[230:250,290:310]!=100) # Detection can remain independently valid.


def test_regions_and_shape_default_off_and_can_be_enabled_independently():
    p,s=sample();baseline=render(p,s,0)
    assert not DEFAULT_LAYERS['Regions'] and not DEFAULT_LAYERS['Shape']
    assert np.all(baseline[195:210,260:390]==100)
    regions=render(p,s,0,dict(DEFAULT_LAYERS,Regions=True))
    shape=render(p,s,0,dict(DEFAULT_LAYERS,Shape=True))
    assert np.any(regions[195:210,260:390]!=100)
    assert np.all(shape[195:210,260:390]==100)
    assert np.any(shape[275:285,310:330]!=100)


def test_live_canvas_has_same_contrast_when_new_raw_frame_arrives(tmp_path):
    from PySide6 import QtWidgets as W,QtGui as G
    from ui import Window
    app=W.QApplication.instance() or W.QApplication([])
    cfg=json.loads((ROOT/'configs/default.json').read_text());cfg['source_mode']='LIVE_CAMERA'
    window=Window(cfg,tmp_path);window.timer.stop();window.show();app.processEvents()
    p,s=sample()
    try:
        assert not window.checks['Regions'].isChecked() and not window.checks['Shape'].isChecked()
        assert window.checks['PIV'].isChecked() and window.checks['Candidates'].isChecked()
        window.canvas.layers=DEFAULT_LAYERS.copy()
        def image(packet,age,fresh=True):
            window.canvas.snap=dict(raw=packet,state=s,valid=fresh,fresh=fresh,age_ms=age,ttl_ms=250.,status='FRESH' if fresh else 'STALE')
            window.canvas.update();app.processEvents()
            q=window.canvas.grab().toImage().convertToFormat(G.QImage.Format_RGBA8888)
            return np.frombuffer(q.bits(),np.uint8,count=q.sizeInBytes()).reshape(q.height(),q.bytesPerLine())[:,:q.width()*4].reshape(q.height(),q.width(),4).copy()
        a=image(p,40.);b=image(replace(p,frame_id=2,arrival_monotonic=10.033),73.)
        assert np.array_equal(a[100:-100],b[100:-100])
        # Cached paint must respond to layer changes and new measurements, and
        # must never resurrect an expired observation.
        window.canvas.layers=dict(DEFAULT_LAYERS,Candidates=False,PIV=False)
        hidden=image(p,73.)
        assert not np.array_equal(a[100:-100],hidden[100:-100])
        window.canvas.layers=DEFAULT_LAYERS.copy()
        expired=image(p,260.,fresh=False)
        assert np.array_equal(hidden[100:-100],expired[100:-100])
        s=copy.deepcopy(s);s.update(frame_id=2,candidates=[],vectors=[])
        updated=image(replace(p,frame_id=2),20.)
        assert np.array_equal(hidden[100:-100],updated[100:-100])
    finally:window.close()
