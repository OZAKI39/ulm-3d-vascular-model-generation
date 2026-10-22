import os
os.environ.setdefault('QT_QPA_PLATFORM','offscreen')
import json
import sys
import threading
import time
from pathlib import Path

import av
import numpy as np
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from experiment_video import ExperimentVideoRecorder,DEFAULT_LAYERS
from live_geometry import new_geometry
from protocol import FramePacket


def packet(i=0,arrival=None,session='synthetic',shape=(192,256)):
    a=np.full(shape,120,np.uint8);a[90:150,110:160]=40
    return FramePacket('LIVE_CAMERA','SYNTHETIC TEST - NO CAMERA',session,i,None,a,
        (shape[1],shape[0]),arrival_monotonic=time.monotonic() if arrival is None else arrival,
        host_decode_timestamp={'clock':'windows_monotonic','value':987654,'meaning':'decode, not exposure'})


def geometry():return new_geometry([[10,60],[240,60],[240,175],[10,175]],[256,192])


def state(p,revision=1):
    return dict(session_id=p.session_id,frame_id=p.frame_id,arrival_monotonic=p.arrival_monotonic,
        completed_monotonic=p.arrival_monotonic+.01,geometry_revision=revision,geometry=geometry(),
        candidates=np.array([[130.,120.]]),vectors=np.array([[130.,120.,20.,0.]]),aggregation=None,
        valid=True,motion_reference='IMAGE_SPACE_UNREGISTERED')


def decode(path):
    with av.open(str(path)) as container:
        frames=list(container.decode(video=0))
        return [f.to_ndarray(format='bgr24') for f in frames],[float(f.pts*f.time_base) for f in frames]


def test_timestamps_exact_past_stale_revision_and_no_future_overlay(tmp_path):
    start=time.monotonic();r=ExperimentVideoRecorder(tmp_path,dict(ttl_ms=80,video_annotation_wait_ms=0),start)
    packets=[packet(i,start+t,session='different' if i==4 else 'synthetic')
             for i,t in enumerate([0.,.037,.100,.251,.310])]
    r.add_analysis(state(packets[0]));r.add_analysis(state(packets[3])) # Future state must not annotate earlier frames.
    for i,p in enumerate(packets):assert r.submit(p,2 if i==3 else 1,geometry(),DEFAULT_LAYERS)
    report=r.close();assert report['video_complete'] and report['written_frames']==5
    assert r.close()['written_frames']==5
    rows=[json.loads(line) for line in (tmp_path/'video_frames.jsonl').read_text().splitlines()]
    assert [row['annotation_match'] for row in rows]==['EXACT_FRAME','RECENT_PAST_FRAME',*['NO_CURRENT_ANALYSIS']*3]
    assert rows[1]['analysis_frame_id']==0 and rows[1]['annotation_lag_ms']==pytest.approx(37)
    assert rows[0]['capture_timestamp'] is None and 'UNKNOWN' in rows[0]['capture_timestamp_status']
    assert rows[0]['host_decode_timestamp']['value']==987654
    raw,pts=decode(tmp_path/'experiment_raw.mp4');overlay,op=decode(tmp_path/'experiment_overlay.mp4')
    assert len(raw)==len(overlay)==len(rows)==5 and pts==op
    assert pts==pytest.approx([0.,.037,.100,.251,.310],abs=2e-6)
    assert all(a.shape==(192,256,3) for a in raw) # Raw camera pixels remain unchanged.
    assert all(a.shape==(720,736,3) for a in overlay)
    camera_overlay=[a[264:456,:256] for a in overlay]
    assert np.abs(raw[0].astype(float)-np.repeat(packets[0].image[...,None],3,axis=2)).mean()<2
    assert np.abs(raw[0][110:130,120:140].astype(float)-camera_overlay[0][110:130,120:140]).mean()>2
    assert np.abs(raw[2][110:130,120:140].astype(float)-camera_overlay[2][110:130,120:140]).mean()<2
    assert report['parts'][0]['overlay_camera_xy']==[0,264]


def test_resolution_changes_split_files_without_losing_timeline(tmp_path):
    start=time.monotonic();r=ExperimentVideoRecorder(tmp_path,dict(video_annotation_wait_ms=0),start)
    for p in [packet(99,start),packet(0,start+.5,session='reconnected',shape=(193,257))]:
        assert r.submit(p,0,None,DEFAULT_LAYERS)
    report=r.close();assert len(report['parts'])==2 and report['video_complete']
    assert report['parts'][1]['encoded_wh']==[258,194]
    assert report['parts'][1]['start_elapsed_s']==pytest.approx(.5)
    for part in report['parts']:
        for name in [part['raw'],part['overlay']]:assert len(decode(tmp_path/name)[0])==1
    rows=[json.loads(line) for line in (tmp_path/'video_frames.jsonl').read_text().splitlines()]
    assert rows[1]['source_frame_id']==0 and rows[1]['source_session_id']=='reconnected'
    assert rows[1]['elapsed_since_first_frame_s']==pytest.approx(.5)
    assert rows[1]['video_pts_us']==0


def test_no_frames_does_not_manufacture_video(tmp_path):
    r=ExperimentVideoRecorder(tmp_path,{},time.monotonic());result=r.close()
    assert result['state']=='NO_FRAMES' and not list(tmp_path.glob('*.mp4'))
    assert not result['video_complete']


def test_source_stop_error_still_finalizes_received_video(tmp_path):
    from pipeline import Pipeline
    class Source:
        def __init__(self,emit,status):self.emit=emit;self.status=status
        def start(self):self.status('STREAMING');self.emit(packet());return self
        def stop(self):raise RuntimeError('simulated camera close error')
    config=dict(source_mode='LIVE_CAMERA',backend='B1',ttl_ms=3000,max_pts_gap_s=.25,
        record_experiment_video=True,video_annotation_wait_ms=0)
    p=Pipeline(config,tmp_path).start(Source)
    result=p.stop()
    assert 'simulated camera close error' in result['error']
    assert result['video_recording']['state']=='SAVED'
    assert len(decode(tmp_path/'experiment_raw.mp4')[0])==1
    assert len(decode(tmp_path/'experiment_overlay.mp4')[0])==1
    assert not p.video_record.thread.is_alive()


def test_encoder_error_is_exposed_and_recording_not_claimed_successful(tmp_path,monkeypatch):
    import experiment_video
    def failed(*args,**kwargs):raise OSError('simulated disk full')
    monkeypatch.setattr(experiment_video,'VideoPair',failed)
    r=ExperimentVideoRecorder(tmp_path,dict(video_annotation_wait_ms=0),time.monotonic())
    r.submit(packet(),0,None,DEFAULT_LAYERS);report=r.close()
    assert report['state']=='ERROR' and 'simulated disk full' in report['error']
    assert not report['video_complete'] and report['written_frames']==0
    assert json.loads((tmp_path/'video_manifest.json').read_text())['state']=='ERROR'


def test_queue_overflow_is_counted_without_blocking_capture(tmp_path,monkeypatch):
    entered=threading.Event();release=threading.Event()
    r=ExperimentVideoRecorder(tmp_path,dict(video_queue=1,video_annotation_wait_ms=0),time.monotonic())
    original=r.choose_state
    def slow(item):entered.set();assert release.wait(3);return original(item)
    monkeypatch.setattr(r,'choose_state',slow)
    try:
        assert r.submit(packet(0),0,None,DEFAULT_LAYERS);assert entered.wait(2)
        assert r.submit(packet(1),0,None,DEFAULT_LAYERS)
        start=time.monotonic();assert not r.submit(packet(2),0,None,DEFAULT_LAYERS)
        assert time.monotonic()-start<.1
    finally:release.set();report=r.close()
    assert report['written_frames']==2 and report['dropped_frames']==1 and not report['video_complete']
    assert report['dropped_frame_examples'][0]['frame_id']==2


def test_start_stop_buttons_save_camera_and_pzt_video_and_open_session(tmp_path,monkeypatch):
    from PySide6.QtWidgets import QApplication
    from ui import Window
    import ui
    from pipeline import Pipeline
    from ycommon import ROOT
    config=json.loads((ROOT/'configs/default.json').read_text());config.update(source_mode='LIVE_CAMERA',backend='B1')
    class Observer:
        def __init__(self,cfg):pass
        def analyze(self,p,reset=False):
            result=state(p,0);result.update(timings={'analysis_ms':1.},reason='OK',valid_windows=1)
            return result
    class Source:
        def __init__(self,emit,status):self.emit=emit;self.status=status;self.done=threading.Event()
        def start(self):self.thread=threading.Thread(target=self.run);self.thread.start();return self
        def run(self):
            self.status('STREAMING');i=0
            while not self.done.is_set():self.emit(packet(i));i+=1;self.done.wait(.04)
        def stop(self):self.done.set();self.thread.join(2);self.status('STOPPED')
    original_start=Pipeline.start
    monkeypatch.setattr(Pipeline,'start',lambda self:original_start(self,Source))
    monkeypatch.setattr(ui,'Pipeline',lambda cfg,path:Pipeline(cfg,path,Observer))
    opened=[];monkeypatch.setattr(ui,'open_directory',lambda path:(opened.append(path) or True,'ok'))
    app=QApplication.instance() or QApplication([]);w=Window(config,tmp_path);w.show()
    def wait_record():
        until=time.monotonic()+4
        while w.pipe.video_record.written<3 and time.monotonic()<until:app.processEvents();time.sleep(.02)
        assert w.pipe.video_record.written>=3
    try:
        w.buttons['Start'].click();wait_record();w.buttons['Stop'].click()
        first=w.pipe.output;summary=json.loads((first/'video_manifest.json').read_text())
        assert summary['state']=='SAVED' and 'Videos saved' in w.video_status.text()
        w.buttons['Open records'].click();assert opened==[first]
        for name in ['experiment_raw.mp4','experiment_overlay.mp4']:
            frames,_=decode(first/name);assert len(frames)==summary['written_frames']
            assert frames[0].shape==((192,256,3) if name=='experiment_raw.mp4' else (720,736,3))
        assert (first/'pzt_readbacks.jsonl').exists()
        w.buttons['Start'].click();wait_record();second=w.pipe.output
        assert first!=second
        w.close() # Closing a running window must also finalize both files.
        assert json.loads((second/'video_manifest.json').read_text())['state']=='SAVED'
        assert all(not t.is_alive() for t in [w.pipe.worker,w.pipe.video_record.thread,w.pipe.record.thread])
    finally:w.close()
