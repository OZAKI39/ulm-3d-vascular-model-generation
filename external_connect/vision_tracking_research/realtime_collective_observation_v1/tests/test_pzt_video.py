"""Timing, provenance and encoding tests with explicitly synthetic AFG values."""
from collections import deque
import copy
import json
import threading
import time

import numpy as np

from test_experiment_video import packet, decode
from test_pzt_monitor import example_status
from pzt_status import StatusReader, experiment_layout, ENDPOINTS
from pzt_video import PZTVideoSource, frame_views
from experiment_video import ExperimentVideoRecorder, DEFAULT_LAYERS


def source():
    reader = StatusReader.__new__(StatusReader)
    reader.lock = threading.Lock()
    first = example_status(1000.)
    second = example_status(1002.)
    second['channels']['PZT1']['output'] = False
    reader.history = deque([dict(data=first,observed_monotonic=99.,observed_unix_s=1000.),
                            dict(data=second,observed_monotonic=101.,observed_unix_s=1002.)])
    layout = experiment_layout()
    for key,target in zip(ENDPOINTS,['Y-','X-','Y+','X+']):layout['assignments'][key]['pzt']=target
    context = PZTVideoSource(reader,layout)
    context.layouts = deque([(90.,0,layout)])
    return context


def test_future_readback_mapping_edits_expiry_and_replay_are_not_misattributed():
    s = source()
    old = s.capture(packet(arrival=100.))
    assert old['sample']['observation']['observed_monotonic']==99.
    assert frame_views(old)['A/CH1']['output']=='ON'
    new = s.capture(packet(arrival=102.))
    assert frame_views(new)['A/CH1']['output']=='OFF'
    expired = s.capture(packet(arrival=111.))
    assert frame_views(expired)['A/CH1']['output']=='UNKNOWN'
    assert not frame_views(expired)['A/CH1']['parameters']
    layout = copy.deepcopy(s.layouts[0][2])
    layout['assignments']['A/CH1']['pzt'],layout['assignments']['A/CH2']['pzt']='X-','Y-'
    s.layouts.append((103.,1,layout))
    assert s.capture(packet(arrival=102.))['sample']['layout_revision']==0
    assert s.capture(packet(arrival=104.))['sample']['layout_revision']==1
    replay = packet(arrival=104.);replay.source_mode='REALTIME_REPLAY'
    context = s.capture(replay)
    assert context['sample']['scope']=='REPLAY_NO_SYNCHRONIZED_PZT'
    assert context['sample']['observation'] is None
    assert not frame_views(context)['A/CH1']['parameters']
    assert s.capture(packet(arrival=98.))['sample']['observation'] is None
    # Snapshots cannot be rewritten by later monitor updates.
    s.reader.history[0]['data']['channels']['PZT1']['output']=False
    assert frame_views(old)['A/CH1']['output']=='ON'


def test_video_has_timestamp_matched_readbacks_and_keeps_raw_pixels(tmp_path):
    s=source();start=time.monotonic()
    r=ExperimentVideoRecorder(tmp_path,dict(video_annotation_wait_ms=0),start)
    contexts=[s.capture(packet(arrival=t)) for t in [100.,100.04,102.,111.]]
    for i,context in enumerate(contexts):
        assert r.submit(packet(i,start+i*.04),0,None,DEFAULT_LAYERS,context)
    report=r.close()
    assert report['video_complete'],report['error']
    rows=[json.loads(line) for line in (tmp_path/'video_frames.jsonl').read_text().splitlines()]
    samples=[json.loads(line) for line in (tmp_path/'pzt_readbacks.jsonl').read_text().splitlines()]
    assert len(samples)==2 # Same readback reused and evaluated for expiry per frame.
    by_id={row['sample_id']:row for row in samples}
    recovered=[frame_views(dict(frame_unix_s=row['pzt_frame_unix_s'],sample=by_id[row['pzt_sample_id']])) for row in rows]
    assert [v['A/CH1']['output'] for v in recovered]==['ON','ON','OFF','UNKNOWN']
    assert by_id[rows[0]['pzt_sample_id']]['layout']['assignments']['A/CH1']['pzt']=='Y-'
    raw,t=decode(tmp_path/'experiment_raw.mp4');over,ot=decode(tmp_path/'experiment_overlay.mp4')
    assert t==ot and len(raw)==len(over)==4
    assert raw[0].shape==(192,256,3) and over[0].shape==(720,736,3)
    # AFG state changes appear in the sidebar, not the untouched raw recording.
    # Inspect the changed ON/OFF text, not the mostly black entire sidebar.
    change=np.abs(over[0][98:118,280:370].astype(float)-over[2][98:118,280:370])
    assert np.count_nonzero(change>20)>50
    original=np.repeat(packet().image[...,None],3,axis=2).astype(float)
    assert all(np.abs(a.astype(float)-original).mean()<2 for a in raw) # H.264 is lossy.
    assert json.loads((tmp_path/'pzt_experiment_model.json').read_text())['block_mm']==[32.5,9,9]
