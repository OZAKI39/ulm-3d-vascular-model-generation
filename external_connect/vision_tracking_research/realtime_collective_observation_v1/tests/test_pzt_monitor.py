import copy
import json
import os
from pathlib import Path
import sys
import time

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'src'))
import pytest
from pzt_status import NAMES, EXPERIMENT_MODEL, ENDPOINTS, experiment_layout, endpoint_view, physical_views, channel_view, read_status, save_layout, validate_layout, StatusReader


def example_status(now=100.):
    return dict(schema=1, source='AFG_GUI_READBACK', controller='CONNECTED', published_unix_s=now,
        channels={name: dict(device='A' if i < 2 else 'B', channel=i%2+1, serial='2147075' if i < 2 else '2439003', state='READ_OK',
            read_unix_s=now, parameters_unix_s=now, output=(i == 0),
            parameters=dict(frequency_hz=2.5e6, amplitude_vpp=.1, phase_deg=90., waveform='SIN', frequency_mode='CW'))
                  for i, name in enumerate(NAMES)})


def test_status_expiry_disconnect_and_clock_error():
    data = example_status()
    assert channel_view(data, 'PZT1', 101.)['output'] == 'ON'
    for now in (109., 90.):
        view = channel_view(data, 'PZT1', now)
        assert view['output'] == 'UNKNOWN' and not view['parameters']
    data['published_unix_s'] = 110.
    assert channel_view(data, 'PZT1', 110.)['connection'] == 'Readback expired'
    data['controller'] = 'CLOSED'
    assert channel_view(data, 'PZT1', 110.)['connection'] == 'Controller closed'


def test_state_and_parameters_have_independent_freshness():
    data = example_status()
    data['published_unix_s'] = data['channels']['PZT1']['read_unix_s'] = 110.
    data['channels']['PZT1']['output'] = False
    view = channel_view(data, 'PZT1', 110.)
    assert view['output'] == 'OFF' and not view['parameters']
    data['channels']['PZT2']['state'] = 'READ_ERROR'
    assert channel_view(data, 'PZT2', 110.)['connection'] == 'Read failed'


@pytest.mark.parametrize('bad', ['{bad', '[]', '{"schema":9}', '{"schema":1,"source":"AFG_GUI_READBACK","channels":[]}'])
def test_missing_or_malformed_file_never_displays_live_on(tmp_path, bad):
    path = tmp_path/'status.json'
    assert channel_view(read_status(path), 'PZT1')['output'] == 'UNKNOWN'
    path.write_text(bad)
    assert channel_view(read_status(path), 'PZT1')['output'] == 'UNKNOWN'


def test_invalid_parameters_are_not_presented(tmp_path):
    data = example_status()
    data['channels']['PZT1']['parameters']['frequency_hz'] = float('nan')
    assert channel_view(data, 'PZT1', 101.)['connection'] == 'Invalid readback'
    data['channels']['PZT1'] = 'ON'
    path = tmp_path/'status.json'
    path.write_text(json.dumps(data))
    assert 'error' in read_status(path)


def test_experiment_geometry_and_channel_assignment_persistence(tmp_path):
    data=experiment_layout()
    assert all(row['pzt'] is None for row in data['assignments'].values())
    assert EXPERIMENT_MODEL['block_mm']==[32.5,9,9]
    assert EXPERIMENT_MODEL['pzts']['X-']['center_mm']==[-16.25,0,0]
    assert EXPERIMENT_MODEL['pzts']['Y+']['center_mm']==[0,4.5,0]
    assert all(row['diameter_mm']==7 for row in EXPERIMENT_MODEL['pzts'].values())
    data['assignments']['A/CH1']['pzt']='Y+'
    path=tmp_path/'layout.json'
    saved=save_layout(path,data)
    assert saved==validate_layout(json.loads(path.read_text()))
    data['assignments']['B/CH2']['pzt']='Y+'
    with pytest.raises(ValueError,match='more than once'):save_layout(path,data)
    assert json.loads(path.read_text())==saved


def test_mapping_routes_the_correct_real_channel_and_rejects_wrong_identity():
    data=example_status()
    layout=experiment_layout()
    layout['assignments']['A/CH1']['pzt']='Y+'
    layout['assignments']['B/CH2']['pzt']='X-'
    views=physical_views(data,layout,101.)
    assert views['Y+']['output']=='ON' and views['Y+']['endpoint']=='A/CH1'
    assert views['X-']['output']=='OFF' and views['X-']['endpoint']=='B/CH2'
    assert views['X+']['output']=='Unassigned' and not views['X+']['parameters']
    # Do not route by legacy labels PZT1/PZT4 or by enumeration order.
    data['channels']['PZT1'],data['channels']['PZT4']=data['channels']['PZT4'],data['channels']['PZT1']
    assert physical_views(data,layout,101.)['Y+']['output']=='ON'
    data['channels']['PZT4']['serial']='WRONG_SERIAL'
    invalid=physical_views(data,layout,101.)['Y+']
    assert invalid['output']=='UNKNOWN' and invalid['connection']=='Identity mismatch'
    assert not invalid['parameters']


def test_mapping_cannot_guess_legacy_geometry_or_change_device_identity():
    with pytest.raises(ValueError):validate_layout({'schema':1,'pzts':{}})
    data=experiment_layout();data['assignments']['A/CH1']['serial']='WRONG'
    with pytest.raises(ValueError,match='identity'):validate_layout(data)
    data=experiment_layout();data['assignments']['A/CH1']['pzt']='PZT1'
    with pytest.raises(ValueError):validate_layout(data)


def test_worker_recovers_from_missing_status_file(tmp_path):
    path = tmp_path/'status.json'
    reader = StatusReader(path)
    try:
        path.write_text(json.dumps(example_status(time.time())))
        deadline = time.monotonic()+2
        while reader.snapshot().get('controller') != 'CONNECTED' and time.monotonic() < deadline:
            time.sleep(.02)
        assert channel_view(reader.snapshot(), 'PZT1')['output'] == 'ON'
        path.unlink()
        deadline = time.monotonic()+2
        while not reader.snapshot().get('error') and time.monotonic() < deadline:
            time.sleep(.02)
        assert channel_view(reader.snapshot(), 'PZT1')['output'] == 'UNKNOWN'
    finally:
        reader.close()
    assert not reader.thread.is_alive()


def test_connection_editor_alignment_duplicate_check_and_monitor(tmp_path):
    from PySide6 import QtWidgets as W
    from pzt_panel import ConnectionEditor,PZTPanel
    app=W.QApplication.instance() or W.QApplication([])
    class Reader:
        def __init__(self,path):pass
        def snapshot(self):return example_status(time.time())
        def close(self):pass
    panel=PZTPanel(dict(pzt_layout_path=str(tmp_path/'layout.json')),reader_factory=Reader)
    panel.resize(440,600);panel.show()
    editor=ConnectionEditor(experiment_layout());editor.show();app.processEvents()
    try:
        combos=list(editor.fields.values())
        assert len({c.geometry().x() for c in combos})==1
        assert len({c.width() for c in combos})==1
        combos[0].setCurrentIndex(combos[0].findData('Y+'))
        assert editor.value()['assignments']['A/CH1']['pzt']=='Y+'
        combos[1].setCurrentIndex(combos[1].findData('Y+'))
        assert not editor.save_button.isEnabled() and 'more than once' in editor.error.text()
        combos[1].setCurrentIndex(0)
        assert editor.save_button.isEnabled()
        panel.layout_data=editor.value();panel.refresh()
        assert 'A/CH1 → Y+'==panel.cards['A/CH1'][0].title()
        assert 'AFG connected | ON' in panel.cards['A/CH1'][1].text()
        assert panel.top.views['Y+']['output']=='ON'
        assert panel.top.views['X-']['output']=='Unassigned'
        assert '0.1 Vpp' in panel.cards['A/CH1'][1].text()
        assert not panel.hint.isVisible()
        assert not panel.findChildren(W.QDoubleSpinBox)
        assert not editor.findChildren(W.QDoubleSpinBox)
        assert not any(b.text() in ('ON','OFF','Apply') for b in panel.findChildren(W.QPushButton))
    finally:
        editor.close();panel.shutdown();panel.close()


def test_last_saved_assignments_are_loaded_on_every_start(tmp_path,monkeypatch):
    from PySide6 import QtWidgets as W
    from pzt_panel import PZTPanel
    app=W.QApplication.instance() or W.QApplication([])
    saved=experiment_layout()
    for key,target in zip(ENDPOINTS,['Y-','X-','Y+','X+']):saved['assignments'][key]['pzt']=target
    path=tmp_path/'saved_layout.json';save_layout(path,saved);original=path.read_bytes()
    class Reader:
        def __init__(self,path):pass
        def snapshot(self):return example_status(time.time())
        def close(self):pass
    monkeypatch.chdir(tmp_path.parent)
    for _ in range(2):
        panel=PZTPanel(dict(pzt_layout_path=str(path)),reader_factory=Reader)
        assert panel.layout_data==saved
        assert panel.cards['A/CH1'][0].title()=='A/CH1 → Y-'
        assert panel.top.views['Y-']['output']=='ON'
        panel.shutdown();panel.close()
    assert path.read_bytes()==original


def test_read_button_sends_new_request_and_waits_for_matching_completion(tmp_path):
    from PySide6 import QtWidgets as W
    from pzt_panel import PZTPanel
    app=W.QApplication.instance() or W.QApplication([])
    class Reader:
        def __init__(self,path):self.data=example_status(time.time());self.calls=0
        def snapshot(self):return self.data
        def close(self):pass
        def request_refresh(self):self.calls+=1;return 'a'*32
    panel=PZTPanel(dict(pzt_layout_path=str(tmp_path/'layout.json')),reader_factory=Reader)
    try:
        panel.read_button.click()
        assert panel.reader.calls==1 and not panel.read_button.isEnabled()
        panel.reader.data['refresh']=dict(request_id='old',state='COMPLETE')
        panel.refresh();assert panel.pending_request
        panel.reader.data['refresh']=dict(request_id='a'*32,state='COMPLETE')
        panel.refresh();assert not panel.pending_request and panel.read_button.isEnabled()
        assert panel.hint.text()=='Status updated.'
    finally:panel.shutdown();panel.close()


def test_status_worker_posts_only_read_command(tmp_path):
    path=tmp_path/'status.json';data=example_status(time.time())
    data.update(session='current-session',refresh_supported=True)
    path.write_text(json.dumps(data))
    reader=StatusReader(path)
    try:
        deadline=time.monotonic()+2
        while not reader.snapshot().get('session') and time.monotonic()<deadline:time.sleep(.01)
        rid=reader.request_refresh();request_path=tmp_path/'observation_read_request.json'
        while not request_path.exists() and time.monotonic()<deadline:time.sleep(.01)
        request=json.loads(request_path.read_text())
        assert request['request_id']==rid and request['session']=='current-session'
        assert request['action']=='READ_STATUS'
    finally:reader.close()


def test_black_background_pixels_and_compact_controls(tmp_path):
    from PySide6 import QtWidgets as W
    from ui import Window
    from ycommon import ROOT
    app=W.QApplication.instance() or W.QApplication([])
    cfg=json.loads((ROOT/'configs/default.json').read_text());cfg.update(source_mode='LIVE_CAMERA',pzt_layout_path=str(tmp_path/'layout.json'),pzt_status_path=str(tmp_path/'status.json'))
    window=Window(cfg,tmp_path/'ui');window.show();app.processEvents()
    try:
        for widget,xy in [(window.canvas,(2,2)),(window.curves,(2,30)),(window.pzt_panel.top,(2,35))]:
            color=widget.grab().toImage().pixelColor(*xy)
            assert (color.red(),color.green(),color.blue())==(0,0,0)
        assert window.pzt_panel.read_button.width()==window.pzt_panel.edit_button.width()
    finally:window.close()
