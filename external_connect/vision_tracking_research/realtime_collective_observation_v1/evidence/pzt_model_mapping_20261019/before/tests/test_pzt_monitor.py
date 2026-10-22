import copy
import json
import os
from pathlib import Path
import sys
import time

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'src'))
import pytest
from pzt_status import NAMES, blank_layout, channel_view, read_status, save_layout, validate_layout, StatusReader


def example_status(now=100.):
    return dict(schema=1, source='AFG_GUI_READBACK', controller='CONNECTED', published_unix_s=now,
        channels={name: dict(device='A' if i < 2 else 'B', channel=i%2+1, state='READ_OK',
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


def test_manual_layout_is_blank_until_enabled_and_persists(tmp_path):
    data = blank_layout()
    assert not any(r['placed'] for r in data['pzts'].values())
    data['pzts']['PZT1'].update(placed=True, x_mm=-16.25, y_mm=0, normal_deg=-360, wiring='Connected')
    path = tmp_path/'layout.json'
    saved = save_layout(path, data)
    assert saved == validate_layout(json.loads(path.read_text()))
    assert saved['pzts']['PZT1']['normal_deg'] == 0
    assert saved['pzts']['PZT1']['x_mm'] == -16.25
    data['pzts']['PZT2']['x_mm'] = float('nan')
    with pytest.raises(ValueError):
        save_layout(path, data)
    assert json.loads(path.read_text()) == saved


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


def test_layout_editor_alignment_and_monitor_read_only(tmp_path):
    from PySide6 import QtWidgets as W
    from pzt_panel import LayoutEditor, PZTPanel
    app = W.QApplication.instance() or W.QApplication([])
    class Reader:
        def __init__(self, path): pass
        def snapshot(self): return example_status(time.time())
        def close(self): pass
    panel = PZTPanel(dict(pzt_layout_path=str(tmp_path/'layout.json')), reader_factory=Reader)
    panel.resize(440, 690)
    panel.show()
    editor = LayoutEditor(blank_layout())
    editor.show()
    app.processEvents()
    try:
        rows = list(editor.fields.values())
        for key in ('x_mm', 'y_mm', 'normal_deg', 'diameter_mm', 'wiring'):
            assert len({r[key].geometry().x() for r in rows}) == 1
            assert len({r[key].width() for r in rows}) == 1
        rows[0]['placed'].setChecked(True)
        rows[0]['x_mm'].setValue(-12.)
        assert editor.value()['pzts']['PZT1']['x_mm'] == -12.
        assert 'AFG connected | ON' in panel.cards['PZT1'][1].text()
        assert '0.1 Vpp' in panel.cards['PZT1'][1].text()
        assert 'Wire: Not checked' in panel.cards['PZT1'][1].text()
        assert not any(b.text() in ('ON', 'OFF', 'Apply') for b in panel.findChildren(W.QPushButton))
        assert panel.top.layout_data == blank_layout()
    finally:
        editor.close()
        panel.shutdown()
        panel.close()
