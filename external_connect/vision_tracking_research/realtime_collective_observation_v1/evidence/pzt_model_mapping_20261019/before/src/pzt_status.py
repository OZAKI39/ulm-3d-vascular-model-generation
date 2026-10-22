"""Read-only Windows telemetry consumer and manually registered lab layout."""
import copy
import json
import math
import os
from pathlib import Path
import threading
import time

NAMES = ('PZT1', 'PZT2', 'PZT3', 'PZT4')
MAX_AGE_S = 8.


def blank_layout():
    return dict(schema=1, view='XY_TOP', units='mm', source='USER_ENTERED',
                pzts={n: dict(placed=False, x_mm=0., y_mm=0., normal_deg=0.,
                    diameter_mm=7., wiring='Not checked') for n in NAMES})


def validate_layout(data):
    if not isinstance(data, dict):
        raise ValueError('Expected a layout object.')
    if data.get('schema') != 1 or data.get('view') != 'XY_TOP' or data.get('units') != 'mm':
        raise ValueError('Expected a top-view layout in millimetres (schema 1).')
    if set(data.get('pzts', {})) != set(NAMES):
        raise ValueError('The layout must contain PZT1 through PZT4.')
    result = blank_layout()
    for name in NAMES:
        row = data['pzts'][name]
        if not isinstance(row, dict):
            raise ValueError(f'{name}: expected position fields.')
        if type(row.get('placed')) is not bool:
            raise ValueError(f'{name}: position enabled must be true or false.')
        clean = dict(placed=row['placed'], wiring=row.get('wiring', 'Not checked'))
        if clean['wiring'] not in ('Not checked', 'Connected', 'Disconnected'):
            raise ValueError(f'{name}: invalid manual wiring status.')
        for key in ('x_mm', 'y_mm', 'normal_deg', 'diameter_mm'):
            value = float(row[key])
            if not math.isfinite(value) or abs(value) > 10000:
                raise ValueError(f'{name}: {key} must be finite and within 10000.')
            clean[key] = value
        if not .01 <= clean['diameter_mm'] <= 1000:
            raise ValueError(f'{name}: diameter must be between 0.01 and 1000 mm.')
        clean['normal_deg'] %= 360.
        result['pzts'][name] = clean
    return result


def save_layout(path, data):
    result = validate_layout(data)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(result, indent=2), encoding='utf-8')
    os.replace(temporary, path)
    return result


def read_status(path):
    try:
        path = Path(path)
        if path.stat().st_size > 262144:
            raise ValueError('Status file is too large.')
        data = json.loads(path.read_text(encoding='utf-8'))
        if not isinstance(data, dict):
            raise ValueError('Expected a status object.')
        if data.get('schema') != 1 or data.get('source') != 'AFG_GUI_READBACK':
            raise ValueError('Unsupported status source.')
        if set(data['channels']) != set(NAMES):
            raise ValueError('Expected four PZT channels.')
        if not all(isinstance(row, dict) for row in data['channels'].values()):
            raise ValueError('Invalid channel status.')
        return data
    except (OSError, ValueError, KeyError, TypeError) as exc:
        return dict(error=str(exc), controller='UNAVAILABLE', channels={})


def channel_view(data, name, now=None):
    now = time.time() if now is None else now
    row = data.get('channels', {}).get(name, {})
    result = dict(name=name, connection='No readback', output='UNKNOWN', parameters={},
                  mapping='', age_s=None, color='#a3afbd', error=data.get('error', ''))
    if row:
        result['mapping'] = f"{row.get('device', '?')} / CH{row.get('channel', '?')}"
    try:
        published = float(data.get('published_unix_s') or 0)
        read = float(row.get('read_unix_s') or 0)
        age = now - read
        if not all(math.isfinite(v) for v in (published, read, age)):
            raise ValueError('Invalid timestamp')
        result['age_s'] = age if read > 0 else None
        if data.get('controller') == 'CLOSED':
            result['connection'] = 'Controller closed'
        elif data.get('controller') == 'DISCONNECTED':
            result['connection'] = 'Not connected'
        elif not -2 <= now - published <= MAX_AGE_S:
            result['connection'] = 'Status expired'
        elif row.get('state') == 'READ_ERROR':
            result.update(connection='Read failed', color='#eda18c', error=row.get('error', ''))
        elif row.get('state') != 'READ_OK' or read <= 0:
            result['connection'] = 'Waiting for read'
        elif not -2 <= age <= MAX_AGE_S:
            result['connection'] = 'Readback expired'
        elif type(row.get('output')) is not bool:
            raise ValueError('Invalid output state')
        else:
            result.update(connection='AFG connected', output='ON' if row['output'] else 'OFF',
                          color='#8ed8b4' if row['output'] else '#aac1d5')
            params = row.get('parameters', {})
            if not isinstance(params, dict):
                raise ValueError('Invalid parameter object')
            ptime = float(row.get('parameters_unix_s') or 0)
            if math.isfinite(ptime) and -2 <= now - ptime <= MAX_AGE_S:
                for key in ('frequency_hz', 'amplitude_vpp', 'phase_deg'):
                    if not math.isfinite(float(params[key])):
                        raise ValueError('Invalid parameter readback')
                result['parameters'] = copy.deepcopy(params)
    except (ValueError, KeyError, TypeError, OverflowError) as exc:
        result.update(connection='Invalid readback', output='UNKNOWN', parameters={}, error=str(exc))
    return result


class StatusReader:
    """Mounted Windows disks are read off the Qt/render thread."""
    def __init__(self, path):
        self.path = Path(path)
        self.lock = threading.Lock()
        self.data = {}
        self.done = threading.Event()
        self.thread = threading.Thread(target=self.run, name='pzt-status-file', daemon=True)
        self.thread.start()

    def run(self):
        while not self.done.is_set():
            data = read_status(self.path)
            with self.lock:
                self.data = data
            self.done.wait(.5)

    def snapshot(self):
        with self.lock:
            return copy.deepcopy(self.data)

    def close(self):
        self.done.set()
        self.thread.join(.6)
