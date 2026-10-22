"""Read-only telemetry, fixed experiment geometry, and user channel assignments."""
import copy
import json
import math
import os
from pathlib import Path
import threading
import time
import uuid
from collections import deque

NAMES = ('PZT1', 'PZT2', 'PZT3', 'PZT4')
MAX_AGE_S = 8.
STATUS_POLL_S = .04


MODEL_PATH = Path(__file__).resolve().parents[1]/'configs/pzt_experiment_model.json'
EXPERIMENT_MODEL = json.loads(MODEL_PATH.read_text(encoding='utf-8'))
PHYSICAL_PZTS = ('X-', 'X+', 'Y-', 'Y+')
ENDPOINTS = {
    'A/CH1': dict(device='A', channel=1, serial='2147075'),
    'A/CH2': dict(device='A', channel=2, serial='2147075'),
    'B/CH1': dict(device='B', channel=1, serial='2439003'),
    'B/CH2': dict(device='B', channel=2, serial='2439003'),
}


def experiment_layout():
    # Geometry is fixed by the experiment model. Only cable assignments are user input.
    return dict(schema=2, model_id=EXPERIMENT_MODEL['model_id'],
                assignments={key: dict(**identity, pzt=None) for key, identity in ENDPOINTS.items()})


def validate_layout(data):
    if not isinstance(data, dict) or data.get('schema') != 2:
        raise ValueError('Select the PZT connected to each AFG channel.')
    if data.get('model_id') != EXPERIMENT_MODEL['model_id']:
        raise ValueError('Connections belong to a different experiment model.')
    if not isinstance(data.get('assignments'), dict) or set(data['assignments']) != set(ENDPOINTS):
        raise ValueError('Expected A/CH1, A/CH2, B/CH1 and B/CH2.')
    result = experiment_layout()
    selected = set()
    for key, identity in ENDPOINTS.items():
        row = data['assignments'][key]
        if not isinstance(row, dict) or any(row.get(k) != v for k, v in identity.items()):
            raise ValueError(f'{key}: generator identity does not match.')
        target = row.get('pzt')
        if target is not None and target not in PHYSICAL_PZTS:
            raise ValueError(f'{key}: select X-, X+, Y- or Y+.')
        if target in selected:
            raise ValueError(f'{target} is selected more than once.')
        if target is not None:
            selected.add(target)
        result['assignments'][key]['pzt'] = target
    return result


def save_layout(path, data):
    result = validate_layout(data)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(result, indent=2), encoding='utf-8')
    os.replace(temporary, path)
    return result


def endpoint_view(data, endpoint, now=None):
    identity = ENDPOINTS[endpoint]
    matches = [name for name, row in data.get('channels', {}).items()
               if all(row.get(key) == value for key, value in identity.items())]
    if len(matches) == 1:
        return channel_view(data, matches[0], now)
    view = channel_view({}, endpoint, now)
    view['connection'] = 'Identity mismatch' if data.get('channels') else 'No readback'
    return view


def physical_views(data, layout, now=None):
    result = {name: dict(connection='Unassigned', output='Unassigned', color='#a3afbd',
                        parameters={}, endpoint=None) for name in PHYSICAL_PZTS}
    for endpoint, row in validate_layout(layout)['assignments'].items():
        if row['pzt'] is not None:
            result[row['pzt']] = dict(endpoint_view(data, endpoint, now), endpoint=endpoint)
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
        self.wake = threading.Event()
        self.outbound = None
        self.request_error = None
        self.file_token = None
        self.history = deque(maxlen=128)
        self.thread = threading.Thread(target=self.run, name='pzt-status-file', daemon=True)
        self.thread.start()

    def run(self):
        while not self.done.is_set():
            with self.lock:
                request = self.outbound
                self.outbound = None
            if request is not None:
                try:
                    path = self.path.with_name('observation_read_request.json')
                    temporary = path.with_name(path.name + '.' + request['request_id'] + '.tmp')
                    temporary.write_text(json.dumps(request), encoding='utf-8')
                    os.replace(temporary, path)
                except OSError as exc:
                    with self.lock:
                        self.request_error = dict(request_id=request['request_id'], error=str(exc))
            # Check metadata often, parse JSON only on a new atomic publication.
            # This keeps slow Windows filesystem work out of Qt and avoids
            # rereading/decoding identical snapshots at every UI refresh.
            try:
                stat = self.path.stat()
                token = (stat.st_mtime_ns, stat.st_size)
            except OSError:
                token = None
            if token != self.file_token or token is None:
                data = read_status(self.path)
                with self.lock:
                    self.data = data
                    self.history.append(dict(data=data, observed_monotonic=time.monotonic(),
                                             observed_unix_s=time.time()))
                self.file_token = token
            self.wake.wait(STATUS_POLL_S)
            self.wake.clear()

    def snapshot(self):
        with self.lock:
            result = copy.deepcopy(self.data)
            if self.request_error:
                result['_request_error'] = copy.deepcopy(self.request_error)
            return result

    def snapshot_at(self, frame_monotonic):
        """Only observations available before this frame arrived, never future readback."""
        with self.lock:
            for sample in reversed(self.history):
                if sample['observed_monotonic'] <= frame_monotonic:
                    return copy.deepcopy(sample)
        return None

    def request_refresh(self):
        with self.lock:
            data = self.data
            now = time.time()
            if data.get('controller') == 'CLOSED' or not -2 <= now-float(data.get('published_unix_s') or 0) <= MAX_AGE_S:
                raise RuntimeError('Open the AFG controller to read status.')
            if not data.get('refresh_supported') or not data.get('session'):
                raise RuntimeError('Restart the updated AFG controller to enable reading.')
            request_id = uuid.uuid4().hex
            self.outbound = dict(schema=1, action='READ_STATUS', request_id=request_id,
                                 session=data['session'], created_unix_s=now)
            self.request_error = None
        self.wake.set()
        return request_id

    def close(self):
        self.done.set()
        self.wake.set()
        self.thread.join(.6)
