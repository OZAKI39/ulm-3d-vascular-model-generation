"""Atomic status publication for the observation UI; never writes SCPI controls.

Called by the existing controller's single USB worker. A separate observer reads
this file without opening another VISA session. No load / PZT cable sensing is
available: instrument output state is not proof of acoustic output.
"""
import copy
import json
import logging
import os
from pathlib import Path
import threading
import time
import uuid


class TelemetryPublisher:
    def __init__(self, path, config):
        self.path = Path(path)
        self.lock = threading.RLock()
        self.data = dict(schema=1, source='AFG_GUI_READBACK', session=uuid.uuid4().hex,
                         pid=os.getpid(), sequence=0, controller='DISCONNECTED',
                         devices={}, channels={}, wiring_sensed=False,
                         refresh_supported=True, refresh=None, read_protocol='sequential')
        self.seen_request = None
        for name, ch in config['channels'].items():
            self.data['channels'][name] = dict(device=ch['device'], channel=ch['channel'],
                serial=config['devices'][ch['device']]['serial'], state='UNKNOWN',
                output=None, read_unix_s=None, parameters_unix_s=None, parameters={})

    def record(self, name, value, now=None):
        with self.lock:
            now = time.time() if now is None else now
            row = self.data['channels'][name]
            row.update(state='READ_OK', output=bool(value['output']), read_unix_s=now, error=None)
            # OFF confirmations carry only the switch state. They must not make
            # a cached frequency/phase look like a new parameter measurement.
            if 'frequency_hz' in value:
                row['parameters'] = copy.deepcopy({k: v for k, v in value.items() if k != 'output'})
                row['parameters_unix_s'] = now

    def invalidate(self, names=None, message='Readback unavailable'):
        with self.lock:
            for name in names or self.data['channels']:
                self.data['channels'][name].update(state='READ_ERROR', output=None, error=str(message))

    def sample(self, rig, name):
        """One channel per worker job, preserving control ordering and USB ownership."""
        with self.lock:
            self.data['devices'] = copy.deepcopy(rig.identities)
            self.data['controller'] = 'CONNECTED'
        try:
            # Use the already verified sequential query path for monitoring.
            # Keep control operations on their existing optimized protocol.
            # The same Rig lock protects this temporary software setting and
            # prevents a control job from observing the override.
            with rig.lock:
                protocol = rig.config.get('query_batch_mode', 'sequential')
                try:
                    rig.config['query_batch_mode'] = 'sequential'
                    value = rig.read_channel(name)
                finally:
                    rig.config['query_batch_mode'] = protocol
            self.record(name, value)
        except Exception as exc:
            device = self.data['channels'][name]['device']
            # A timed-out sequential query can also leave a late reply queued.
            # Reuse the driver's fail-closed guard; never retry a desynchronized
            # session or let that reply be used as an output confirmation.
            with rig.lock:
                session = rig.sessions.get(device)
                if session is not None:
                    rig._unreliable_sessions.add(id(session))
                    rig.prepared.clear()
            self.invalidate([n for n, c in self.data['channels'].items() if c['device'] == device], exc)
        self.publish()

    def take_request(self, now=None):
        """Accept only a fresh read request addressed to this controller session."""
        path = self.path.with_name('observation_read_request.json')
        try:
            if path.stat().st_size > 8192:
                return None
            req = json.loads(path.read_text(encoding='utf-8'))
            if not isinstance(req, dict):
                return None
            now = time.time() if now is None else now
            request_id = req.get('request_id')
            if (req.get('schema') != 1 or req.get('action') != 'READ_STATUS'
                    or req.get('session') != self.data['session']
                    or not isinstance(request_id, str) or len(request_id) != 32
                    or request_id == self.seen_request
                    or not -2 <= now-float(req.get('created_unix_s', 0)) <= 30):
                return None
            int(request_id, 16)
            self.seen_request = request_id
            return req
        except (OSError, ValueError, TypeError):
            return None

    def begin_request(self, req):
        with self.lock:
            self.data['refresh'] = dict(request_id=req['request_id'], state='READING',
                                        started_unix_s=time.time())
        self.publish()

    def read_requested(self, rig, req):
        self.begin_request(req)
        for name in self.data['channels']:
            self.sample(rig, name)
        with self.lock:
            failed = [name for name, row in self.data['channels'].items() if row['state'] != 'READ_OK']
            self.data['refresh'].update(state='FAILED' if len(failed) == len(self.data['channels']) else 'PARTIAL' if failed else 'COMPLETE',
                                        failed_channels=failed, completed_unix_s=time.time())
        self.publish()

    def fail_request(self, req, error):
        self.invalidate(message=error)
        with self.lock:
            self.data['refresh'] = dict(request_id=req['request_id'], state='FAILED',
                                        error=str(error), completed_unix_s=time.time())
        self.publish()

    def publish(self, controller=None):
        # A filesystem failure must not interfere with the existing hardware
        # shutdown path or be reported as a generator control failure.
        with self.lock:
            if controller is not None:
                self.data['controller'] = controller
            self.data['sequence'] += 1
            self.data['published_unix_s'] = time.time()
            temporary = self.path.with_name(self.path.name + f'.{self.data["session"]}.tmp')
            try:
                self.path.parent.mkdir(parents=True, exist_ok=True)
                temporary.write_text(json.dumps(self.data, indent=2, allow_nan=False), encoding='utf-8')
                os.replace(temporary, self.path)
                return True
            except (OSError, ValueError):
                logging.exception('Cannot publish observation status')
                return False
