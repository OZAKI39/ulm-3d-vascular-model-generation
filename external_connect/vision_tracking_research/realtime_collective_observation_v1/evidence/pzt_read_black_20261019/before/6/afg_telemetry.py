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
                         devices={}, channels={}, wiring_sensed=False)
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
            value = rig.read_channel(name)  # Queries only; no error-queue clearing.
            self.record(name, value)
        except Exception as exc:
            device = self.data['channels'][name]['device']
            self.invalidate([n for n, c in self.data['channels'].items() if c['device'] == device], exc)
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
