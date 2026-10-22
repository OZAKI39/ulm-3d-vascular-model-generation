import json
from pathlib import Path
import tempfile
import time
import unittest
from unittest.mock import patch

from afg_control import Rig, load_config
from afg_telemetry import TelemetryPublisher, read_monitor_channel
from test_afg_control import FakeManager


class TelemetryTests(unittest.TestCase):
    def test_seven_query_display_values_match_full_readback(self):
        manager=FakeManager();rig=Rig(load_config(),manager=manager);rig.connect()
        for name in rig.config['channels']:
            full=rig.read_channel(name)
            session,_=rig._target(name)
            session.queries.clear();session.writes.clear()
            display=read_monitor_channel(rig,name)
            self.assertEqual(len(session.queries),7)
            self.assertEqual(session.writes,[])
            self.assertEqual(display,{k:full[k] for k in display})
            self.assertEqual(rig.config['query_batch_mode'],'legacy_lines')
        self.assertFalse(rig.touched);rig.close()

    def test_nonfinite_monitor_values_invalidate_the_device(self):
        rig=Rig(load_config(),manager=FakeManager());rig.connect()
        rig.sessions['A'].override['SOURce1:FREQuency:FIXed?']='nan'
        with tempfile.TemporaryDirectory() as tmp:
            p=TelemetryPublisher(Path(tmp)/'status.json',rig.config)
            p.sample(rig,'PZT1')
            self.assertEqual(p.data['channels']['PZT1']['state'],'READ_ERROR')
            self.assertIsNone(p.data['channels']['PZT1']['output'])
            self.assertIn(id(rig.sessions['A']),rig._unreliable_sessions)
        rig.close()

    def test_requested_read_is_fresh_query_only_and_keeps_control_protocol(self):
        rig = Rig(load_config(), manager=FakeManager())
        rig.connect()
        rig.sessions['A'].values['OUTPut1:STATe'] = '1'
        with tempfile.TemporaryDirectory() as tmp:
            p = TelemetryPublisher(Path(tmp)/'status.json', rig.config)
            request = dict(schema=1, action='READ_STATUS', request_id='1'*32,
                           session=p.data['session'], created_unix_s=time.time())
            path = p.path.with_name('observation_read_request.json')
            path.write_text(json.dumps(request))
            accepted = p.take_request()
            self.assertEqual(accepted, request)
            self.assertIsNone(p.take_request())
            protocol = rig.config['query_batch_mode']
            p.read_requested(rig, accepted)
            result = json.loads(p.path.read_text())
            self.assertEqual(result['refresh']['request_id'], request['request_id'])
            self.assertEqual(result['refresh']['state'], 'COMPLETE')
            self.assertTrue(result['channels']['PZT1']['output'])
            self.assertEqual(rig.config['query_batch_mode'], protocol)
            self.assertFalse(rig.touched)
            self.assertTrue(all(command.endswith('?') for s in rig.sessions.values() for command in s.queries))
            self.assertTrue(all(not s.writes for s in rig.sessions.values()))
        rig.close()

    def test_requests_cannot_control_outputs_or_replay_other_sessions(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = TelemetryPublisher(Path(tmp)/'status.json', load_config())
            path = p.path.with_name('observation_read_request.json')
            req = dict(schema=1, action='READ_STATUS', request_id='2'*32,
                       session=p.data['session'], created_unix_s=100.)
            for override in (dict(action='OUTPUT_ON'), dict(session='old'), dict(created_unix_s=0), dict(request_id='invalid')):
                path.write_text(json.dumps(dict(req, **override)))
                self.assertIsNone(p.take_request(now=101.))

    def test_failed_query_blocks_that_session_and_returns_partial(self):
        rig = Rig(load_config(), manager=FakeManager()); rig.connect()
        original = rig.sessions['A'].query
        rig.sessions['A'].query = lambda c: (_ for _ in ()).throw(TimeoutError('late response'))
        with tempfile.TemporaryDirectory() as tmp:
            p = TelemetryPublisher(Path(tmp)/'status.json', rig.config)
            p.read_requested(rig, dict(request_id='3'*32))
            self.assertEqual(p.data['refresh']['state'], 'PARTIAL')
            self.assertEqual(p.data['refresh']['failed_channels'], ['PZT1','PZT2'])
            self.assertIn(id(rig.sessions['A']), rig._unreliable_sessions)
            self.assertEqual(rig.config['query_batch_mode'], 'legacy_lines')
        rig.sessions['A'].query = original; rig.close()

    def test_read_only_publication_preserves_on_outputs(self):
        manager = FakeManager()
        manager.sessions[0].values['OUTPut1:STATe'] = '1'
        rig = Rig(load_config(), manager=manager)
        rig.connect()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'status.json'
            publisher = TelemetryPublisher(path, rig.config)
            for name in rig.config['channels']:
                publisher.sample(rig, name)
            data = json.loads(path.read_text())
            self.assertTrue(data['channels']['PZT1']['output'])
            self.assertEqual(data['channels']['PZT4']['serial'], '2439003')
            self.assertEqual(data['channels']['PZT1']['parameters']['frequency_hz'], 2.5e6)
            self.assertFalse(rig.touched)
            rig.close()
            self.assertEqual(manager.sessions[0].values['OUTPut1:STATe'], '1')
            self.assertTrue(all(command.endswith('?') for session in manager.sessions
                                for write in session.writes for command in write.split(';')))

    def test_partial_off_does_not_refresh_parameters(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = TelemetryPublisher(Path(tmp)/'status.json', load_config())
            p.record('PZT1', dict(output=True, frequency_hz=2.5e6, amplitude_vpp=.1, phase_deg=90), now=100.)
            p.record('PZT1', dict(output=False), now=108.)
            row = p.data['channels']['PZT1']
            self.assertEqual(row['read_unix_s'], 108.)
            self.assertEqual(row['parameters_unix_s'], 100.)
            self.assertFalse(row['output'])

    def test_read_failure_invalidates_only_the_affected_device(self):
        rig = Rig(load_config(), manager=FakeManager())
        rig.connect()
        with tempfile.TemporaryDirectory() as tmp:
            p = TelemetryPublisher(Path(tmp)/'status.json', rig.config)
            p.sample(rig, 'PZT3')
            with patch('afg_telemetry.read_monitor_channel', side_effect=TimeoutError('USB lost')):
                p.sample(rig, 'PZT1')
            self.assertEqual(p.data['channels']['PZT1']['state'], 'READ_ERROR')
            self.assertEqual(p.data['channels']['PZT2']['state'], 'READ_ERROR')
            self.assertEqual(p.data['channels']['PZT3']['state'], 'READ_OK')
            self.assertIsNone(p.data['channels']['PZT1']['output'])
        rig.close()

    def test_publish_failure_does_not_interrupt_control_callbacks(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = TelemetryPublisher(Path(tmp)/'status.json', load_config())
            with patch('afg_telemetry.os.replace', side_effect=PermissionError('file busy')):
                self.assertFalse(p.publish())


if __name__ == '__main__':
    unittest.main()
