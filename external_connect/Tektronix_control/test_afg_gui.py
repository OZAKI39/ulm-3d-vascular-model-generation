"""GUI callback regression tests with no window and no physical hardware."""
from unittest import TestCase
from unittest.mock import Mock, patch

from afg_control import ChannelSettings
from afg_gui import App
from test_afg_control import FakeManager
from afg_control import Rig, load_config


class OutputCallbackTests(TestCase):
    def test_gui_on_uses_single_readback_and_off_preserves_last_read_parameters(self):
        manager = FakeManager()
        rig = Rig(load_config(), manager=manager)
        rig.connect()
        rig.apply({"PZT1": ChannelSettings(phase_deg=90)})
        app = App.__new__(App)
        app.rig = rig
        app.states = {"PZT1": Mock()}
        app.last_readback = {}
        app.submit = lambda title, action, callback: callback(action())
        session = manager.sessions[0]
        session.queries.clear()
        app.output("PZT1", True)
        self.assertEqual(len(session.queries), 21)
        text = app.states["PZT1"].set.call_args.args[0]
        self.assertIn("Output: ON", text)
        self.assertIn("Phase 90.00", text)
        session.queries.clear()
        app.output("PZT1", False)
        self.assertEqual(len(session.queries), 2)
        text = app.states["PZT1"].set.call_args.args[0]
        self.assertIn("Output: OFF", text)
        self.assertIn("Last read:", text)
        self.assertIn("Phase 90.00", text)
        rig.close()

    def test_gui_off_before_any_read_is_supported(self):
        app = App.__new__(App)
        app.states = {"PZT1": Mock()}
        app.last_readback = {}
        app.show_states({"PZT1": {"output": False}})
        app.states["PZT1"].set.assert_called_once_with("Output: OFF")

    def test_live_callback_keeps_output_on_without_extra_status_query(self):
        manager = FakeManager()
        rig = Rig(load_config(), manager=manager)
        rig.connect()
        rig.apply({"PZT1": ChannelSettings()})
        rig.set_output("PZT1", True)
        app = App.__new__(App)
        app.rig, app.config = rig, rig.config
        app.busy = app.closing = False
        app.states, app.last_readback = {"PZT1": Mock()}, {}
        app.fields = {"PZT1": {key: Mock(get=Mock(return_value=value)) for key, value in
                             {"mode": "Continuous Sine", "frequency": "3", "amplitude": "0.2", "phase": "90"}.items()}}
        app.submit = lambda title, action, callback: callback(action())
        session = manager.sessions[0]
        session.writes.clear()
        session.queries.clear()
        self.assertEqual(app.live_one("PZT1"), "break")
        self.assertEqual(len(session.queries), 25)
        self.assertTrue(all(not w.startswith("OUTPut1:STATe ") for w in session.writes))
        self.assertIn("Output: ON", app.states["PZT1"].set.call_args.args[0])
        self.assertIn("3 MHz / 0.2 Vpp", app.states["PZT1"].set.call_args.args[0])
        rig.close()

    def test_live_enter_while_busy_is_not_queued(self):
        app = App.__new__(App)
        app.busy, app.closing = True, False
        app.note, app.submit = Mock(), Mock()
        self.assertEqual(app.live_one("PZT1"), "break")
        app.submit.assert_not_called()
        app.note.assert_called_once()
