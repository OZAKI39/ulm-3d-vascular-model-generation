"""GUI callback regression tests with no window and no physical hardware."""
from unittest import TestCase
from unittest.mock import Mock

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
