"""Offline regression tests: no VISA driver or hardware is used."""
from copy import deepcopy
from dataclasses import replace
import math
import re
import unittest

from afg_control import AFGError, ChannelSettings, MODULATIONS, Rig, load_config


class FakeSession:
    def __init__(self, serial):
        self.serial = serial
        self.resource_name = f"USB0::0x0699::0x0353::{serial}::INSTR"
        self.writes = []
        self.queries = []
        self.closed = False
        self.fail_write = None
        self.errors = []
        self.override = {}
        self.pending = []
        self.fail_read_after = None
        self.read_count = 0
        self.values = {}
        for c in (1, 2):
            p = f"SOURce{c}:"
            self.values.update({f"OUTPut{c}:STATe": "0", f"OUTPut{c}:IMPedance": "50",
                p+"FUNCtion:SHAPe": "SIN", p+"FREQuency:FIXed": "2500000", p+"FREQuency:MODE": "CW",
                p+"VOLTage:LEVel:IMMediate:AMPLitude": "0.1", p+"VOLTage:LEVel:IMMediate:OFFSet": "0",
                p+"PHASe:ADJust": "0", p+"BURSt:STATe": "0", p+"FREQuency:CONCurrent": "0",
                p+"BURSt:MODE": "TRIG", p+"BURSt:NCYCles": "10", p+"BURSt:SOURce": "MAN",
                p+"FREQuency:STARt": "2000000", p+"FREQuency:STOP": "3000000",
                p+"SWEep:TIME": "0.1", p+"SWEep:SPACing": "LIN", p+"SWEep:SOURce": "MAN"})
            self.values.update({p+m+":STATe": "0" for m in MODULATIONS})

    def write(self, command):
        self.writes.append(command)
        if self.fail_write and self.fail_write(command):
            raise TimeoutError("Injected USB failure")
        if command.endswith("?"):
            self.pending.extend(self.query(c) for c in command.split(";"))
            return
        if " " not in command:
            return
        key, value = command.split(" ", 1)
        if value in ("ON", "OFF"):
            value = "1" if value == "ON" else "0"
        elif value.endswith("DEG"):
            value = str(math.radians(float(value[:-3])))
        elif re.match(r"^[+\-0-9.]", value):
            value = re.sub(r"(?:Vpp|Hz|V|s)$", "", value)
        self.values[key] = value

    def query(self, command):
        self.queries.append(command)
        if command in self.override:
            return self.override[command]
        if command == "*IDN?":
            return f"TEKTRONIX,AFG1062,{self.serial},1.2.0"
        if command == "SYSTem:ERRor:NEXT?":
            return self.errors.pop(0) if self.errors else '0,"No error"'
        return self.values[command[:-1]]

    def close(self):
        self.closed = True

    def read(self):
        if self.fail_read_after is not None and self.read_count >= self.fail_read_after:
            raise TimeoutError("Injected incomplete batch")
        self.read_count += 1
        if not self.pending:
            raise TimeoutError("No queued response")
        return self.pending.pop(0)


class FakeManager:
    def __init__(self):
        self.sessions = [FakeSession("2147075"), FakeSession("2439003")]
        self.closed = False

    def list_resources(self, _):
        return tuple(s.resource_name for s in reversed(self.sessions))

    def open_resource(self, resource):
        return next(s for s in self.sessions if s.resource_name == resource)

    def close(self):
        self.closed = True


class ControllerTests(unittest.TestCase):
    batch_mode = "sequential"

    def setUp(self):
        self.cfg = load_config()
        self.cfg["query_batch_mode"] = self.batch_mode
        self.manager = FakeManager()
        self.rig = Rig(self.cfg, manager=self.manager)
        self.rig.connect()
        self.a, self.b = self.manager.sessions

    def tearDown(self):
        self.rig.touched = False
        self.rig.close()

    def assert_no_enable(self):
        self.assertFalse(any(re.match(r"OUTPut[12]:STATe ON$", w) for s in (self.a, self.b) for w in s.writes))

    def test_read_only_connect_status_close(self):
        self.rig.status()
        self.rig.close()
        self.assertTrue(all(all(c.endswith("?") for c in w.split(";"))
                            for w in self.a.writes + self.b.writes))

    def test_serial_binding_ignores_enumeration_order(self):
        self.assertIs(self.rig.sessions["A"], self.a)
        self.assertIs(self.rig.sessions["B"], self.b)

    def test_missing_device_closes_partial_connection_without_output_changes(self):
        manager = FakeManager()
        manager.sessions.pop()
        with self.assertRaises(AFGError):
            Rig(self.cfg, manager=manager).connect()
        self.assertTrue(manager.closed)
        self.assertTrue(manager.sessions[0].closed)
        self.assertEqual(manager.sessions[0].writes, [])

    def test_wrong_idn_rejected(self):
        manager = FakeManager()
        manager.sessions[0].override["*IDN?"] = "TEKTRONIX,AFG1022,2147075,1.0"
        with self.assertRaises(AFGError):
            Rig(self.cfg, manager=manager).connect()
        self.assertTrue(manager.closed)

    def test_validate_entire_batch_before_any_write(self):
        with self.assertRaises(ValueError):
            self.rig.apply({"PZT1": ChannelSettings(), "PZT4": ChannelSettings(amplitude_vpp=0.501)})
        self.assertEqual(self.a.writes + self.b.writes, [])

    def test_nonfinite_and_invalid_parameters_rejected(self):
        for kw in ({"amplitude_vpp": float('nan')}, {"frequency_hz": float('inf')},
                   {"frequency_hz": 10_000_001}, {"amplitude_vpp": 0},
                   {"amplitude_vpp": 0.1001}, {"phase_deg": -1},
                   {"mode": "sweep", "sweep_time_s": 0},
                   {"mode": "burst", "burst_cycles": 1.1},
                   {"mode": "gated", "trigger_source": "MAN"}):
            with self.subTest(kw=kw), self.assertRaises(ValueError):
                self.rig.apply({"PZT1": ChannelSettings(**kw)})
        self.assertEqual(self.a.writes + self.b.writes, [])

    def test_hard_amplifier_limit_cannot_be_relaxed_by_configuration(self):
        self.rig.config["channels"]["PZT1"]["max_input_vpp"] = 4.1
        with self.assertRaises(ValueError):
            self.rig.apply({"PZT1": ChannelSettings()})
        self.assertEqual(self.a.writes, [])

    def test_apply_is_off_and_correct_phase_conversion(self):
        values = self.rig.apply({"PZT1": ChannelSettings(phase_deg=90)})
        self.assertAlmostEqual(values["PZT1"]["phase_deg"], 90)
        self.assertEqual(self.a.writes[0], "OUTPut1:STATe OFF")
        self.assertIn("SOURce1:PHASe:ADJust 90DEG", self.a.writes)
        self.assertIn("OUTPut1:IMPedance 50", self.a.writes)
        self.assert_no_enable()

    def test_independent_channel_disables_other_frequency_coupling_first(self):
        self.rig.apply({"PZT1": ChannelSettings(frequency_hz=490000)})
        self.assertLess(self.a.writes.index("SOURce2:FREQuency:CONCurrent OFF"),
                        self.a.writes.index("SOURce1:FREQuency:FIXed 490000Hz"))

    def test_phase_correction_wraps(self):
        self.cfg["channels"]["PZT1"]["phase_correction_deg"] = 30
        self.rig.apply({"PZT1": ChannelSettings(phase_deg=350)})
        self.assertIn("SOURce1:PHASe:ADJust 20DEG", self.a.writes)

    def test_explicit_enable_requires_successful_apply(self):
        with self.assertRaises(AFGError):
            self.rig.set_output("PZT1", True)
        self.assert_no_enable()
        self.rig.apply({"PZT1": ChannelSettings()})
        self.rig.set_output("PZT1", True)
        self.assertEqual(self.a.values["OUTPut1:STATe"], "1")

    def test_panel_change_prevents_enabling(self):
        self.rig.apply({"PZT1": ChannelSettings()})
        self.a.values["SOURce1:VOLTage:LEVel:IMMediate:AMPLitude"] = "0.6"
        with self.assertRaises(AFGError):
            self.rig.set_output("PZT1", True)
        self.assert_no_enable()
        self.assertIn("OUTPut2:STATe OFF", self.b.writes)

    def test_error_queue_blocks_apply_and_shuts_down(self):
        self.a.errors = ['-113,"Undefined header"']
        with self.assertRaises(AFGError):
            self.rig.apply({"PZT1": ChannelSettings()})
        self.assert_no_enable()
        self.assertIn("OUTPut2:STATe OFF", self.b.writes)

    def test_usb_failure_still_attempts_remaining_outputs(self):
        self.a.fail_write = lambda c: c == "OUTPut1:STATe OFF"
        failures = self.rig.all_off()
        self.assertEqual(len(failures), 1)
        self.assertIn("OUTPut2:STATe OFF", self.a.writes)
        self.assertIn("OUTPut1:STATe OFF", self.b.writes)
        self.assertIn("OUTPut2:STATe OFF", self.b.writes)

    def test_sweep_and_triggered_burst_transitions(self):
        self.rig.apply({"PZT1": ChannelSettings(mode="sweep", sweep_start_hz=3e6, sweep_stop_hz=2e6, sweep_spacing="LOG")})
        self.assertTrue(self.rig.read_channel("PZT1")["frequency_mode"].startswith("SWE"))
        self.rig.apply({"PZT1": ChannelSettings(mode="burst", burst_cycles=123, trigger_source="EXT")})
        value = self.rig.read_channel("PZT1")
        self.assertEqual(value["frequency_mode"], "CW")
        self.assertEqual(value["burst_cycles"], 123)
        self.rig.apply({"PZT1": ChannelSettings()})
        self.assertFalse(self.rig.read_channel("PZT1")["burst_enabled"])

    def test_gated_does_not_send_invalid_burst_source_command(self):
        self.rig.apply({"PZT4": ChannelSettings(mode="gated", trigger_source="EXT")})
        self.assertIn("SOURce2:BURSt:MODE GATed", self.b.writes)
        self.assertFalse(any("BURSt:SOURce" in w for w in self.b.writes))
        self.assert_no_enable()

    def test_stop_latch_prevents_on(self):
        self.rig.apply({"PZT1": ChannelSettings()})
        self.rig.request_stop()
        with self.assertRaises(AFGError):
            self.rig.set_output("PZT1", True)
        self.assert_no_enable()

    def test_manual_trigger_targets_one_device(self):
        self.rig.apply({"PZT1": ChannelSettings(mode="burst")})
        self.rig.set_output("PZT1", True)
        self.rig.trigger("A")
        self.assertIn("*TRG", self.a.writes)
        self.assertNotIn("*TRG", self.b.writes)

    def test_unknown_active_channel_blocks_trigger(self):
        self.rig.apply({"PZT1": ChannelSettings(mode="burst")})
        self.rig.set_output("PZT1", True)
        self.a.values["OUTPut2:STATe"] = "1"
        with self.assertRaises(AFGError):
            self.rig.trigger("A")
        self.assertNotIn("*TRG", self.a.writes)

    def test_phase_align_is_local_and_requires_outputs_off(self):
        self.a.values["OUTPut1:STATe"] = "1"
        with self.assertRaises(AFGError):
            self.rig.align_phase("A")
        self.a.values["OUTPut1:STATe"] = "0"
        self.rig.align_phase("A")
        self.assertIn("SOURce1:PHASe:INITiate", self.a.writes)
        self.assertEqual(self.b.writes, [])

    def test_close_after_control_shuts_down(self):
        self.rig.apply({"PZT4": ChannelSettings()})
        self.rig.set_output("PZT4", True)
        self.rig.close()
        self.assertEqual(self.b.values["OUTPut2:STATe"], "0")
        self.assertTrue(self.a.closed and self.b.closed)

    def test_switch_returns_verified_state_without_second_full_read(self):
        self.rig.apply({"PZT1": ChannelSettings(phase_deg=37)})
        self.a.queries.clear()
        state = self.rig.set_output("PZT1", True, return_state=True)
        self.assertTrue(state["output"])
        self.assertAlmostEqual(state["phase_deg"], 37)
        self.assertEqual(len(self.a.queries), 21)
        self.assertEqual(self.a.queries.count("SOURce1:VOLTage:LEVel:IMMediate:AMPLitude?"), 1)
        self.a.queries.clear()
        self.assertEqual(self.rig.set_output("PZT1", False, return_state=True), {"output": False})
        self.assertEqual(self.a.queries, ["OUTPut1:STATe?", "SYSTem:ERRor:NEXT?"])

    def test_error_after_precheck_prevents_on(self):
        self.rig.apply({"PZT1": ChannelSettings()})
        self.a.errors = ['0,"No error"', '-102,"Syntax error"']
        with self.assertRaises(AFGError):
            self.rig.set_output("PZT1", True)
        self.assert_no_enable()

    def test_both_precheck_errors_are_reported_and_block_on(self):
        self.rig.apply({"PZT1": ChannelSettings()})
        self.a.errors = ['-102,"Syntax error"', '-113,"Undefined header"', '-222,"Data out of range"']
        with self.assertRaises(AFGError) as raised:
            self.rig.set_output("PZT1", True)
        message = str(raised.exception)
        for code in ("-102", "-113", "-222"):
            self.assertIn(code, message)
        self.assertFalse(self.a.errors)
        self.assert_no_enable()

    def test_error_in_switch_confirmation_triggers_shutdown(self):
        self.rig.apply({"PZT1": ChannelSettings()})
        self.a.errors = ['0,"No error"', '0,"No error"', '-200,"Execution error"']
        with self.assertRaisesRegex(AFGError, "-200"):
            self.rig.set_output("PZT1", True)
        self.assertEqual(self.a.values["OUTPut1:STATe"], "0")
        self.assertFalse(self.rig.prepared)

    def test_error_from_mode_specific_query_blocks_on(self):
        self.rig.apply({"PZT1": ChannelSettings(mode="sweep")})
        self.a.errors = ['0,"No error"', '0,"No error"', '-102,"Syntax error"']
        with self.assertRaisesRegex(AFGError, "-102"):
            self.rig.set_output("PZT1", True)
        self.assert_no_enable()


class BatchControllerTests(ControllerTests):
    """Repeat all control/safety regression cases with the optimized transport."""
    batch_mode = "legacy_lines"

    def test_batch_matches_individual_queries_and_consumes_all_responses(self):
        self.rig.config["query_batch_mode"] = "sequential"
        before = self.rig.status()
        self.rig.config["query_batch_mode"] = "legacy_lines"
        self.assertEqual(self.rig.status(), before)
        self.assertFalse(self.a.pending or self.b.pending)
        self.assertEqual(len(self.a.writes), 2)
        self.assertEqual(len(self.b.writes), 2)
        for command in self.a.writes + self.b.writes:
            self.assertEqual(len(command.split(";")), 17)
            self.assertNotIn(";:", command)

    def test_partial_batch_blocks_on_and_still_sends_all_off(self):
        self.rig.apply({"PZT1": ChannelSettings()})
        self.a.fail_read_after = self.a.read_count + 4
        with self.assertRaises(TimeoutError):
            self.rig.set_output("PZT1", True)
        self.assert_no_enable()
        self.assertFalse(self.rig.prepared)
        self.assertIn("OUTPut2:STATe OFF", self.a.writes)
        self.assertIn("OUTPut2:STATe OFF", self.b.writes)
        query_count = len(self.a.queries)
        with self.assertRaisesRegex(AFGError, "Incomplete batch"):
            self.rig.read_channel("PZT1")
        self.assertEqual(len(self.a.queries), query_count)

    def test_unexpected_response_framing_blocks_on(self):
        self.rig.apply({"PZT1": ChannelSettings()})
        self.a.override["OUTPut1:STATe?"] = "0;50"
        with self.assertRaisesRegex(AFGError, "framing"):
            self.rig.set_output("PZT1", True)
        self.assert_no_enable()

    def test_read_only_batch_failure_does_not_change_outputs(self):
        self.a.fail_read_after = 0
        with self.assertRaises(TimeoutError):
            self.rig.read_channel("PZT1")
        self.rig.close()
        self.assertTrue(all(w.endswith("?") for w in self.a.writes))
        self.assertEqual(self.b.writes, [])


if __name__ == "__main__":
    unittest.main(verbosity=2)
