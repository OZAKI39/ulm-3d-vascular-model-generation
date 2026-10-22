"""Two Tektronix AFG1062 units / four High Wave 3.2 bridge amplifiers.

Commands follow Tektronix 077112902, with query framing tested on local firmware.
Opening a connection
and reading status never reset an instrument or enable/disable its outputs.
All frequencies are Hz, generator amplitudes Vpp into 50 ohms, phases degrees.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
import logging
import math
from pathlib import Path
import threading
from time import perf_counter
from typing import Any

BASE = Path(__file__).resolve().parent
LOG = logging.getLogger("afg")
MODULATIONS = ("AM", "FM", "PM", "PWM", "FSKey", "ASKey", "PSKey")
HARD_MAX_VPP = 4.0  # High Wave 3.2 INPUT limit, not PZT voltage.
HARD_MAX_HZ = 10_000_000.0  # Amplifier nominal bandwidth; AFG itself is 60 MHz.


class AFGError(RuntimeError):
    pass


def finite(value: Any, name: str, low: float, high: float) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{name}: boolean is not a number")
    v = float(value)
    if not math.isfinite(v) or not low <= v <= high:
        raise ValueError(f"{name} must be finite and in [{low:g}, {high:g}], got {value!r}")
    return v


@dataclass(frozen=True)
class ChannelSettings:
    frequency_hz: float = 2_500_000.0
    amplitude_vpp: float = 0.1
    phase_deg: float = 0.0
    mode: str = "continuous"  # continuous / sweep / burst / gated
    sweep_start_hz: float = 2_000_000.0
    sweep_stop_hz: float = 3_000_000.0
    sweep_time_s: float = 0.1
    sweep_spacing: str = "LIN"
    trigger_source: str = "MAN"  # TIM / MAN / EXT; gated requires EXT
    burst_cycles: int = 10

    def validate(self, limit: dict) -> None:
        cap = finite(limit["max_input_vpp"], "max_input_vpp", 0.001, HARD_MAX_VPP)
        freq_cap = finite(limit["max_frequency_hz"], "max_frequency_hz", 0.002, HARD_MAX_HZ)
        if self.mode not in ("continuous", "sweep", "burst", "gated"):
            raise ValueError("mode must be continuous / sweep / burst / gated")
        finite(self.amplitude_vpp, "amplitude_vpp", 0.001, cap)
        # Disallow sub-mV requests which the AFG would round unpredictably.
        if not math.isclose(self.amplitude_vpp * 1000, round(self.amplitude_vpp * 1000), abs_tol=1e-7):
            raise ValueError("Use 0.001 Vpp increments")
        finite(self.frequency_hz, "frequency_hz", 0.002 if self.mode in ("burst", "gated") else 1e-6, freq_cap)
        finite(self.phase_deg, "phase_deg", 0, 360)
        if self.trigger_source not in ("TIM", "MAN", "EXT"):
            raise ValueError("trigger_source must be TIM / MAN / EXT")
        if self.mode == "gated" and self.trigger_source != "EXT":
            raise ValueError("Gated burst uses the external trigger input; select EXT")
        if self.mode == "sweep":
            finite(self.sweep_start_hz, "sweep_start_hz", 1e-6, freq_cap)
            finite(self.sweep_stop_hz, "sweep_stop_hz", 1e-6, freq_cap)
            if self.sweep_start_hz == self.sweep_stop_hz:
                raise ValueError("Sweep endpoints must differ; descending sweeps are allowed")
            finite(self.sweep_time_s, "sweep_time_s", 0.001, 500)
            if self.sweep_spacing not in ("LIN", "LOG"):
                raise ValueError("sweep_spacing must be LIN / LOG")
        if self.mode == "burst":
            if type(self.burst_cycles) is not int or not 1 <= self.burst_cycles <= 1_000_000:
                raise ValueError("burst_cycles must be an integer in [1, 1000000]")


def load_config(path: str | Path = BASE / "config.json") -> dict:
    cfg = json.loads(Path(path).read_text(encoding="utf-8-sig"))
    if set(cfg["devices"]) != {"A", "B"}:
        raise ValueError("devices must contain A and B")
    serials = [str(d["serial"]) for d in cfg["devices"].values()]
    if len(set(serials)) != 2 or any(not s.strip() for s in serials):
        raise ValueError("A and B must have distinct nonempty serial numbers")
    if set(cfg["channels"]) != {"PZT1", "PZT2", "PZT3", "PZT4"}:
        raise ValueError("channels must contain PZT1 through PZT4")
    pairs = []
    for name, ch in cfg["channels"].items():
        if ch["device"] not in cfg["devices"] or type(ch["channel"]) is not int or ch["channel"] not in (1, 2):
            raise ValueError(f"Invalid channel mapping: {name}")
        pairs.append((ch["device"], ch["channel"]))
        finite(ch["max_input_vpp"], name + " voltage limit", 0.001, HARD_MAX_VPP)
        finite(ch["max_frequency_hz"], name + " frequency limit", 0.002, HARD_MAX_HZ)
        finite(ch.get("phase_correction_deg", 0), "phase correction", -360, 360)
    if len(set(pairs)) != 4:
        raise ValueError("Each PZT must map to a different physical AFG channel")
    finite(cfg.get("timeout_ms", 3000), "timeout_ms", 100, 10000)
    if cfg.get("query_batch_mode", "sequential") not in ("sequential", "legacy_lines"):
        raise ValueError("query_batch_mode must be sequential or legacy_lines")
    return cfg


def resource_manager(cfg: dict):
    try:
        import pyvisa
    except ImportError as e:
        raise AFGError("PyVISA is missing. Run setup_python.ps1 or pip install -r requirements.txt") from e
    try:
        library = cfg.get("visa_library", "")
        return pyvisa.ResourceManager(library) if library else pyvisa.ResourceManager()
    except Exception as e:
        raise AFGError("No usable VISA library. Install TekVISA or NI-VISA (including USBTMC). "
                       "PyVISA alone is not a Windows USB driver. " + str(e)) from e


def parse_idn(text: str) -> dict:
    fields = [s.strip() for s in text.split(",")]
    if len(fields) < 4 or "TEKTRONIX" not in fields[0].upper() or fields[1].upper() != "AFG1062":
        raise AFGError(f"Not a Tektronix AFG1062: {text!r}")
    return dict(manufacturer=fields[0], model=fields[1], serial=fields[2], firmware=",".join(fields[3:]))


def prepare_session(session, cfg):
    session.timeout = int(cfg.get("timeout_ms", 3000))
    session.write_termination = "\n"
    session.read_termination = "\n"


def discover(cfg: dict) -> list[dict]:
    """Read-only discovery. Only probes Tektronix USB instruments, not COM ports."""
    rm = resource_manager(cfg)
    found = []
    try:
        for address in rm.list_resources("USB?*::INSTR"):
            parts = address.split("::")
            try:
                if int(parts[1], 0) != 0x0699:
                    continue
            except (ValueError, IndexError):
                continue
            session = None
            try:
                session = rm.open_resource(address)
                prepare_session(session, cfg)
                found.append(dict(resource=address, **parse_idn(session.query("*IDN?").strip())))
            except Exception as e:
                found.append(dict(resource=address, error=str(e)))
            finally:
                if session is not None:
                    session.close()
    finally:
        rm.close()
    return found


class Rig:
    """Thread-safe four-channel controller. Use as a context manager.

    close() shuts down channels after this instance has performed any control
    operation. A read-only connection closes without changing existing outputs.
    USB loss/process termination cannot guarantee physical shutdown.
    """
    def __init__(self, config: dict | None = None, *, manager=None):
        self.config = config or load_config()
        self.rm = manager
        self.sessions = {}
        self.identities = {}
        self.prepared: dict[str, ChannelSettings] = {}
        self.lock = threading.RLock()
        self.stop_requested = threading.Event()
        self.touched = False
        self._unreliable_sessions = set()

    def connect(self):
        with self.lock:
            if self.sessions:
                raise AFGError("Already connected")
            self.rm = self.rm or resource_manager(self.config)
            try:
                resources = self.rm.list_resources("USB?*::INSTR")
                for alias, info in self.config["devices"].items():
                    serial = str(info["serial"])
                    matches = []
                    for address in resources:
                        parts = address.split("::")
                        try:
                            if len(parts) >= 5 and int(parts[1], 0) == 0x0699 and int(parts[2], 0) == 0x0353 and parts[3] == serial:
                                matches.append(address)
                        except ValueError:
                            pass
                    # Resolve from enumeration; never assign A/B by enumeration order.
                    if len(matches) != 1:
                        raise AFGError(f"{alias}: expected one USB AFG1062 with serial {serial}, found {len(matches)}. "
                                       "Check driver, power, USB cable and config.json.")
                    session = self.rm.open_resource(matches[0])
                    self.sessions[alias] = session
                    prepare_session(session, self.config)
                    identity = parse_idn(session.query("*IDN?").strip())
                    if identity["serial"] != serial:
                        raise AFGError(f"Serial mismatch for {alias}: {identity}")
                    self.identities[alias] = dict(resource=matches[0], **identity)
                return dict(self.identities)
            except BaseException:
                self._close_sessions()
                raise

    def _close_sessions(self):
        for session in self.sessions.values():
            try:
                session.close()
            except Exception:
                LOG.exception("Session close failed")
        self.sessions.clear()
        self.identities.clear()
        self.prepared.clear()
        self._unreliable_sessions.clear()
        if self.rm is not None:
            try:
                self.rm.close()
            finally:
                self.rm = None

    def close(self):
        with self.lock:
            failures = self.all_off() if self.touched and self.sessions else []
            self._close_sessions()
            if failures:
                raise AFGError("Shutdown could not be verified: " + "; ".join(failures))

    def __enter__(self):
        self.connect()
        return self

    def __exit__(self, exc_type, exc, tb):
        try:
            self.close()
        except Exception:
            if exc is None:
                raise
            LOG.exception("Shutdown also failed; check physical outputs")

    def _target(self, name):
        if name not in self.config["channels"]:
            raise ValueError(f"Unknown PZT name: {name}")
        ch = self.config["channels"][name]
        if ch["device"] not in self.sessions:
            raise AFGError("Not connected")
        return self.sessions[ch["device"]], ch

    def _write(self, session, command, *, shutdown=False):
        if not shutdown and self.stop_requested.is_set():
            raise AFGError("Stop requested. Apply settings again to resume.")
        if not shutdown:
            self._check_session(session)
        LOG.info("WRITE %s %s", getattr(session, "resource_name", ""), command)
        session.write(command)

    def _query(self, session, command):
        self._check_session(session)
        result = session.query(command).strip()
        LOG.debug("QUERY %s -> %s", command, result)
        return result

    def _check_session(self, session):
        if id(session) in self._unreliable_sessions:
            raise AFGError("Incomplete batch response; further queries and enabling are blocked. "
                           "Check physical outputs, restart the controller and reconnect USB.")

    def _query_many(self, session, commands):
        """Read the same scalar queries using fewer USB writes.

        AFG1062 firmware V1.0.9 and V1.1.0 were tested with complete headers
        separated by ';', WITHOUT leading colons, and one read per response.
        Their framing differs from the programmer manual's SCPI examples.
        Do not silently retry after a partial reply: unread data could be
        mistaken for a later output-state confirmation.
        """
        if self.config.get("query_batch_mode", "sequential") == "sequential":
            return [self._query(session, command) for command in commands]
        self._check_session(session)
        if not commands or any(not c.endswith("?") or any(x in c for x in ";\r\n") for c in commands):
            raise ValueError("Batch reads require separate scalar queries")
        try:
            session.write(";".join(commands))
            values = []
            for command in commands:
                value = session.read().strip()
                if not value or any(x in value for x in ";\r\n"):
                    raise AFGError(f"Unexpected batch response framing: {value!r}")
                LOG.debug("BATCH QUERY %s -> %s", command, value)
                values.append(value)
            return values
        except BaseException:
            self._unreliable_sessions.add(id(session))
            self.prepared.clear()
            raise

    def _errors(self, session, initial=()):
        """Check already-read error replies, then drain remaining errors if needed.

        A checked batch may have consumed two error entries. Preserve BOTH in
        the exception rather than losing the second entry or clearing the queue.
        """
        errors = []
        def code_of(value):
            try:
                return int(value.split(",", 1)[0])
            except ValueError as e:
                raise AFGError(f"Invalid error response: {value!r}") from e
        for value in initial:
            if code_of(value) != 0:
                errors.append(value)
        if initial and not errors:
            return
        for _ in range(65):
            value = self._query(session, "SYSTem:ERRor:NEXT?")
            code = code_of(value)
            if code == 0:
                if errors:
                    raise AFGError("Instrument error queue: " + "; ".join(errors))
                return
            errors.append(value)
        raise AFGError("Instrument error queue did not empty: " + "; ".join(errors))

    @staticmethod
    def _state(value):
        v = value.upper().strip()
        if v in ("1", "ON"):
            return True
        if v in ("0", "OFF"):
            return False
        raise AFGError(f"Invalid ON/OFF response: {value!r}")

    def read_channel(self, name, *, check_errors=False) -> dict:
        with self.lock:
            s, ch = self._target(name)
            prefix = f"SOURce{ch['channel']}:"
            q = lambda tail: self._query(s, prefix + tail + "?")
            tails = ["FUNCtion:SHAPe", "FREQuency:FIXed",
                     "VOLTage:LEVel:IMMediate:AMPLitude", "VOLTage:LEVel:IMMediate:OFFSet",
                     "PHASe:ADJust", "FREQuency:MODE", "BURSt:STATe", "FREQuency:CONCurrent"]
            commands = [f"OUTPut{ch['channel']}:STATe?", f"OUTPut{ch['channel']}:IMPedance?"]
            commands += [prefix + tail + "?" for tail in tails]
            commands += [prefix + m + ":STATe?" for m in MODULATIONS]
            if check_errors:
                # All replies are consumed before validation or any ON write.
                values = self._query_many(s, ["SYSTem:ERRor:NEXT?", *commands, "SYSTem:ERRor:NEXT?"])
                self._errors(s, initial=(values[0], values[-1]))
                values = values[1:-1]
            else:
                values = self._query_many(s, commands)
            result = dict(
                output=self._state(values[0]), load_ohm=float(values[1]), waveform=values[2],
                frequency_hz=float(values[3]), amplitude_vpp=float(values[4]), offset_v=float(values[5]),
                phase_deg=math.degrees(float(values[6])) % 360,
                frequency_mode=values[7], burst_enabled=self._state(values[8]),
                concurrent=self._state(values[9]),
                modulation={m: self._state(v) for m, v in zip(MODULATIONS, values[10:])},
            )
            if result["frequency_mode"].upper().startswith("SWE"):
                v = self._query_many(s, [prefix + tail + "?" for tail in (
                    "FREQuency:STARt", "FREQuency:STOP", "SWEep:TIME", "SWEep:SPACing", "SWEep:SOURce")])
                result.update(sweep_start_hz=float(v[0]), sweep_stop_hz=float(v[1]),
                              sweep_time_s=float(v[2]), sweep_spacing=v[3], trigger_source=v[4])
            if result["burst_enabled"]:
                result["burst_mode"] = q("BURSt:MODE")
                if result["burst_mode"].upper().startswith("TRIG"):
                    v = self._query_many(s, [prefix + "BURSt:NCYCles?", prefix + "BURSt:SOURce?"])
                    result.update(burst_cycles=float(v[0]), trigger_source=v[1])
            if check_errors and (result["frequency_mode"].upper().startswith("SWE") or result["burst_enabled"]):
                self._errors(s)  # Include any errors from the extra mode queries.
            return result

    def status(self):
        with self.lock:
            return {name: self.read_channel(name) for name in self.config["channels"]}

    def _verify(self, name, target, actual):
        _, ch = self._target(name)
        target.validate(ch)
        checks = {
            "load_ohm": (50, 0.1), "offset_v": (0, 0.00001),
            "amplitude_vpp": (target.amplitude_vpp, 0.00051),
        }
        if target.mode != "sweep":
            checks["frequency_hz"] = (target.frequency_hz, max(1e-6, target.frequency_hz * 1e-9))
        for key, (expected, tolerance) in checks.items():
            if not math.isfinite(actual[key]) or abs(actual[key] - expected) > tolerance:
                raise AFGError(f"{name}: {key} readback {actual[key]!r} != requested {expected!r}")
        if actual["amplitude_vpp"] > min(ch["max_input_vpp"], HARD_MAX_VPP) + 1e-9:
            raise AFGError(f"{name}: amplitude readback exceeds limit")
        phase = (target.phase_deg + ch.get("phase_correction_deg", 0)) % 360
        if not math.isfinite(actual["phase_deg"]) or abs((actual["phase_deg"] - phase + 180) % 360 - 180) > 0.02:
            raise AFGError(f"{name}: phase readback mismatch (query unit is RAD)")
        if not actual["waveform"].upper().startswith("SIN") or actual["concurrent"] or any(actual["modulation"].values()):
            raise AFGError(f"{name}: waveform/modulation/coupling readback is not the requested sine configuration")
        sweep = actual["frequency_mode"].upper().startswith("SWE")
        if sweep != (target.mode == "sweep") or actual["burst_enabled"] != (target.mode in ("burst", "gated")):
            raise AFGError(f"{name}: mode readback mismatch")
        if target.mode == "sweep":
            for key in ("sweep_start_hz", "sweep_stop_hz", "sweep_time_s"):
                if not math.isclose(actual[key], getattr(target, key), rel_tol=1e-6, abs_tol=1e-6):
                    raise AFGError(f"{name}: {key} readback mismatch")
            if not actual["sweep_spacing"].upper().startswith(target.sweep_spacing):
                raise AFGError(f"{name}: sweep spacing mismatch")
        if target.mode == "burst":
            if not actual["burst_mode"].upper().startswith("TRIG") or actual["burst_cycles"] != target.burst_cycles:
                raise AFGError(f"{name}: burst mode/count mismatch")
        if target.mode == "gated" and not actual["burst_mode"].upper().startswith("GAT"):
            raise AFGError(f"{name}: gated mode mismatch")
        if target.mode in ("sweep", "burst") and not actual["trigger_source"].upper().startswith(target.trigger_source):
            raise AFGError(f"{name}: trigger source mismatch")

    def _configure(self, name, settings):
        s, ch = self._target(name)
        c = ch["channel"]
        p = f"SOURce{c}:"
        self.prepared.pop(name, None)
        self.touched = True
        self._write(s, f"OUTPut{c}:STATe OFF", shutdown=True)
        if self._state(self._query(s, f"OUTPut{c}:STATe?")):
            raise AFGError(f"{name}: failed to turn output OFF")
        self._errors(s)  # Report prior errors; never hide them with *CLS.
        commands = [p + "FREQuency:CONCurrent OFF", p + "BURSt:STATe OFF", p + "FREQuency:MODE CW"]
        commands += [p + m + ":STATe OFF" for m in MODULATIONS]
        commands += [p + "FUNCtion:SHAPe SINusoid", f"OUTPut{c}:IMPedance 50",
                     p + "VOLTage:LEVel:IMMediate:OFFSet 0V",
                     p + f"VOLTage:LEVel:IMMediate:AMPLitude {settings.amplitude_vpp:.12g}Vpp",
                     p + f"FREQuency:FIXed {settings.frequency_hz:.12g}Hz",
                     p + f"PHASe:ADJust {(settings.phase_deg + ch.get('phase_correction_deg', 0)) % 360:.12g}DEG"]
        if settings.mode == "sweep":
            commands += [p + f"FREQuency:STARt {settings.sweep_start_hz:.12g}Hz",
                         p + f"FREQuency:STOP {settings.sweep_stop_hz:.12g}Hz",
                         p + f"SWEep:TIME {settings.sweep_time_s:.12g}s",
                         p + f"SWEep:SPACing {settings.sweep_spacing}",
                         p + f"SWEep:SOURce {settings.trigger_source}", p + "FREQuency:MODE SWEep"]
        elif settings.mode == "burst":
            commands += [p + "BURSt:MODE TRIGgered", p + f"BURSt:NCYCles {settings.burst_cycles}",
                         p + f"BURSt:SOURce {settings.trigger_source}", p + "BURSt:STATe ON"]
        elif settings.mode == "gated":
            # BURSt:SOURce is documented ONLY for triggered burst; do not send it here.
            commands += [p + "BURSt:MODE GATed", p + "BURSt:STATe ON"]
        for command in commands:
            self._write(s, command)
        self._write(s, "*WAI")
        self._errors(s)
        actual = self.read_channel(name)
        self._verify(name, settings, actual)
        if actual["output"]:
            raise AFGError(f"{name}: output unexpectedly enabled during configuration")
        self.prepared[name] = settings
        return actual

    def apply(self, settings_by_name: dict[str, ChannelSettings]) -> dict:
        """Validate the WHOLE request, configure sequentially, leave outputs OFF."""
        self.stop_requested.clear()  # An explicit new apply resumes after a prior stop.
        with self.lock:
            if not settings_by_name:
                raise ValueError("No channels supplied")
            for name, setting in settings_by_name.items():
                _, ch = self._target(name)
                setting.validate(ch)
            try:
                # Turn off ALL selected channels before reconfiguration begins.
                for name in settings_by_name:
                    s, ch = self._target(name)
                    self.touched = True
                    self._write(s, f"OUTPut{ch['channel']}:STATe OFF", shutdown=True)
                    if self._state(self._query(s, f"OUTPut{ch['channel']}:STATe?")):
                        raise AFGError(f"{name}: failed to turn output OFF")
                # Either channel may have frequency coupling left ON at the panel.
                for alias in {self.config['channels'][n]['device'] for n in settings_by_name}:
                    for c in (1, 2):
                        self._write(self.sessions[alias], f"SOURce{c}:FREQuency:CONCurrent OFF")
                return {name: self._configure(name, setting) for name, setting in settings_by_name.items()}
            except BaseException:
                failures = self.all_off()
                if failures:
                    LOG.error("Physical OFF not verified: %s", failures)
                raise

    def set_channel(self, name, *, freq_hz, vpp, phase_deg=0, on=False):
        setting = ChannelSettings(frequency_hz=freq_hz, amplitude_vpp=vpp, phase_deg=phase_deg)
        self.apply({name: setting})
        if on:
            self.set_output(name, True)

    def set_output(self, name, on: bool, *, return_state=False):
        """Switch and confirm output, checking fresh parameters before ON.

        Normally returns a bool. With return_state=True, returns the verified
        ON parameters plus confirmed output state, or just output state for OFF.
        This lets the GUI avoid a second full parameter read after switching.
        Timings in the log measure software/USB work, not analog output onset.
        """
        with self.lock:
            started = perf_counter()
            if type(on) is not bool:
                raise ValueError("on must be a bool")
            s, ch = self._target(name)
            if on and name not in self.prepared:
                raise AFGError(f"Apply {name} settings successfully before enabling output")
            try:
                state = {}
                if on:
                    state = self.read_channel(name, check_errors=True)
                    self._verify(name, self.prepared[name], state)
                checked = perf_counter()
                self.touched = True
                self._write(s, f"OUTPut{ch['channel']}:STATe {'ON' if on else 'OFF'}", shutdown=not on)
                written = perf_counter()
                replies = self._query_many(s, [f"OUTPut{ch['channel']}:STATe?", "SYSTem:ERRor:NEXT?"])
                self._errors(s, initial=(replies[1],))
                actual = self._state(replies[0])
                if actual != on:
                    raise AFGError(f"{name}: output switch readback mismatch")
                state["output"] = actual
                LOG.info("SWITCH %s %s precheck_ms=%.1f write_ms=%.1f confirmation_ms=%.1f total_ms=%.1f",
                         name, "ON" if on else "OFF", (checked-started)*1000, (written-checked)*1000,
                         (perf_counter()-written)*1000, (perf_counter()-started)*1000)
                return state if return_state else actual
            except BaseException:
                self.all_off()
                raise

    def request_stop(self):
        """Can be called immediately by a GUI while another operation holds the lock."""
        self.stop_requested.set()

    def all_off(self) -> list[str]:
        """Best effort across ALL connected outputs, even if an earlier one fails.

        Returns failures rather than hiding them. This is software shutdown,
        not a hardware emergency-stop or a guarantee during USB disconnection.
        """
        self.request_stop()
        with self.lock:
            self.prepared.clear()
            failures = []
            for alias, s in self.sessions.items():
                for c in (1, 2):
                    try:
                        self._write(s, f"OUTPut{c}:STATe OFF", shutdown=True)
                        if self._state(self._query(s, f"OUTPut{c}:STATe?")):
                            raise AFGError("output is still ON")
                    except Exception as e:
                        failures.append(f"{alias}/CH{c}: {e}")
            return failures

    def align_phase(self, alias: str):
        """Align the two channels WITHIN one AFG; not cross-instrument sync."""
        with self.lock:
            if alias not in self.sessions:
                raise AFGError("Unknown or disconnected device")
            s = self.sessions[alias]
            if any(self._state(self._query(s, f"OUTPut{c}:STATe?")) for c in (1, 2)):
                raise AFGError("Turn both outputs OFF before phase alignment")
            self.touched = True
            self._write(s, "SOURce1:PHASe:INITiate", shutdown=True)
            self._errors(s)

    def trigger(self, alias: str):
        """Device-wide *TRG; may trigger both armed MAN channels on that AFG."""
        with self.lock:
            if alias not in self.sessions:
                raise AFGError("Unknown or disconnected device")
            eligible = False
            try:
                for name, ch in self.config["channels"].items():
                    if ch["device"] != alias:
                        continue
                    actual = self.read_channel(name)
                    if actual["output"]:
                        if name not in self.prepared:
                            raise AFGError(f"{name} is ON but was not prepared by this program")
                        setting = self.prepared[name]
                        self._verify(name, setting, actual)
                        eligible |= setting.mode in ("sweep", "burst") and setting.trigger_source == "MAN"
                if not eligible:
                    raise AFGError("No enabled MAN sweep/burst channel on this device")
                self._errors(self.sessions[alias])
                self._write(self.sessions[alias], "*TRG")
                self._errors(self.sessions[alias])
            except BaseException:
                self.all_off()
                raise


def main():
    parser = argparse.ArgumentParser(description="AFG1062 four-channel USB control; default commands are read-only")
    parser.add_argument("--config", type=Path, default=BASE / "config.json")
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("discover")
    sub.add_parser("status")
    sub.add_parser("off")
    a = sub.add_parser("apply", help="Read settings JSON, configure all listed channels and keep outputs OFF")
    a.add_argument("settings", type=Path)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    try:
        cfg = load_config(args.config)
        if args.command == "discover":
            print(json.dumps(discover(cfg), indent=2, ensure_ascii=False))
        else:
            with Rig(cfg) as rig:
                if args.command == "status":
                    result = dict(devices=rig.identities, channels=rig.status())
                elif args.command == "off":
                    result = rig.all_off()
                    if result:
                        raise AFGError("OFF not verified: " + "; ".join(result))
                else:
                    raw = json.loads(args.settings.read_text(encoding="utf-8-sig"))
                    result = rig.apply({name: ChannelSettings(**value) for name, value in raw.items()})
                print(json.dumps(result, indent=2, ensure_ascii=False))
    except (Exception, KeyboardInterrupt) as e:
        parser.exit(1, f"ERROR: {e}\n")


if __name__ == "__main__":
    main()
