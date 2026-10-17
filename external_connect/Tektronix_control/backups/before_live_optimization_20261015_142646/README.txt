AFG1062 Four-Channel PZT Controller
Version: 2026-10-15 (English interface; reduced output-switch latency)

1. Amplifiers and wiring

The confirmed setup uses four High Wave 3.2 amplifiers. On each amplifier,
the two red 4 mm output terminals connect to the two electrodes of one PZT.
This is bridge operation. The black cable is connected to the red Out 2
terminal: it is an active output, not ground.
Each amplifier has one input. Out 1 and Out 2 are fixed 180 degrees apart;
they are not independent signal sources.

Amplifier manual: printed pages 4-7 and 14-15 (PDF pages 5-8 and 15-16):
  - Input impedance: 50 ohms when powered on. Maximum input: 4 Vpp.
  - Approximate single-ended voltage gain: 7; differential bridge gain: 14.
    Actual gain depends on frequency, load and waveform.
  - Stated bridge capability is approximately 50-60 Vpp and 0.5 A. These
    values are not guaranteed simultaneously at every frequency and load.
  - Each 4 mm output has approximately 1 ohm output impedance. The 50 ohm
    BNC outputs are distinct from the amplifier's 50 ohm input.
  - Do not connect the bridge load to output ground or the enclosure.
    Measure across the PZT with a suitable differential probe. Do not attach
    a conventional grounded oscilloscope clip to either active output.
  - A PZT is a capacitive load. Staying within the input voltage limit does
    not by itself guarantee acceptable current, heating or distortion.
  - The manual requires resistive termination of unused outputs. The BNC
    outputs appear unused in the supplied photograph. The manual does not
    specify their termination resistance or power rating for this bridge
    arrangement; confirm with the manufacturer before adding a load.

The controller uses a sine carrier, a 50 ohm load setting and zero DC offset.
The default software input limit is 0.5 Vpp per channel; the initial setting
is 0.1 Vpp. Ideal bridge output estimates are about 7 Vpp and 1.4 Vpp,
respectively. These estimates and initial settings have not been validated
as safe operating points for the specific PZT, impedance or temperature rise.
Each max_input_vpp value can be adjusted in config.json. The program always
enforces a 4 Vpp ceiling. Frequency is limited to the amplifier's nominal
10 MHz bandwidth, rather than the generator's 60 MHz capability.

2. Instrument and channel mapping

The following USB identities were detected on this computer:
  A = AFG1062, USB serial 2147075, VID 0699, PID 0353
  B = AFG1062, USB serial 2439003, VID 0699, PID 0353

  PZT1 = A / CH1 -> Amplifier 1 Input -> two red output terminals -> PZT1
  PZT2 = A / CH2 -> Amplifier 2 Input -> two red output terminals -> PZT2
  PZT3 = B / CH1 -> Amplifier 3 Input -> two red output terminals -> PZT3
  PZT4 = B / CH2 -> Amplifier 4 Input -> two red output terminals -> PZT4

Instruments are selected by USB serial number, not connection order.
PZT1-PZT4 are software labels; verify their physical locations using the cables.
Communication uses USBTMC/VISA, not a COM serial port. No baud rate, serial
data format or COM port number is required. After enumerating VISA resources,
the program checks the model and serial number using *IDN?. Firmware versions
are read from the instruments rather than assumed.

3. Installation and startup

Last read-only hardware check (2026-10-15): both USB instruments responded
successfully. A reports firmware V1.0.9 and B reports V1.1.0. All four channels
were read successfully. The existing installation can run start_gui.cmd directly.
The following driver instructions apply to a new or broken installation.

The official TekVISA 5.11.1 installer is included. Its Tektronix, Incorporated
digital signature was verified by Windows.

Right-click drivers\TekVISA_5.11.1.exe, choose Run as administrator and complete
the installation wizard. Save your work before any requested restart; the
controller does not restart the computer automatically. The installer release
notes list Windows 10; this computer reported Windows build 26200. If that
package does not support this system, use an NI-VISA release supporting your
Windows version and include USB/USBTMC support:
https://www.ni.com/en/support/downloads/drivers/download.ni-visa.html

The dedicated Python environment is in .venv. To recreate it, open PowerShell
in the project directory and run:
  powershell -NoProfile -ExecutionPolicy Bypass -File .\setup_python.ps1

After installing the driver:
  1. Double-click check_connection.cmd. It should list both AFG1062 instruments,
     their serial numbers and firmware versions.
  2. Double-click start_gui.cmd and select Connect Devices (Read Only).
  3. Enter frequency in MHz, amplifier input in Vpp, and phase in degrees.
     Choose Continuous Sine, Sweep, N-Cycle Burst or External Gate.
  4. Select Apply (OFF). After successful readback, select ON.
  5. Apply again after editing any field. ON uses the last successfully applied
     settings, not uncommitted edits in the interface.
  6. All Outputs OFF clears the prepared settings. Apply again before enabling.

Connecting and reading status preserve existing instrument output states.
An output already enabled on the front panel stays enabled when connecting.
Apply turns the selected outputs off, configures the mode and verifies readback;
it does not enable output. ON rechecks actual instrument settings. If front-panel
changes cause a mismatch, the program refuses to enable output and attempts to
turn all outputs off. Run only one controller at a time and avoid simultaneous
changes from another program or operator.

Log file: control.log.
Read-only diagnosis: .venv\Scripts\python.exe -X utf8 diagnose.py
connection_report.json is a snapshot of a past check, not a live status display.

4. Sweep, burst, gating and phase

Continuous Sine: independent frequency, Vpp, phase and output state per channel.
Sweep: independent sine-wave endpoints; ascending or descending; LIN or LOG;
sweep time from 1 ms to 500 s.
N-Cycle Burst: 1-1,000,000 cycles, with MAN, TIM or EXT trigger source.
External Gate: uses the instrument's Fsk/Ext Trig In input; EXT is selected
automatically. BURSt:SOURce is only documented for triggered burst, so the
controller does not send it in gated mode.

MAN: enable output to arm the channel, then select A: Trigger MAN or B: Trigger MAN.
*TRG is a device-wide command. It can trigger all enabled MAN sweep/burst channels
on that instrument.
TIM: internal repeated triggering. Burst interval uses the current front-panel
setting; the programmer manual does not document a command to set that interval.
EXT: external TTL triggering. Each instrument shares one external trigger input;
this does not provide four fully independent hardware gates. For independent
gating, characterize instrument behavior and design the additional hardware.

A: Align CH1/CH2 and B: Align CH1/CH2 align channels within the selected instrument
only. First turn off both outputs on that instrument. Phase writes explicitly
use DEG. The phase query returns RAD, which the program converts to degrees.
Set phase_correction_deg in config.json to each channel's measured compensation:
  transmitted phase = (requested phase + phase correction) modulo 360.
Changing frequency, restarting, realigning or changing wiring may require
recalibration.

Clock wiring: A Ref Clk Out -> B Ref Clk/Counter In.
Set A to Internal. On B, open the Utility system settings, find CLK Ref and
select External. See the AFG1000 user manual, printed page 51 (PDF page 61),
and confirm the settings on the actual front panel. The reference input is
shared with the frequency-counter input, so these input functions cannot be
used simultaneously.

The controller does not send undocumented ROSC commands and cannot query or
prove external clock lock. A shared clock avoids drift from separate free-running
clocks but does not establish a repeatable starting phase between instruments.
Cross-instrument coherent control requires a common timebase, a common hardware
trigger where needed, and measurement/calibration of the actual PZT phase.
USB commands are sequential, and this program performs multiple readback checks.
It does not guarantee simultaneous switching or a 100 ms real-time deadline.

5. Python API and command line

This example configures the channels and leaves outputs off:

  from afg_control import Rig, ChannelSettings
  with Rig() as rig:
      rig.apply({
          "PZT1": ChannelSettings(frequency_hz=2.5e6, amplitude_vpp=0.1, phase_deg=0),
          "PZT2": ChannelSettings(frequency_hz=2.5e6, amplitude_vpp=0.1, phase_deg=90),
          "PZT3": ChannelSettings(frequency_hz=490e3, amplitude_vpp=0.1, phase_deg=180),
          "PZT4": ChannelSettings(frequency_hz=490e3, amplitude_vpp=0.1, phase_deg=270)
      })
      print(rig.status())
      # After checking the experiment, explicitly enable an output:
      # rig.set_output("PZT1", True)
      # rig.set_output("PZT1", False)

Sweep:
  ChannelSettings(mode="sweep", sweep_start_hz=2e6, sweep_stop_hz=3e6,
                  sweep_time_s=0.1, sweep_spacing="LIN", trigger_source="MAN")
Burst:
  ChannelSettings(mode="burst", burst_cycles=100, trigger_source="MAN")
Gate:
  ChannelSettings(mode="gated", trigger_source="EXT")

Read-only discovery:
  .venv\Scripts\python.exe -X utf8 afg_control.py discover
Read-only status:
  .venv\Scripts\python.exe -X utf8 afg_control.py status
Turn all outputs off:
  .venv\Scripts\python.exe -X utf8 afg_control.py off
Apply a profile saved from the interface, keeping outputs off:
  .venv\Scripts\python.exe -X utf8 afg_control.py apply settings.json

example_control.py configures a channel but leaves it off by default. Adding
--enable explicitly enables PZT1 at 0.1 Vpp for two seconds. Leaving the Rig
context attempts to turn off all four outputs, including when an exception occurs.
Outputs are not intended to remain on after the example exits.
USB disconnection, forced process termination or a system crash can prevent
software shutdown; use physical switches to stop the experiment when necessary.
all_off() returns a list of outputs whose shutdown could not be verified.
Check that return value rather than treating a failed shutdown as successful.

6. Output-switch latency

The previous version issued 19 separate queries BEFORE sending ON in continuous
mode. Each query took about 35 ms on these instruments, delaying the ON command
by about 0.67 s. It then read the whole channel again to update the interface.

config.json now selects query_batch_mode = "legacy_lines". This groups the
same 17 base parameter queries into one USB write and reads all 17 responses.
It keeps the two separate error-queue checks before ON, all parameter validation,
and output-state/error confirmation after switching. No voltage, load, phase or
mode check was removed. The GUI reuses the fresh pre-ON parameter readback and
confirmed output state; OFF requires only output-state and error confirmation.
"Last read" parameters remain a snapshot, not a live monitor of front-panel edits.
The GUI checks for completed background work every 20 ms instead of 100 ms.

Measured on both instruments, all four channels, three samples per mode:
  Full channel read:        about 600-603 ms -> about 112-113 ms.
  Pre-ON query stage:       671-673 ms -> 183-184 ms (about 73% shorter).
These are measurements of the actual read-only query sequences. Outputs were
not enabled during benchmarking. They do not measure click-to-analog-output
delay, relay settling, PZT response or precise switching synchronization.
Sweep and burst require additional mode-specific queries and can take longer.

Firmware compatibility: on the tested V1.0.9 and V1.1.0 units, batch messages
must contain full headers separated by semicolons WITHOUT leading colons, and
each reply must be read separately. This differs from the manual's examples.
If changing firmware or hardware, select "sequential" and restart the program
until the batch behavior is validated. An absent query_batch_mode also selects
sequential reads. An incomplete/malformed batch blocks further queries and ON;
the controller still attempts OFF commands during shutdown but cannot verify
those states on a desynchronized connection. Check the physical outputs, close
the controller and reconnect USB before retrying. No automatic retry is used.

Each successful ON/OFF operation writes precheck_ms, write_ms, confirmation_ms
and total_ms to control.log. These are software/USB timings, not measured analog
transition times. Lowering timeout_ms will not speed up successful responses;
it only changes how long a missing response is allowed to take.

Repeat the read-only benchmark with other control programs closed:
  .venv\Scripts\python.exe -X utf8 benchmark_reads.py
Results are written to latency_report.json. This sends identity/status/parameter
queries and reads the error queue (consuming errors), without changing settings
or outputs. It stops on errors instead of clearing them. For a fixed timing
deadline or synchronized onset, use and characterize hardware triggering/gating.

7. Validation scope

Offline tests cover serial mapping, whole-request validation, voltage/frequency
and NaN rejection, output-off sequencing, phase conversion, sweep/burst mode
changes, gated command restrictions, instrument errors, front-panel changes,
and attempts to turn off remaining channels after a USB failure. The same control
cases run in sequential and batched modes. Additional tests cover partial batch
responses, unexpected response framing, retained pre-ON checks, and GUI updates
without redundant hardware queries. All 52 offline tests passed for this update.
Run the tests with:
  .venv\Scripts\python.exe -X utf8 -m unittest -v test_afg_control.py test_afg_gui.py

Offline tests do not verify physical instruments or analog outputs. Read-only
batch query acceptance and timing were also checked on both connected instruments.
Analog voltage, phase, clock lock, loaded operation, physical switching/camera
timing and PZT acoustic pressure still require experiment-specific measurements.

References:
  User-supplied High Wave 3.2 German manual, revision 2020-02-24.
  Tektronix AFG1000 Programmer Manual 077112902:
  https://download.tek.com/manual/AFG1000-Programmer-Manual-EN-077112902.pdf
  Tektronix AFG1000 Quick Start User Manual 077113001:
  https://download.tek.com/manual/AFG1000-Quick-Start-User-Manual-EN-077113001.pdf
  Official TekVISA download:
  https://www.tek.com/en/support/software/driver/tekvisa-connectivity-software-v5111
