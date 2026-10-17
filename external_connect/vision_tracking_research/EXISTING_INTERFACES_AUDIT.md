# Existing interfaces audit — 2026-10-16

Scope: read-only source review and sample-video replay. No AFG module was imported, no VISA session opened, no instrument command sent, no camera setting changed. New files and dependencies are under `vision_tracking_research/`. No applicable AGENTS.md was found under `/home/lzy/projects/external_connect` or its ancestors. Both existing directories are not Git checkouts; their identity is the path and SHA-256 manifest, not an invented Git commit.

## Inventory and runtime

`evidence/existing_tree.json` records all 4,657 existing files with paths and sizes, including environments, vendor tools, backups and previous data. `evidence/protected_sources_before.json` hashes 46 current source/config/launcher files. Backups, environments and recordings were inspected as metadata only and were not edited. The final verification compares the manifest again.

- `CHRONOS_control/`: `chronos_launcher.py`, `test_chronos.py`, `chronos_common.py`, `chronos_rtsp_control.py`, `chronos_http_preview.py`, connection/latency/stability/calibration tools, tests, JSON config, launch/setup scripts; `camera_rtsp/` (existing sender patch and evidence); `tools/vlc-3.0.24/`; `data/` (sessions/images/logs); `backups/`; Windows `.venv/`.
- `Tektronix_control/`: `afg_control.py`, `afg_gui.py`, examples/diagnosis/tests/benchmark, config and previous latency/readback reports, Windows launch/setup scripts; `drivers/`, `backups/`, Windows `.venv/`.
- Both `.venv/pyvenv.cfg` identify **Windows CPython 3.14.2** (`C:\Python314`), not a WSL runtime. Chronos pins OpenCV 5.0.0.93 and NumPy 2.5.4; camera-maintenance extra is Paramiko 4.0.0. Tektronix pins PyVISA 1.16.2. These have not been upgraded or replaced.
- Research uses a separate WSL Python 3.12 venv. Norfair 2.3 requires NumPy <2; therefore do not install it into the current camera environment. The compatibility choice here is NumPy 1.26.4/OpenCV 4.11.0.86, recorded in `evidence/dependencies_classical.txt`.

All source paths below are relative to `/home/lzy/projects/external_connect/` and refer to audited unchanged files.

## Chronos: reuse the real receiver

| Actual source and line | Interface / finding |
|---|---|
| `CHRONOS_control/chronos_launcher.py:15` | `main()` creates the existing GUI; `run()` at line 54 launches scripts and prevents a second preview launched by this GUI. Launching the GUI is not a read-only audit: `save_settings()` line 43 writes config. Not launched. |
| `CHRONOS_control/test_chronos.py:14` | `main()` is the continuous preview/measurement consumer. Lines 101–115 call `stream.snapshot()`, reject stale/repeated frames and count skipped sequence numbers. Display is throttled separately. |
| `CHRONOS_control/chronos_common.py:30` | `load_config()` validates IP/transport/timeouts; current IP 192.168.0.10, UDP ports 5000/5001, read timeout 1.5 s, reconnect interval 0.5 s. |
| `CHRONOS_control/chronos_common.py:66` | `video_options()` supplies FFmpeg UDP/low-latency flags. It deliberately retains codec probing data: comments document earlier failed `nobuffer`/tiny-probe attempts. Do not blindly replace these tested settings. |
| `CHRONOS_control/chronos_common.py:80` | `capture_parameters()` opens FFmpeg with timeouts and one decoder thread in low-latency mode. |
| `CHRONOS_control/chronos_common.py:113` | **`LatestFrameStream`**, with `start()` line 147, `_run()` line 155, **`snapshot()` line 240**, `refresh()` line 244, `stats()` line 251, `stop()` line 269. No `read_latest_frame()` exists. |
| `CHRONOS_control/chronos_common.py:219` | Actual tuple: `(frame_sequence, receive_monotonic, receive_utc_string, frame_ndarray)`. OpenCV color frames are H×W×3 BGR uint8. The monotonic timestamp is **after decoding**, not exposure or network arrival. |
| `CHRONOS_control/chronos_common.py:120` | One background producer and lock-protected latest slot. Slow consumers overwrite older decoded frames. This bounds Python storage, **not** camera encoder, UDP, jitter or decoder buffering. `snapshot()` returns the same array reference; consumers must not annotate it in place. |
| `CHRONOS_control/chronos_common.py:222` | Receive exception clears latest, records a bounded failure history, releases backend on its own thread, then reconnects. `refresh()` clears latest immediately; close waits for configured backend timeouts. |
| `CHRONOS_control/chronos_common.py:277` | `SnapshotWriter` is a separate bounded queue (default 4) for PNG saving; it is not the frame-delivery queue. |
| `CHRONOS_control/chronos_rtsp_control.py:35` | `ChronosRtspControl` repairs request framing with fresh control TCP connections. RTP video goes directly to the receiver; this relay is not a video buffer. Preserve this device-specific workaround. |
| `CHRONOS_control/chronos_http_preview.py:11` | `HttpJpegCapture.read()` line 28 performs sequential GET `/cgi-bin/screenCap` at default 4 Hz, with no pipelining and size checks. This mode is an occasional JPEG preview, not a 60 fps acquisition path. |
| `CHRONOS_control/analyze_stability.py:10` | Existing `measure_frame(frame, roi, polarity, threshold, min_area)` already uses OpenCV threshold/contours/moments and selects the largest complete dark region. Reusable for an isolated known object; it is not a multi-bubble detector. |
| `CHRONOS_control/chronos_launcher.py:93` | Saved-image measurement ROI entry exists. There is no live sensor-ROI setter in the main receiver; research ROI is a downstream array crop. |
| `CHRONOS_control/calibrate_scale.py:74` | Existing scale calibration assumes unchanged optics/ROI/resolution/scaling. No matching calibration was supplied for the cropped sample; all research output remains px and px per replay second. |
| `CHRONOS_control/measure_latency.py:9` | `latency_summary()` and the displayed optical clock support approximate scene-to-display checks. This does not establish camera-to-controller latency. |

`evidence/live_probe.json`: bounded read-only TCP attempts to ports 80 and 554 both timed out. **LIVE_CAMERA_NOT_AVAILABLE** in that original WSL session. A subsequent [2026-10-16 follow-up](LIVE_CAMERA_FOLLOWUP_20261016.md) verified reachable ports, approximately 60 fps actual Windows reception and WSL HTTP snapshots. Native WSL RTSP/UDP remains unresolved. No camera firmware or sender patch was applied. Previous saved session `data/capture_20261015_185340_460247/summary.json` reports 185,171 received frames and 60.0386 fps over the receiving interval, then 583 failures and a long outage. Those are historical logs, not current live performance and not a guarantee of uninterrupted streaming.

The research replay subclass changes only the producer to PyAV file decoding; it inherits the original `start()`, `snapshot()`, lock/slot and `stop()`. It returns the actual tuple and stores replay PTS separately. RTSP delivery was **not** simulated or reimplemented. H.264 decoder reference/reorder buffering remains inherent; no tracking algorithm receives future presentation frames. File preroll before each selected clip is excluded from steady-state decode latency; cold camera connection time is not measured.

## Tektronix: architecture only

| Actual source and line | Interface / constraints |
|---|---|
| `Tektronix_control/afg_control.py:40` | `ChannelSettings`: frequency Hz, **generator/amplifier-input Vpp into 50 Ω**, phase degrees, continuous/sweep/burst/gated modes. These are not pressure Pa or measured PZT voltage. |
| `Tektronix_control/afg_control.py:53` | `ChannelSettings.validate()` checks finite values, modes, integer mV steps, configured 0.5 Vpp limit, 4 Vpp hard ceiling and 10 MHz maximum. These software limits do not calibrate experimental safety. |
| `Tektronix_control/afg_control.py:81` | `load_config()` validates unique A/B serials and four unique channel mappings. A=2147075, B=2439003; PZT1=A/CH1, PZT2=A/CH2, PZT3=B/CH1, PZT4=B/CH2. Spatial X/Y locations cannot be inferred from software labels. |
| `Tektronix_control/afg_control.py:132` | `discover()` enumerates USBTMC/VISA and queries `*IDN?`. `Rig.connect()` line 177 selects serial/model/VID/PID, independent of enumeration order. Not called in this task. |
| `Tektronix_control/afg_control.py:337` | `read_channel()` / `status()` line 381 read output and parameters; phase reply radians is converted to degrees. No calls made. |
| `Tektronix_control/afg_control.py:461` | `apply()` configures with outputs OFF and verifies. `apply_live()` line 489 preserves current ON/OFF but performs actual sequential writes and readback. Neither is safe to invoke merely for architecture testing. |
| `Tektronix_control/afg_control.py:568` | `set_channel()`, `set_output()` line 574, `align_phase()` line 636 and `trigger()` line 648 change hardware. Not invoked. |
| `Tektronix_control/afg_control.py:612` | `request_stop()` sets a thread event; `all_off()` line 616 attempts and verifies shutdown. `close()` line 226 performs shutdown after a control operation (`touched=True`). Physical USB loss/process death cannot guarantee OFF. |
| `Tektronix_control/afg_control.py:253` | `_write()` logs SCPI; `_query()` line 261 logs readback. `apply_live()` line 558 and `set_output()` line 604 log software/USB phase durations. These timings are not analog output onset timestamps. |
| `Tektronix_control/afg_gui.py:32` | Single `ThreadPoolExecutor` for hardware; line 33 has an unbounded completion-event Queue. `submit()` line 164 rejects work while busy, so it is not presently a per-frame visual-action queue. `poll()` line 174 runs nominally every 20 ms. |
| `Tektronix_control/afg_gui.py:364` | GUI `main()` writes `control.log` with wall-clock timestamps. `finish_close()` line 346 invokes Rig close. Not launched. |

The task identifies High Wave **3.3**; current code/config/README describe **3.2**. This discrepancy is recorded, not silently corrected or used to change limits. Actual task geometry is approximately **1.3 mm × 30 mm**, distinct from F10S simulation geometry (1.6 mm × 32.5 mm). Neither nominal diameter is a pixel-scale calibration.

## Minimal future integration (DESIGN, not hardware implementation)

1. Keep the camera producer and AFG driver intact. Attach a consumer to `LatestFrameStream.snapshot()`, copy only cropped arrays when annotation/modification is required, and reject duplicate IDs/stale packets before processing.
2. Deploy matching Python dependencies in an independent vision process. On Windows the existing receiver is already working; a Windows vision worker avoids making WSL networking the first integration obstacle. If vision remains WSL, transport frames or compact observations with explicit clock/session metadata and bounded latest-value buffers.
3. Separate observation messages from action proposals. The tested JSON is emitted to a **mock-only local ZeroMQ receiver**; no AFG API adapter is created here. Per-frame contours should eventually be on a slower diagnostics channel; compact state/quality belongs on the control channel.
4. Cross-OS monotonic clocks cannot be subtracted directly. Add clock-offset calibration/round-trip uncertainty or stamp everything within one host clock domain, a stream session ID, frame ID and action correlation ID. Current WSL loopback timings do not measure Windows↔WSL transport.
5. Future command admission must require explicit user control enable, fresh reliable target state, phase/voltage limits already in the driver, emergency OFF priority and a watchdog. A proposed safe stop is not implemented by this research mock. Do not queue one AFG command per camera frame; current synchronous USB transactions and readbacks have a separate rate limit.
