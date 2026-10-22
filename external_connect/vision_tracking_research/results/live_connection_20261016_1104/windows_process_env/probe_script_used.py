"""Bounded read-only test of the existing Chronos latest-frame receiver.

Run under either the working Windows camera venv or the research WSL venv.
Imports the receiver without modifying its source/configuration. Opens only its
own RTSP preview session; never sends camera setting or recording commands.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import socket
import subprocess
import sys
import time


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=True, allow_nan=False), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--seconds", type=float, default=30)
    args = parser.parse_args()
    if not 1 <= args.seconds <= 60:
        parser.error("Use a bounded duration between 1 and 60 seconds")
    source = Path(args.source).resolve()
    output = Path(args.output).resolve()
    sys.path.insert(0, str(source))
    from chronos_common import (LatestFrameStream, configure_video_environment,
                                load_config, save_image, video_options)
    cfg = load_config(source / "chronos_config.json")
    # Preserve the working test_chronos.py launcher contract: Windows FFmpeg
    # needs these options in the process environment, not only a later putenv.
    desired = video_options(cfg)
    if cfg["latency_mode"] != "fresh" and os.environ.get("OPENCV_FFMPEG_CAPTURE_OPTIONS") != desired:
        environment = os.environ.copy()
        environment["OPENCV_FFMPEG_CAPTURE_OPTIONS"] = desired
        return subprocess.call([sys.executable, "-B", str(Path(__file__).resolve()), *sys.argv[1:]],
                               env=environment)
    configure_video_environment(cfg)
    output.mkdir(parents=True, exist_ok=False)
    from check_connection import diagnose
    import numpy as np

    names = ("chronos_common.py", "chronos_rtsp_control.py", "chronos_http_preview.py",
             "check_connection.py", "chronos_config.json")
    def source_hashes():
        return {n: hashlib.sha256((source / n).read_bytes()).hexdigest() for n in names}

    before = source_hashes()
    write_json(output / "identity.json", {
        "utc": datetime.now(timezone.utc).isoformat(), "platform": platform.platform(),
        "python": sys.version, "executable": sys.executable, "source": str(source),
        "config_in_memory": cfg, "source_sha256_before": before,
        "probe_script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "process_capture_options": os.environ.get("OPENCV_FFMPEG_CAPTURE_OPTIONS"),
        "camera_setting_changes": False, "recording_started": False,
    })
    diagnostic = diagnose(cfg)
    write_json(output / "diagnostic.json", diagnostic)
    if not diagnostic["rtsp_port_554"]["reachable"]:
        raise RuntimeError("RTSP port is unreachable; no preview session opened")
    # Refuse to compete for the configured RTP/RTCP ports with another receiver.
    if cfg["rtsp_transport"] == "udp":
        for port in (cfg["udp_port"], cfg["udp_port"] + 1):
            with socket.socket(socket.AF_INET, socket.SOCK_DGRAM) as sock:
                if hasattr(socket, "SO_EXCLUSIVEADDRUSE"):
                    sock.setsockopt(socket.SOL_SOCKET, socket.SO_EXCLUSIVEADDRUSE, 1)
                sock.bind(("0.0.0.0", port))

    stream = LatestFrameStream(cfg)
    samples, frames, checkpoints = [], [], []
    selected_times = [0, args.seconds / 2, args.seconds - 2]
    next_selection, previous_seq = 0, 0
    previous_small, changed_pairs = None, 0
    started = time.monotonic()
    stream.start()
    print(json.dumps({"event": "started", "output": str(output), "seconds": args.seconds}), flush=True)
    next_progress = 5
    error = None
    try:
        while time.monotonic() - started < args.seconds:
            now = time.monotonic()
            elapsed = now - started
            packet = stream.snapshot()
            if packet is not None and packet[0] != previous_seq:
                now = time.monotonic()
                seq, received, stamp, frame = packet
                frames.append({"seq": seq, "receive_monotonic": received, "utc": stamp,
                               "observed_elapsed_s": elapsed,
                               "post_decode_age_ms": (now - received) * 1000,
                               "consumer_skipped_since_previous": seq - previous_seq - 1})
                previous_seq = seq
                small = frame[::16, ::16].copy()
                if previous_small is not None and not np.array_equal(small, previous_small):
                    changed_pairs += 1
                previous_small = small
                if next_selection < len(selected_times) and elapsed >= selected_times[next_selection]:
                    samples.append((seq, frame.copy()))
                    next_selection += 1
            if elapsed >= next_progress:
                stats = stream.stats()
                checkpoints.append({"elapsed_s": elapsed, **stats})
                print(json.dumps({"event": "progress", "elapsed_s": elapsed,
                                  "frames": stats["received_frames"], "failures": stats["failures"],
                                  "state": stats["state"], "last_error": stats["last_error"]}), flush=True)
                next_progress += 5
            time.sleep(0.001)
        before_stop = stream.stats()
    except Exception as exc:
        error = repr(exc)
        before_stop = stream.stats()
        raise
    finally:
        stream.stop()
        after_stop = stream.stats()
        after = source_hashes()
        for seq, frame in samples:
            save_image(output / f"frame_{seq:06d}.png", frame)
        ages = [f["post_decode_age_ms"] for f in frames]
        report = {
            "status": "LIVE_RTSP_DECODED" if frames else "NO_DECODED_FRAMES",
            "requested_seconds": args.seconds, "wall_seconds_including_stop": time.monotonic() - started,
            "before_stop": before_stop, "after_stop": after_stop,
            "consumer_unique_frames": len(frames), "changed_subsample_pairs": changed_pairs,
            "frame_shape_bgr": list(samples[0][1].shape) if samples else None,
            "post_decode_age_p95_ms": float(np.percentile(ages, 95)) if ages else None,
            "post_decode_age_note": "Consumer scheduling after decode; not sensor-to-display latency",
            "sensor_dropped_frames": "NOT_MEASURED; no sensor sequence or timestamp in this interface",
            "source_and_config_unchanged": before == after,
            "source_sha256_after": after, "checkpoints": checkpoints, "error": error,
            "limits": "Short live preview test only; no bubble detection accuracy or control validation",
        }
        write_json(output / "summary.json", report)
        with (output / "observed_frames.jsonl").open("w", encoding="utf-8") as f:
            for row in frames:
                f.write(json.dumps(row) + "\n")
        print(json.dumps({k: report[k] for k in ("status", "consumer_unique_frames", "frame_shape_bgr",
                                               "source_and_config_unchanged")}), flush=True)
    return 0 if frames and before == after else 2


if __name__ == "__main__":
    sys.dont_write_bytecode = True
    raise SystemExit(main())
