"""Chronos live preview, raw PNG snapshots, receive-rate logging and reconnect."""
import argparse
import csv
import os
import subprocess
from pathlib import Path
import sys
import time

from chronos_common import (DEFAULT_CONFIG, LatestFrameStream, configure_video_environment, load_config, video_options,
                            save_image, session_dir, utc_now, write_json)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=str(DEFAULT_CONFIG))
    parser.add_argument("--ip")
    parser.add_argument("--transport", choices=["tcp", "udp"])
    parser.add_argument("--headless", action="store_true", help="Run without a display window")
    parser.add_argument("--seconds", type=float, default=0, help="Session duration; 0 = until Q/Ctrl+C")
    parser.add_argument("--save-every", type=float, default=0, help="Auto-save a PNG every N seconds; 0 disables")
    parser.add_argument("--max-images", type=int, default=20, help="Limit automatic snapshots")
    parser.add_argument("--connect-wait", type=float, default=20, help="Fail if no frame arrives before this deadline")
    args = parser.parse_args()
    if args.seconds < 0 or args.save_every < 0 or args.max_images < 0 or args.connect_wait <= 0:
        parser.error("Time values and image limits must be nonnegative; connect-wait must be positive")
    config = load_config(args.config, args.ip)
    if args.transport:
        config["rtsp_transport"] = args.transport
    # The FFmpeg DLL on Windows reads its environment at process startup.
    # Passing env to a child also works when launching this .py directly.
    desired_options = video_options(config)
    if os.environ.get("OPENCV_FFMPEG_CAPTURE_OPTIONS") != desired_options:
        environment = os.environ.copy()
        environment["OPENCV_FFMPEG_CAPTURE_OPTIONS"] = desired_options
        try:
            return subprocess.call([sys.executable, str(Path(__file__).resolve()), *sys.argv[1:]], env=environment)
        except KeyboardInterrupt:
            return 130
    configure_video_environment(config)
    import cv2
    import numpy as np
    folder = session_dir(config, "capture")
    write_json(folder / "session.json", {"started_utc": utc_now(), "config": config, "arguments": vars(args),
               "timestamp_meaning": "PC time immediately after decoding, not camera exposure time",
               "latency_note": "Local frame age is NOT end-to-end scene-to-screen latency."})
    print(f"Chronos: rtsp://{config['camera_ip']} | {config['rtsp_transport']}", flush=True)
    print(f"Output: {folder}\nS: save raw PNG | Q/Esc: quit | Ctrl+C: stop", flush=True)
    stream = LatestFrameStream(config).start()
    started = time.monotonic()
    last_log = last_console = started
    next_save = started
    consumed = skipped = automatic = saved = 0
    logged_failures = 0
    last_sequence = 0
    first_shape = None
    result = 0
    reason = "normal"
    window = "Chronos 1.4 | S: save | Q: quit"
    csv_file = (folder / "stream_metrics.csv").open("w", newline="", encoding="utf-8-sig")
    fields = ["utc", "elapsed_s", "received_frames", "receive_fps_last_5s", "consumed_frames",
              "skipped_by_consumer", "seconds_since_last_frame", "connections", "failures", "state"]
    writer = csv.DictWriter(csv_file, fieldnames=fields)
    writer.writeheader()

    def snapshot(packet):
        nonlocal saved
        sequence, received, stamp, raw = packet
        path = folder / "images" / f"frame_{sequence:09d}_{saved:04d}.png"
        save_image(path, raw)
        write_json(path.with_suffix(".json"), {"sequence": sequence, "received_utc": stamp,
                   "received_monotonic_s": received, "width": raw.shape[1], "height": raw.shape[0],
                   "camera_ip": config["camera_ip"], "saved_utc": utc_now()})
        saved += 1
        print(f"Saved: {path.name}", flush=True)

    try:
        if not args.headless:
            cv2.namedWindow(window, cv2.WINDOW_NORMAL)
            cv2.resizeWindow(window, 1000, 800)
        while True:
            now = time.monotonic()
            if args.seconds and now - started >= args.seconds:
                reason = "duration_reached"
                break
            packet = stream.snapshot()
            stats = stream.stats()
            if stats["failures"] != logged_failures:
                logged_failures = stats["failures"]
                print(f"Receive failure #{logged_failures}: {stats['last_error']}; reconnecting.", flush=True)
            if not stats["received_frames"] and now - started >= args.connect_wait:
                result, reason = 2, "no_frames_before_deadline"
                print("No decoded frames. Run check_connection.py; verify UDP firewall rule and transport.", flush=True)
                break
            fresh = packet is not None and now - packet[1] < 1.0
            if fresh and packet[0] != last_sequence:
                skipped += max(0, packet[0] - last_sequence - 1)
                consumed += 1
                last_sequence = packet[0]
                raw = packet[3]
                if first_shape is None:
                    first_shape = [raw.shape[1], raw.shape[0]]
                    print(f"First frame decoded: {first_shape[0]} x {first_shape[1]}", flush=True)
                if args.save_every and now >= next_save and automatic < args.max_images:
                    snapshot(packet)
                    automatic += 1
                    next_save = now + args.save_every
            if not args.headless:
                if fresh:
                    display = packet[3].copy()
                    label = f"Receive {stats['receive_fps_last_5s']:.1f} fps | #{packet[0]} | S save  Q quit"
                else:
                    display = np.zeros((600, 900, 3), dtype=np.uint8)
                    label = f"{stats['state']} | waiting for a fresh frame | Q quit"
                cv2.rectangle(display, (0, 0), (display.shape[1], 42), (0, 0, 0), -1)
                cv2.putText(display, label, (12, 28), cv2.FONT_HERSHEY_SIMPLEX, 0.65, (80, 255, 80), 1)
                cv2.imshow(window, display)
                key = cv2.waitKey(1) & 0xFF
                if key in (ord("q"), ord("Q"), 27) or cv2.getWindowProperty(window, cv2.WND_PROP_VISIBLE) < 1:
                    reason = "user_quit"
                    break
                if key in (ord("s"), ord("S")) and fresh:
                    snapshot(packet)
            if now - last_log >= 1:
                row = {name: stats[name] for name in fields if name in stats}
                row.update(utc=utc_now(), elapsed_s=now - started,
                           consumed_frames=consumed, skipped_by_consumer=skipped)
                writer.writerow(row)
                csv_file.flush()
                last_log = now
            if now - last_console >= 5:
                print(f"{now-started:.0f}s | {stats['state']} | receive {stats['receive_fps_last_5s']:.1f} fps"
                      f" | frames {stats['received_frames']} | failures {stats['failures']}", flush=True)
                last_console = now
            time.sleep(0.001 if fresh else 0.02)
    except KeyboardInterrupt:
        reason = "keyboard_interrupt"
    except Exception as exc:
        result, reason = 1, f"error: {exc}"
        print(reason, file=sys.stderr, flush=True)
    finally:
        ended = time.monotonic()
        final_stats = stream.stats()
        try:
            stream.stop()
        except RuntimeError as exc:
            result, reason = 1, str(exc)
        csv_file.close()
        if not args.headless:
            cv2.destroyAllWindows()
        if final_stats["received_frames"] == 0:
            result = 2
        write_json(folder / "summary.json", {**final_stats, "ended_utc": utc_now(),
                   "elapsed_s": ended - started, "image_size": first_shape,
                   "consumed_frames": consumed, "skipped_by_consumer": skipped,
                   "saved_images": saved, "exit_code": result, "stop_reason": reason,
                   "end_to_end_latency_ms": None,
                   "note": "Skipped consumer frames are deliberate latest-frame selection, not measured network loss."})
        print(f"Finished: {reason}; received {final_stats['received_frames']} frames; saved {saved} images.", flush=True)
    return result


if __name__ == "__main__":
    sys.exit(main())
