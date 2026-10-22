"""Chronos live preview, raw PNG snapshots, receive-rate logging and reconnect."""
import argparse
import csv
import os
import subprocess
from pathlib import Path
import sys
import time

from chronos_common import (DEFAULT_CONFIG, LatestFrameStream, SnapshotWriter, configure_video_environment, load_config, video_options,
                            save_image, session_dir, utc_now, write_json)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=str(DEFAULT_CONFIG))
    parser.add_argument("--ip")
    parser.add_argument("--transport", choices=["tcp", "udp"])
    parser.add_argument("--latency-mode", choices=["fresh", "low", "compatibility"], help="Fresh HTTP snapshots, tuned RTSP, or original RTSP buffering")
    parser.add_argument("--jpeg-fps", type=float, help="Maximum sequential screenshot requests per second, 1 to 10")
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
    if args.latency_mode:
        config["latency_mode"] = args.latency_mode
    if args.jpeg_fps is not None:
        if not 1 <= args.jpeg_fps <= 10:
            parser.error('--jpeg-fps must be between 1 and 10')
        config['jpeg_max_fps'] = args.jpeg_fps
    # The FFmpeg DLL on Windows reads its environment at process startup.
    # Passing env to a child also works when launching this .py directly.
    desired_options = video_options(config)
    if config['latency_mode'] != 'fresh' and os.environ.get("OPENCV_FFMPEG_CAPTURE_OPTIONS") != desired_options:
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
    source_url = f"http://{config['camera_ip']}/cgi-bin/screenCap" if config['latency_mode'] == 'fresh' else f"rtsp://{config['camera_ip']}"
    print(f"Chronos: {source_url} | {config['latency_mode']}", flush=True)
    print(f"Output: {folder}\nS: save raw PNG | R: reconnect to refresh | Q/Esc: quit | Ctrl+C: stop", flush=True)
    stream = LatestFrameStream(config).start()
    snapshots = SnapshotWriter(folder, config['camera_ip'])
    started = time.monotonic()
    last_log = last_console = started
    next_save = started
    consumed = skipped = automatic = 0
    logged_failures = 0
    last_sequence = 0
    first_shape = None
    result = 0
    reason = "normal"
    window = "Chronos 1.4 | S save | R refresh | Q quit"
    display_sequence = -1
    display_count = 0
    next_display = started
    display_period = 1.0 / config['preview_max_fps']
    stale_seconds = config['stale_frame_ms'] / 1000
    if config['latency_mode'] == 'fresh':
        stale_seconds = max(stale_seconds, 1.0, 2.0/config['jpeg_max_fps'])
    from collections import deque
    display_ages_ms = deque(maxlen=6000)
    csv_file = (folder / "stream_metrics.csv").open("w", newline="", encoding="utf-8-sig")
    fields = ["utc", "elapsed_s", "received_frames", "receive_fps_last_5s", "consumed_frames",
              "skipped_by_consumer", "seconds_since_last_frame", "connections", "failures", "state",
              "decoder_threads", "displayed_frames", "latest_decode_to_display_ms",
              "preview_source", "latest_http_request_ms"]
    writer = csv.DictWriter(csv_file, fieldnames=fields)
    writer.writeheader()

    def snapshot(packet):
        accepted = snapshots.submit(packet)
        print(f"Snapshot {'queued' if accepted else 'skipped: writer busy'}: frame {packet[0]}", flush=True)
        return accepted

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
                print("No decoded frames. Run check_connection.py; verify the camera webpage, selected preview mode and UDP rule for RTSP.", flush=True)
                break
            fresh = packet is not None and now - packet[1] < stale_seconds
            if fresh and packet[0] != last_sequence:
                skipped += max(0, packet[0] - last_sequence - 1)
                consumed += 1
                last_sequence = packet[0]
                raw = packet[3]
                if first_shape is None:
                    first_shape = [raw.shape[1], raw.shape[0]]
                    print(f"First frame decoded: {first_shape[0]} x {first_shape[1]} | decoder threads={stats['decoder_threads']}", flush=True)
                if args.save_every and now >= next_save and automatic < args.max_images:
                    if snapshot(packet):
                        automatic += 1
                    next_save = now + args.save_every
            if not args.headless:
                needs_display = (fresh and packet[0] != display_sequence) or not fresh
                if needs_display and now >= next_display:
                    if fresh:
                        raw = packet[3]
                        width = min(raw.shape[1], int(config['preview_width']))
                        display = cv2.resize(raw, (width, round(raw.shape[0]*width/raw.shape[1])),
                                             interpolation=cv2.INTER_LINEAR) if width < raw.shape[1] else raw.copy()
                        mode_label = 'FRESH JPEG' if config['latency_mode'] == 'fresh' else config['latency_mode'].upper()
                        label = f"{mode_label} | {stats['receive_fps_last_5s']:.1f} fps | S save  R refresh  Q quit"
                    else:
                        display = np.zeros((600, 900, 3), dtype=np.uint8)
                        label = f"{stats['state']} | waiting for a fresh frame | Q quit"
                    cv2.rectangle(display, (0, 0), (display.shape[1], 42), (0, 0, 0), -1)
                    cv2.putText(display, label, (12, 28), cv2.FONT_HERSHEY_SIMPLEX, 0.60, (80, 255, 80), 1)
                    if fresh:
                        display_ages_ms.append((time.monotonic()-packet[1])*1000)
                        display_sequence = packet[0]
                        display_count += 1
                    cv2.imshow(window, display)
                    next_display = max(next_display + display_period, now)
                # pollKey drains the Win32 message queue without waitKey's extra wait.
                key = cv2.pollKey() & 0xFF
                try:
                    visible = cv2.getWindowProperty(window, cv2.WND_PROP_VISIBLE) >= 1
                except cv2.error:
                    visible = False
                if key in (ord("q"), ord("Q"), 27) or not visible:
                    reason = "user_quit"
                    break
                if key in (ord("s"), ord("S")) and fresh:
                    snapshot(packet)
                if key in (ord("r"), ord("R")):
                    stream.refresh()
                    display_sequence = -1
                    print('Manual refresh requested.', flush=True)
            if now - last_log >= 1:
                row = {name: stats[name] for name in fields if name in stats}
                row.update(utc=utc_now(), elapsed_s=now - started,
                           consumed_frames=consumed, skipped_by_consumer=skipped,
                           displayed_frames=display_count,
                           latest_decode_to_display_ms=display_ages_ms[-1] if display_ages_ms else None)
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
        try:
            snapshots.close()
        except RuntimeError as exc:
            result, reason = 1, str(exc)
        saved_stats = snapshots.stats()
        if saved_stats['errors']:
            result, reason = 1, 'snapshot_write_failed'
        csv_file.close()
        if not args.headless:
            cv2.destroyAllWindows()
        if final_stats["received_frames"] == 0:
            result = 2
        write_json(folder / "summary.json", {**final_stats, "ended_utc": utc_now(),
                   "elapsed_s": ended - started, "image_size": first_shape,
                   "consumed_frames": consumed, "skipped_by_consumer": skipped,
                   "saved_images": saved_stats['saved'], "snapshot_writer": saved_stats,
                   "displayed_frames": display_count,
                   "decode_to_display_ms": {'samples':len(display_ages_ms),
                       'median':float(np.median(display_ages_ms)) if display_ages_ms else None,
                       'p95':float(np.percentile(display_ages_ms,95)) if display_ages_ms else None,
                       'meaning':'PC decode completion to imshow submission, last 6000 displayed frames; excludes camera, network/decoder buffers and screen scanout.'},
                   "exit_code": result, "stop_reason": reason,
                   "end_to_end_latency_ms": None,
                   "note": "Skipped consumer frames are deliberate latest-frame selection, not measured network loss."})
        print(f"Finished: {reason}; received {final_stats['received_frames']} frames; saved {saved_stats['saved']} images.", flush=True)
    return result


if __name__ == "__main__":
    sys.exit(main())
