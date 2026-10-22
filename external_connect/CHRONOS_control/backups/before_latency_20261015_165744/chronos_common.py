"""Shared configuration and latest-frame RTSP receiver for Chronos 1.4."""
from __future__ import annotations

import ipaddress
import json
import os
from collections import deque
from datetime import datetime, timezone
from pathlib import Path
import threading
import time

BASE = Path(__file__).resolve().parent
DEFAULT_CONFIG = BASE / "chronos_config.json"


def utc_now():
    return datetime.now(timezone.utc).isoformat(timespec="milliseconds")


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")
    temporary.replace(path)


def load_config(path=DEFAULT_CONFIG, ip=None):
    config = json.loads(Path(path).read_text(encoding="utf-8-sig"))
    config["camera_ip"] = str(ipaddress.IPv4Address(ip or config["camera_ip"]))
    if config.get("rtsp_transport") not in ("tcp", "udp"):
        raise ValueError("rtsp_transport must be tcp or udp")
    for key in ("open_timeout_ms", "read_timeout_ms"):
        if not isinstance(config[key], int) or config[key] <= 0:
            raise ValueError(f"{key} must be a positive integer")
    if not 0 < float(config["reconnect_seconds"]) <= 60:
        raise ValueError("reconnect_seconds must be between 0 and 60")
    output = Path(config.get("output_dir", "data"))
    config["output_dir"] = str(output if output.is_absolute() else Path(path).resolve().parent / output)
    return config


def session_dir(config, prefix):
    name = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    path = Path(config["output_dir"]) / f"{prefix}_{name}"
    path.mkdir(parents=True, exist_ok=False)
    return path


def video_options(config):
    options = f"rtsp_transport;{config['rtsp_transport']}"
    if config['rtsp_transport'] == 'udp' and config.get('udp_port'):
        port = int(config['udp_port'])
        if port < 1024 or port > 65534 or port % 2:
            raise ValueError('udp_port must be an even port from 1024 to 65534')
        options += f"|min_port;{port}|max_port;{port+1}"
    return options


def configure_video_environment(config):
    os.environ['OPENCV_FFMPEG_CAPTURE_OPTIONS'] = video_options(config)


def load_image(path):
    import cv2
    import numpy as np
    image = cv2.imdecode(np.fromfile(str(path), dtype=np.uint8), cv2.IMREAD_COLOR)
    if image is None:
        raise ValueError(f"Cannot read image: {path}")
    return image


def save_image(path, frame):
    import cv2
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    ok, encoded = cv2.imencode(path.suffix, frame)
    if not ok:
        raise OSError(f"Image encoding failed: {path}")
    encoded.tofile(str(path))


class LatestFrameStream:
    """One reader drains FFmpeg; consumers see at most one latest decoded frame.

    receive_monotonic marks decode completion on this PC, not camera exposure.
    This bounds the Python queue only; it cannot remove camera/decoder buffers.
    """

    def __init__(self, config):
        self.config = config
        self.url = f"rtsp://{config['camera_ip']}"
        self.lock = threading.Lock()
        self.stop_event = threading.Event()
        self.thread = None
        self.latest = None
        self.frames = 0
        self.failures = 0
        self.connections = 0
        self.attempts = 0
        self.state = "not_started"
        self.last_error = ""
        self.first_received = None
        self.last_received = None
        self.times = deque(maxlen=1200)
        self.reported_fps = None
        self.failure_events = deque(maxlen=1000)

    def start(self):
        if self.thread is not None:
            raise RuntimeError("A stream can be started only once")
        self.thread = threading.Thread(target=self._run, name="chronos-receiver", daemon=True)
        self.thread.start()
        return self

    def _run(self):
        # The attached camera uses UDP; set options before loading cv2 on Windows.
        configure_video_environment(self.config)
        import cv2
        params = [cv2.CAP_PROP_OPEN_TIMEOUT_MSEC, self.config["open_timeout_ms"],
                  cv2.CAP_PROP_READ_TIMEOUT_MSEC, self.config["read_timeout_ms"]]
        while not self.stop_event.is_set():
            cap = None
            try:
                with self.lock:
                    self.state = "connecting"
                    self.attempts += 1
                    self.latest = None
                cap = cv2.VideoCapture(self.url, cv2.CAP_FFMPEG, params)
                if not cap.isOpened():
                    raise RuntimeError("Cannot open RTSP; check IP, Web Control and transport")
                with self.lock:
                    self.connections += 1
                    reported = cap.get(cv2.CAP_PROP_FPS)
                    self.reported_fps = reported if 0 < reported < 100000 else None
                    self.state = "waiting_for_frame"
                while not self.stop_event.is_set():
                    ok, frame = cap.read()
                    if not ok or frame is None:
                        raise RuntimeError("RTSP read failed or timed out")
                    now = time.monotonic()
                    stamp = utc_now()
                    with self.lock:
                        self.frames += 1
                        self.first_received = now if self.first_received is None else self.first_received
                        self.last_received = now
                        self.times.append(now)
                        while self.times and now - self.times[0] > 5:
                            self.times.popleft()
                        self.latest = (self.frames, now, stamp, frame)
                        self.state = "streaming"
                        self.last_error = ""
            except Exception as exc:
                with self.lock:
                    self.failures += 1
                    self.state = "reconnecting"
                    self.last_error = str(exc)
                    self.failure_events.append({"utc": utc_now(), "received_frames": self.frames,
                                                "error": str(exc)})
                    self.latest = None  # Do not offer an old frame as a new image.
            finally:
                if cap is not None:
                    cap.release()  # Only the reader thread touches the capture object.
            self.stop_event.wait(self.config["reconnect_seconds"])
        with self.lock:
            self.state = "stopped"

    def snapshot(self):
        with self.lock:
            return self.latest

    def stats(self):
        now = time.monotonic()
        with self.lock:
            recent = [t for t in self.times if now - t <= 5]
            fps = (len(recent) - 1) / (recent[-1] - recent[0]) if len(recent) > 1 else 0.0
            span = (self.last_received - self.first_received) if self.first_received is not None else 0
            return {"state": self.state, "received_frames": self.frames,
                    "receive_fps_last_5s": fps,
                    "receive_fps_first_to_last": (self.frames - 1) / span if span > 0 else 0.0,
                    "seconds_since_last_frame": now - self.last_received if self.last_received is not None else None,
                    "connections": self.connections, "connection_attempts": self.attempts,
                    "failures": self.failures, "reported_stream_fps": self.reported_fps,
                    "last_error": self.last_error, "recent_failure_events": list(self.failure_events)}

    def stop(self):
        self.stop_event.set()
        if self.thread is not None:
            self.thread.join((self.config["open_timeout_ms"] + self.config["read_timeout_ms"]) / 1000 + 2)
            if self.thread.is_alive():
                raise RuntimeError("Video backend did not stop within its configured timeout")
