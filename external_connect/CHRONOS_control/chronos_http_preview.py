"""Fresh, sequential JPEG requests using the Chronos camera's own web endpoint.

This bypasses the camera's H.264/RTSP path. It is an occasional-image preview,
not a 60 fps stream and not raw sensor data. Never pipeline image requests.
"""
import http.client
import threading
import time


class HttpJpegCapture:
    MAX_IMAGE_BYTES = 4_000_000

    def __init__(self, config, stop_event=None, refresh_event=None):
        self.host = config['camera_ip']
        self.period = 1.0 / config.get('jpeg_max_fps', 4)
        self.timeout = config['read_timeout_ms'] / 1000
        self.stop_event = stop_event or threading.Event()
        self.refresh_event = refresh_event or threading.Event()
        self.next_request = 0
        self.latest_request_ms = None
        self.closed = False
        self.connection = None

    def isOpened(self):
        return not self.closed

    def read(self):
        import cv2
        import numpy as np
        while time.monotonic() < self.next_request:
            if self.stop_event.wait(min(.02, self.next_request-time.monotonic())):
                return False, None
            if self.refresh_event.is_set():
                return False, None
        if self.closed or self.stop_event.is_set() or self.refresh_event.is_set():
            return False, None
        started = time.monotonic()
        self.next_request = started + self.period
        # Direct LAN connection: do not route the camera through an HTTP proxy.
        connection = http.client.HTTPConnection(self.host, timeout=self.timeout)
        self.connection = connection
        try:
            connection.request('GET', f'/cgi-bin/screenCap?t={time.monotonic_ns()}',
                               headers={'Cache-Control':'no-cache, no-store', 'Connection':'close'})
            response = connection.getresponse()
            if response.status != 200:
                raise RuntimeError(f'Camera screenshot HTTP {response.status}')
            data = response.read(self.MAX_IMAGE_BYTES + 1)
            if len(data) > self.MAX_IMAGE_BYTES:
                raise RuntimeError('Camera screenshot exceeded 4 MB')
            if not data.startswith(b'\xff\xd8'):
                raise RuntimeError('Camera screenshot response is not JPEG')
            frame = cv2.imdecode(np.frombuffer(data, dtype=np.uint8), cv2.IMREAD_COLOR)
            if frame is None:
                raise RuntimeError('Cannot decode camera screenshot')
            self.latest_request_ms = (time.monotonic() - started) * 1000
            return True, frame
        finally:
            connection.close()
            self.connection = None

    def release(self):
        # Called only by the receiver thread, after read() has returned.
        self.closed = True
        if self.connection is not None:
            self.connection.close()
