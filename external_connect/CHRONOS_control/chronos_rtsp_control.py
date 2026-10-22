"""Loopback-only control relay for the old Chronos RTSP request-buffer bug.

Each control request gets a fresh camera TCP connection. The camera keeps its
RTSP session across connections. RTP video continues directly camera -> PC;
this relay never receives, decodes or queues video packets.
"""
import socket
import threading
import time
from urllib.parse import urlsplit

MAX_MESSAGE = 65536
ALLOWED = {'OPTIONS', 'DESCRIBE', 'SETUP', 'PLAY', 'PAUSE', 'TEARDOWN', 'GET_PARAMETER'}


def split_message(buffer):
    end = buffer.find(b'\r\n\r\n')
    if end < 0:
        if len(buffer) > MAX_MESSAGE:
            raise ValueError('RTSP headers too large')
        return None
    headers = buffer[:end].split(b'\r\n')
    length = 0
    for line in headers[1:]:
        if line.lower().startswith(b'content-length:'):
            length = int(line.split(b':', 1)[1].strip())
    total = end + 4 + length
    if length < 0 or total > MAX_MESSAGE:
        raise ValueError('RTSP message too large')
    if len(buffer) < total:
        return None
    return buffer[:total], buffer[total:]


class ChronosRtspControl:
    def __init__(self, camera_ip, camera_port=554):
        self.camera_ip = camera_ip
        self.camera_port = camera_port
        self.stop_event = threading.Event()
        self.listener = None
        self.thread = None
        self.client = None
        self.requests = 0
        self.errors = []

    def start(self):
        self.listener = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.listener.bind(('127.0.0.1', 0))
        self.listener.listen(1)
        self.listener.settimeout(.2)
        self.url = f'rtsp://127.0.0.1:{self.listener.getsockname()[1]}/'
        self.thread = threading.Thread(target=self._run, name='chronos-rtsp-control', daemon=True)
        self.thread.start()
        return self

    def _camera_request(self, request):
        first, rest = request.split(b'\r\n', 1)
        method, uri, protocol = first.decode('ascii').split(' ', 2)
        if method not in ALLOWED or protocol != 'RTSP/1.0':
            raise ValueError('Unsupported RTSP control method')
        parsed = urlsplit(uri)
        target = f'rtsp://{self.camera_ip}:{self.camera_port}{parsed.path or "/"}'
        if parsed.query:
            target += '?' + parsed.query
        forwarded = f'{method} {target} {protocol}\r\n'.encode() + rest
        with socket.create_connection((self.camera_ip, self.camera_port), timeout=2) as upstream:
            upstream.settimeout(2)
            upstream.sendall(forwarded)
            received = b''
            deadline = time.monotonic()+3
            while time.monotonic() < deadline:
                chunk = upstream.recv(8192)
                if not chunk:
                    raise ConnectionError('Camera closed RTSP before completing its response')
                received += chunk
                parts = split_message(received)
                if parts:
                    self.requests += 1
                    return parts[0]
            raise TimeoutError('Camera RTSP response deadline exceeded')

    def _run(self):
        while not self.stop_event.is_set():
            try:
                client, _ = self.listener.accept()
            except socket.timeout:
                continue
            except OSError:
                break
            self.client = client
            client.settimeout(.2)
            buffer = b''
            try:
                with client:
                    while not self.stop_event.is_set():
                        try:
                            chunk = client.recv(8192)
                        except socket.timeout:
                            continue
                        if not chunk:
                            break
                        buffer += chunk
                        while True:
                            parts = split_message(buffer)
                            if not parts:
                                break
                            request, buffer = parts
                            client.sendall(self._camera_request(request))
            except (OSError, ValueError) as exc:
                if not self.stop_event.is_set():
                    self.errors.append(str(exc))
                    self.errors = self.errors[-20:]
            finally:
                self.client = None

    def close(self):
        self.stop_event.set()
        if self.listener:
            self.listener.close()
        if self.client:
            try:
                self.client.shutdown(socket.SHUT_RDWR)
            except OSError:
                pass
        if self.thread:
            self.thread.join(5)
            if self.thread.is_alive():
                raise RuntimeError('RTSP control relay did not stop')
