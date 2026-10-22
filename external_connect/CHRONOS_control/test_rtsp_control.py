"""Protocol regression tests; the camera's real 5-minute test is separate."""
import socket
import threading
import unittest

from chronos_rtsp_control import ChronosRtspControl, split_message


def receive(sock):
    data=b''
    while True:
        data+=sock.recv(4096)
        result=split_message(data)
        if result:return result[0]


class RelayTests(unittest.TestCase):
    def test_fragmented_and_coalesced_messages_keep_body_boundaries(self):
        first=b'RTSP/1.0 200 OK\r\nContent-Length: 5\r\n\r\nhello'
        second=b'RTSP/1.0 200 OK\r\nCSeq: 2\r\n\r\n'
        for end in range(len(first)):
            self.assertIsNone(split_message(first[:end]))
        self.assertEqual(split_message(first+second),(first,second))
        with self.assertRaises(ValueError):
            split_message(b'RTSP/1.0 200 OK\r\nContent-Length: 9999999\r\n\r\n')

    def test_many_keepalives_get_new_camera_connections_and_preserve_session(self):
        listener=socket.socket()
        listener.bind(('127.0.0.1',0));listener.listen(3);listener.settimeout(3)
        received=[];errors=[]
        def camera():
            try:
                for i in range(24):
                    conn,_=listener.accept()
                    with conn:
                        conn.settimeout(3)
                        request=receive(conn);received.append(request)
                        response=f'RTSP/1.0 200 OK\r\nCSeq: {i+1}\r\nContent-Length: 3\r\n\r\nabc'.encode()
                        conn.sendall(response[:13]);conn.sendall(response[13:])
            except Exception as exc:errors.append(exc)
        worker=threading.Thread(target=camera,daemon=True);worker.start()
        proxy=ChronosRtspControl('127.0.0.1',listener.getsockname()[1]).start()
        try:
            with socket.create_connection(proxy.listener.getsockname(),timeout=3) as client:
                for i in range(24):
                    request=f'GET_PARAMETER {proxy.url} RTSP/1.0\r\nCSeq: {i+1}\r\nSession: 123456\r\n\r\n'.encode()
                    client.sendall(request)
                    response=receive(client)
                    self.assertIn(f'CSeq: {i+1}\r\n'.encode(),response)
                    self.assertTrue(response.endswith(b'abc'))
            worker.join(3)
            self.assertFalse(worker.is_alive())
            self.assertEqual(errors,[])
            self.assertEqual(len(received),24)
            self.assertEqual(proxy.requests,24)
            for request in received:
                self.assertIn(b'Session: 123456\r\n',request)
                self.assertIn(f'rtsp://127.0.0.1:{listener.getsockname()[1]}/'.encode(),request)
        finally:
            proxy.close();listener.close()
        self.assertFalse(proxy.thread.is_alive())


if __name__=='__main__':unittest.main()
