"""Offline math/segmentation tests. Synthetic images are NOT camera validation."""
import math
from pathlib import Path
import tempfile
import threading
import time
import unittest
from unittest.mock import patch

import cv2
import numpy as np

from analyze_stability import measure_frame
from calibrate_scale import scale_from_points
from chronos_common import LatestFrameStream, SnapshotWriter, capture_parameters, load_image, save_image, video_options
from measure_latency import latency_summary


class ToolTests(unittest.TestCase):
    def test_calibration_geometry(self):
        pixels, scale = scale_from_points([10, 20, 40, 60], 100)
        self.assertEqual(pixels, 50)
        self.assertEqual(scale, 2)

    def test_calibration_rejects_invalid_inputs(self):
        for points, distance in [([1, 1, 1, 1], 100), ([0, 0, 10, 0], 0), ([0, 0, 10, 0], math.nan)]:
            with self.assertRaises(ValueError):
                scale_from_points(points, distance)

    def test_dark_circle_and_roi_offset(self):
        frame = np.full((200, 240, 3), 255, np.uint8)
        cv2.circle(frame, (120, 90), 25, (0, 0, 0), -1)
        result, _ = measure_frame(frame, (70, 40, 100, 100), threshold=127)
        self.assertTrue(result["detected"])
        self.assertAlmostEqual(result["centroid_x_px"], 120)
        self.assertAlmostEqual(result["centroid_y_px"], 90)
        self.assertLess(abs(result["area_px2"] - math.pi*25**2), 120)
        self.assertAlmostEqual(result["aspect_ratio"], 1, places=5)

    def test_bright_circle(self):
        frame = np.zeros((100, 100, 3), np.uint8)
        cv2.circle(frame, (50, 50), 12, (255, 255, 255), -1)
        result, _ = measure_frame(frame, polarity="bright", threshold=127)
        self.assertTrue(result["detected"])

    def test_empty_and_clipped_target_rejected(self):
        frame = np.full((100, 100, 3), 255, np.uint8)
        self.assertFalse(measure_frame(frame, threshold=127)[0]["detected"])
        cv2.circle(frame, (0, 50), 20, (0, 0, 0), -1)
        self.assertFalse(measure_frame(frame, threshold=127)[0]["detected"])

    def test_roi_bounds(self):
        with self.assertRaises(ValueError):
            measure_frame(np.zeros((20, 20, 3), np.uint8), (10, 10, 20, 20))

    def test_unicode_image_roundtrip(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "标定图.png"
            frame = np.full((30, 40, 3), 123, np.uint8)
            save_image(path, frame)
            np.testing.assert_array_equal(frame, load_image(path))

    def test_latency_math_and_rejection(self):
        result = latency_summary([(1000, 700), (2000, 1600)])
        self.assertEqual(result["mean_ms"], 350)
        self.assertEqual(result["n"], 2)
        with self.assertRaises(ValueError):
            latency_summary([(100, 101)])

    def test_camera_udp_options(self):
        self.assertEqual(video_options({"rtsp_transport": "udp", "udp_port": 5000}),
                         "rtsp_transport;udp|min_port;5000|max_port;5001")

    def test_low_latency_preserves_camera_headers_and_sets_one_thread(self):
        config = {'rtsp_transport':'udp', 'udp_port':5000, 'latency_mode':'low',
                  'open_timeout_ms':5000, 'read_timeout_ms':1500}
        self.assertIn('max_delay;20000', video_options(config))
        self.assertIn('reorder_queue_size;32', video_options(config))
        self.assertNotIn('nobuffer', video_options(config))
        params = capture_parameters(config, cv2)
        self.assertEqual(dict(zip(params[::2], params[1::2]))[cv2.CAP_PROP_N_THREADS], 1)

    def test_background_snapshots_are_bounded_and_preserve_raw_pixels(self):
        entered, release = threading.Event(), threading.Event()
        images = []
        def slow_write(path, frame):
            entered.set()
            release.wait(2)
            images.append(frame.copy())
        raw = np.full((20, 20, 3), 55, np.uint8)
        packet = (1, time.monotonic(), '2026-10-15T00:00:00+00:00', raw)
        with tempfile.TemporaryDirectory() as folder, patch('chronos_common.save_image', side_effect=slow_write):
            writer = SnapshotWriter(folder, '192.168.0.10', queue_size=1)
            try:
                self.assertTrue(writer.submit(packet))
                self.assertTrue(entered.wait(1))
                self.assertTrue(writer.submit(packet))
                self.assertFalse(writer.submit(packet))
                self.assertEqual(writer.stats()['dropped_queue_full'], 1)
            finally:
                release.set()
                writer.close()
            self.assertEqual(writer.stats()['saved'], 2)
            self.assertEqual(writer.stats()['errors'], [])
            self.assertFalse(writer.submit(packet))
            for image in images:
                np.testing.assert_array_equal(image, raw)

    def test_manual_refresh_reconnects_on_reader_thread(self):
        captures = []
        class FakeCapture:
            def __init__(self, *args):
                self.released = False
                captures.append(self)
            def isOpened(self):
                return True
            def get(self, key):
                return 60
            def read(self):
                time.sleep(0.005)
                return True, np.zeros((10, 10, 3), np.uint8)
            def release(self):
                self.released = True
        config = {'camera_ip':'192.168.0.10', 'rtsp_transport':'udp', 'udp_port':5000,
                  'latency_mode':'low', 'open_timeout_ms':1000, 'read_timeout_ms':1000, 'reconnect_seconds':5}
        def wait_for(predicate):
            end = time.monotonic() + 2
            while not predicate() and time.monotonic() < end:
                time.sleep(0.01)
            self.assertTrue(predicate())
        with patch('cv2.VideoCapture', side_effect=FakeCapture):
            stream = LatestFrameStream(config).start()
            try:
                wait_for(lambda: stream.snapshot() is not None)
                stream.refresh()
                wait_for(lambda: stream.stats()['manual_refreshes'] == 1 and stream.stats()['connections'] == 2)
                self.assertEqual(stream.stats()['failures'], 0)
            finally:
                stream.stop()
        self.assertTrue(all(cap.released for cap in captures))

    def test_disconnect_removes_stale_frame_and_releases_capture(self):
        ready, fail = threading.Event(), threading.Event()
        class FakeCapture:
            def __init__(self):
                self.reads = 0
                self.released = False
            def isOpened(self):
                return True
            def get(self, key):
                return 60.0
            def read(self):
                self.reads += 1
                if self.reads == 1:
                    return True, np.zeros((20, 20, 3), np.uint8)
                ready.set()
                fail.wait(2)
                return False, None
            def release(self):
                self.released = True
        capture = FakeCapture()
        config = {"camera_ip": "192.168.0.10", "rtsp_transport": "udp", "udp_port": 5000,
                  "open_timeout_ms": 1000, "read_timeout_ms": 1000, "reconnect_seconds": 5}
        with patch("cv2.VideoCapture", return_value=capture):
            stream = LatestFrameStream(config).start()
            try:
                self.assertTrue(ready.wait(2))
                self.assertIsNotNone(stream.snapshot())
                fail.set()
                deadline = time.monotonic() + 2
                while stream.stats()["failures"] == 0 and time.monotonic() < deadline:
                    time.sleep(0.01)
                self.assertEqual(stream.stats()["failures"], 1)
                self.assertIsNone(stream.snapshot())
            finally:
                fail.set()
                stream.stop()
            self.assertTrue(capture.released)


if __name__ == "__main__":
    unittest.main()
