"""Explicit microscope profiles for the existing dark-component detector.

No channel geometry or physical calibration is inferred by this module.
"""
import json
import time
from pathlib import Path

import cv2
import numpy as np

LEGACY = dict(sigma_px=4., threshold=10., min_area_px=8, max_area_px=450,
              max_aspect=2.8, border_px=2)


def validate_config(config):
    result = {**LEGACY, **config}
    if set(result) != set(LEGACY):
        raise ValueError('Unknown detector parameter')
    if not all(np.isfinite(v) for v in result.values()):
        raise ValueError('Nonfinite detector parameter')
    if not (0 < result['sigma_px'] <= 32 and 0 < result['threshold'] <= 255
            and 1 <= result['min_area_px'] <= result['max_area_px']
            and 1 <= result['max_aspect'] <= 20 and result['border_px'] >= 0):
        raise ValueError('Invalid detector bounds')
    return result


def load_profile(path):
    p = json.loads(Path(path).read_text())
    if p.get('schema') != 'MICROSCOPE_OPENCV_PROFILE_V1' or not p.get('admitted'):
        raise ValueError('Detector profile has not passed its admission checks')
    wh=p.get('camera_image_wh',[])
    if len(wh)!=2 or any(not isinstance(v,int) or v<48 for v in wh):
        raise ValueError('Detector profile requires the original camera resolution')
    p['parameters'] = validate_config(p['parameters'])
    return p


def detect(gray, config=None):
    c = validate_config(config or {})
    if gray.ndim != 2 or gray.dtype != np.uint8:
        raise ValueError('Expected original uint8 grayscale pixels')
    h, w = gray.shape
    signal = cv2.subtract(cv2.GaussianBlur(gray, (0, 0), c['sigma_px']), gray)
    _, labels, stats, centers = cv2.connectedComponentsWithStats(
        np.uint8(signal >= c['threshold']), 8)
    s = stats[1:]; q = centers[1:]; border = c['border_px']
    keep = ((s[:, 4] >= c['min_area_px']) & (s[:, 4] <= c['max_area_px'])
            & (s[:, 2] / np.maximum(s[:, 3], 1) <= c['max_aspect'])
            & (s[:, 3] / np.maximum(s[:, 2], 1) <= c['max_aspect'])
            & (s[:, 0] > border) & (s[:, 1] > border)
            & (s[:, 0] + s[:, 2] < w-border) & (s[:, 1] + s[:, 3] < h-border))
    ids = np.flatnonzero(keep) + 1
    objects = [dict(label=int(n), xy_local=q[n-1].tolist(), area=int(s[n-1, 4])) for n in ids]
    return signal, labels, objects


def profiled_channel(config, profile):
    """Keep the existing PIV geometry/registration and replace candidate scanning."""
    from channel import Channel

    class ProfiledChannel(Channel):
        def scan(self, gray):
            start=time.perf_counter()
            result = super().scan(gray)
            signal, _, objects = detect(self.cut(gray), profile['parameters'])
            q = np.asarray([o['xy_local'] for o in objects], float).reshape(-1, 2)
            areas = np.asarray([o['area'] for o in objects], int)
            inside = self.local[np.rint(q[:, 1]).astype(int), np.rint(q[:, 0]).astype(int)].astype(bool)
            outside = q[~inside]; q, areas = q[inside], areas[inside]
            pointmap = np.zeros(self.local.shape, np.float32)
            if len(q):
                np.add.at(pointmap, (np.rint(q[:, 1]).astype(int), np.rint(q[:, 0]).astype(int)), 1)
            counts = self.window_mean(pointmap) * self.ws**2
            # Preserve existing texture/motion activation and PIV preprocessing.
            selected = self.safe & ((counts >= 1) | ((result['texture'] >= 2) & (result['activity'] >= .4)))
            result.update(candidates=q, areas=areas, excluded_candidates=outside,
                          count=counts, selected=selected, detector_profile_id=profile['id'],
                          scan_ms=(time.perf_counter()-start)*1000)
            return result

    return ProfiledChannel(config)
