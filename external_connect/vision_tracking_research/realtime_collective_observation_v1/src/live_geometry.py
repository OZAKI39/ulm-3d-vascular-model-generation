"""Validate an interactive LIVE mask; no implicit camera-motion calibration."""
from copy import deepcopy

import cv2
import numpy as np


def validate_geometry(value, image_wh):
    g = deepcopy(value)
    w, h = map(int, image_wh)
    if g.get('source_mode') != 'LIVE_CAMERA' or not g.get('human_verified'):
        raise ValueError('Please use a mask drawn on a LIVE image.')
    if list(g.get('image_wh', [])) != [w, h]:
        raise ValueError('Mask resolution differs from the current camera image. Draw a new mask.')
    if min(w, h) < 48:
        raise ValueError('The camera image must be at least 48 pixels in each direction.')
    points = np.asarray(g.get('channel_polygon', []), dtype=float)
    if points.ndim != 2 or points.shape[1] != 2 or len(points) < 3:
        raise ValueError('Click at least three different points around the channel.')
    if np.linalg.norm(points[-1] - points[0]) < 1e-6:
        points = points[:-1]
    if len(points) < 3 or not np.isfinite(points).all():
        raise ValueError('The polygon needs at least three finite, different points.')
    if np.any(points < 0) or np.any(points >= [w, h]):
        raise ValueError('All mask points must be inside the current image.')
    # Reject crossed polygons instead of silently turning them into another ROI.
    def orient(a, b, c):
        return float((b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0]))
    def on_segment(a, b, c):
        return abs(orient(a, b, c)) < 1e-7 and np.all(c >= np.minimum(a, b)-1e-7) and np.all(c <= np.maximum(a, b)+1e-7)
    for i, a in enumerate(points):
        b = points[(i+1) % len(points)]
        if np.linalg.norm(b-a) < 1e-6:
            raise ValueError('Remove duplicate neighbouring points.')
        for j in range(i+1, len(points)):
            if j == i+1 or (i == 0 and j == len(points)-1):
                continue
            c, d = points[j], points[(j+1) % len(points)]
            crosses = orient(a,b,c)*orient(a,b,d) < 0 and orient(c,d,a)*orient(c,d,b) < 0
            if crosses or any((on_segment(a,b,c), on_segment(a,b,d), on_segment(c,d,a), on_segment(c,d,b))):
                raise ValueError('The mask crosses itself. Click points in order around the boundary.')
    if abs(cv2.contourArea(points.astype(np.float32))) < 16:
        raise ValueError('The selected area is too small. Draw around the channel interior.')
    crop = np.asarray(g.get('crop_xywh', [0,0,w,h]), dtype=float)
    if crop.shape != (4,) or not np.isfinite(crop).all() or np.any(crop != np.round(crop)):
        raise ValueError('Invalid image crop.')
    x, y, cw, ch = crop.astype(int)
    if min(x,y) < 0 or min(cw,ch) < 48 or x+cw > w or y+ch > h:
        raise ValueError('The image crop must be inside the camera image and at least 48 pixels wide/high.')
    roi = g.get('registration_roi')
    if roi is not None:
        box = np.asarray(roi, dtype=float)
        if box.shape != (4,) or not np.isfinite(box).all() or np.any(box != np.round(box)):
            raise ValueError('Invalid fixed-wall reference rectangle.')
        rx, ry, rw, rh = box.astype(int)
        if min(rx,ry) < 0 or min(rw,rh) < 8 or rx+rw > w or ry+rh > h:
            raise ValueError('The fixed-wall reference must lie inside the camera image.')
        g['registration_roi'] = box.astype(int).tolist()
    g.update(channel_polygon=points.tolist(), crop_xywh=crop.astype(int).tolist())
    g.setdefault('excluded_polygons', [])
    g.setdefault('minimum_valid_window_fraction', 1.0)
    g.setdefault('registration_min_response', .08)
    g.setdefault('registration_max_translation_px', 6)
    # This flag describes the coordinate frame, not a passed drift check.
    g['motion_reference'] = 'FIXED_WALL_CHECK' if roi is not None else 'IMAGE_SPACE_UNREGISTERED'
    return g


def new_geometry(points, image_wh):
    w, h = image_wh
    return validate_geometry(dict(image_wh=[w,h], crop_xywh=[0,0,w,h],
        channel_polygon=points, excluded_polygons=[], human_verified=True,
        source_mode='LIVE_CAMERA', registration_roi=None), image_wh)
