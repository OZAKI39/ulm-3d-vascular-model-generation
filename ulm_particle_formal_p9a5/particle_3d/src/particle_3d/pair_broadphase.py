"""Deterministic AABB sweep; shape support plus swept motion, no neighbor skin."""
from itertools import combinations
import numpy as np
from .particle_shapes import roundoff_length


def all_pairs(shapes):
    return list(combinations(sorted(shapes), 2))


def broadphase(shapes, *, end_shapes=None, surface_motion_bounds=None):
    records = []
    for particle_id, shape in shapes.items():
        axes = np.eye(3)
        low = np.array([shape.support(-a)[k] for k, a in enumerate(axes)])
        high = np.array([shape.support(a)[k] for k, a in enumerate(axes)])
        pad = roundoff_length(shape.center_m, shape.bounding_radius_m)
        motion = 0. if surface_motion_bounds is None else surface_motion_bounds[particle_id]
        if end_shapes is not None:
            end = end_shapes[particle_id]
            low = np.minimum(low, [end.support(-a)[k] for k, a in enumerate(axes)])
            high = np.maximum(high, [end.support(a)[k] for k, a in enumerate(axes)])
            pad = max(pad, roundoff_length(end.center_m, end.bounding_radius_m))
        records.append((low[0]-pad-motion, particle_id, low-pad-motion, high+pad+motion))
    records.sort(key=lambda r: (r[0], r[1])); active = []; pairs = []
    for x, particle_id, low, high in records:
        active = [r for r in active if r[2][0] >= x]
        for other_id, other_low, other_high in active:
            if np.all(other_high >= low) and np.all(high >= other_low):
                pairs.append(tuple(sorted([particle_id, other_id])))
        active.append((particle_id, low, high))
    return sorted(pairs)
