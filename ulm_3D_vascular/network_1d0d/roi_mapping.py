"""Exact interpolation and virtual-node insertion; no nearest-neighbor fallback."""
import numpy as np


def edge_position(x0, x1, r0, r1, t):
    if not np.isfinite(t) or t < 0 or t > 1:
        raise ValueError('Cut fraction outside edge')
    return np.asarray(x0) + t*(np.asarray(x1)-x0), r0+t*(r1-r0)


def match_edge_position(x0, x1, r0, r1, point, radius, tolerance_um=1e-8):
    d = np.asarray(x1)-x0
    if np.dot(d, d) == 0:
        raise ValueError('Zero-length edge')
    t = float(np.dot(np.asarray(point)-x0, d)/np.dot(d, d))
    if t < -1e-12 or t > 1+1e-12:
        raise ValueError('Point outside the saved original edge')
    xyz, r = edge_position(x0, x1, r0, r1, np.clip(t, 0, 1))
    err, rerr = float(np.linalg.norm(xyz-point)), float(abs(r-radius))
    if err > tolerance_um or rerr > tolerance_um:
        raise ValueError('Saved edge identity does not match geometry/radius')
    return dict(fraction_parent_to_child=float(np.clip(t, 0, 1)), position_error_um=err,
                radius_error_um=rerr, candidate_count=1, method='SAVED_EDGE_ID_EXACT_INTERPOLATION',
                tangent_parent_to_child=(d/np.linalg.norm(d)).tolist())


def insert_cuts(ids, xyz_m, radii_m, edges, cuts):
    """cuts={edge_index: [(fraction, unique_virtual_id), ...]}.

    Returns geometry plus original-edge index and fraction interval per segment.
    Repeated cuts at the same point are rejected rather than making zero edges.
    """
    new_ids = list(map(int, ids)); xyz = list(np.asarray(xyz_m)); rad = list(radii_m)
    result, origin, interval = [], [], []
    if set(cuts)-set(range(len(edges))):
        raise ValueError('Unknown cut edge')
    for eid, (u, v) in enumerate(edges):
        sequence = [(0., int(u))]
        prior = 0.
        for t, virtual_id in sorted(cuts.get(eid, [])):
            if not prior < t < 1 or virtual_id in new_ids:
                raise ValueError('Cuts must be strictly interior, distinct, with unique IDs')
            point, radius = edge_position(xyz_m[u], xyz_m[v], radii_m[u], radii_m[v], t)
            sequence.append((t, len(new_ids))); new_ids.append(virtual_id)
            xyz.append(point); rad.append(radius); prior = t
        sequence.append((1., int(v)))
        for (a, i), (b, j) in zip(sequence, sequence[1:]):
            result.append((i, j)); origin.append(eid); interval.append((a, b))
    return dict(ids=np.asarray(new_ids), xyz_m=np.asarray(xyz), radius_m=np.asarray(rad),
                edges=np.asarray(result, dtype=int), original_edge=np.asarray(origin), fractions=np.asarray(interval))


def clip_interval(x0, x1, lower, upper):
    enter, leave = 0., 1.
    for x, dx, lo, hi in zip(x0, np.asarray(x1)-x0, lower, upper):
        if abs(dx) < 1e-14:
            if x < lo or x > hi:
                return None
        else:
            a, b = sorted(((lo-x)/dx, (hi-x)/dx))
            enter, leave = max(enter, a), min(leave, b)
    return (enter, leave) if leave-enter > 1e-12 else None
