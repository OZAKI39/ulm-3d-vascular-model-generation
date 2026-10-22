"""Diagnostic arithmetic only; production samplers/solvers supply all physics."""
import numpy as np


def held_alignment(t, x, saved_velocity):
    dt = np.diff(t)
    if not (np.isfinite(t).all() and np.isfinite(x).all() and
            np.isfinite(saved_velocity).all() and np.all(dt > 0)):
        raise ValueError('NONFINITE_OR_NONMONOTONE_ACCEPTED_STATES')
    held = saved_velocity[1:]
    error = np.diff(x, axis=0) - dt[:, None]*held
    return dt, held, error


def increment_loads(R0, R1, b0, b1, U0, U1, q0, q1):
    r0 = R0 @ U0 - b0 - q0
    r1 = R1 @ U1 - b1 - q1
    loads = np.array([b1-b0, -(R1-R0)@U0, q1-q0, r1-r0])
    return loads, R1@(U1-U0), r0, r1


def segment_cells(field, a, b):
    """Clip the actual segment against every candidate tetra, not endpoint IDs.

    No field polynomial is evaluated outside its cell. Tiny gaps/overlaps at
    face roots are allowed only within the production geometry roundoff bound.
    Midpoints must independently be strictly inside the enumerated cell.
    """
    import vtk
    d = b-a
    length = np.linalg.norm(d)
    if length == 0:
        return [], []
    pad = field.geometry.candidate_padding_m
    ids = vtk.vtkIdList()
    bounds = np.column_stack((np.minimum(a, b)-pad, np.maximum(a, b)+pad)).ravel()
    field._locator.FindCellsWithinBounds(bounds, ids)
    ids = np.array(sorted(ids.GetId(i) for i in range(ids.GetNumberOfIds())), dtype=int)
    w0 = field.geometry.weights(a, ids)
    wd = np.einsum('nij,j->ni', field.geometry.shape_gradients_m_inv[ids], d)
    spans = []
    for cid, w, dw in zip(ids, w0, wd):
        lo, hi = 0., 1.
        for z, slope in zip(w, dw):
            if slope > 0:
                lo = max(lo, -z/slope)
            elif slope < 0:
                hi = min(hi, -z/slope)
            elif z < -field.geometry.weight_tolerance[cid]:
                hi = -1
        if hi > lo and (hi-lo)*length > 4*pad:
            mid = a + ((lo+hi)/2)*d
            weights = field.geometry.weights(mid, np.array([cid]))[0]
            if np.min(weights) < 0:
                raise ValueError(('SEGMENT_MIDPOINT_NOT_STRICTLY_CONTAINED', int(cid)))
            spans.append((float(lo), float(hi), int(cid)))
    spans.sort()
    if not spans:
        # An entire segment shorter than roundoff resolution is not an interface test.
        if length <= 8*pad:
            return [], []
        raise ValueError('SEGMENT_WITHOUT_VALID_CELL')
    gaps = []
    covered = 0.
    for lo, hi, cid in spans:
        if (lo-covered)*length > 8*pad:
            raise ValueError(('UNCOVERED_REAL_SEGMENT', lo, covered))
        covered = max(covered, hi)
    if (1-covered)*length > 8*pad:
        raise ValueError('UNCOVERED_REAL_SEGMENT_END')
    crossings = []
    for left, right in zip(spans, spans[1:]):
        if left[2] == right[2]:
            continue
        discrepancy = abs(left[1]-right[0])*length
        if discrepancy > 8*pad:
            raise ValueError(('OVERLAPPING_OR_GAPPED_CELL_PARTITIONS', discrepancy))
        crossings.append(dict(fraction=(left[1]+right[0])/2, left=left[2], right=right[2],
                              left_fraction=left[0], right_fraction=right[1],
                              root_roundoff_discrepancy_m=discrepancy))
    return spans, crossings


def strict_probe(field, position, expected_cell=None):
    cid, weights = field.locate(position)
    if cid < 0 or not np.isfinite(weights).all():
        raise ValueError('INTERFACE_PROBE_OUTSIDE')
    if np.min(weights) < 0 or np.max(weights) > 1:
        raise ValueError('INTERFACE_PROBE_NOT_STRICTLY_CONTAINED')
    if expected_cell is not None and cid != expected_cell:
        raise ValueError(('INTERFACE_PROBE_WRONG_SIDE', cid, expected_cell))
    return cid, weights


def direction_angle(a, b):
    aa, bb = np.linalg.norm(a), np.linalg.norm(b)
    if aa == 0 or bb == 0:
        return None
    return float(np.degrees(np.arccos(np.clip(np.dot(a, b)/(aa*bb), -1, 1))))
