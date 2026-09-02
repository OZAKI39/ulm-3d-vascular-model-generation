"""Kinematics derived from geometry-accepted, saved particle positions."""

from __future__ import annotations

import numpy as np


def realized_center_velocities_um_s(
    frame_offsets: np.ndarray,
    bubble_id: np.ndarray,
    positions_um: np.ndarray,
    output_dt_s: float,
) -> np.ndarray:
    """Return same-ID frame displacement velocities for ragged trajectory records.

    The mobility right-hand side describes the instantaneous velocity requested
    by the physical model.  A wall constraint may accept a different center
    displacement.  This function measures that realized displacement without
    ever differencing two different permanent bubble IDs.
    """

    offsets = np.asarray(frame_offsets, dtype=np.int64)
    ids = np.asarray(bubble_id, dtype=np.int64)
    positions = np.asarray(positions_um, dtype=np.float64)
    if offsets.ndim != 1 or offsets.size < 1:
        raise ValueError("frame_offsets must be a non-empty one-dimensional array.")
    if offsets[0] != 0 or offsets[-1] != ids.size:
        raise ValueError("frame_offsets must span every ragged trajectory record.")
    if positions.shape != (ids.size, 3):
        raise ValueError("positions_um must have shape (record_count, 3).")
    dt_s = float(output_dt_s)
    if not np.isfinite(dt_s) or dt_s <= 0.0:
        raise ValueError("output_dt_s must be finite and positive.")

    realized = np.zeros_like(positions, dtype=np.float64)
    for frame_index in range(max(0, offsets.size - 2)):
        start0 = int(offsets[frame_index])
        stop0 = int(offsets[frame_index + 1])
        start1 = int(offsets[frame_index + 1])
        stop1 = int(offsets[frame_index + 2])
        if start0 == stop0 or start1 == stop1:
            continue
        _, local0, local1 = np.intersect1d(
            ids[start0:stop0],
            ids[start1:stop1],
            assume_unique=True,
            return_indices=True,
        )
        if local0.size == 0:
            continue
        record0 = start0 + local0
        record1 = start1 + local1
        velocity = (positions[record1] - positions[record0]) / dt_s
        realized[record0] = velocity
        realized[record1] = velocity
    return realized
