"""Internal-step molecular target exposure along accepted particle paths.

This module deliberately owns no persistent trajectory state.  A caller gives
it one ragged polyline per permanent particle ID for one physical time step;
the routine integrates the target-positive capture exposure over that path and
returns one compact row per particle.  Deterministic motion and any subsequent
RBC stochastic reflection/topology pieces therefore share *one* ``dt`` rather
than each consuming a separate physical time interval.

All path coordinates are physical X-Z coordinates in micrometres.  The wall
gap is queried from the authoritative continuous geometry and the target area
is delegated to :meth:`MolecularTargetField.reaction_area_um2`; no raster wall
or vessel labels are consulted here.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from ..geometry.continuous_vessel_geometry import ContinuousVesselGeometry
from ..molecular.molecular_binding import reaction_disk_radius_um
from .particle_geometry_tolerance import (
    particle_geometry_roundoff_tolerance_um,
)


# Eight fixed Gauss-Legendre nodes give a deterministic, allocation-bounded
# integration rule *within each exposed sub-interval*.  All accepted path
# vertices and the segment midpoint are separate discovery probes.  Whenever
# adjacent probes disagree about exposure, the transition is located by a
# deterministic batched bisection before quadrature.  These are numerical
# implementation constants, not configurable physical parameters.
_GAUSS_NODES, _GAUSS_WEIGHTS = np.polynomial.legendre.leggauss(8)
_UNIT_GAUSS_NODES = 0.5 * (_GAUSS_NODES + 1.0)
_UNIT_GAUSS_WEIGHTS = 0.5 * _GAUSS_WEIGHTS
_DISCOVERY_FRACTIONS = np.unique(
    np.concatenate(
        (
            np.asarray([0.0, 0.5, 1.0], dtype=np.float64),
            _UNIT_GAUSS_NODES,
        )
    )
)
_DISCOVERY_GAUSS_INDICES = np.searchsorted(
    _DISCOVERY_FRACTIONS, _UNIT_GAUSS_NODES
)
_ROOT_FRACTION_TOLERANCE = 64.0 * np.finfo(np.float64).eps
_ROOT_MAX_ITERATIONS = 52


@dataclass(frozen=True, slots=True)
class InternalTargetExposureBatch:
    """Compact per-particle exposure increments for one physical step.

    ``exposure_open_at_end`` is the state to feed into the next step.
    ``right_censored_at_end`` is set only when the caller explicitly marks the
    end of its observation window; merely reaching an internal step boundary
    does not right-censor an open event.
    """

    permanent_ids: np.ndarray
    exposure_time_s: np.ndarray
    reaction_area_time_um2_s: np.ndarray
    quantitatively_applicable_exposure_time_s: np.ndarray
    event_start_count: np.ndarray
    event_end_count: np.ndarray
    event_started: np.ndarray
    event_ended: np.ndarray
    exposure_open_at_end: np.ndarray
    right_censored_at_end: np.ndarray


def integrate_internal_target_exposure(
    *,
    permanent_ids: np.ndarray,
    path_points_xz_um: np.ndarray,
    path_point_offsets: np.ndarray,
    dt_s: float | np.ndarray,
    radius_um: float | np.ndarray,
    capture_distance_um: float,
    boundary_geometry: ContinuousVesselGeometry,
    molecular_target_field: Any,
    initially_exposed: bool | np.ndarray = False,
    quantitatively_applicable: bool | np.ndarray = True,
    terminated_at_end: bool | np.ndarray = False,
    observation_ends_at_step_end: bool | np.ndarray = False,
    path_time_fraction: np.ndarray | None = None,
) -> InternalTargetExposureBatch:
    """Integrate target exposure for a ragged batch of accepted path polylines.

    Parameters
    ----------
    permanent_ids:
        Unique, stable particle IDs.  The result preserves this order.
    path_points_xz_um, path_point_offsets:
        CSR-style ragged polylines.  Lane ``i`` uses
        ``points[offsets[i]:offsets[i + 1]]`` and must contain at least one
        point.  Reflection contacts and topology crossings should be retained
        as intermediate vertices so transient capture-layer visits are visible.
    dt_s:
        The single physical-step duration, scalar or one value per lane.
    path_time_fraction:
        Optional CSR-aligned normalized time coordinate.  Each lane must start
        at zero, end at one, and be nondecreasing.  When omitted, the one ``dt``
        is apportioned by cumulative polyline arclength.  This is the stable
        operator-split convention used when a stochastic geometry sweep has no
        separate physical time.  A stationary path is evaluated at its sole
        position for the full ``dt``.

    Notes
    -----
    Geometry and target queries are issued in batches, never once per particle.
    Every accepted path vertex is sampled, including stochastic reflection and
    topology-event vertices.  Exposure transitions bracketed by the vertices,
    segment midpoint, or quadrature discovery probes are refined to a fixed
    machine-precision-scaled tolerance.  Exposure time is then accumulated from
    the resulting positive sub-interval lengths, while area-time uses an
    eight-point Gauss-Legendre rule independently on each positive interval.
    Thus a very short enter-and-exit visit around an accepted path vertex cannot
    disappear merely because no fixed whole-segment quadrature node lands in it.
    """

    ids = np.asarray(permanent_ids, dtype=np.int64).reshape(-1)
    lane_count = int(ids.size)
    if np.unique(ids).size != lane_count:
        raise ValueError("permanent_ids must be unique within a batch.")

    points = np.asarray(path_points_xz_um, dtype=np.float64)
    if points.ndim != 2 or points.shape[1] != 2:
        raise ValueError("path_points_xz_um must have shape (point_count, 2).")
    if np.any(~np.isfinite(points)):
        raise ValueError("path_points_xz_um must contain only finite values.")

    offsets = np.asarray(path_point_offsets, dtype=np.int64).reshape(-1)
    if offsets.shape != (lane_count + 1,):
        raise ValueError("path_point_offsets must have length particle_count + 1.")
    if (
        offsets.size == 0
        or int(offsets[0]) != 0
        or int(offsets[-1]) != int(points.shape[0])
        or np.any(np.diff(offsets) < 1)
    ):
        raise ValueError(
            "path_point_offsets must start at zero, end at point_count, and "
            "assign at least one point to every particle."
        )

    durations = _broadcast_lane_values(dt_s, lane_count, "dt_s", np.float64)
    if np.any(~np.isfinite(durations)) or np.any(durations < 0.0):
        raise ValueError("dt_s must contain only finite non-negative values.")
    radii = _broadcast_lane_values(radius_um, lane_count, "radius_um", np.float64)
    if np.any(~np.isfinite(radii)) or np.any(radii <= 0.0):
        raise ValueError("radius_um must contain only finite positive values.")
    capture = float(capture_distance_um)
    if not np.isfinite(capture) or capture < 0.0:
        raise ValueError("capture_distance_um must be finite and non-negative.")
    if np.any(capture > 2.0 * radii):
        raise ValueError("capture_distance_um must not exceed twice radius_um.")

    initially_open = _broadcast_lane_values(
        initially_exposed, lane_count, "initially_exposed", np.bool_
    )
    applicable = _broadcast_lane_values(
        quantitatively_applicable,
        lane_count,
        "quantitatively_applicable",
        np.bool_,
    )
    terminated = _broadcast_lane_values(
        terminated_at_end, lane_count, "terminated_at_end", np.bool_
    )
    observation_ends = _broadcast_lane_values(
        observation_ends_at_step_end,
        lane_count,
        "observation_ends_at_step_end",
        np.bool_,
    )

    if lane_count == 0:
        return _empty_result(ids)

    point_counts = np.diff(offsets)
    segment_counts = point_counts - 1
    segment_lane = np.repeat(np.arange(lane_count, dtype=np.int64), segment_counts)
    segment_start_index = _segment_start_indices(offsets)
    segment_start = points[segment_start_index]
    segment_end = points[segment_start_index + 1]
    segment_delta = segment_end - segment_start
    segment_length = np.linalg.norm(segment_delta, axis=1)

    segment_fraction, _segment_fraction_start = _segment_time_fractions(
        points=points,
        offsets=offsets,
        segment_lane=segment_lane,
        segment_length=segment_length,
        path_time_fraction=path_time_fraction,
    )

    positive_time_segment = segment_fraction > 0.0
    active_segment_indices = np.flatnonzero(positive_time_segment)
    active_segment_count = int(active_segment_indices.size)
    active_lanes = segment_lane[active_segment_indices]
    active_start = segment_start[active_segment_indices]
    active_delta = segment_delta[active_segment_indices]

    # Query every supplied vertex exactly once.  Interior discovery probes are
    # added only for positive-time segments.  This explicit vertex sampling is
    # what makes an arbitrarily short visit centred on a reflection/topology
    # vertex observable even when every whole-segment Gauss node is outside.
    if active_segment_count:
        discovery_points = (
            active_start[:, None, :]
            + _DISCOVERY_FRACTIONS[None, :, None] * active_delta[:, None, :]
        )
        interior_discovery_points = discovery_points[:, 1:-1, :].reshape(-1, 2)
        interior_discovery_lanes = np.repeat(
            active_lanes, _DISCOVERY_FRACTIONS.size - 2
        )
    else:
        interior_discovery_points = np.empty((0, 2), dtype=np.float64)
        interior_discovery_lanes = np.empty(0, dtype=np.int64)

    point_lane = np.repeat(np.arange(lane_count, dtype=np.int64), point_counts)
    initial_query_points = np.vstack((points, interior_discovery_points))
    initial_query_lanes = np.concatenate((point_lane, interior_discovery_lanes))
    initial_area = _target_reaction_area_at_points(
        points_xz_um=initial_query_points,
        lane_index=initial_query_lanes,
        lane_radii_um=radii,
        capture_distance_um=capture,
        boundary_geometry=boundary_geometry,
        molecular_target_field=molecular_target_field,
    )
    if np.count_nonzero(initial_area) == 0:
        zeros_float = np.zeros(lane_count, dtype=np.float64)
        zeros_int = np.zeros(lane_count, dtype=np.int64)
        zeros_bool = np.zeros(lane_count, dtype=bool)
        event_end_count = initially_open.astype(np.int64, copy=True)
        return InternalTargetExposureBatch(
            permanent_ids=ids.copy(),
            exposure_time_s=zeros_float,
            reaction_area_time_um2_s=zeros_float.copy(),
            quantitatively_applicable_exposure_time_s=zeros_float.copy(),
            event_start_count=zeros_int,
            event_end_count=event_end_count,
            event_started=zeros_bool,
            event_ended=initially_open.copy(),
            exposure_open_at_end=zeros_bool.copy(),
            right_censored_at_end=zeros_bool.copy(),
        )
    vertex_area = initial_area[: points.shape[0]]
    vertex_inside = vertex_area > 0.0

    discovery_area = np.empty(
        (active_segment_count, _DISCOVERY_FRACTIONS.size), dtype=np.float64
    )
    if active_segment_count:
        discovery_area[:, 0] = vertex_area[
            segment_start_index[active_segment_indices]
        ]
        discovery_area[:, -1] = vertex_area[
            segment_start_index[active_segment_indices] + 1
        ]
        discovery_area[:, 1:-1] = initial_area[points.shape[0] :].reshape(
            active_segment_count, _DISCOVERY_FRACTIONS.size - 2
        )
    discovery_inside = discovery_area > 0.0

    left_fraction = np.broadcast_to(
        _DISCOVERY_FRACTIONS[:-1],
        (active_segment_count, _DISCOVERY_FRACTIONS.size - 1),
    )
    right_fraction = np.broadcast_to(
        _DISCOVERY_FRACTIONS[1:],
        (active_segment_count, _DISCOVERY_FRACTIONS.size - 1),
    )
    left_inside = discovery_inside[:, :-1]
    right_inside = discovery_inside[:, 1:]
    transition = left_inside != right_inside
    transition_segment_row, transition_cell = np.nonzero(transition)
    transition_root = np.full(transition.shape, np.nan, dtype=np.float64)
    if transition_segment_row.size:
        roots = _bisect_exposure_transitions(
            segment_start=active_start[transition_segment_row],
            segment_delta=active_delta[transition_segment_row],
            lane_index=active_lanes[transition_segment_row],
            lane_radii_um=radii,
            left_fraction=left_fraction[transition_segment_row, transition_cell],
            right_fraction=right_fraction[transition_segment_row, transition_cell],
            left_inside=left_inside[transition_segment_row, transition_cell],
            capture_distance_um=capture,
            boundary_geometry=boundary_geometry,
            molecular_target_field=molecular_target_field,
        )
        transition_root[transition_segment_row, transition_cell] = roots

    # Each discovery cell is either wholly positive/negative or has exactly
    # one bracketed, refined transition.  Convert the positive portion into a
    # disjoint sub-interval.  Exposure time follows directly from its length.
    exposed_cell = left_inside | transition
    exposed_start = np.where(
        transition & ~left_inside,
        transition_root,
        left_fraction,
    )
    exposed_end = np.where(
        transition & left_inside,
        transition_root,
        right_fraction,
    )
    exposed_width = np.where(
        exposed_cell, np.maximum(exposed_end - exposed_start, 0.0), 0.0
    )
    segment_exposed_fraction = np.sum(exposed_width, axis=1)
    exposure_time = np.bincount(
        active_lanes,
        weights=(
            durations[active_lanes]
            * segment_fraction[active_segment_indices]
            * segment_exposed_fraction
        ),
        minlength=lane_count,
    ).astype(np.float64, copy=False)

    # Reuse the original whole-segment Gauss samples when the entire segment is
    # exposed.  Only partially exposed segments need a second batched query on
    # their root-delimited positive intervals.
    fully_exposed_segment = np.all(discovery_inside, axis=1)
    area_time_segment = np.zeros(active_segment_count, dtype=np.float64)
    if np.any(fully_exposed_segment):
        whole_rows = np.flatnonzero(fully_exposed_segment)
        whole_gauss_area = discovery_area[
            whole_rows[:, None], _DISCOVERY_GAUSS_INDICES[None, :]
        ]
        area_time_segment[whole_rows] = np.sum(
            whole_gauss_area * _UNIT_GAUSS_WEIGHTS[None, :], axis=1
        )

    partial_exposed = exposed_cell & ~fully_exposed_segment[:, None]
    partial_row, partial_cell = np.nonzero(partial_exposed)
    if partial_row.size:
        partial_start = exposed_start[partial_row, partial_cell]
        partial_width = exposed_width[partial_row, partial_cell]
        integration_fraction = (
            partial_start[:, None]
            + partial_width[:, None] * _UNIT_GAUSS_NODES[None, :]
        )
        integration_points = (
            active_start[partial_row, None, :]
            + integration_fraction[:, :, None]
            * active_delta[partial_row, None, :]
        ).reshape(-1, 2)
        integration_lanes = np.repeat(
            active_lanes[partial_row], _UNIT_GAUSS_NODES.size
        )
        integration_area = _target_reaction_area_at_points(
            points_xz_um=integration_points,
            lane_index=integration_lanes,
            lane_radii_um=radii,
            capture_distance_um=capture,
            boundary_geometry=boundary_geometry,
            molecular_target_field=molecular_target_field,
        ).reshape(partial_row.size, _UNIT_GAUSS_NODES.size)
        partial_area_integral = partial_width * np.sum(
            integration_area * _UNIT_GAUSS_WEIGHTS[None, :], axis=1
        )
        np.add.at(area_time_segment, partial_row, partial_area_integral)

    area_time = np.bincount(
        active_lanes,
        weights=(
            durations[active_lanes]
            * segment_fraction[active_segment_indices]
            * area_time_segment
        ),
        minlength=lane_count,
    ).astype(np.float64, copy=False)

    # A lane without a positive-time segment is stationary and consumes its
    # complete dt at the first (and necessarily only physical) position.
    moving_lane = np.bincount(
        segment_lane,
        weights=positive_time_segment.astype(np.int64),
        minlength=lane_count,
    ) > 0
    stationary_lanes = np.flatnonzero(~moving_lane)
    stationary_area = vertex_area[offsets[stationary_lanes]]
    stationary_inside = stationary_area > 0.0
    exposure_time[stationary_lanes] = (
        durations[stationary_lanes] * stationary_inside
    )
    area_time[stationary_lanes] = durations[stationary_lanes] * stationary_area
    applicable_exposure_time = exposure_time * applicable.astype(np.float64)

    # Event counting follows the same ordered discovery states.  Each segment
    # contributes its interior transitions once; shared vertices are not
    # double-counted.  Zero-time geometric pieces still preserve event state,
    # while contributing no exposure duration.
    first_inside = vertex_inside[offsets[:-1]]
    event_start_count = ((~initially_open) & first_inside).astype(np.int64)
    event_end_count = (initially_open & ~first_inside).astype(np.int64)
    if segment_lane.size:
        inactive_segment_indices = np.flatnonzero(~positive_time_segment)
        segment_start_count = np.zeros(segment_lane.size, dtype=np.int64)
        segment_end_count = np.zeros(segment_lane.size, dtype=np.int64)
        if active_segment_count:
            segment_start_count[active_segment_indices] = np.sum(
                (~left_inside) & right_inside, axis=1, dtype=np.int64
            )
            segment_end_count[active_segment_indices] = np.sum(
                left_inside & (~right_inside), axis=1, dtype=np.int64
            )
        if inactive_segment_indices.size:
            inactive_start_point = segment_start_index[inactive_segment_indices]
            inactive_left = vertex_inside[inactive_start_point]
            inactive_right = vertex_inside[inactive_start_point + 1]
            segment_start_count[inactive_segment_indices] = (
                (~inactive_left) & inactive_right
            )
            segment_end_count[inactive_segment_indices] = (
                inactive_left & (~inactive_right)
            )
        event_start_count += np.bincount(
            segment_lane,
            weights=segment_start_count,
            minlength=lane_count,
        ).astype(np.int64, copy=False)
        event_end_count += np.bincount(
            segment_lane,
            weights=segment_end_count,
            minlength=lane_count,
        ).astype(np.int64, copy=False)

    final_inside = vertex_inside[offsets[1:] - 1]
    open_at_end = final_inside.copy()
    terminated_while_open = terminated & open_at_end
    event_end_count = event_end_count + terminated_while_open.astype(np.int64)
    open_at_end = open_at_end & ~terminated
    right_censored = open_at_end & observation_ends

    return InternalTargetExposureBatch(
        permanent_ids=ids.copy(),
        exposure_time_s=exposure_time,
        reaction_area_time_um2_s=area_time,
        quantitatively_applicable_exposure_time_s=applicable_exposure_time,
        event_start_count=event_start_count,
        event_end_count=event_end_count,
        event_started=event_start_count > 0,
        event_ended=event_end_count > 0,
        exposure_open_at_end=open_at_end,
        right_censored_at_end=right_censored,
    )


def _bisect_exposure_transitions(
    *,
    segment_start: np.ndarray,
    segment_delta: np.ndarray,
    lane_index: np.ndarray,
    lane_radii_um: np.ndarray,
    left_fraction: np.ndarray,
    right_fraction: np.ndarray,
    left_inside: np.ndarray,
    capture_distance_um: float,
    boundary_geometry: ContinuousVesselGeometry,
    molecular_target_field: Any,
) -> np.ndarray:
    """Locate bracketed ``area > 0`` transitions in deterministic batches.

    The predicate is intentionally the same one used for event counting and
    exposure integration.  Bisection therefore remains valid at a capture-layer
    boundary where the continuous reaction area is exactly zero on the outside
    rather than changing algebraic sign.  The iteration cap and tolerance are
    fixed numerical constants and do not alter model physics.
    """

    start = np.asarray(segment_start, dtype=np.float64).reshape(-1, 2)
    delta = np.asarray(segment_delta, dtype=np.float64).reshape(-1, 2)
    lanes = np.asarray(lane_index, dtype=np.int64).reshape(-1)
    left = np.asarray(left_fraction, dtype=np.float64).reshape(-1).copy()
    right = np.asarray(right_fraction, dtype=np.float64).reshape(-1).copy()
    state_left = np.asarray(left_inside, dtype=bool).reshape(-1).copy()
    if not (
        start.shape == delta.shape == (lanes.size, 2)
        and left.shape == right.shape == state_left.shape == (lanes.size,)
    ):
        raise ValueError("Exposure-transition bracket arrays are inconsistent.")
    if lanes.size == 0:
        return np.empty(0, dtype=np.float64)

    for _ in range(_ROOT_MAX_ITERATIONS):
        unresolved = (right - left) > _ROOT_FRACTION_TOLERANCE
        if not np.any(unresolved):
            break
        rows = np.flatnonzero(unresolved)
        midpoint = 0.5 * (left[rows] + right[rows])
        midpoint_points = start[rows] + midpoint[:, None] * delta[rows]
        midpoint_area = _target_reaction_area_at_points(
            points_xz_um=midpoint_points,
            lane_index=lanes[rows],
            lane_radii_um=lane_radii_um,
            capture_distance_um=capture_distance_um,
            boundary_geometry=boundary_geometry,
            molecular_target_field=molecular_target_field,
        )
        midpoint_inside = midpoint_area > 0.0
        same_as_left = midpoint_inside == state_left[rows]
        left_rows = rows[same_as_left]
        right_rows = rows[~same_as_left]
        left[left_rows] = midpoint[same_as_left]
        right[right_rows] = midpoint[~same_as_left]
        state_left[left_rows] = midpoint_inside[same_as_left]

    return 0.5 * (left + right)


def _target_reaction_area_at_points(
    *,
    points_xz_um: np.ndarray,
    lane_index: np.ndarray,
    lane_radii_um: np.ndarray,
    capture_distance_um: float,
    boundary_geometry: ContinuousVesselGeometry,
    molecular_target_field: Any,
) -> np.ndarray:
    """Evaluate continuous gap and target reaction area in one batch."""

    if points_xz_um.shape[0] == 0:
        return np.empty(0, dtype=np.float64)
    exact = boundary_geometry.exact_solid_wall_state_xz_um_accelerated(
        points_xz_um
    )
    distance = np.asarray(exact.distance_um, dtype=np.float64).reshape(-1)
    inward_normal = np.asarray(exact.inward_normal_xz, dtype=np.float64).reshape(-1, 2)
    if distance.shape != (points_xz_um.shape[0],) or inward_normal.shape != points_xz_um.shape:
        raise ValueError("The continuous wall query returned incompatible array shapes.")
    normal_norm = np.linalg.norm(inward_normal, axis=1)
    if np.any(~np.isfinite(normal_norm)) or np.any(normal_norm <= 0.0):
        raise ValueError("The continuous wall query returned invalid inward normals.")
    sample_radius = lane_radii_um[lane_index]
    physical_gap = distance - sample_radius
    geometry_roundoff = particle_geometry_roundoff_tolerance_um(
        points_xz_um,
        sample_radius,
        capture_distance_um,
    )
    if np.any(physical_gap < -geometry_roundoff):
        raise ValueError(
            "gap_um contains a materially negative true wall gap; molecular "
            "capture cannot repair an invalid particle position."
        )
    physical_gap = np.maximum(physical_gap, 0.0)
    reaction_radius = reaction_disk_radius_um(
        sample_radius, physical_gap, capture_distance_um
    )
    positive_reaction = reaction_radius > 0.0
    positive_count = int(np.count_nonzero(positive_reaction))
    area = np.zeros(points_xz_um.shape[0], dtype=np.float64)
    if positive_count:
        positive_unit_normal = (
            inward_normal[positive_reaction]
            / normal_norm[positive_reaction, None]
        )
        positive_tangent = np.column_stack(
            (-positive_unit_normal[:, 1], positive_unit_normal[:, 0])
        )
        positive_radius = reaction_radius[positive_reaction]
        positive_area = np.asarray(
            molecular_target_field.reaction_area_um2(
                points_xz_um[positive_reaction],
                positive_tangent,
                positive_radius,
            ),
            dtype=np.float64,
        ).reshape(-1)
        if positive_area.shape != (positive_count,):
            raise ValueError("reaction_area_um2 must return one value per query point.")
        if np.any(~np.isfinite(positive_area)) or np.any(positive_area < 0.0):
            raise ValueError("reaction_area_um2 returned invalid reaction areas.")
        positive_maximum = np.pi * positive_radius * positive_radius
        tolerance = 256.0 * np.finfo(np.float64).eps * np.maximum(
            1.0, positive_maximum
        )
        if np.any(positive_area > positive_maximum + tolerance):
            raise ValueError(
                "reaction_area_um2 exceeded the physical reaction disk area."
            )
        area[positive_reaction] = np.minimum(positive_area, positive_maximum)
    return area


def _segment_start_indices(offsets: np.ndarray) -> np.ndarray:
    point_count = int(offsets[-1])
    if point_count <= 1:
        return np.empty(0, dtype=np.int64)
    candidates = np.arange(point_count - 1, dtype=np.int64)
    lane_last_points = offsets[1:] - 1
    return candidates[~np.isin(candidates, lane_last_points)]


def _segment_time_fractions(
    *,
    points: np.ndarray,
    offsets: np.ndarray,
    segment_lane: np.ndarray,
    segment_length: np.ndarray,
    path_time_fraction: np.ndarray | None,
) -> tuple[np.ndarray, np.ndarray]:
    lane_count = offsets.size - 1
    segment_counts = np.diff(offsets) - 1
    segment_offsets = np.concatenate(
        (np.zeros(1, dtype=np.int64), np.cumsum(segment_counts, dtype=np.int64))
    )
    if path_time_fraction is not None:
        time = np.asarray(path_time_fraction, dtype=np.float64).reshape(-1)
        if time.shape != (points.shape[0],) or np.any(~np.isfinite(time)):
            raise ValueError(
                "path_time_fraction must be finite and aligned with path points."
            )
        tolerance = 128.0 * np.finfo(np.float64).eps
        point_count = np.diff(offsets)
        first = offsets[:-1]
        last = offsets[1:] - 1
        stationary = point_count == 1
        if np.any(stationary & (np.abs(time[first]) > tolerance)):
            raise ValueError(
                "A stationary one-point path_time_fraction lane must equal 0."
            )
        moving = ~stationary
        starts = _segment_start_indices(offsets)
        if (
            np.any(moving & (np.abs(time[first]) > tolerance))
            or np.any(moving & (np.abs(time[last] - 1.0) > tolerance))
            or np.any(time[starts + 1] - time[starts] < -tolerance)
        ):
            raise ValueError(
                "Each path_time_fraction lane must run monotonically from 0 to 1."
            )
        if segment_length.size == 0:
            return np.empty(0, dtype=np.float64), np.empty(0, dtype=np.float64)
        fraction_start = time[starts]
        fraction = np.maximum(time[starts + 1] - fraction_start, 0.0)
        return fraction, fraction_start

    total_length = np.bincount(
        segment_lane, weights=segment_length, minlength=lane_count
    ).astype(np.float64, copy=False)
    fraction = np.divide(
        segment_length,
        total_length[segment_lane],
        out=np.zeros_like(segment_length),
        where=total_length[segment_lane] > 0.0,
    )
    if fraction.size == 0:
        return fraction, np.empty(0, dtype=np.float64)
    exclusive_cumulative = np.cumsum(fraction, dtype=np.float64) - fraction
    lane_base = np.zeros(lane_count, dtype=np.float64)
    nonempty = segment_counts > 0
    lane_base[nonempty] = exclusive_cumulative[
        segment_offsets[:-1][nonempty]
    ]
    fraction_start = exclusive_cumulative - lane_base[segment_lane]
    return fraction, fraction_start


def _broadcast_lane_values(
    value: Any,
    lane_count: int,
    name: str,
    dtype: Any,
) -> np.ndarray:
    array = np.asarray(value, dtype=dtype)
    if array.ndim == 0:
        return np.full(lane_count, array.item(), dtype=dtype)
    array = array.reshape(-1)
    if array.shape != (lane_count,):
        raise ValueError(f"{name} must be scalar or have one value per particle.")
    return array.copy()


def _empty_result(ids: np.ndarray) -> InternalTargetExposureBatch:
    zeros_float = np.zeros(0, dtype=np.float64)
    zeros_int = np.zeros(0, dtype=np.int64)
    zeros_bool = np.zeros(0, dtype=bool)
    return InternalTargetExposureBatch(
        permanent_ids=ids.copy(),
        exposure_time_s=zeros_float,
        reaction_area_time_um2_s=zeros_float.copy(),
        quantitatively_applicable_exposure_time_s=zeros_float.copy(),
        event_start_count=zeros_int,
        event_end_count=zeros_int.copy(),
        event_started=zeros_bool,
        event_ended=zeros_bool.copy(),
        exposure_open_at_end=zeros_bool.copy(),
        right_censored_at_end=zeros_bool.copy(),
    )


__all__ = [
    "InternalTargetExposureBatch",
    "integrate_internal_target_exposure",
]
