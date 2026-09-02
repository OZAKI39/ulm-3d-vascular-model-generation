"""Compiled broad-phase queries for Revised-v16 continuous vessel walls.

The geometry owner constructs conservative CSR candidate bins.  Kernels in
this module only evaluate exact point/segment and swept-disc/segment formulas;
they never infer a wall from the CFD mask.
"""

from __future__ import annotations

import math

import numpy as np

try:
    from numba import njit, prange
except ImportError:  # pragma: no cover
    njit = None
    prange = range


NUMBA_AVAILABLE = njit is not None

# Four-worker, pre-warmed measurements on resampled batches from the checked-in
# 20260721 formal snapshot put the stable crossover at 28 particles.  Keep
# smaller batches on the lower-overhead serial dispatcher while sending the
# sustained production batches (typically 70+) through the parallel kernel.
_EXACT_STATE_PARALLEL_MIN_BATCH = 28


def exact_continuous_wall_states(
    positions_xz_um: np.ndarray,
    segment_start_xz_um: np.ndarray,
    segment_end_xz_um: np.ndarray,
    segment_inward_normal_xz: np.ndarray,
    segment_ring_index: np.ndarray,
    ring_size: int,
    bin_edge_origin_xz_um: np.ndarray,
    bin_size_um: float,
    bin_shape: np.ndarray,
    bin_offsets: np.ndarray,
    bin_segment_indices: np.ndarray,
    *,
    use_numba: bool,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    if use_numba and njit is not None:
        kernel = (
            _exact_continuous_wall_states_parallel_numba
            if int(positions_xz_um.shape[0]) >= _EXACT_STATE_PARALLEL_MIN_BATCH
            else _exact_continuous_wall_states_numba
        )
    else:
        kernel = _exact_continuous_wall_states_kernel
    return kernel(
        positions_xz_um,
        segment_start_xz_um,
        segment_end_xz_um,
        segment_inward_normal_xz,
        segment_ring_index,
        int(ring_size),
        bin_edge_origin_xz_um,
        float(bin_size_um),
        bin_shape,
        bin_offsets,
        bin_segment_indices,
    )


def inspect_swept_continuous_wall_paths(
    starts_xz_um: np.ndarray,
    ends_xz_um: np.ndarray,
    radii_um: np.ndarray,
    broad_phase_radius_um: np.ndarray,
    segment_start_xz_um: np.ndarray,
    segment_end_xz_um: np.ndarray,
    segment_inward_normal_xz: np.ndarray,
    segment_ring_index: np.ndarray,
    ring_size: int,
    segment_arclength_start_um: np.ndarray,
    segment_arclength_end_um: np.ndarray,
    ring_length_um: float,
    bin_edge_origin_xz_um: np.ndarray,
    bin_size_um: float,
    bin_shape: np.ndarray,
    bin_offsets: np.ndarray,
    bin_segment_indices: np.ndarray,
    tolerance_um: float,
    *,
    use_numba: bool,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    if use_numba and njit is not None:
        kernel = (
            _inspect_swept_continuous_wall_paths_parallel_numba
            if int(starts_xz_um.shape[0]) >= 64
            else _inspect_swept_continuous_wall_paths_numba
        )
    else:
        kernel = _inspect_swept_continuous_wall_paths_kernel
    return kernel(
        starts_xz_um,
        ends_xz_um,
        radii_um,
        broad_phase_radius_um,
        segment_start_xz_um,
        segment_end_xz_um,
        segment_inward_normal_xz,
        segment_ring_index,
        int(ring_size),
        segment_arclength_start_um,
        segment_arclength_end_um,
        float(ring_length_um),
        bin_edge_origin_xz_um,
        float(bin_size_um),
        bin_shape,
        bin_offsets,
        bin_segment_indices,
        float(tolerance_um),
    )


def inspect_swept_continuous_wall_paths_with_end_state(
    starts_xz_um: np.ndarray,
    ends_xz_um: np.ndarray,
    radii_um: np.ndarray,
    start_distance_um: np.ndarray,
    segment_start_xz_um: np.ndarray,
    segment_end_xz_um: np.ndarray,
    segment_inward_normal_xz: np.ndarray,
    segment_ring_index: np.ndarray,
    ring_size: int,
    segment_arclength_start_um: np.ndarray,
    segment_arclength_end_um: np.ndarray,
    ring_length_um: float,
    exact_bin_edge_origin_xz_um: np.ndarray,
    exact_bin_size_um: float,
    exact_bin_shape: np.ndarray,
    exact_bin_offsets: np.ndarray,
    exact_bin_segment_indices: np.ndarray,
    sweep_bin_edge_origin_xz_um: np.ndarray,
    sweep_bin_size_um: float,
    sweep_bin_shape: np.ndarray,
    sweep_bin_offsets: np.ndarray,
    sweep_bin_segment_indices: np.ndarray,
    tolerance_um: float,
    *,
    use_numba: bool,
) -> tuple[np.ndarray, ...]:
    """Evaluate the shared path endpoint and sweep in one compiled dispatch."""

    kernel = (
        _inspect_swept_continuous_wall_paths_with_end_state_numba
        if use_numba and njit is not None
        else _inspect_swept_continuous_wall_paths_with_end_state_kernel
    )
    return kernel(
        starts_xz_um,
        ends_xz_um,
        radii_um,
        start_distance_um,
        segment_start_xz_um,
        segment_end_xz_um,
        segment_inward_normal_xz,
        segment_ring_index,
        int(ring_size),
        segment_arclength_start_um,
        segment_arclength_end_um,
        float(ring_length_um),
        exact_bin_edge_origin_xz_um,
        float(exact_bin_size_um),
        exact_bin_shape,
        exact_bin_offsets,
        exact_bin_segment_indices,
        sweep_bin_edge_origin_xz_um,
        float(sweep_bin_size_um),
        sweep_bin_shape,
        sweep_bin_offsets,
        sweep_bin_segment_indices,
        float(tolerance_um),
    )


def inspect_swept_continuous_wall_paths_with_end_distance(
    starts_xz_um: np.ndarray,
    ends_xz_um: np.ndarray,
    radii_um: np.ndarray,
    start_distance_um: np.ndarray,
    segment_start_xz_um: np.ndarray,
    segment_end_xz_um: np.ndarray,
    segment_inward_normal_xz: np.ndarray,
    segment_ring_index: np.ndarray,
    ring_size: int,
    segment_arclength_start_um: np.ndarray,
    segment_arclength_end_um: np.ndarray,
    ring_length_um: float,
    exact_bin_edge_origin_xz_um: np.ndarray,
    exact_bin_size_um: float,
    exact_bin_shape: np.ndarray,
    exact_bin_offsets: np.ndarray,
    exact_bin_segment_indices: np.ndarray,
    sweep_bin_edge_origin_xz_um: np.ndarray,
    sweep_bin_size_um: float,
    sweep_bin_shape: np.ndarray,
    sweep_bin_offsets: np.ndarray,
    sweep_bin_segment_indices: np.ndarray,
    tolerance_um: float,
    *,
    use_numba: bool,
) -> tuple[np.ndarray, ...]:
    """Evaluate the needed endpoint distance and sweep in one dispatch."""

    if use_numba and njit is not None:
        kernel = (
            _inspect_swept_continuous_wall_paths_with_end_distance_parallel_numba
            if int(starts_xz_um.shape[0]) >= 64
            else _inspect_swept_continuous_wall_paths_with_end_distance_numba
        )
    else:
        kernel = _inspect_swept_continuous_wall_paths_with_end_distance_kernel
    return kernel(
        starts_xz_um,
        ends_xz_um,
        radii_um,
        start_distance_um,
        segment_start_xz_um,
        segment_end_xz_um,
        segment_inward_normal_xz,
        segment_ring_index,
        int(ring_size),
        segment_arclength_start_um,
        segment_arclength_end_um,
        float(ring_length_um),
        exact_bin_edge_origin_xz_um,
        float(exact_bin_size_um),
        exact_bin_shape,
        exact_bin_offsets,
        exact_bin_segment_indices,
        sweep_bin_edge_origin_xz_um,
        float(sweep_bin_size_um),
        sweep_bin_shape,
        sweep_bin_offsets,
        sweep_bin_segment_indices,
        float(tolerance_um),
    )


def _bin_for_point(x, z, origin, bin_size, shape):
    bx = int(math.floor((x - origin[0]) / bin_size))
    bz = int(math.floor((z - origin[1]) / bin_size))
    bx = min(max(bx, 0), int(shape[0]) - 1)
    bz = min(max(bz, 0), int(shape[1]) - 1)
    return bx, bz


def _point_segment_projection(x, z, ax, az, bx, bz):
    dx = bx - ax
    dz = bz - az
    denominator = dx * dx + dz * dz
    if denominator <= 0.0:
        px, pz = ax, az
    else:
        fraction = ((x - ax) * dx + (z - az) * dz) / denominator
        fraction = min(max(fraction, 0.0), 1.0)
        px = ax + fraction * dx
        pz = az + fraction * dz
    return math.hypot(x - px, z - pz), px, pz


def _ring_adjacent_kernel(first, second, ring_size):
    difference = abs(int(first) - int(second))
    return difference <= 1 or (ring_size > 1 and difference == ring_size - 1)


def _exact_continuous_wall_states_kernel(
    positions,
    starts,
    ends,
    inward,
    ring_indices,
    ring_size,
    bin_origin,
    bin_size,
    bin_shape,
    bin_offsets,
    bin_segments,
):
    count = positions.shape[0]
    distance = np.full(count, math.inf, dtype=np.float64)
    normal = np.zeros((count, 2), dtype=np.float64)
    primary = np.full(count, -1, dtype=np.int64)
    projected = np.full((count, 2), np.nan, dtype=np.float64)
    unique = np.ones(count, dtype=np.bool_)
    for lane in range(count):
        x, z = positions[lane, 0], positions[lane, 1]
        bin_x, bin_z = _bin_for_point(x, z, bin_origin, bin_size, bin_shape)
        flat = bin_x * int(bin_shape[1]) + bin_z
        first = int(bin_offsets[flat])
        last = int(bin_offsets[flat + 1])
        for cursor in range(first, last):
            segment = int(bin_segments[cursor])
            candidate, px, pz = _point_segment_projection(
                x,
                z,
                starts[segment, 0],
                starts[segment, 1],
                ends[segment, 0],
                ends[segment, 1],
            )
            if candidate < distance[lane] or (
                candidate == distance[lane]
                and (primary[lane] < 0 or segment < primary[lane])
            ):
                distance[lane] = candidate
                primary[lane] = segment
                projected[lane, 0] = px
                projected[lane, 1] = pz

        selected = int(primary[lane])
        if selected < 0:
            continue
        normal_x = inward[selected, 0]
        normal_z = inward[selected, 1]
        coordinate_scale = max(abs(x), abs(z), distance[lane], 1.0)
        tie_tolerance = 256.0 * np.finfo(np.float64).eps * coordinate_scale
        # Average normals only at a shared vertex of adjacent tessellation
        # elements.  Equal-distance non-adjacent walls remain distinct.
        for cursor in range(first, last):
            segment = int(bin_segments[cursor])
            if segment == selected:
                continue
            candidate, _, _ = _point_segment_projection(
                x,
                z,
                starts[segment, 0],
                starts[segment, 1],
                ends[segment, 0],
                ends[segment, 1],
            )
            if abs(candidate - distance[lane]) <= tie_tolerance:
                if _ring_adjacent_kernel(
                    ring_indices[selected], ring_indices[segment], ring_size
                ):
                    normal_x += inward[segment, 0]
                    normal_z += inward[segment, 1]
                else:
                    unique[lane] = False
        norm = math.hypot(normal_x, normal_z)
        if norm > 0.0:
            normal[lane, 0] = normal_x / norm
            normal[lane, 1] = normal_z / norm
    return distance, normal, primary, projected, unique


def _exact_continuous_wall_states_parallel_kernel(
    positions,
    starts,
    ends,
    inward,
    ring_indices,
    ring_size,
    bin_origin,
    bin_size,
    bin_shape,
    bin_offsets,
    bin_segments,
):
    """Independent particle-parallel form of the exact-state query.

    Each lane owns all of its output rows and visits candidate segments in the
    same CSR order as the serial reference.  The separate function name is
    intentional: it prevents one Python function from being compiled into two
    dispatchers with different parallel semantics and makes threshold routing
    directly testable.
    """

    count = positions.shape[0]
    distance = np.full(count, math.inf, dtype=np.float64)
    normal = np.zeros((count, 2), dtype=np.float64)
    primary = np.full(count, -1, dtype=np.int64)
    projected = np.full((count, 2), np.nan, dtype=np.float64)
    unique = np.ones(count, dtype=np.bool_)
    for lane in prange(count):
        x, z = positions[lane, 0], positions[lane, 1]
        bin_x, bin_z = _bin_for_point(x, z, bin_origin, bin_size, bin_shape)
        flat = bin_x * int(bin_shape[1]) + bin_z
        first = int(bin_offsets[flat])
        last = int(bin_offsets[flat + 1])
        for cursor in range(first, last):
            segment = int(bin_segments[cursor])
            candidate, px, pz = _point_segment_projection(
                x,
                z,
                starts[segment, 0],
                starts[segment, 1],
                ends[segment, 0],
                ends[segment, 1],
            )
            if candidate < distance[lane] or (
                candidate == distance[lane]
                and (primary[lane] < 0 or segment < primary[lane])
            ):
                distance[lane] = candidate
                primary[lane] = segment
                projected[lane, 0] = px
                projected[lane, 1] = pz

        selected = int(primary[lane])
        if selected < 0:
            continue
        normal_x = inward[selected, 0]
        normal_z = inward[selected, 1]
        coordinate_scale = max(abs(x), abs(z), distance[lane], 1.0)
        tie_tolerance = 256.0 * np.finfo(np.float64).eps * coordinate_scale
        # Preserve the serial cursor order so vertex-normal accumulation and
        # distinct-wall classification remain bitwise deterministic per lane.
        for cursor in range(first, last):
            segment = int(bin_segments[cursor])
            if segment == selected:
                continue
            candidate, _, _ = _point_segment_projection(
                x,
                z,
                starts[segment, 0],
                starts[segment, 1],
                ends[segment, 0],
                ends[segment, 1],
            )
            if abs(candidate - distance[lane]) <= tie_tolerance:
                if _ring_adjacent_kernel(
                    ring_indices[selected], ring_indices[segment], ring_size
                ):
                    normal_x += inward[segment, 0]
                    normal_z += inward[segment, 1]
                else:
                    unique[lane] = False
        norm = math.hypot(normal_x, normal_z)
        if norm > 0.0:
            normal[lane, 0] = normal_x / norm
            normal[lane, 1] = normal_z / norm
    return distance, normal, primary, projected, unique


def _exact_continuous_wall_distances_kernel(
    positions,
    starts,
    ends,
    bin_origin,
    bin_size,
    bin_shape,
    bin_offsets,
    bin_segments,
):
    count = positions.shape[0]
    distance = np.full(count, math.inf, dtype=np.float64)
    primary = np.full(count, -1, dtype=np.int64)
    for lane in range(count):
        x, z = positions[lane, 0], positions[lane, 1]
        bin_x, bin_z = _bin_for_point(x, z, bin_origin, bin_size, bin_shape)
        flat = bin_x * int(bin_shape[1]) + bin_z
        first = int(bin_offsets[flat])
        last = int(bin_offsets[flat + 1])
        for cursor in range(first, last):
            segment = int(bin_segments[cursor])
            candidate, _, _ = _point_segment_projection(
                x,
                z,
                starts[segment, 0],
                starts[segment, 1],
                ends[segment, 0],
                ends[segment, 1],
            )
            if candidate < distance[lane] or (
                candidate == distance[lane]
                and (primary[lane] < 0 or segment < primary[lane])
            ):
                distance[lane] = candidate
                primary[lane] = segment
    return distance


def _cross(ax, az, bx, bz):
    return ax * bz - az * bx


def _segments_intersect(ax, az, bx, bz, cx, cz, dx, dz):
    rx, rz = bx - ax, bz - az
    sx, sz = dx - cx, dz - cz
    first_length_squared = rx * rx + rz * rz
    second_length_squared = sx * sx + sz * sz
    coordinate_scale = max(
        abs(ax),
        abs(az),
        abs(bx),
        abs(bz),
        abs(cx),
        abs(cz),
        abs(dx),
        abs(dz),
        1.0,
    )
    coordinate_tolerance = 128.0 * np.finfo(np.float64).eps * coordinate_scale
    if first_length_squared <= 0.0:
        distance, _, _ = _point_segment_projection(ax, az, cx, cz, dx, dz)
        return distance <= coordinate_tolerance
    if second_length_squared <= 0.0:
        distance, _, _ = _point_segment_projection(cx, cz, ax, az, bx, bz)
        return distance <= coordinate_tolerance

    denominator = _cross(rx, rz, sx, sz)
    denominator_scale = abs(rx * sz) + abs(rz * sx)
    denominator_tolerance = 128.0 * np.finfo(np.float64).eps * denominator_scale
    if abs(denominator) > denominator_tolerance:
        offset_x, offset_z = cx - ax, cz - az
        first = _cross(offset_x, offset_z, sx, sz) / denominator
        second = _cross(offset_x, offset_z, rx, rz) / denominator
        parameter_tolerance = 128.0 * np.finfo(np.float64).eps
        return (
            -parameter_tolerance <= first <= 1.0 + parameter_tolerance
            and -parameter_tolerance <= second <= 1.0 + parameter_tolerance
        )
    offset_x, offset_z = cx - ax, cz - az
    collinearity = _cross(offset_x, offset_z, rx, rz)
    collinearity_scale = abs(offset_x * rz) + abs(offset_z * rx)
    collinearity_tolerance = 128.0 * np.finfo(np.float64).eps * collinearity_scale
    if abs(collinearity) > collinearity_tolerance:
        return False
    if abs(rx) >= abs(rz):
        return max(min(ax, bx), min(cx, dx)) <= (
            min(max(ax, bx), max(cx, dx)) + coordinate_tolerance
        )
    return max(min(az, bz), min(cz, dz)) <= (
        min(max(az, bz), max(cz, dz)) + coordinate_tolerance
    )


def _segment_segment_distance(ax, az, bx, bz, cx, cz, dx, dz):
    first_dx, first_dz = bx - ax, bz - az
    second_dx, second_dz = dx - cx, dz - cz
    first_length_squared = first_dx * first_dx + first_dz * first_dz
    second_length_squared = second_dx * second_dx + second_dz * second_dz
    if first_length_squared <= 0.0:
        distance, _, _ = _point_segment_projection(ax, az, cx, cz, dx, dz)
        return distance
    if second_length_squared <= 0.0:
        distance, _, _ = _point_segment_projection(cx, cz, ax, az, bx, bz)
        return distance
    if _segments_intersect(ax, az, bx, bz, cx, cz, dx, dz):
        return 0.0
    first, _, _ = _point_segment_projection(ax, az, cx, cz, dx, dz)
    second, _, _ = _point_segment_projection(bx, bz, cx, cz, dx, dz)
    third, _, _ = _point_segment_projection(cx, cz, ax, az, bx, bz)
    fourth, _, _ = _point_segment_projection(dx, dz, ax, az, bx, bz)
    return min(first, second, third, fourth)


def _first_capsule_contact_fraction(sx, sz, ex, ez, ax, az, bx, bz, radius, tolerance):
    start_distance, _, _ = _point_segment_projection(sx, sz, ax, az, bx, bz)
    if start_distance <= radius + tolerance:
        return 0.0
    direction_x, direction_z = ex - sx, ez - sz
    quadratic_a = direction_x * direction_x + direction_z * direction_z
    if quadratic_a == 0.0:
        return math.inf
    tangent_x, tangent_z = bx - ax, bz - az
    length = math.hypot(tangent_x, tangent_z)
    if length <= 0.0:
        return math.inf
    tangent_x /= length
    tangent_z /= length
    normal_x, normal_z = -tangent_z, tangent_x
    center_x, center_z = 0.5 * (ax + bx), 0.5 * (az + bz)
    relative_x, relative_z = sx - center_x, sz - center_z
    u0 = relative_x * tangent_x + relative_z * tangent_z
    v0 = relative_x * normal_x + relative_z * normal_z
    du = direction_x * tangent_x + direction_z * tangent_z
    dv = direction_x * normal_x + direction_z * normal_z
    half = 0.5 * length
    best = math.inf
    if dv != 0.0:
        for sign in (-1.0, 1.0):
            fraction = (sign * radius - v0) / dv
            if -tolerance <= fraction <= 1.0 + tolerance:
                clipped = min(max(fraction, 0.0), 1.0)
                u = u0 + clipped * du
                if -half - tolerance <= u <= half + tolerance:
                    best = min(best, clipped)

    for endpoint_x, endpoint_z in ((ax, az), (bx, bz)):
        offset_x, offset_z = sx - endpoint_x, sz - endpoint_z
        quadratic_b = 2.0 * (
            offset_x * direction_x + offset_z * direction_z
        )
        quadratic_c = offset_x * offset_x + offset_z * offset_z - radius * radius
        discriminant = quadratic_b * quadratic_b - 4.0 * quadratic_a * quadratic_c
        discriminant_tolerance = 256.0 * np.finfo(np.float64).eps * max(
            quadratic_b * quadratic_b,
            abs(4.0 * quadratic_a * quadratic_c),
            1.0,
        )
        if discriminant < -discriminant_tolerance:
            continue
        root = math.sqrt(max(discriminant, 0.0))
        first = (-quadratic_b - root) / (2.0 * quadratic_a)
        second = (-quadratic_b + root) / (2.0 * quadratic_a)
        if -tolerance <= first <= 1.0 + tolerance:
            best = min(best, min(max(first, 0.0), 1.0))
        if -tolerance <= second <= 1.0 + tolerance:
            best = min(best, min(max(second, 0.0), 1.0))
    return best


def _inspect_swept_continuous_wall_paths_kernel(
    starts,
    ends,
    radii,
    broad_phase_radii,
    wall_starts,
    wall_ends,
    inward,
    ring_indices,
    ring_size,
    arc_start,
    arc_end,
    ring_length,
    bin_origin,
    bin_size,
    bin_shape,
    bin_offsets,
    bin_segments,
    tolerance,
):
    count = starts.shape[0]
    minimum_gap = np.full(count, math.inf, dtype=np.float64)
    first_fraction = np.full(count, np.nan, dtype=np.float64)
    first_segment = np.full(count, -1, dtype=np.int64)
    multiple = np.zeros(count, dtype=np.bool_)
    fraction_tolerance = 256.0 * np.finfo(np.float64).eps
    for lane in prange(count):
        sx, sz = starts[lane, 0], starts[lane, 1]
        ex, ez = ends[lane, 0], ends[lane, 1]
        radius = radii[lane]
        search_radius = broad_phase_radii[lane]
        first_x, first_z = _bin_for_point(
            min(sx, ex) - search_radius - tolerance,
            min(sz, ez) - search_radius - tolerance,
            bin_origin,
            bin_size,
            bin_shape,
        )
        last_x, last_z = _bin_for_point(
            max(sx, ex) + search_radius + tolerance,
            max(sz, ez) + search_radius + tolerance,
            bin_origin,
            bin_size,
            bin_shape,
        )
        best_distance = math.inf
        best_fraction = math.inf
        selected = -1
        selected_normal_x = 0.0
        selected_normal_z = 0.0
        for bin_x in range(first_x, last_x + 1):
            for bin_z in range(first_z, last_z + 1):
                flat = bin_x * int(bin_shape[1]) + bin_z
                for cursor in range(int(bin_offsets[flat]), int(bin_offsets[flat + 1])):
                    segment = int(bin_segments[cursor])
                    ax, az = wall_starts[segment, 0], wall_starts[segment, 1]
                    bx, bz = wall_ends[segment, 0], wall_ends[segment, 1]
                    candidate_distance = _segment_segment_distance(
                        sx, sz, ex, ez, ax, az, bx, bz
                    )
                    best_distance = min(best_distance, candidate_distance)
                    fraction = _first_capsule_contact_fraction(
                        sx, sz, ex, ez, ax, az, bx, bz, radius, tolerance
                    )
                    if not math.isfinite(fraction):
                        continue
                    if fraction < best_fraction - fraction_tolerance:
                        best_fraction = fraction
                        selected = segment
                        selected_normal_x = inward[segment, 0]
                        selected_normal_z = inward[segment, 1]
                        multiple[lane] = False
                    elif (
                        abs(fraction - best_fraction) <= fraction_tolerance
                        and segment != selected
                    ):
                        first_arc = 0.5 * (arc_start[selected] + arc_end[selected])
                        second_arc = 0.5 * (arc_start[segment] + arc_end[segment])
                        arc_distance = abs(first_arc - second_arc)
                        arc_distance = min(arc_distance, ring_length - arc_distance)
                        same_local_wall = arc_distance <= (
                            2.0 * radius
                            + (arc_end[selected] - arc_start[selected])
                            + (arc_end[segment] - arc_start[segment])
                            + tolerance
                        )
                        if not same_local_wall:
                            dot = (
                                inward[segment, 0] * selected_normal_x
                                + inward[segment, 1] * selected_normal_z
                            )
                            if dot < 1.0 - 1.0e-12:
                                multiple[lane] = True
                        if segment < selected:
                            selected = segment
                            selected_normal_x = inward[segment, 0]
                            selected_normal_z = inward[segment, 1]
        minimum_gap[lane] = best_distance - radius
        if math.isfinite(best_fraction):
            first_fraction[lane] = min(max(best_fraction, 0.0), 1.0)
            first_segment[lane] = selected
    return minimum_gap, first_fraction, first_segment, multiple


def _inspect_swept_continuous_wall_paths_with_end_state_kernel(
    path_starts,
    path_ends,
    radii,
    start_distance,
    wall_starts,
    wall_ends,
    inward,
    ring_indices,
    ring_size,
    arc_start,
    arc_end,
    ring_length,
    exact_bin_origin,
    exact_bin_size,
    exact_bin_shape,
    exact_bin_offsets,
    exact_bin_segments,
    sweep_bin_origin,
    sweep_bin_size,
    sweep_bin_shape,
    sweep_bin_offsets,
    sweep_bin_segments,
    tolerance,
):
    end_distance, end_normal, end_primary, end_projected, _ = (
        _exact_continuous_wall_states_numba(
            path_ends,
            wall_starts,
            wall_ends,
            inward,
            ring_indices,
            ring_size,
            exact_bin_origin,
            exact_bin_size,
            exact_bin_shape,
            exact_bin_offsets,
            exact_bin_segments,
        )
    )
    broad_phase_radius = np.maximum(
        radii, np.minimum(start_distance, end_distance)
    )
    minimum_gap, first_fraction, first_segment, multiple = (
        _inspect_swept_continuous_wall_paths_numba(
            path_starts,
            path_ends,
            radii,
            broad_phase_radius,
            wall_starts,
            wall_ends,
            inward,
            ring_indices,
            ring_size,
            arc_start,
            arc_end,
            ring_length,
            sweep_bin_origin,
            sweep_bin_size,
            sweep_bin_shape,
            sweep_bin_offsets,
            sweep_bin_segments,
            tolerance,
        )
    )
    return (
        end_distance,
        end_normal,
        end_primary,
        end_projected,
        minimum_gap,
        first_fraction,
        first_segment,
        multiple,
    )


def _inspect_swept_continuous_wall_paths_with_end_distance_kernel(
    path_starts,
    path_ends,
    radii,
    start_distance,
    wall_starts,
    wall_ends,
    inward,
    ring_indices,
    ring_size,
    arc_start,
    arc_end,
    ring_length,
    exact_bin_origin,
    exact_bin_size,
    exact_bin_shape,
    exact_bin_offsets,
    exact_bin_segments,
    sweep_bin_origin,
    sweep_bin_size,
    sweep_bin_shape,
    sweep_bin_offsets,
    sweep_bin_segments,
    tolerance,
):
    end_distance = _exact_continuous_wall_distances_numba(
        path_ends,
        wall_starts,
        wall_ends,
        exact_bin_origin,
        exact_bin_size,
        exact_bin_shape,
        exact_bin_offsets,
        exact_bin_segments,
    )
    broad_phase_radius = np.maximum(
        radii, np.minimum(start_distance, end_distance)
    )
    minimum_gap, first_fraction, first_segment, multiple = (
        _inspect_swept_continuous_wall_paths_numba(
            path_starts,
            path_ends,
            radii,
            broad_phase_radius,
            wall_starts,
            wall_ends,
            inward,
            ring_indices,
            ring_size,
            arc_start,
            arc_end,
            ring_length,
            sweep_bin_origin,
            sweep_bin_size,
            sweep_bin_shape,
            sweep_bin_offsets,
            sweep_bin_segments,
            tolerance,
        )
    )
    return (
        end_distance,
        minimum_gap,
        first_fraction,
        first_segment,
        multiple,
    )


def _inspect_swept_continuous_wall_paths_with_end_distance_parallel_kernel(
    path_starts,
    path_ends,
    radii,
    start_distance,
    wall_starts,
    wall_ends,
    inward,
    ring_indices,
    ring_size,
    arc_start,
    arc_end,
    ring_length,
    exact_bin_origin,
    exact_bin_size,
    exact_bin_shape,
    exact_bin_offsets,
    exact_bin_segments,
    sweep_bin_origin,
    sweep_bin_size,
    sweep_bin_shape,
    sweep_bin_offsets,
    sweep_bin_segments,
    tolerance,
):
    end_distance = _exact_continuous_wall_distances_numba(
        path_ends,
        wall_starts,
        wall_ends,
        exact_bin_origin,
        exact_bin_size,
        exact_bin_shape,
        exact_bin_offsets,
        exact_bin_segments,
    )
    broad_phase_radius = np.maximum(
        radii, np.minimum(start_distance, end_distance)
    )
    minimum_gap, first_fraction, first_segment, multiple = (
        _inspect_swept_continuous_wall_paths_parallel_numba(
            path_starts,
            path_ends,
            radii,
            broad_phase_radius,
            wall_starts,
            wall_ends,
            inward,
            ring_indices,
            ring_size,
            arc_start,
            arc_end,
            ring_length,
            sweep_bin_origin,
            sweep_bin_size,
            sweep_bin_shape,
            sweep_bin_offsets,
            sweep_bin_segments,
            tolerance,
        )
    )
    return (
        end_distance,
        minimum_gap,
        first_fraction,
        first_segment,
        multiple,
    )


if njit is not None:
    _bin_for_point = njit(cache=True, nogil=True)(_bin_for_point)
    _point_segment_projection = njit(cache=True, nogil=True)(
        _point_segment_projection
    )
    _ring_adjacent_kernel = njit(cache=True, nogil=True)(_ring_adjacent_kernel)
    _cross = njit(cache=True, nogil=True)(_cross)
    _segments_intersect = njit(cache=True, nogil=True)(_segments_intersect)
    _segment_segment_distance = njit(cache=True, nogil=True)(
        _segment_segment_distance
    )
    _first_capsule_contact_fraction = njit(cache=True, nogil=True)(
        _first_capsule_contact_fraction
    )
    _exact_continuous_wall_states_numba = njit(cache=True, nogil=True)(
        _exact_continuous_wall_states_kernel
    )
    _exact_continuous_wall_states_parallel_numba = njit(
        cache=True, nogil=True, parallel=True
    )(_exact_continuous_wall_states_parallel_kernel)
    _exact_continuous_wall_distances_numba = njit(cache=True, nogil=True)(
        _exact_continuous_wall_distances_kernel
    )
    _inspect_swept_continuous_wall_paths_numba = njit(cache=True, nogil=True)(
        _inspect_swept_continuous_wall_paths_kernel
    )
    _inspect_swept_continuous_wall_paths_parallel_numba = njit(
        cache=True, nogil=True, parallel=True
    )(_inspect_swept_continuous_wall_paths_kernel)
    _inspect_swept_continuous_wall_paths_with_end_state_numba = njit(
        cache=True, nogil=True
    )(_inspect_swept_continuous_wall_paths_with_end_state_kernel)
    _inspect_swept_continuous_wall_paths_with_end_distance_numba = njit(
        cache=True, nogil=True
    )(_inspect_swept_continuous_wall_paths_with_end_distance_kernel)
    _inspect_swept_continuous_wall_paths_with_end_distance_parallel_numba = njit(
        cache=True, nogil=True
    )(_inspect_swept_continuous_wall_paths_with_end_distance_parallel_kernel)
else:  # pragma: no cover
    _exact_continuous_wall_states_numba = _exact_continuous_wall_states_kernel
    _exact_continuous_wall_states_parallel_numba = (
        _exact_continuous_wall_states_parallel_kernel
    )
    _exact_continuous_wall_distances_numba = _exact_continuous_wall_distances_kernel
    _inspect_swept_continuous_wall_paths_numba = (
        _inspect_swept_continuous_wall_paths_kernel
    )
    _inspect_swept_continuous_wall_paths_parallel_numba = (
        _inspect_swept_continuous_wall_paths_kernel
    )
    _inspect_swept_continuous_wall_paths_with_end_state_numba = (
        _inspect_swept_continuous_wall_paths_with_end_state_kernel
    )
    _inspect_swept_continuous_wall_paths_with_end_distance_numba = (
        _inspect_swept_continuous_wall_paths_with_end_distance_kernel
    )
    _inspect_swept_continuous_wall_paths_with_end_distance_parallel_numba = (
        _inspect_swept_continuous_wall_paths_with_end_distance_kernel
    )
