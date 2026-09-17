"""Deterministic short-range microbubble collision forces."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

try:
    from numba import njit
    from numba.extending import register_jitable
except ImportError:  # pragma: no cover - exercised only in deliberately minimal environments
    njit = None

    def register_jitable(function):
        return function


@dataclass(frozen=True)
class CollisionResult:
    """Per-particle forces and aggregate overlap diagnostics for one state snapshot."""

    force_xz_pn: np.ndarray
    neighbor_count: np.ndarray
    maximum_physical_overlap_um: float
    maximum_collision_compression_um: float
    interacting_pair_count: int
    search_strategy: str


_AUTO_CELL_LIST_THRESHOLD = 1024


def compute_collision_forces(
    positions_xz_um: np.ndarray,
    radii_um: np.ndarray,
    translational_mobility_xz_um_pn_s: np.ndarray,
    active: np.ndarray,
    bubble_ids: np.ndarray,
    *,
    collision_layer_um: float,
    relaxation_time_s: float,
    strategy: str = "all_pairs",
    use_numba: bool = True,
) -> CollisionResult:
    """Compute equal-and-opposite central repulsion from an immutable particle snapshot."""

    positions = np.ascontiguousarray(positions_xz_um, dtype=np.float64)
    radii = np.ascontiguousarray(radii_um, dtype=np.float64)
    mobility = np.ascontiguousarray(translational_mobility_xz_um_pn_s, dtype=np.float64)
    active_mask = np.ascontiguousarray(active, dtype=np.bool_)
    ids = np.ascontiguousarray(bubble_ids, dtype=np.int64)
    n = int(positions.shape[0])
    if positions.shape != (n, 2):
        raise ValueError("positions_xz_um must have shape (n, 2).")
    if radii.shape != (n,) or active_mask.shape != (n,) or ids.shape != (n,):
        raise ValueError("radii_um, active, and bubble_ids must each have shape (n,).")
    if mobility.shape != (n, 2, 2):
        raise ValueError("translational_mobility_xz_um_pn_s must have shape (n, 2, 2).")
    if float(collision_layer_um) < 0.0:
        raise ValueError("collision_layer_um must be non-negative.")
    if float(relaxation_time_s) <= 0.0:
        raise ValueError("relaxation_time_s must be positive.")
    requested_strategy = str(strategy).lower()
    resolved_strategy = resolve_collision_search_strategy(
        requested_strategy,
        int(np.count_nonzero(active_mask)),
    )

    if resolved_strategy == "cell_list":
        kernel = (
            _collision_cell_list_numba
            if use_numba and njit is not None
            else _collision_cell_list_kernel
        )
    else:
        kernel = _collision_numba if use_numba and njit is not None else _collision_kernel
    force, neighbors, max_overlap, max_compression, pair_count, invalid_mobility = kernel(
        positions,
        radii,
        mobility,
        active_mask,
        ids,
        float(collision_layer_um),
        float(relaxation_time_s),
    )
    if invalid_mobility:
        raise ValueError("Collision relative mobility must be finite and strictly positive.")
    return CollisionResult(
        force_xz_pn=np.asarray(force, dtype=np.float64),
        neighbor_count=np.asarray(neighbors, dtype=np.int32),
        maximum_physical_overlap_um=float(max_overlap),
        maximum_collision_compression_um=float(max_compression),
        interacting_pair_count=int(pair_count),
        search_strategy=resolved_strategy,
    )


def resolve_collision_search_strategy(strategy: str, active_particle_count: int) -> str:
    """Resolve ``auto`` without changing the collision force law.

    The threshold is an implementation choice, not a physical parameter.  Small
    populations keep the compact all-pairs loop; larger populations use the
    sparse cell list so distant particles are never presented to the force law.
    """

    value = str(strategy).lower()
    if value == "auto":
        return (
            "cell_list"
            if int(active_particle_count) >= _AUTO_CELL_LIST_THRESHOLD
            else "all_pairs"
        )
    if value in {"all_pairs", "cell_list"}:
        return value
    raise ValueError("Collision search must be one of: all_pairs, cell_list, auto.")


def compute_collision_forces_prevalidated_numba(
    positions_xz_um: np.ndarray,
    radii_um: np.ndarray,
    translational_mobility_xz_um_pn_s: np.ndarray,
    active: np.ndarray,
    bubble_ids: np.ndarray,
    *,
    collision_layer_um: float,
    relaxation_time_s: float,
    strategy: str,
) -> tuple[np.ndarray, np.ndarray, float, float, int, str]:
    """Run the checked production arrays without repeating public API guards."""

    resolved_strategy = resolve_collision_search_strategy(
        strategy,
        int(np.count_nonzero(active)),
    )
    kernel = (
        _collision_cell_list_numba
        if resolved_strategy == "cell_list"
        else _collision_numba
    )
    result = kernel(
        positions_xz_um,
        radii_um,
        translational_mobility_xz_um_pn_s,
        active,
        bubble_ids,
        float(collision_layer_um),
        float(relaxation_time_s),
    )
    if bool(result[5]):
        raise ValueError(
            "Collision relative mobility must be finite and strictly positive."
        )
    return (
        result[0],
        result[1],
        float(result[2]),
        float(result[3]),
        int(result[4]),
        resolved_strategy,
    )


def _collision_kernel(
    positions,
    radii,
    mobility,
    active,
    bubble_ids,
    collision_layer,
    relaxation_time,
):
    n_particles = positions.shape[0]
    force = np.zeros((n_particles, 2), dtype=np.float64)
    neighbor_count = np.zeros(n_particles, dtype=np.int32)
    maximum_physical_overlap = 0.0
    maximum_compression = 0.0
    pair_count = 0
    invalid_mobility = False
    for i in range(n_particles - 1):
        if not active[i]:
            continue
        for j in range(i + 1, n_particles):
            if not active[j]:
                continue
            dx = positions[i, 0] - positions[j, 0]
            dz = positions[i, 1] - positions[j, 1]
            distance_squared = dx * dx + dz * dz
            distance = np.sqrt(max(distance_squared, 0.0))
            physical_overlap = radii[i] + radii[j] - distance
            compression = physical_overlap + collision_layer
            if physical_overlap > maximum_physical_overlap:
                maximum_physical_overlap = physical_overlap
            if compression <= 0.0:
                continue
            if compression > maximum_compression:
                maximum_compression = compression
            if distance > 1.0e-14:
                nx = dx / distance
                nz = dz / distance
            else:
                low_id = min(bubble_ids[i], bubble_ids[j])
                high_id = max(bubble_ids[i], bubble_ids[j])
                code = (low_id * 73856093 + high_id * 19349663) % 104729
                angle = 2.0 * np.pi * (code + 0.5) / 104729.0
                nx = np.cos(angle)
                nz = np.sin(angle)
                if bubble_ids[i] != low_id:
                    nx = -nx
                    nz = -nz
            relative_mobility = (
                nx * (mobility[i, 0, 0] + mobility[j, 0, 0]) * nx
                + nx * (mobility[i, 0, 1] + mobility[j, 0, 1]) * nz
                + nz * (mobility[i, 1, 0] + mobility[j, 1, 0]) * nx
                + nz * (mobility[i, 1, 1] + mobility[j, 1, 1]) * nz
            )
            if not np.isfinite(relative_mobility) or relative_mobility <= 0.0:
                invalid_mobility = True
                continue
            magnitude = compression / (relative_mobility * relaxation_time)
            fx = magnitude * nx
            fz = magnitude * nz
            force[i, 0] += fx
            force[i, 1] += fz
            force[j, 0] -= fx
            force[j, 1] -= fz
            neighbor_count[i] += 1
            neighbor_count[j] += 1
            pair_count += 1
    return (
        force,
        neighbor_count,
        maximum_physical_overlap,
        maximum_compression,
        pair_count,
        invalid_mobility,
    )


@register_jitable
def _lower_bound(sorted_values, target):
    low = 0
    high = sorted_values.shape[0]
    while low < high:
        middle = (low + high) // 2
        if sorted_values[middle] < target:
            low = middle + 1
        else:
            high = middle
    return low


@register_jitable
def _upper_bound(sorted_values, target):
    low = 0
    high = sorted_values.shape[0]
    while low < high:
        middle = (low + high) // 2
        if sorted_values[middle] <= target:
            low = middle + 1
        else:
            high = middle
    return low


def _collision_cell_list_kernel(
    positions,
    radii,
    mobility,
    active,
    bubble_ids,
    collision_layer,
    relaxation_time,
):
    """Evaluate the same pair law after a deterministic sparse spatial sort.

    Each particle accumulates only its own force.  A physical pair is therefore
    visited once from each endpoint, which avoids shared writes and keeps this
    kernel safe for a future parallel particle loop.  Equal-and-opposite forces
    still use the same canonical permanent-ID direction for coincident centres.
    """

    n_particles = positions.shape[0]
    force = np.zeros((n_particles, 2), dtype=np.float64)
    neighbor_count = np.zeros(n_particles, dtype=np.int32)
    maximum_overlap_by_particle = np.zeros(n_particles, dtype=np.float64)
    maximum_compression_by_particle = np.zeros(n_particles, dtype=np.float64)
    invalid_by_particle = np.zeros(n_particles, dtype=np.uint8)
    cell_x = np.zeros(n_particles, dtype=np.int64)
    cell_z = np.zeros(n_particles, dtype=np.int64)

    active_count = 0
    maximum_radius = 0.0
    minimum_x = 0.0
    minimum_z = 0.0
    first_active = True
    for index in range(n_particles):
        if not active[index]:
            continue
        active_count += 1
        if radii[index] > maximum_radius:
            maximum_radius = radii[index]
        if first_active:
            minimum_x = positions[index, 0]
            minimum_z = positions[index, 1]
            first_active = False
        else:
            minimum_x = min(minimum_x, positions[index, 0])
            minimum_z = min(minimum_z, positions[index, 1])

    if active_count < 2:
        return force, neighbor_count, 0.0, 0.0, 0, False

    cell_size = max(2.0 * maximum_radius + collision_layer, 1.0e-12)
    maximum_cell_z = 0
    for index in range(n_particles):
        if not active[index]:
            continue
        cell_x[index] = int(np.floor((positions[index, 0] - minimum_x) / cell_size))
        cell_z[index] = int(np.floor((positions[index, 1] - minimum_z) / cell_size))
        if cell_z[index] > maximum_cell_z:
            maximum_cell_z = cell_z[index]

    stride_z = maximum_cell_z + 3
    inactive_key = np.iinfo(np.int64).max
    keys = np.full(n_particles, inactive_key, dtype=np.int64)
    for index in range(n_particles):
        if active[index]:
            keys[index] = cell_x[index] * stride_z + cell_z[index]

    order = np.argsort(keys)
    sorted_keys = keys[order]

    # Keep permanent IDs ordered inside each cell.  This makes force summation
    # reproducible even when a caller presents active lanes in another order.
    cell_start = 0
    while cell_start < active_count:
        cell_end = cell_start + 1
        while cell_end < active_count and sorted_keys[cell_end] == sorted_keys[cell_start]:
            cell_end += 1
        for current in range(cell_start + 1, cell_end):
            saved_index = order[current]
            saved_id = bubble_ids[saved_index]
            insertion = current - 1
            while insertion >= cell_start and bubble_ids[order[insertion]] > saved_id:
                order[insertion + 1] = order[insertion]
                insertion -= 1
            order[insertion + 1] = saved_index
        cell_start = cell_end

    # The cell keys do not change when indices inside one equal-key group are
    # reordered, so ``sorted_keys`` remains valid for binary range lookup.
    for i in range(n_particles):
        if not active[i]:
            continue
        for offset_x in range(-1, 2):
            neighbour_cell_x = cell_x[i] + offset_x
            if neighbour_cell_x < 0:
                continue
            for offset_z in range(-1, 2):
                neighbour_cell_z = cell_z[i] + offset_z
                if neighbour_cell_z < 0:
                    continue
                neighbour_key = neighbour_cell_x * stride_z + neighbour_cell_z
                begin = _lower_bound(sorted_keys, neighbour_key)
                finish = _upper_bound(sorted_keys, neighbour_key)
                for ordered_position in range(begin, finish):
                    j = order[ordered_position]
                    if j == i or not active[j]:
                        continue
                    dx = positions[i, 0] - positions[j, 0]
                    dz = positions[i, 1] - positions[j, 1]
                    distance_squared = dx * dx + dz * dz
                    distance = np.sqrt(max(distance_squared, 0.0))
                    physical_overlap = radii[i] + radii[j] - distance
                    compression = physical_overlap + collision_layer
                    if physical_overlap > maximum_overlap_by_particle[i]:
                        maximum_overlap_by_particle[i] = physical_overlap
                    if compression <= 0.0:
                        continue
                    if compression > maximum_compression_by_particle[i]:
                        maximum_compression_by_particle[i] = compression
                    if distance > 1.0e-14:
                        nx = dx / distance
                        nz = dz / distance
                    else:
                        low_id = min(bubble_ids[i], bubble_ids[j])
                        high_id = max(bubble_ids[i], bubble_ids[j])
                        code = (low_id * 73856093 + high_id * 19349663) % 104729
                        angle = 2.0 * np.pi * (code + 0.5) / 104729.0
                        nx = np.cos(angle)
                        nz = np.sin(angle)
                        if bubble_ids[i] != low_id:
                            nx = -nx
                            nz = -nz
                    relative_mobility = (
                        nx * (mobility[i, 0, 0] + mobility[j, 0, 0]) * nx
                        + nx * (mobility[i, 0, 1] + mobility[j, 0, 1]) * nz
                        + nz * (mobility[i, 1, 0] + mobility[j, 1, 0]) * nx
                        + nz * (mobility[i, 1, 1] + mobility[j, 1, 1]) * nz
                    )
                    if not np.isfinite(relative_mobility) or relative_mobility <= 0.0:
                        invalid_by_particle[i] = 1
                        continue
                    magnitude = compression / (relative_mobility * relaxation_time)
                    force[i, 0] += magnitude * nx
                    force[i, 1] += magnitude * nz
                    neighbor_count[i] += 1

    maximum_overlap = 0.0
    maximum_compression = 0.0
    invalid_mobility = False
    neighbour_sum = 0
    for index in range(n_particles):
        maximum_overlap = max(maximum_overlap, maximum_overlap_by_particle[index])
        maximum_compression = max(
            maximum_compression,
            maximum_compression_by_particle[index],
        )
        neighbour_sum += neighbor_count[index]
        invalid_mobility = invalid_mobility or invalid_by_particle[index] != 0
    return (
        force,
        neighbor_count,
        maximum_overlap,
        maximum_compression,
        neighbour_sum // 2,
        invalid_mobility,
    )


if njit is not None:
    _collision_numba = njit(cache=True)(_collision_kernel)
    _collision_cell_list_numba = njit(cache=True)(_collision_cell_list_kernel)
else:  # pragma: no cover - the project environment includes Numba
    _collision_numba = _collision_kernel
    _collision_cell_list_numba = _collision_cell_list_kernel
