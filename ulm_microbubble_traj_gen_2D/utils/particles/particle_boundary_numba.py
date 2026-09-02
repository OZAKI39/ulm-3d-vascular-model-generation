"""Compiled kernels for authoritative particle boundary queries.

The public geometry object owns validation and physical semantics.  These
kernels only replace allocation-heavy NumPy expressions inside the internal
particle loop.  Their serial lane/outlet order deliberately preserves the
reference implementation's deterministic first-index tie break.
"""

from __future__ import annotations

import math

import numpy as np

try:
    from numba import njit
except ImportError:  # pragma: no cover - used only in minimal environments
    njit = None


NUMBA_AVAILABLE = njit is not None


def first_outlet_crossing_arrays(
    starts_xz_um: np.ndarray,
    ends_xz_um: np.ndarray,
    section_points_xz_um: np.ndarray,
    section_normals_xz: np.ndarray,
    section_tangents_xz: np.ndarray,
    section_half_width_um: np.ndarray,
    bin_edge_origin_xz_um: np.ndarray,
    bin_size_um: float,
    bin_shape: np.ndarray,
    bin_offsets: np.ndarray,
    bin_section_indices: np.ndarray,
    tolerance_um: float,
    *,
    use_numba: bool,
) -> tuple[np.ndarray, ...]:
    """Return compact earliest-crossing arrays for a prevalidated batch."""

    kernel = (
        _first_outlet_crossing_arrays_numba
        if use_numba and njit is not None
        else _first_outlet_crossing_arrays_kernel
    )
    return kernel(
        starts_xz_um,
        ends_xz_um,
        section_points_xz_um,
        section_normals_xz,
        section_tangents_xz,
        section_half_width_um,
        bin_edge_origin_xz_um,
        float(bin_size_um),
        bin_shape,
        bin_offsets,
        bin_section_indices,
        float(tolerance_um),
    )


def directed_open_section_crossing_mask(
    starts_xz_um: np.ndarray,
    ends_xz_um: np.ndarray,
    section_points_xz_um: np.ndarray,
    section_normals_xz: np.ndarray,
    section_tangents_xz: np.ndarray,
    section_half_width_um: np.ndarray,
    tolerance_um: float,
    *,
    use_numba: bool,
) -> np.ndarray:
    """Return lanes crossing any supplied open section from inside to outside.

    This dense kernel is intended for the one (or very few) inlet caps.  Outlet
    queries keep their CSR spatial index.  Checking the inlet explicitly lets a
    solid-wall-certified chord prove that its endpoint is still in the lumen,
    without invoking a Shapely point-in-polygon query on every internal step.
    """

    kernel = (
        _directed_open_section_crossing_mask_numba
        if use_numba and njit is not None
        else _directed_open_section_crossing_mask_kernel
    )
    return kernel(
        starts_xz_um,
        ends_xz_um,
        section_points_xz_um,
        section_normals_xz,
        section_tangents_xz,
        section_half_width_um,
        float(tolerance_um),
    )


def first_directed_open_section_crossing_arrays(
    starts_xz_um: np.ndarray,
    ends_xz_um: np.ndarray,
    section_points_xz_um: np.ndarray,
    section_normals_xz: np.ndarray,
    section_tangents_xz: np.ndarray,
    section_half_width_um: np.ndarray,
    tolerance_um: float,
    *,
    use_numba: bool,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return the earliest directed crossing of a small open-section set.

    Unlike the outlet CSR query this dense query is intended for the normally
    one-element inlet catalog.  Returning the crossing fraction, rather than
    only a mask, lets callers compare a reverse-inlet escape chronologically
    with an earlier wall, outlet, or topology event on the same chord.
    """

    kernel = (
        _first_directed_open_section_crossing_arrays_numba
        if use_numba and njit is not None
        else _first_directed_open_section_crossing_arrays_kernel
    )
    return kernel(
        starts_xz_um,
        ends_xz_um,
        section_points_xz_um,
        section_normals_xz,
        section_tangents_xz,
        section_half_width_um,
        float(tolerance_um),
    )


def _first_outlet_crossing_arrays_kernel(
    starts,
    ends,
    points,
    normals,
    tangents,
    widths,
    bin_origin,
    bin_size,
    bin_shape,
    bin_offsets,
    bin_section_indices,
    tolerance,
):
    particle_count = starts.shape[0]
    fractions = np.full(particle_count, np.nan, dtype=np.float64)
    outlet_indices = np.full(particle_count, -1, dtype=np.int32)
    crossing_positions = np.full((particle_count, 2), np.nan, dtype=np.float64)
    signed_start = np.full(particle_count, np.nan, dtype=np.float64)
    signed_end = np.full(particle_count, np.nan, dtype=np.float64)

    for lane in range(particle_count):
        start_x = starts[lane, 0]
        start_z = starts[lane, 1]
        delta_x = ends[lane, 0] - start_x
        delta_z = ends[lane, 1] - start_z
        best_fraction = math.inf
        best_index = -1
        best_start = math.nan
        best_end = math.nan
        path_min_x = min(start_x, ends[lane, 0]) - tolerance
        path_max_x = max(start_x, ends[lane, 0]) + tolerance
        path_min_z = min(start_z, ends[lane, 1]) - tolerance
        path_max_z = max(start_z, ends[lane, 1]) + tolerance
        min_bin_x = int(math.floor((path_min_x - bin_origin[0]) / bin_size))
        max_bin_x = int(math.floor((path_max_x - bin_origin[0]) / bin_size))
        min_bin_z = int(math.floor((path_min_z - bin_origin[1]) / bin_size))
        max_bin_z = int(math.floor((path_max_z - bin_origin[1]) / bin_size))
        if (
            max_bin_x >= 0
            and max_bin_z >= 0
            and min_bin_x < bin_shape[0]
            and min_bin_z < bin_shape[1]
        ):
            min_bin_x = max(min_bin_x, 0)
            max_bin_x = min(max_bin_x, bin_shape[0] - 1)
            min_bin_z = max(min_bin_z, 0)
            max_bin_z = min(max_bin_z, bin_shape[1] - 1)
            for bin_x in range(min_bin_x, max_bin_x + 1):
                for bin_z in range(min_bin_z, max_bin_z + 1):
                    bin_index = bin_x * bin_shape[1] + bin_z
                    first = bin_offsets[bin_index]
                    last = bin_offsets[bin_index + 1]
                    for candidate in range(first, last):
                        outlet = bin_section_indices[candidate]
                        normal_x = normals[outlet, 0]
                        normal_z = normals[outlet, 1]
                        psi_start = (
                            (points[outlet, 0] - start_x) * normal_x
                            + (points[outlet, 1] - start_z) * normal_z
                        )
                        psi_end = (
                            (points[outlet, 0] - ends[lane, 0]) * normal_x
                            + (points[outlet, 1] - ends[lane, 1]) * normal_z
                        )
                        effective_start = (
                            0.0 if abs(psi_start) <= tolerance else psi_start
                        )
                        effective_end = (
                            0.0 if abs(psi_end) <= tolerance else psi_end
                        )
                        denominator = effective_start - effective_end
                        if (
                            psi_start < -tolerance
                            or psi_end > tolerance
                            or psi_end >= psi_start - tolerance
                            or denominator <= 0.0
                        ):
                            continue
                        fraction = effective_start / denominator
                        if fraction < -tolerance or fraction > 1.0 + tolerance:
                            continue
                        fraction = min(max(fraction, 0.0), 1.0)
                        crossing_x = start_x + fraction * delta_x
                        crossing_z = start_z + fraction * delta_z
                        lateral = abs(
                            (crossing_x - points[outlet, 0]) * tangents[outlet, 0]
                            + (crossing_z - points[outlet, 1]) * tangents[outlet, 1]
                        )
                        if lateral > widths[outlet] + tolerance:
                            continue
                        if (
                            fraction < best_fraction
                            or (
                                fraction == best_fraction
                                and (best_index < 0 or outlet < best_index)
                            )
                        ):
                            best_fraction = fraction
                            best_index = outlet
                            best_start = psi_start
                            best_end = psi_end

        if best_index >= 0:
            fractions[lane] = best_fraction
            outlet_indices[lane] = best_index
            crossing_positions[lane, 0] = start_x + best_fraction * delta_x
            crossing_positions[lane, 1] = start_z + best_fraction * delta_z
            signed_start[lane] = best_start
            signed_end[lane] = best_end

    return (
        fractions,
        outlet_indices,
        crossing_positions,
        signed_start,
        signed_end,
    )


def _directed_open_section_crossing_mask_kernel(
    starts,
    ends,
    points,
    normals,
    tangents,
    widths,
    tolerance,
):
    particle_count = starts.shape[0]
    section_count = points.shape[0]
    crossed = np.zeros(particle_count, dtype=np.bool_)
    for lane in range(particle_count):
        start_x = starts[lane, 0]
        start_z = starts[lane, 1]
        end_x = ends[lane, 0]
        end_z = ends[lane, 1]
        delta_x = end_x - start_x
        delta_z = end_z - start_z
        for section in range(section_count):
            normal_x = normals[section, 0]
            normal_z = normals[section, 1]
            psi_start = (
                (points[section, 0] - start_x) * normal_x
                + (points[section, 1] - start_z) * normal_z
            )
            psi_end = (
                (points[section, 0] - end_x) * normal_x
                + (points[section, 1] - end_z) * normal_z
            )
            effective_start = 0.0 if abs(psi_start) <= tolerance else psi_start
            effective_end = 0.0 if abs(psi_end) <= tolerance else psi_end
            denominator = effective_start - effective_end
            if (
                psi_start < -tolerance
                or psi_end > tolerance
                or psi_end >= psi_start - tolerance
                or denominator <= 0.0
            ):
                continue
            fraction = effective_start / denominator
            if fraction < -tolerance or fraction > 1.0 + tolerance:
                continue
            fraction = min(max(fraction, 0.0), 1.0)
            crossing_x = start_x + fraction * delta_x
            crossing_z = start_z + fraction * delta_z
            lateral = abs(
                (crossing_x - points[section, 0]) * tangents[section, 0]
                + (crossing_z - points[section, 1]) * tangents[section, 1]
            )
            if lateral <= widths[section] + tolerance:
                crossed[lane] = True
                break
    return crossed


def _first_directed_open_section_crossing_arrays_kernel(
    starts,
    ends,
    points,
    normals,
    tangents,
    widths,
    tolerance,
):
    particle_count = starts.shape[0]
    section_count = points.shape[0]
    fractions = np.full(particle_count, np.nan, dtype=np.float64)
    section_indices = np.full(particle_count, -1, dtype=np.int32)
    positions = np.full((particle_count, 2), np.nan, dtype=np.float64)
    for lane in range(particle_count):
        start_x = starts[lane, 0]
        start_z = starts[lane, 1]
        end_x = ends[lane, 0]
        end_z = ends[lane, 1]
        delta_x = end_x - start_x
        delta_z = end_z - start_z
        best_fraction = math.inf
        best_section = -1
        for section in range(section_count):
            normal_x = normals[section, 0]
            normal_z = normals[section, 1]
            psi_start = (
                (points[section, 0] - start_x) * normal_x
                + (points[section, 1] - start_z) * normal_z
            )
            psi_end = (
                (points[section, 0] - end_x) * normal_x
                + (points[section, 1] - end_z) * normal_z
            )
            effective_start = 0.0 if abs(psi_start) <= tolerance else psi_start
            effective_end = 0.0 if abs(psi_end) <= tolerance else psi_end
            denominator = effective_start - effective_end
            if (
                psi_start < -tolerance
                or psi_end > tolerance
                or psi_end >= psi_start - tolerance
                or denominator <= 0.0
            ):
                continue
            fraction = effective_start / denominator
            if fraction < -tolerance or fraction > 1.0 + tolerance:
                continue
            fraction = min(max(fraction, 0.0), 1.0)
            crossing_x = start_x + fraction * delta_x
            crossing_z = start_z + fraction * delta_z
            lateral = abs(
                (crossing_x - points[section, 0]) * tangents[section, 0]
                + (crossing_z - points[section, 1]) * tangents[section, 1]
            )
            if lateral > widths[section] + tolerance:
                continue
            if (
                fraction < best_fraction
                or (
                    fraction == best_fraction
                    and (best_section < 0 or section < best_section)
                )
            ):
                best_fraction = fraction
                best_section = section
        if best_section >= 0:
            fractions[lane] = best_fraction
            section_indices[lane] = best_section
            positions[lane, 0] = start_x + best_fraction * delta_x
            positions[lane, 1] = start_z + best_fraction * delta_z
    return fractions, section_indices, positions


if njit is not None:
    _first_outlet_crossing_arrays_numba = njit(cache=True, nogil=True)(
        _first_outlet_crossing_arrays_kernel
    )
    _directed_open_section_crossing_mask_numba = njit(cache=True, nogil=True)(
        _directed_open_section_crossing_mask_kernel
    )
    _first_directed_open_section_crossing_arrays_numba = njit(
        cache=True, nogil=True
    )(_first_directed_open_section_crossing_arrays_kernel)
else:  # pragma: no cover - exercised only without Numba
    _first_outlet_crossing_arrays_numba = _first_outlet_crossing_arrays_kernel
    _directed_open_section_crossing_mask_numba = (
        _directed_open_section_crossing_mask_kernel
    )
    _first_directed_open_section_crossing_arrays_numba = (
        _first_directed_open_section_crossing_arrays_kernel
    )
