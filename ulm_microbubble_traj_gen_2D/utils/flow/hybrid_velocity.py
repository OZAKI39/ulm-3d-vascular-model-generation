"""Hybrid Cartesian/finite-element velocity evaluation.

The Cartesian field is a fast far-field cache.  Close to the authoritative
continuous wall, velocity and its gradient come from a cell-local polynomial
that exactly reconstructs the affine-triangle DOLFINx finite-element function.
The two values are blended smoothly in the transition band.
"""

from __future__ import annotations

import math

import numpy as np

from ..core.types import FiniteElementVelocityField, HybridVelocityField

try:
    from numba import njit
except ImportError:  # pragma: no cover
    njit = None


FINITE_ELEMENT_REGION = np.uint8(1)
TRANSITION_REGION = np.uint8(2)
REGULAR_GRID_REGION = np.uint8(3)


def polynomial_exponents(degree: int) -> np.ndarray:
    """Return ``r**a * s**b`` powers for a complete triangle polynomial."""

    if int(degree) < 1:
        raise ValueError("Finite-element velocity degree must be positive.")
    powers = []
    for total in range(int(degree) + 1):
        for a in range(total, -1, -1):
            powers.append((a, total - a))
    return np.ascontiguousarray(powers, dtype=np.int16)


def interpolation_nodes(degree: int) -> np.ndarray:
    """Return a unisolvent equispaced node set on the reference triangle."""

    degree = int(degree)
    if degree < 1:
        raise ValueError("Finite-element velocity degree must be positive.")
    nodes = []
    for i in range(degree + 1):
        for j in range(degree + 1 - i):
            nodes.append((i / degree, j / degree))
    return np.ascontiguousarray(nodes, dtype=np.float64)


def fit_velocity_polynomials(
    sampled_velocity_um_s: np.ndarray,
    degree: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Fit reference-triangle polynomial coefficients for every cell."""

    values = np.asarray(sampled_velocity_um_s, dtype=np.float64)
    powers = polynomial_exponents(degree)
    nodes = interpolation_nodes(degree)
    expected = (nodes.shape[0], 2)
    if values.ndim != 3 or values.shape[1:] != expected:
        raise ValueError(
            "sampled_velocity_um_s must have shape "
            f"(cell_count, {expected[0]}, 2)."
        )
    vandermonde = np.empty((nodes.shape[0], powers.shape[0]), dtype=np.float64)
    for column, (a, b) in enumerate(powers):
        vandermonde[:, column] = nodes[:, 0] ** int(a) * nodes[:, 1] ** int(b)
    coefficients = np.linalg.solve(vandermonde, values)
    return powers, np.ascontiguousarray(coefficients, dtype=np.float64)


def build_finite_element_velocity_field(
    cell_vertices_xz_um: np.ndarray,
    sampled_velocity_um_s: np.ndarray,
    degree: int,
    *,
    preferred_bin_size_um: float,
) -> FiniteElementVelocityField:
    """Build the serializable finite-element evaluator and its spatial bins."""

    vertices = np.ascontiguousarray(cell_vertices_xz_um, dtype=np.float64)
    if vertices.ndim != 3 or vertices.shape[1:] != (3, 2):
        raise ValueError("cell_vertices_xz_um must have shape (cell_count, 3, 2).")
    if vertices.shape[0] == 0:
        raise ValueError("The finite-element velocity field must contain cells.")
    determinants = (
        (vertices[:, 1, 0] - vertices[:, 0, 0])
        * (vertices[:, 2, 1] - vertices[:, 0, 1])
        - (vertices[:, 2, 0] - vertices[:, 0, 0])
        * (vertices[:, 1, 1] - vertices[:, 0, 1])
    )
    if np.any(~np.isfinite(determinants)) or np.any(np.abs(determinants) <= 1.0e-14):
        raise ValueError("Finite-element mesh contains a degenerate triangle.")

    powers, coefficients = fit_velocity_polynomials(
        sampled_velocity_um_s, int(degree)
    )
    bin_size = float(preferred_bin_size_um)
    if not math.isfinite(bin_size) or bin_size <= 0.0:
        raise ValueError("preferred_bin_size_um must be finite and positive.")
    origin, shape, offsets, cell_indices = _build_cell_bins(vertices, bin_size)
    return FiniteElementVelocityField(
        degree=int(degree),
        cell_vertices_xz_um=vertices,
        polynomial_exponents=powers,
        velocity_coefficients_um_s=coefficients,
        bin_origin_xz_um=origin,
        bin_size_um=bin_size,
        bin_shape=shape,
        bin_offsets=offsets,
        bin_cell_indices=cell_indices,
    )


def _build_cell_bins(
    vertices: np.ndarray,
    bin_size_um: float,
) -> tuple[np.ndarray, tuple[int, int], np.ndarray, np.ndarray]:
    minimum = np.min(vertices, axis=(0, 1))
    maximum = np.max(vertices, axis=(0, 1))
    origin = np.floor(minimum / bin_size_um) * bin_size_um
    shape_array = np.maximum(
        np.floor((maximum - origin) / bin_size_um).astype(np.int64) + 1,
        1,
    )
    shape = (int(shape_array[0]), int(shape_array[1]))
    bins = [[] for _ in range(shape[0] * shape[1])]
    for cell, triangle in enumerate(vertices):
        low = np.floor((np.min(triangle, axis=0) - origin) / bin_size_um).astype(
            np.int64
        )
        high = np.floor((np.max(triangle, axis=0) - origin) / bin_size_um).astype(
            np.int64
        )
        low = np.maximum(low, 0)
        high = np.minimum(high, shape_array - 1)
        for i in range(int(low[0]), int(high[0]) + 1):
            for j in range(int(low[1]), int(high[1]) + 1):
                bins[i * shape[1] + j].append(cell)
    offsets = np.zeros(len(bins) + 1, dtype=np.int64)
    for index, members in enumerate(bins):
        offsets[index + 1] = offsets[index] + len(members)
    cell_indices = np.empty(int(offsets[-1]), dtype=np.int32)
    for index, members in enumerate(bins):
        cell_indices[offsets[index] : offsets[index + 1]] = members
    return (
        np.ascontiguousarray(origin, dtype=np.float64),
        shape,
        offsets,
        cell_indices,
    )


def validate_hybrid_velocity_field(field: HybridVelocityField) -> None:
    """Reject incomplete fields instead of falling back to grid-only sampling."""

    near = float(field.finite_element_distance_um)
    far = float(field.regular_grid_distance_um)
    if not math.isfinite(near) or near <= 0.0:
        raise ValueError("finite_element_distance_um must be finite and positive.")
    if not math.isfinite(far) or far <= near:
        raise ValueError(
            "regular_grid_distance_um must be finite and greater than "
            "finite_element_distance_um."
        )
    finite_element = field.finite_element
    vertices = np.asarray(finite_element.cell_vertices_xz_um)
    powers = np.asarray(finite_element.polynomial_exponents)
    coefficients = np.asarray(finite_element.velocity_coefficients_um_s)
    if vertices.ndim != 3 or vertices.shape[1:] != (3, 2):
        raise ValueError("Finite-element cell vertices have an invalid shape.")
    if powers.ndim != 2 or powers.shape[1] != 2:
        raise ValueError("Finite-element polynomial exponents have an invalid shape.")
    if coefficients.shape != (vertices.shape[0], powers.shape[0], 2):
        raise ValueError("Finite-element velocity coefficients have an invalid shape.")
    expected_terms = (int(finite_element.degree) + 1) * (
        int(finite_element.degree) + 2
    ) // 2
    if powers.shape[0] != expected_terms:
        raise ValueError(
            "Finite-element polynomial term count does not match its degree."
        )
    bin_origin = np.asarray(finite_element.bin_origin_xz_um)
    offsets = np.asarray(finite_element.bin_offsets)
    bin_cells = np.asarray(finite_element.bin_cell_indices)
    bin_count = int(finite_element.bin_shape[0]) * int(
        finite_element.bin_shape[1]
    )
    if (
        bin_origin.shape != (2,)
        or min(finite_element.bin_shape) < 1
        or offsets.shape != (bin_count + 1,)
        or offsets[0] != 0
        or np.any(np.diff(offsets) < 0)
        or offsets[-1] != bin_cells.size
        or np.any((bin_cells < 0) | (bin_cells >= vertices.shape[0]))
    ):
        raise ValueError("Finite-element cell-bin index is inconsistent.")
    if not (
        np.isfinite(vertices).all()
        and np.isfinite(coefficients).all()
        and np.isfinite(bin_origin).all()
        and np.isfinite(float(finite_element.bin_size_um))
        and float(finite_element.bin_size_um) > 0.0
    ):
        raise ValueError("Finite-element velocity field contains non-finite values.")


def finite_element_weight(
    wall_distance_um: np.ndarray,
    field: HybridVelocityField,
) -> np.ndarray:
    """Return one near-wall weight per point."""

    distance = np.asarray(wall_distance_um, dtype=np.float64)
    near = float(field.finite_element_distance_um)
    far = float(field.regular_grid_distance_um)
    coordinate = np.clip((far - distance) / (far - near), 0.0, 1.0)
    return coordinate * coordinate * (3.0 - 2.0 * coordinate)


def hybrid_region_map(
    wall_distance_um: np.ndarray,
    lumen_mask: np.ndarray,
    field: HybridVelocityField,
) -> tuple[np.ndarray, np.ndarray]:
    """Build categorical regions and smooth finite-element weights on a grid."""

    distance = np.asarray(wall_distance_um, dtype=np.float64)
    lumen = np.asarray(lumen_mask, dtype=bool)
    if distance.shape != lumen.shape:
        raise ValueError("wall_distance_um and lumen_mask must have the same shape.")
    weight = finite_element_weight(distance, field)
    region = np.zeros(lumen.shape, dtype=np.uint8)
    region[lumen & (distance <= field.finite_element_distance_um)] = (
        FINITE_ELEMENT_REGION
    )
    region[
        lumen
        & (distance > field.finite_element_distance_um)
        & (distance < field.regular_grid_distance_um)
    ] = TRANSITION_REGION
    region[lumen & (distance >= field.regular_grid_distance_um)] = (
        REGULAR_GRID_REGION
    )
    weight[~lumen] = np.nan
    return region, np.ascontiguousarray(weight, dtype=np.float32)


def sample_finite_element_velocity(
    field: FiniteElementVelocityField,
    points_xz_um: np.ndarray,
    required: np.ndarray,
    *,
    use_numba: bool,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Evaluate velocity and physical gradient in containing triangles."""

    points = np.ascontiguousarray(points_xz_um, dtype=np.float64)
    required_mask = np.ascontiguousarray(required, dtype=np.bool_)
    if points.ndim != 2 or points.shape[1] != 2:
        raise ValueError("points_xz_um must have shape (N, 2).")
    if required_mask.shape != (points.shape[0],):
        raise ValueError("required must contain one Boolean per point.")
    kernel = (
        _sample_finite_element_velocity_numba
        if use_numba and njit is not None
        else _sample_finite_element_velocity_kernel
    )
    velocity, gradient, cells = kernel(
        points,
        required_mask,
        field.cell_vertices_xz_um,
        field.polynomial_exponents,
        field.velocity_coefficients_um_s,
        field.bin_origin_xz_um,
        float(field.bin_size_um),
        int(field.bin_shape[0]),
        int(field.bin_shape[1]),
        field.bin_offsets,
        field.bin_cell_indices,
    )
    missing = required_mask & (cells < 0)
    if np.any(missing):
        first = int(np.flatnonzero(missing)[0])
        raise ValueError(
            "A point requiring finite-element velocity is outside the DOLFINx "
            f"mesh: index={first}, xz_um={points[first].tolist()}."
        )
    return velocity, gradient, cells


def blend_velocity_and_gradient(
    grid_velocity: np.ndarray,
    grid_gradient: np.ndarray,
    finite_element_velocity: np.ndarray,
    finite_element_gradient: np.ndarray,
    wall_distance_um: np.ndarray,
    inward_normal_xz: np.ndarray,
    field: HybridVelocityField,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Blend values and include the spatial derivative of the blend weight."""

    distance = np.asarray(wall_distance_um, dtype=np.float64)
    normal = np.asarray(inward_normal_xz, dtype=np.float64)
    weight = finite_element_weight(distance, field)
    velocity_difference = finite_element_velocity - grid_velocity
    velocity = grid_velocity + weight[:, None] * velocity_difference
    gradient = grid_gradient + weight[:, None, None] * (
        finite_element_gradient - grid_gradient
    )

    near = float(field.finite_element_distance_um)
    far = float(field.regular_grid_distance_um)
    coordinate = np.clip((far - distance) / (far - near), 0.0, 1.0)
    derivative_distance = (
        -6.0 * coordinate * (1.0 - coordinate) / (far - near)
    )
    weight_gradient = derivative_distance[:, None] * normal
    gradient += velocity_difference[:, :, None] * weight_gradient[:, None, :]
    return velocity, gradient, weight


def _sample_finite_element_velocity_kernel(
    points,
    required,
    vertices,
    powers,
    coefficients,
    bin_origin,
    bin_size,
    bins_x,
    bins_z,
    bin_offsets,
    bin_cells,
):
    count = points.shape[0]
    velocity = np.zeros((count, 2), dtype=np.float64)
    gradient = np.zeros((count, 2, 2), dtype=np.float64)
    located_cells = np.full(count, -1, dtype=np.int32)
    tolerance = 1.0e-10

    for lane in range(count):
        if not required[lane]:
            continue
        x = points[lane, 0]
        z = points[lane, 1]
        i = int(math.floor((x - bin_origin[0]) / bin_size))
        j = int(math.floor((z - bin_origin[1]) / bin_size))
        if i < 0 or i >= bins_x or j < 0 or j >= bins_z:
            continue
        flat_bin = i * bins_z + j
        for cursor in range(bin_offsets[flat_bin], bin_offsets[flat_bin + 1]):
            cell = bin_cells[cursor]
            x0 = vertices[cell, 0, 0]
            z0 = vertices[cell, 0, 1]
            e1x = vertices[cell, 1, 0] - x0
            e1z = vertices[cell, 1, 1] - z0
            e2x = vertices[cell, 2, 0] - x0
            e2z = vertices[cell, 2, 1] - z0
            determinant = e1x * e2z - e2x * e1z
            dx = x - x0
            dz = z - z0
            r = (e2z * dx - e2x * dz) / determinant
            s = (-e1z * dx + e1x * dz) / determinant
            if r < -tolerance or s < -tolerance or r + s > 1.0 + tolerance:
                continue

            dr_dx = e2z / determinant
            dr_dz = -e2x / determinant
            ds_dx = -e1z / determinant
            ds_dz = e1x / determinant
            for term in range(powers.shape[0]):
                a = int(powers[term, 0])
                b = int(powers[term, 1])
                basis = r**a * s**b
                derivative_r = 0.0 if a == 0 else a * r ** (a - 1) * s**b
                derivative_s = 0.0 if b == 0 else b * r**a * s ** (b - 1)
                for component in range(2):
                    coefficient = coefficients[cell, term, component]
                    velocity[lane, component] += coefficient * basis
                    gradient[lane, component, 0] += coefficient * (
                        derivative_r * dr_dx + derivative_s * ds_dx
                    )
                    gradient[lane, component, 1] += coefficient * (
                        derivative_r * dr_dz + derivative_s * ds_dz
                    )
            located_cells[lane] = cell
            break
    return velocity, gradient, located_cells


if njit is not None:
    _sample_finite_element_velocity_numba = njit(cache=True)(
        _sample_finite_element_velocity_kernel
    )
else:  # pragma: no cover
    _sample_finite_element_velocity_numba = _sample_finite_element_velocity_kernel
