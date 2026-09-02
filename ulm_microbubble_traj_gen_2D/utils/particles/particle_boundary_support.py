"""Grid overlays sampled from the continuous vessel boundary for indexing."""

from __future__ import annotations

import math

import numpy as np
from scipy import ndimage

from ..geometry.continuous_vessel_geometry import ContinuousVesselGeometry
from ..core.types import GridDomain


def sample_continuous_boundary_support_to_grid(
    domain: GridDomain,
    geometry: ContinuousVesselGeometry,
) -> tuple[np.ndarray, np.ndarray]:
    """Return solid-wall and open-section grid overlays for indexing only."""

    shape = tuple(int(value) for value in domain.shape)
    if geometry.shape != shape:
        raise ValueError("Continuous vessel geometry must match domain.shape.")
    count = shape[0] * shape[1]
    inside = np.empty(count, dtype=bool)
    chunk = 100_000
    for start in range(0, count, chunk):
        stop = min(start + chunk, count)
        flat = np.arange(start, stop, dtype=np.int64)
        ix = flat // shape[1]
        iz = flat - ix * shape[1]
        points = np.column_stack(
            (domain.x_coordinates_um[ix], domain.z_coordinates_um[iz])
        )
        inside[start:stop] = np.asarray(
            geometry.contains_xz_um(points), dtype=bool
        )

    inside_grid = inside.reshape(shape)
    open_grid = np.zeros(shape, dtype=bool)
    half_cell = 0.75 * float(domain.spacing_um)
    for point, normal, tangent, width in zip(
        geometry.open_section_point_xz_um,
        geometry.open_section_outward_normal_xz,
        geometry.open_section_tangent_xz,
        geometry.open_section_half_width_um,
    ):
        reach = float(width) + 2.0 * float(domain.spacing_um)
        minimum = point - reach
        maximum = point + reach
        ix0 = max(
            0,
            int(
                math.floor(
                    (minimum[0] - domain.origin_um[0]) / domain.spacing_um
                )
            ),
        )
        ix1 = min(
            shape[0],
            int(
                math.ceil(
                    (maximum[0] - domain.origin_um[0]) / domain.spacing_um
                )
            )
            + 1,
        )
        iz0 = max(
            0,
            int(
                math.floor(
                    (minimum[1] - domain.origin_um[2]) / domain.spacing_um
                )
            ),
        )
        iz1 = min(
            shape[1],
            int(
                math.ceil(
                    (maximum[1] - domain.origin_um[2]) / domain.spacing_um
                )
            )
            + 1,
        )
        if ix0 >= ix1 or iz0 >= iz1:
            continue
        x, z = np.meshgrid(
            domain.x_coordinates_um[ix0:ix1],
            domain.z_coordinates_um[iz0:iz1],
            indexing="ij",
        )
        relative_x = x - float(point[0])
        relative_z = z - float(point[1])
        normal_coordinate = (
            relative_x * float(normal[0]) + relative_z * float(normal[1])
        )
        tangent_coordinate = (
            relative_x * float(tangent[0]) + relative_z * float(tangent[1])
        )
        open_grid[ix0:ix1, iz0:iz1] |= (
            np.abs(normal_coordinate) <= half_cell
        ) & (np.abs(tangent_coordinate) <= float(width) + half_cell)

    boundary_sites = inside_grid & ~ndimage.binary_erosion(
        inside_grid,
        structure=ndimage.generate_binary_structure(2, 1),
        border_value=0,
    )
    solid_sites = boundary_sites & ~open_grid
    if not np.any(solid_sites):
        raise ValueError(
            "The sampled continuous lumen contains no solid boundary cells."
        )
    return (
        np.ascontiguousarray(solid_sites),
        np.ascontiguousarray(open_grid),
    )
