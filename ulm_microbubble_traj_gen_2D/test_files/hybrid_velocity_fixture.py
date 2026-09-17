"""Small affine-triangle hybrid fields used by unit tests."""

from __future__ import annotations

import numpy as np

from ulm_microbubble_traj_gen_2D.utils.core.types import HybridVelocityField
from ulm_microbubble_traj_gen_2D.utils.flow.hybrid_velocity import (
    build_finite_element_velocity_field,
    interpolation_nodes,
)


def rectangular_hybrid_velocity(
    domain,
    velocity_xz=(0.0, 0.0),
    *,
    finite_element_distance_um=3.0,
    regular_grid_distance_um=4.0,
    bounds_xz=None,
):
    if bounds_xz is None:
        x = np.asarray(domain.x_coordinates_um, dtype=np.float64)
        z = np.asarray(domain.z_coordinates_um, dtype=np.float64)
    else:
        xmin, xmax, zmin, zmax = (float(value) for value in bounds_xz)
        spacing = float(domain.spacing_um)
        x = np.arange(xmin, xmax + 0.5 * spacing, spacing)
        z = np.arange(zmin, zmax + 0.5 * spacing, spacing)
    triangles = []
    for i in range(x.size - 1):
        for j in range(z.size - 1):
            lower_left = (x[i], z[j])
            lower_right = (x[i + 1], z[j])
            upper_left = (x[i], z[j + 1])
            upper_right = (x[i + 1], z[j + 1])
            triangles.append((lower_left, lower_right, upper_left))
            triangles.append((upper_right, upper_left, lower_right))
    vertices = np.asarray(triangles, dtype=np.float64)
    node_count = interpolation_nodes(1).shape[0]
    values = np.broadcast_to(
        np.asarray(velocity_xz, dtype=np.float64),
        (vertices.shape[0], node_count, 2),
    ).copy()
    finite_element = build_finite_element_velocity_field(
        vertices,
        values,
        1,
        preferred_bin_size_um=float(domain.spacing_um),
    )
    return HybridVelocityField(
        finite_element=finite_element,
        finite_element_distance_um=float(finite_element_distance_um),
        regular_grid_distance_um=float(regular_grid_distance_um),
    )
