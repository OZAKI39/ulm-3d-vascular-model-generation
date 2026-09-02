"""Grid-domain construction."""

import numpy as np
from ..core.types import GridDomain


def build_domain_from_vessels(vessels, cfg):
    """
    Build an X-Z grid covering all vessels.
    """

    if not vessels:
        raise ValueError("Cannot build a grid domain without vessels.")

    spacing       = cfg.grid_spacing_um
    padding       = cfg.padding_um
    max_grid_size = cfg.max_grid_cells
    if not np.isfinite(spacing) or spacing <= 0:
        raise ValueError("domain.grid_spacing_um must be finite and positive.")
    if not np.isfinite(padding) or padding < 0:
        raise ValueError("domain.padding_um must be finite and non-negative.")
    if max_grid_size <= 0:
        raise ValueError("domain.max_grid_cells must be positive.")

    points = np.asarray(
        [point for vessel in vessels for point in (vessel.x_p, vessel.x_d)],
        dtype=float,
    )
    if not np.isfinite(points).all():
        raise ValueError("Vessel endpoints must contain finite coordinates.")

    min_xyz = points.min(axis=0) - padding
    max_xyz = points.max(axis=0) + padding

    nx = int(np.ceil((max_xyz[0] - min_xyz[0]) / spacing)) + 1
    nz = int(np.ceil((max_xyz[2] - min_xyz[2]) / spacing)) + 1
    if nx <= 1 or nz <= 1:
        raise ValueError("The generated X-Z grid is too small for field simulation.")
    if nx * nz > max_grid_size:
        raise ValueError(
            f"The requested grid has {nx * nz} cells, exceeding "
            f"domain.max_grid_cells={max_grid_size}. "
            "Increase grid_spacing_um or max_grid_cells."
        )

    x_coordinates = min_xyz[0] + spacing * np.arange(nx)
    z_coordinates = min_xyz[2] + spacing * np.arange(nz)
    fixed_y        = points[:, 1].mean()

    return GridDomain(
        origin_um=np.asarray([min_xyz[0], fixed_y, min_xyz[2]]),
        spacing_um=spacing,
        shape=(nx, nz),
        fixed_y_um=fixed_y,
        x_coordinates_um=x_coordinates,
        z_coordinates_um=z_coordinates,
    )
