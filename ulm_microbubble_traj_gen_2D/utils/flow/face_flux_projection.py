"""Face-flux diagnostics for Cartesian samples of DOLFINx flow."""

from __future__ import annotations

import numpy as np

from ..core.types import GridDomain
from .flow_boundaries import BoundaryFluxFields


def cell_velocity_to_face_flux(
    velocity_xz_um_s: np.ndarray,
    lumen_mask: np.ndarray,
    spacing_um: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Convert cell-centered velocity to internal face fluxes with zero wall flux."""

    velocity = np.asarray(velocity_xz_um_s, dtype=float)
    lumen = np.asarray(lumen_mask, dtype=bool)
    nx, nz = lumen.shape
    spacing = float(spacing_um)
    flux_x = np.zeros((nx + 1, nz), dtype=np.float64)
    flux_z = np.zeros((nx, nz + 1), dtype=np.float64)
    if nx > 1:
        internal_x = lumen[:-1, :] & lumen[1:, :]
        flux_x[1:nx, :][internal_x] = (
            0.5
            * (
                velocity[:-1, :, 0][internal_x]
                + velocity[1:, :, 0][internal_x]
            )
            * spacing
        )
    if nz > 1:
        internal_z = lumen[:, :-1] & lumen[:, 1:]
        flux_z[:, 1:nz][internal_z] = (
            0.5
            * (
                velocity[:, :-1, 1][internal_z]
                + velocity[:, 1:, 1][internal_z]
            )
            * spacing
        )
    return flux_x, flux_z


def net_outflow_from_face_flux(
    flux_x_um2_s: np.ndarray,
    flux_z_um2_s: np.ndarray,
    lumen_mask: np.ndarray,
) -> np.ndarray:
    """Return per-cell net outflow from the face flux arrays."""

    lumen = np.asarray(lumen_mask, dtype=bool)
    flux_x = np.asarray(flux_x_um2_s, dtype=float)
    flux_z = np.asarray(flux_z_um2_s, dtype=float)
    net = flux_x[1:, :] - flux_x[:-1, :] + flux_z[:, 1:] - flux_z[:, :-1]
    net = np.asarray(net, dtype=np.float64)
    net[~lumen] = 0.0
    return net


def divergence_from_face_flux(
    flux_x_um2_s: np.ndarray,
    flux_z_um2_s: np.ndarray,
    domain: GridDomain,
    lumen_mask: np.ndarray,
) -> np.ndarray:
    """Return finite-volume divergence of the sampled DOLFINx velocity."""

    spacing = max(float(domain.spacing_um), np.finfo(float).eps)
    return net_outflow_from_face_flux(
        flux_x_um2_s, flux_z_um2_s, lumen_mask
    ) / (spacing * spacing)


def normalized_divergence_error_from_flux(
    velocity_xz_um_s: np.ndarray,
    divergence_s_inv: np.ndarray,
    domain: GridDomain,
    mask: np.ndarray,
) -> float:
    """Normalize finite-volume divergence by the sampled velocity scale."""

    selected = np.asarray(mask, dtype=bool)
    if not np.any(selected):
        return float("inf")
    divergence = np.asarray(divergence_s_inv, dtype=float)
    div_norm = float(np.linalg.norm(divergence[selected]))
    if div_norm <= np.finfo(float).eps:
        return 0.0
    velocity = np.asarray(velocity_xz_um_s, dtype=float)
    speed_norm = float(np.linalg.norm(np.linalg.norm(velocity[selected], axis=-1)))
    spacing = max(float(domain.spacing_um), np.finfo(float).eps)
    return div_norm / max(speed_norm / spacing, np.finfo(float).eps)


def solid_wall_penetration_field(
    flux_x_um2_s: np.ndarray,
    flux_z_um2_s: np.ndarray,
    lumen_mask: np.ndarray,
    boundaries: BoundaryFluxFields,
    spacing_um: float,
) -> np.ndarray:
    """Measure sampled normal velocity on solid wall faces."""

    lumen = np.asarray(lumen_mask, dtype=bool)
    spacing = max(float(spacing_um), np.finfo(float).eps)
    open_x, open_z = _open_face_masks(
        boundaries, flux_x_um2_s.shape, flux_z_um2_s.shape
    )
    wall = np.zeros(lumen.shape, dtype=np.float64)
    flux_x = np.asarray(flux_x_um2_s, dtype=float)
    flux_z = np.asarray(flux_z_um2_s, dtype=float)
    nx, nz = lumen.shape
    for i, j in np.argwhere(lumen):
        values = []
        if (i == 0 or not lumen[i - 1, j]) and not open_x[i, j]:
            values.append(abs(float(flux_x[i, j])) / spacing)
        if (i == nx - 1 or not lumen[i + 1, j]) and not open_x[i + 1, j]:
            values.append(abs(float(flux_x[i + 1, j])) / spacing)
        if (j == 0 or not lumen[i, j - 1]) and not open_z[i, j]:
            values.append(abs(float(flux_z[i, j])) / spacing)
        if (j == nz - 1 or not lumen[i, j + 1]) and not open_z[i, j + 1]:
            values.append(abs(float(flux_z[i, j + 1])) / spacing)
        if values:
            wall[i, j] = max(values)
    return wall.astype(np.float32)


def _open_face_masks(
    boundaries: BoundaryFluxFields,
    flux_x_shape: tuple[int, int],
    flux_z_shape: tuple[int, int],
) -> tuple[np.ndarray, np.ndarray]:
    open_x = np.zeros(flux_x_shape, dtype=bool)
    open_z = np.zeros(flux_z_shape, dtype=bool)
    for row in range(boundaries.open_face_axis.size):
        i, j = (int(x) for x in boundaries.open_face_index_ij[row])
        if int(boundaries.open_face_axis[row]) == 0:
            open_x[i, j] = True
        else:
            open_z[i, j] = True
    return open_x, open_z
