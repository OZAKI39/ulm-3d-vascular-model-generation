"""Clean-room numerical reference for the McClure D3Q19 flux boundary.

This module implements the equations from McClure et al. (2020), Eq. (9)
and Appendix B, independently of the GPL LBPM implementation.  Its direction
ordering is stated explicitly because it is a mathematical reference, not the
ordering used by Musubi.
"""

from __future__ import annotations

import numpy as np


# McClure/LBPM D3Q19 ordering used only by this reference implementation.
LBPM_D3Q19_DIRECTIONS = np.asarray(
    [
        (0, 0, 0),
        (1, 0, 0),
        (-1, 0, 0),
        (0, 1, 0),
        (0, -1, 0),
        (0, 0, 1),
        (0, 0, -1),
        (1, 1, 0),
        (-1, -1, 0),
        (1, -1, 0),
        (-1, 1, 0),
        (1, 0, 1),
        (-1, 0, -1),
        (1, 0, -1),
        (-1, 0, 1),
        (0, 1, 1),
        (0, -1, -1),
        (0, 1, -1),
        (0, -1, 1),
    ],
    dtype=np.float64,
)


def _boundary_rows(pdfs: np.ndarray) -> tuple[np.ndarray, bool]:
    values = np.asarray(pdfs, dtype=np.float64)
    scalar = values.ndim == 1
    if scalar:
        values = values.reshape(1, -1)
    if values.ndim != 2 or values.shape[1] != 19:
        raise ValueError("pdfs must have shape (19,) or (n, 19)")
    if not np.all(np.isfinite(values)):
        raise ValueError("pdfs must be finite")
    return values, scalar


def d3q19_consistency_k(pdfs: np.ndarray) -> np.ndarray | float:
    """Evaluate the known-population term in McClure Eq. (9)."""

    values, scalar = _boundary_rows(pdfs)
    result = (
        values[:, 0]
        + values[:, 1]
        + values[:, 2]
        + values[:, 3]
        + values[:, 4]
        + values[:, 7]
        + values[:, 8]
        + values[:, 9]
        + values[:, 10]
        + 2.0
        * (
            values[:, 6]
            + values[:, 12]
            + values[:, 13]
            + values[:, 16]
            + values[:, 17]
        )
    )
    return float(result[0]) if scalar else result


def solve_uniform_boundary_density(
    pdfs: np.ndarray,
    target_total_flux_lattice: float,
    *,
    reference_density: float = 1.0,
) -> float:
    """Solve McClure Eq. (9) for one density shared by all inlet nodes."""

    values, _ = _boundary_rows(pdfs)
    if reference_density <= 0.0 or not np.isfinite(reference_density):
        raise ValueError("reference_density must be finite and positive")
    target = float(target_total_flux_lattice)
    if not np.isfinite(target):
        raise ValueError("target_total_flux_lattice must be finite")
    return float(
        reference_density * target / len(values)
        + np.mean(d3q19_consistency_k(values))
    )


def reconstruct_positive_z_pressure_boundary(
    pdfs: np.ndarray,
    boundary_density: float,
    *,
    reference_density: float = 1.0,
    tangential_velocity_x: float = 0.0,
    tangential_velocity_y: float = 0.0,
) -> np.ndarray:
    """Apply the Appendix-B D3Q19 pressure closure at a positive-z inlet."""

    values, scalar = _boundary_rows(pdfs)
    rho = float(boundary_density)
    rho0 = float(reference_density)
    ux = float(tangential_velocity_x)
    uy = float(tangential_velocity_y)
    if not all(np.isfinite(item) for item in (rho, rho0, ux, uy)) or rho0 <= 0:
        raise ValueError("density and tangential velocities must be finite")

    reconstructed = values.copy()
    uz = (rho - np.asarray(d3q19_consistency_k(values))) / rho0
    nzx = (
        0.5
        * (
            values[:, 1]
            + values[:, 7]
            + values[:, 9]
            - values[:, 2]
            - values[:, 10]
            - values[:, 8]
        )
        - rho0 * ux / 3.0
    )
    nzy = (
        0.5
        * (
            values[:, 3]
            + values[:, 7]
            + values[:, 10]
            - values[:, 4]
            - values[:, 9]
            - values[:, 8]
        )
        - rho0 * uy / 3.0
    )

    reconstructed[:, 5] = values[:, 6] + rho0 * uz / 3.0
    reconstructed[:, 11] = (
        values[:, 12] + rho0 * (ux + uz) / 6.0 - nzx
    )
    reconstructed[:, 14] = (
        values[:, 13] + rho0 * (-ux + uz) / 6.0 + nzx
    )
    reconstructed[:, 15] = (
        values[:, 16] + rho0 * (uy + uz) / 6.0 - nzy
    )
    reconstructed[:, 18] = (
        values[:, 17] + rho0 * (-uy + uz) / 6.0 + nzy
    )
    return reconstructed[0] if scalar else reconstructed


def integrated_positive_z_flux(
    pdfs: np.ndarray, *, reference_density: float = 1.0
) -> float:
    """Integrate the signed z-volume flux represented by reconstructed PDFs."""

    values, _ = _boundary_rows(pdfs)
    return float(
        np.sum(values @ LBPM_D3Q19_DIRECTIONS[:, 2]) / reference_density
    )


def physical_volume_flux_to_lattice(
    volume_flux_m3_s: float, *, dx_m: float, dt_s: float
) -> float:
    """Convert total physical volume flux directly to total lattice flux."""

    if dx_m <= 0.0 or dt_s <= 0.0:
        raise ValueError("dx_m and dt_s must be positive")
    return float(volume_flux_m3_s) * float(dt_s) / float(dx_m) ** 3

