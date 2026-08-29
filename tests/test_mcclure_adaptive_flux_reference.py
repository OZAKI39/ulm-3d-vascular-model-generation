from __future__ import annotations

import numpy as np
import pytest

from utils.cfd_flow.mcclure_adaptive_flux_reference import (
    integrated_positive_z_flux,
    physical_volume_flux_to_lattice,
    reconstruct_positive_z_pressure_boundary,
    solve_uniform_boundary_density,
)


@pytest.mark.parametrize(
    ("node_count", "target", "perturbation"),
    [
        (1, 2.5e-4, 0.0),
        (7, -1.7e-3, 2.0e-5),
        (41, 3.1e-3, 7.0e-5),
    ],
)
def test_eq9_and_appendix_b_reconstruct_total_flux(
    node_count: int, target: float, perturbation: float
) -> None:
    rng = np.random.default_rng(8931 + node_count)
    pdfs = np.full((node_count, 19), 1.0 / 19.0)
    if perturbation:
        pdfs += rng.uniform(-perturbation, perturbation, size=pdfs.shape)
    rho = solve_uniform_boundary_density(pdfs, target)
    reconstructed = reconstruct_positive_z_pressure_boundary(pdfs, rho)
    actual = integrated_positive_z_flux(reconstructed)
    assert abs(actual - target) / abs(target) <= 1.0e-12


def test_physical_flux_conversion_is_total_flux_and_has_no_area_argument() -> None:
    q_phys = 7.693508475538942e-16
    dx = 2.0e-7
    dt = 2.44140625e-8
    expected = q_phys * dt / dx**3
    assert physical_volume_flux_to_lattice(q_phys, dx_m=dx, dt_s=dt) == pytest.approx(
        expected, rel=0.0, abs=1.0e-18
    )

