from __future__ import annotations

import inspect

import numpy as np
import pytest

from utils.cfd_flow.adaptive_flux_pressure_audit import (
    musubi_pressure_flux_affine_coefficients,
    musubi_pressure_flux_lattice,
    physical_mass_factor,
    solve_boundary_density,
)
from utils.cfd_flow.mcclure_adaptive_flux_reference import (
    physical_volume_flux_to_lattice,
)


def _synthetic_boundary() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rng = np.random.default_rng(28114)
    stored = np.full((9, 19), 1.0 / 19.0) + rng.uniform(-1e-4, 1e-4, (9, 19))
    masks = np.zeros((9, 18), dtype=bool)
    masks[:, [2, 6, 8, 10, 11]] = True
    masks[::2, 0] = True
    velocity = rng.uniform(-8e-4, 8e-4, (9, 3))
    return stored, masks, velocity


def test_musubi_native_pressure_flux_is_affine_and_solve_is_exact() -> None:
    stored, masks, velocity = _synthetic_boundary()
    alpha, beta = musubi_pressure_flux_affine_coefficients(
        stored_boundary_pdfs=stored,
        incoming_masks=masks,
        extrapolated_velocity=velocity,
    )
    rho1, rho2 = 0.91, 1.08
    f1 = musubi_pressure_flux_lattice(
        rho1,
        stored_boundary_pdfs=stored,
        incoming_masks=masks,
        extrapolated_velocity=velocity,
    )
    f2 = musubi_pressure_flux_lattice(
        rho2,
        stored_boundary_pdfs=stored,
        incoming_masks=masks,
        extrapolated_velocity=velocity,
    )
    assert (f2 - f1) / (rho2 - rho1) == pytest.approx(alpha, rel=1e-14)
    assert f1 - alpha * rho1 == pytest.approx(beta, abs=1e-14)
    target = 2.7e-3
    rho = solve_boundary_density(target, alpha, beta)
    actual = musubi_pressure_flux_lattice(
        rho,
        stored_boundary_pdfs=stored,
        incoming_masks=masks,
        extrapolated_velocity=velocity,
    )
    assert abs(actual - target) / target <= 1.0e-12


def test_signed_flux_and_physical_conversion_are_consistent() -> None:
    dx, dt, rho0 = 2.0e-7, 2.44140625e-8, 1056.0
    q_phys = 7.693508475538942e-16
    q_lat = physical_volume_flux_to_lattice(q_phys, dx_m=dx, dt_s=dt)
    mdot = q_lat * physical_mass_factor(density_kg_m3=rho0, dx_m=dx, dt_s=dt)
    assert mdot == pytest.approx(rho0 * q_phys, rel=1e-15)


def test_adaptive_solve_has_no_mfr_area_proxy_dependency() -> None:
    source = inspect.getsource(solve_boundary_density)
    assert "area" not in source.lower()

