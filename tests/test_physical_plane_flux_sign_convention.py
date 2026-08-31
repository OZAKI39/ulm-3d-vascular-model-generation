from __future__ import annotations

import numpy as np
import pytest

from utils.cfd_flow.tau1_grid_convergence import PhysicalPortPlane, build_plane_quadrature, integrate_plane_flux


def test_zero_velocity_and_normal_reversal() -> None:
    aperture = np.array([[-0.1, -0.1], [0.1, -0.1], [0.1, 0.1], [-0.1, 0.1]])
    forward = PhysicalPortPlane("outlet_01", np.zeros(3), np.array([0.0, 0.0, 1.0]), np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0]), aperture, "a")
    reverse = PhysicalPortPlane("outlet_01", np.zeros(3), np.array([0.0, 0.0, -1.0]), np.array([1.0, 0.0, 0.0]), np.array([0.0, -1.0, 0.0]), aperture, "b")
    centers = np.zeros((1, 3))
    velocity = np.array([[0.0, 0.0, 3.0]])
    qf = integrate_plane_flux(velocity, plane=forward, quadrature=build_plane_quadrature(centers, dx_m=1.0, plane=forward))["physical_q_m3_s"]
    qr = integrate_plane_flux(velocity, plane=reverse, quadrature=build_plane_quadrature(centers, dx_m=1.0, plane=reverse))["physical_q_m3_s"]
    q0 = integrate_plane_flux(np.zeros((1, 3)), plane=forward, quadrature=build_plane_quadrature(centers, dx_m=1.0, plane=forward))["physical_q_m3_s"]
    assert qr == pytest.approx(-qf)
    assert q0 == pytest.approx(0.0)
