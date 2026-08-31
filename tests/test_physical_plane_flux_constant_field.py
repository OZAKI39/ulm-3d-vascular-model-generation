from __future__ import annotations

import numpy as np
import pytest

from utils.cfd_flow.tau1_grid_convergence import (
    PhysicalPortPlane,
    build_plane_quadrature,
    integrate_plane_flux,
)


def test_constant_velocity_equals_normal_velocity_times_aperture_area() -> None:
    plane = PhysicalPortPlane(
        "outlet_01",
        np.zeros(3),
        np.array([0.0, 0.0, 1.0]),
        np.array([1.0, 0.0, 0.0]),
        np.array([0.0, 1.0, 0.0]),
        np.array([[-0.25, -0.25], [0.25, -0.25], [0.25, 0.25], [-0.25, 0.25]]),
        "synthetic",
    )
    centers = np.zeros((1, 3))
    quadrature = build_plane_quadrature(centers, dx_m=1.0, plane=plane)
    result = integrate_plane_flux(np.array([[0.0, 0.0, 2.0]]), plane=plane, quadrature=quadrature)
    assert result["area_coverage_relative_error"] < 1e-14
    assert result["physical_q_m3_s"] == pytest.approx(0.5)
