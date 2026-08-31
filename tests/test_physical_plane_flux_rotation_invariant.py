from __future__ import annotations

import numpy as np
import pytest

from utils.cfd_flow.port_grid_sensitivity import orthonormal_plane_basis
from utils.cfd_flow.tau1_grid_convergence import PhysicalPortPlane, build_plane_quadrature, integrate_plane_flux


def _result(normal: np.ndarray) -> float:
    u, v = orthonormal_plane_basis(normal)
    plane = PhysicalPortPlane(
        "outlet_01", np.zeros(3), normal, u, v,
        np.array([[-0.1, -0.1], [0.1, -0.1], [0.1, 0.1], [-0.1, 0.1]]),
        "synthetic",
    )
    quadrature = build_plane_quadrature(np.zeros((1, 3)), dx_m=1.0, plane=plane)
    return float(integrate_plane_flux((2.0 * normal)[None, :], plane=plane, quadrature=quadrature)["physical_q_m3_s"])


def test_rotated_plane_preserves_constant_normal_flux() -> None:
    axis = np.array([0.0, 0.0, 1.0])
    rotated = np.ones(3) / np.sqrt(3.0)
    assert _result(rotated) == pytest.approx(_result(axis), rel=1e-13)
