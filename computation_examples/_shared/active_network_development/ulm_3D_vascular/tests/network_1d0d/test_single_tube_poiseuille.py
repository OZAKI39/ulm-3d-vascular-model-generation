"""Permanent analytic validation for steady Newtonian networks."""
import numpy as np
import pytest
from network_1d0d.hydraulic_resistance import linear_radius_resistance
from network_1d0d.network_solver import solve_network


def test_single_tube_poiseuille():
    L, r, mu, dp = 100e-6, 2e-6, .00345312, 123.
    R = linear_radius_resistance(L, r, r, mu)
    s = solve_network(2, [[0, 1]], np.array([R]), {0: dp, 1: 0.})
    assert s.flow[0] == pytest.approx(np.pi*r**4*dp/(8*mu*L), rel=1e-13)
