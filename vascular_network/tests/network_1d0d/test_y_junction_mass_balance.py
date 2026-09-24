"""Permanent analytic validation for steady Newtonian networks."""
import numpy as np
import pytest
from network_1d0d.network_solver import solve_network


def test_y_junction_mass_balance():
    s = solve_network(4, [[0, 1], [1, 2], [1, 3]], np.array([1., 2., 3.])*1e17, {0: 20., 2: 0., 3: 0.})
    assert s.flow[0] == pytest.approx(s.flow[1]+s.flow[2], rel=1e-13, abs=0)
    assert s.flow[1]/s.flow[2] == pytest.approx(1.5, rel=1e-13)
