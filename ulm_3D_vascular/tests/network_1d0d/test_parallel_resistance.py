"""Permanent analytic validation for steady Newtonian networks."""
import numpy as np
import pytest
from network_1d0d.network_solver import solve_network


def test_parallel_resistance():
    s = solve_network(2, [[0, 1], [0, 1]], np.array([2e17, 5e17]), {0: 10., 1: 0.})
    assert s.node_outflow[0] == pytest.approx(10*(1/2e17+1/5e17), rel=1e-13, abs=0)
