"""Permanent analytic validation for steady Newtonian networks."""
import numpy as np
import pytest
from network_1d0d.network_solver import solve_network
from network_1d0d.boundary_conditions import require_source, SourceAmbiguous, cap_pressures


def test_disconnected_ungrounded_network_rejected():
    with pytest.raises(ValueError, match='Unreferenced'):
        solve_network(4, [[0, 1], [2, 3]], np.array([1., 1.]), {0: 1., 1: 0.})


def test_source_gate_does_not_promote_serialization_root():
    with pytest.raises(SourceAmbiguous): require_source({'structural_root_id': 2410})
    assert require_source({'hydraulic_source_id': 2410, 'identification_kind': 'EXPLICIT_MODEL_ASSUMPTION', 'evidence_reference': 'explicit source declaration'}) == 2410


def test_extension_pressure_sign():
    np.testing.assert_allclose(cap_pressures([100., 100.], [2e-16, -2e-16], 1e17), [80., 120.])
