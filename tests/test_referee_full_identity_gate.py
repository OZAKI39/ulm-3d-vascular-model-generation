from __future__ import annotations

from utils.cfd_flow.full_timestep_mass_referee import (
    FULL_IDENTITY_GATE,
    full_identity_pass,
)


def test_full_identity_gate_is_not_loosened() -> None:
    assert FULL_IDENTITY_GATE == 1.0e-8
    assert full_identity_pass([1.0e-12, 1.0e-8], 5.0e-9)
    assert not full_identity_pass([1.0e-12, 1.0000001e-8], 5.0e-9)
    assert not full_identity_pass([1.0e-12], 1.0000001e-8)
    assert not full_identity_pass([], 0.0)
