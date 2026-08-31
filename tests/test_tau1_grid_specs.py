from __future__ import annotations

import pytest

from utils.cfd_flow.tau1_grid_convergence import GRID_SPECS, REFINEMENT_RATIO, validate_grid_spec


def test_tau1_grid_specs() -> None:
    assert REFINEMENT_RATIO == pytest.approx(1.3)
    assert GRID_SPECS["coarse"].dx_m == pytest.approx(2.6e-7)
    assert GRID_SPECS["base"].dx_m == pytest.approx(2.0e-7)
    assert GRID_SPECS["fine"].dx_m == pytest.approx(1.5384615384615385e-7)
    for spec in GRID_SPECS.values():
        assert all(validate_grid_spec(spec).values())
        assert spec.tau == pytest.approx(1.0, abs=1e-12)
        assert spec.omega == pytest.approx(1.0, abs=1e-12)
