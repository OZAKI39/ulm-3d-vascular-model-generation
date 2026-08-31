from __future__ import annotations

import pytest

from utils.cfd_flow.full_timestep_mass_referee import PHYSICAL_FLUX_DEFINITION
from utils.cfd_flow.tau1_grid_convergence import HISTORICAL_CLASSIFICATION, validate_physical_flux_observables


def test_historical_fallback_q_grid_values_are_rejected() -> None:
    records = {label: {"flux_definition": PHYSICAL_FLUX_DEFINITION} for label in ("coarse", "base", "fine")}
    records["coarse"]["historical_classification"] = HISTORICAL_CLASSIFICATION
    with pytest.raises(ValueError):
        validate_physical_flux_observables(records)
