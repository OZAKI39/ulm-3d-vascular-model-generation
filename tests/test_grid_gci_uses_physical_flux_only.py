from __future__ import annotations

import pytest

from utils.cfd_flow.full_timestep_mass_referee import BOUNDARY_FLUX_DEFINITION, PHYSICAL_FLUX_DEFINITION
from utils.cfd_flow.tau1_grid_convergence import validate_physical_flux_observables


def test_boundary_pdf_accounting_flux_cannot_enter_primary_gci() -> None:
    records = {label: {"flux_definition": PHYSICAL_FLUX_DEFINITION} for label in ("coarse", "base", "fine")}
    validate_physical_flux_observables(records)
    records["base"]["flux_definition"] = BOUNDARY_FLUX_DEFINITION
    with pytest.raises(ValueError):
        validate_physical_flux_observables(records)
