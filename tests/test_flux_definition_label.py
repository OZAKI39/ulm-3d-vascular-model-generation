from __future__ import annotations

from utils.cfd_flow.full_timestep_mass_referee import (
    BOUNDARY_FLUX_DEFINITION,
    DEFERRED_PHYSICAL_FLUX,
    PHYSICAL_FLUX_DEFINITION,
)


def test_flux_definitions_are_explicit_and_non_aliasing() -> None:
    assert BOUNDARY_FLUX_DEFINITION == "BOUNDARY_PDF_ACCOUNTING_FLUX"
    assert PHYSICAL_FLUX_DEFINITION == "PHYSICAL_CROSS_SECTION_FLUX"
    assert BOUNDARY_FLUX_DEFINITION != PHYSICAL_FLUX_DEFINITION
    assert DEFERRED_PHYSICAL_FLUX == "DEFERRED_TO_GRID_CONVERGENCE_STAGE"
