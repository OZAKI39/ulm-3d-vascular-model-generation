import numpy as np
import pytest

from utils.cfd_flow.musubi_wall_force_diagnostics import independent_cross_section_flux


def test_flux_is_direct_axial_plane_quadrature() -> None:
    velocity = np.asarray([[1.0, 2.0, 3.0], [2.0, -1.0, 4.0]])
    flux = independent_cross_section_flux(velocity, axis=(0.0, 0.0, 1.0), dx_m=0.25)
    assert flux == pytest.approx((3.0 + 4.0) * 0.25**2)
