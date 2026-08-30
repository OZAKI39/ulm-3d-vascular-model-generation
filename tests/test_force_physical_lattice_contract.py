from __future__ import annotations

import numpy as np

from utils.cfd_flow.musubi_wall_force_diagnostics import body_force_conversion


def test_force_physical_lattice_contract() -> None:
    force = np.array([0.0, 0.0, 3_884_423.6432636492])
    rho0 = 1056.0
    dx = 1.1686589466964423e-7
    dt = 8.3359602886574595e-9
    result = body_force_conversion(force, rho0_kg_m3=rho0, dx_m=dx, dt_s=dt)
    expected = force * dt**2 / (rho0 * dx)
    assert result["formula"] == "F_lat = F_phy * dt^2 / (rho0 * dx)"
    np.testing.assert_allclose(result["lattice_force_density"], expected, rtol=0, atol=0)
