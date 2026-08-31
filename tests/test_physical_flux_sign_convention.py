from __future__ import annotations

import numpy as np

from utils.cfd_flow.physical_port_flux import (
    InteriorPlane,
    build_polygon_quadrature,
    integrate_reconstructed_flux,
)


def _plane(label: str) -> InteriorPlane:
    return InteriorPlane(
        label,
        "synthetic",
        np.zeros(3),
        np.asarray((0.0, 0.0, 1.0)),
        np.asarray((1.0, 0.0, 0.0)),
        np.asarray((0.0, 1.0, 0.0)),
        np.asarray(((-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0))),
        "synthetic",
    )


def test_inlet_and_outlet_outward_sign_convention() -> None:
    outlet = _plane("outlet_01")
    inlet = _plane("inlet")
    quadrature = build_polygon_quadrature(outlet.aperture)
    outward = np.repeat([[0.0, 0.0, 2.0]], len(quadrature.high_points_uv_m), axis=0)
    outward_low = np.repeat([[0.0, 0.0, 2.0]], len(quadrature.low_points_uv_m), axis=0)
    qout = integrate_reconstructed_flux(
        plane=outlet,
        quadrature=quadrature,
        high_velocity_m_s=outward,
        low_velocity_m_s=outward_low,
    )["physical_q_m3_s"]
    qin = integrate_reconstructed_flux(
        plane=inlet,
        quadrature=quadrature,
        high_velocity_m_s=-outward,
        low_velocity_m_s=-outward_low,
    )["physical_q_m3_s"]
    assert qout > 0.0
    assert np.isclose(qin, qout)
