from __future__ import annotations

import numpy as np
from shapely.geometry import Polygon

from utils.cfd_flow.physical_port_flux import build_polygon_quadrature


def test_constrained_triangle_area_equals_aperture_area() -> None:
    polygon = Polygon(((0.0, 0.0), (2.0, 0.0), (1.6, 0.7), (0.8, 0.4), (0.0, 1.0)))
    quadrature = build_polygon_quadrature(polygon)
    assert quadrature.area_relative_error <= 1.0e-12
    assert np.isclose(
        quadrature.high_weights_m2.sum(),
        quadrature.aperture_area_m2,
        rtol=0.0,
        atol=2.0e-15,
    )
