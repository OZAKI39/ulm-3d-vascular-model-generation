from __future__ import annotations

import numpy as np

from utils.cfd_flow.physical_port_flux import (
    MLS_CONDITION_GATE,
    build_mls_stencil_map,
    evaluate_mls_stencil_map,
)


def test_mls_reproduces_affine_field_with_valid_stencils() -> None:
    axis = np.arange(-3, 4, dtype=np.float64)
    x, y, z = np.meshgrid(axis, axis, axis, indexing="ij")
    points = np.column_stack((x.ravel(), y.ravel(), z.ravel()))
    queries = np.asarray(((0.2, -0.1, 0.3), (-0.4, 0.25, -0.15)))
    values = 2.0 + 0.3 * points[:, 0] - 0.2 * points[:, 1] + 0.7 * points[:, 2]
    stencil_map = build_mls_stencil_map(
        points, queries, dx_m=1.0, polynomial_order=1
    )
    actual = evaluate_mls_stencil_map(stencil_map, values)
    expected = 2.0 + 0.3 * queries[:, 0] - 0.2 * queries[:, 1] + 0.7 * queries[:, 2]
    assert stencil_map.invalid_count == 0
    assert stencil_map.max_condition_number <= MLS_CONDITION_GATE
    assert np.allclose(actual, expected, rtol=0.0, atol=1.0e-13)
