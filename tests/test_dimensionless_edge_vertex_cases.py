from __future__ import annotations

import numpy as np

from utils.cfd_flow.dimensionless_geometry_kernel import dimensionless_ray_triangle


TRIANGLE = np.asarray(((0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (0.0, 1.0, 1.0)))


def test_vertex_and_edge_classification_is_deterministic_across_units() -> None:
    for xy, expected in (((0.0, 0.0), "HIT_VERTEX"), ((0.5, 0.0), "HIT_EDGE")):
        classifications = set()
        for scale in 10.0 ** np.arange(-9, 10, 3):
            origin = np.asarray((xy[0], xy[1], 0.0)) * scale
            result = dimensionless_ray_triangle(
                origin,
                np.asarray((0.0, 0.0, 2.0)) * scale,
                TRIANGLE * scale,
            )
            assert result.hit is True
            classifications.add(result.classification)
        assert classifications == {expected}


def test_grazing_ray_has_deterministic_no_hit_classification() -> None:
    classifications = {
        dimensionless_ray_triangle(
            np.zeros(3),
            np.asarray((1.0, 0.0, 1.0e-16)) * scale,
            TRIANGLE * scale,
        ).classification
        for scale in 10.0 ** np.arange(-9, 10, 3)
    }

    assert classifications == {"PARALLEL_OR_GRAZING"}
