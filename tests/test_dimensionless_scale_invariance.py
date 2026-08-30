from __future__ import annotations

import numpy as np

from utils.cfd_flow.dimensionless_geometry_kernel import dimensionless_ray_triangle


SCALES = 10.0 ** np.arange(-9, 10, 3)


def test_well_conditioned_hit_is_scale_invariant_to_one_e_minus_twelve() -> None:
    origin = np.asarray((3.0, -2.0, 4.0))
    link = np.asarray((1.0, 0.5, -0.25))
    point = origin + 0.37 * link
    triangle = np.asarray((point + (0, -1, -2), point + (0, 2, -1), point + (0, -1, 2)))
    results = [
        dimensionless_ray_triangle(origin * scale, link * scale, triangle * scale)
        for scale in SCALES
    ]

    assert {result.hit for result in results} == {True}
    assert max(abs(result.q - results[0].q) for result in results) <= 1.0e-12


def test_no_hit_parallel_and_degenerate_classifications_are_scale_invariant() -> None:
    cases = (
        (
            np.zeros(3),
            np.asarray((1.0, 0.0, 0.0)),
            np.asarray(((0, 0, 1), (1, 0, 1), (0, 1, 1))),
        ),
        (
            np.zeros(3),
            np.asarray((0.0, 0.0, 1.0)),
            np.asarray(((0, 0, 1), (1, 0, 1), (2, 0, 1))),
        ),
    )
    for origin, link, triangle in cases:
        classifications = {
            dimensionless_ray_triangle(
                origin * scale, link * scale, triangle * scale
            ).classification
            for scale in SCALES
        }
        assert len(classifications) == 1
