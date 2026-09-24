from __future__ import annotations

import numpy as np

from utils.cfd_flow.wall_qvalue_oracle import ray_cylinder_fraction


def test_cylinder_qvalue_oracle_known_radial_links() -> None:
    origins = np.asarray([[0.0, 0.0, 0.0], [0.5, 0.0, 1.0]])
    directions = np.asarray([[2.0, 0.0, 0.0], [1.0, 0.0, 0.0]])

    qvalues = ray_cylinder_fraction(
        origins,
        directions,
        axis=(0.0, 0.0, 1.0),
        radius_m=1.0,
    )

    assert np.allclose(qvalues, [0.5, 0.5])
