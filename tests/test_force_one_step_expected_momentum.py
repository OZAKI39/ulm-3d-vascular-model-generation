from __future__ import annotations

import numpy as np

from utils.cfd_flow.musubi_wall_force_diagnostics import (
    expected_force_momentum_increment,
)


def test_force_one_step_expected_momentum() -> None:
    force = np.array([1.25, -3.5, 7.0])
    result = expected_force_momentum_increment(force, volume_m3=2.5, dt_s=0.125)
    np.testing.assert_allclose(result, [0.390625, -1.09375, 2.1875], rtol=0, atol=0)
