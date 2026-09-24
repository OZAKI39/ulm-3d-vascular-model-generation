from __future__ import annotations

import numpy as np

from utils.cfd_flow.full_timestep_mass_referee import (
    stable_delta,
    tau1_bgk_collision,
)


def test_tau1_bgk_is_numerically_constructed_and_mass_conserving() -> None:
    generator = np.random.default_rng(274)
    fetched = generator.uniform(1.0e-5, 8.0e-4, size=(31, 19))

    collided = tau1_bgk_collision(fetched)

    assert collided.shape == fetched.shape
    assert not np.array_equal(collided, fetched)
    assert abs(stable_delta(collided, fetched)) <= 2.0e-15
    assert np.allclose(
        np.sum(collided, axis=1), np.sum(fetched, axis=1), rtol=0.0, atol=2.0e-18
    )
