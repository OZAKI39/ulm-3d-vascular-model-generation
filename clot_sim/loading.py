"""Frozen uniform +X tangential slip surrogate, no spatial window or controller."""

import numpy as np
from .config import BETA, RISE, HOLD, OFF, LABEL, MU, H_EFF, RHO, NOMINAL_SHEAR


class ShearField:
    def __init__(self, mode, theta=None, surface_z=0.0):
        assert mode == "uniform" and theta is None
        self.mode = mode
        self.beta = BETA
        self.U = NOMINAL_SHEAR / BETA
        self.theta = dict(
            nominal_shear_Pa=NOMINAL_SHEAR,
            direction=[1.0, 0.0, 0.0],
            drive_duration_s=OFF,
        )

    def envelope(self, t):
        t = np.asarray(t)
        u = np.clip(t / RISE, 0.0, 1.0)
        v = np.clip((t - HOLD) / (OFF - HOLD), 0.0, 1.0)
        return np.where(
            t < 0,
            0.0,
            np.where(
                t < RISE,
                3 * u * u - 2 * u * u * u,
                np.where(t < HOLD, 1.0, 1 - 3 * v * v + 2 * v * v * v),
            ),
        )

    def spatial(self, x, n):
        raw = np.zeros_like(x)
        raw[:, 0] = 1.0
        return raw - n * n[:, 0, None]

    def __call__(self, x, n, t):
        return self.U * float(self.envelope(t)) * self.spatial(x, n)

    def amplitude(self, t):
        return NOMINAL_SHEAR

    def applicability(self):
        return dict(
            label=LABEL,
            mode="UNIFORM_SURFACE_SHEAR_X",
            nominal_stationary_shear_cap_Pa=NOMINAL_SHEAR,
            outer_speed_bound_m_s=self.U,
            mu_Pa_s=MU,
            h_eff_m=H_EFF,
            Re_h_bound=RHO * self.U * H_EFF / MU,
            acoustic_realizability="UNESTABLISHED",
            meaning="Prescribed current-surface tangential slip proxy. No solved acoustic streaming or calibrated bubble power. Passive drag remains after drive-off.",
        )
