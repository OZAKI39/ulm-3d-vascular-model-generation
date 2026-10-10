"""Slow-time drag interface and cycle average, separated from radial solvers."""

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class QuasiSteadyStokesTranslation:
    viscosity_Pa_s: float

    def velocity(self, force_N, radii_m):
        if self.viscosity_Pa_s <= 0 or np.any(np.asarray(radii_m) <= 0):
            raise ValueError("Positive viscosity and radii required")
        return np.asarray(force_N) / (6 * np.pi * self.viscosity_Pa_s * np.asarray(radii_m)[..., None])


def primary_cycle_average(time_s, volume_m3, gradient_Pa_m):
    """Full numerical integral over caller-supplied complete-cycle window.

    Endpoints included. Gradients may include pressure harmonics or transients.
    This function does not assume a sinusoidal volume response.
    """
    t, v, g = np.asarray(time_s), np.asarray(volume_m3), np.asarray(gradient_Pa_m)
    if t.ndim != 1 or len(t) < 3 or np.any(np.diff(t) <= 0):
        raise ValueError("Increasing integration times required")
    if v.shape != t.shape or g.shape != (len(t), 3) or not np.isfinite(g).all() or not np.isfinite(v).all():
        raise ValueError("Finite volume and 3D gradient samples required")
    return -np.trapezoid(v[:, None] * g, t, axis=0) / (t[-1] - t[0])
