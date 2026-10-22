"""Existing point-particle translation/rotation in a resolved analytical tube field.

RBCs use the project's oblate Jeffery mobility. The normal biconcave surface
is a rigid display geometry, not a newly solved membrane or mobility model.
MB orientation is a passive body-frame tracer of the existing half-curl spin.
"""
from types import SimpleNamespace
import numpy as np
from .rbc_orientation import angular_velocity, advance_orientation, short_axis


class PoiseuilleField:
    def __init__(self, radius_m=7e-6, mean_m_s=.002):
        self.radius_m = float(radius_m)
        self.mean_m_s = float(mean_m_s)
        if self.radius_m <= 0 or self.mean_m_s <= 0:
            raise ValueError('Positive radius and mean flow required')

    def sample(self, position_m):
        p = np.asarray(position_m, dtype=float)
        inside = bool(p[1:] @ p[1:] <= self.radius_m**2)
        g = np.zeros((3, 3))
        g[0, 1:] = -4 * self.mean_m_s * p[1:] / self.radius_m**2
        u = np.array([2*self.mean_m_s*(1-p[1:]@p[1:]/self.radius_m**2), 0., 0.])
        curl = np.array([0., g[0, 2], -g[0, 1]])
        return SimpleNamespace(inside_lumen=inside, tetra_id=0 if inside else -1,
            velocity_m_s=u, velocity_gradient_s_inv=g, vorticity_s_inv=curl)


def integrate_orientations(times_s, initial_q, rbc_count, positions0_m,
                           jeffery_lambda, field, max_dt_s):
    """RBC: existing frozen-old-Omega quaternion steps, MB: exact constant spin.

    All RBC centers share one transverse lane (and therefore one gradient).
    Substeps end exactly at every video sample; no orientation interpolation.
    """
    times = np.asarray(times_s)
    if times[0] != 0 or np.any(np.diff(times) <= 0) or max_dt_s <= 0:
        raise ValueError('Increasing times starting at zero and positive dt required')
    q0 = np.asarray(initial_q, dtype=float)
    samples = [field.sample(p) for p in positions0_m]
    g = samples[0].velocity_gradient_s_inv
    assert all(np.array_equal(s.velocity_gradient_s_inv, g) for s in samples[:rbc_count])
    assert all(s.inside_lumen for s in samples)
    q = np.empty((len(times), len(q0), 4)); q[0] = q0
    omega = np.empty((len(times), len(q0), 3))
    step_count = 0
    for k, t in enumerate(times):
        if k:
            htotal = t-times[k-1]
            n = int(np.ceil(htotal/max_dt_s)); h = htotal/n
            current = q[k-1, :rbc_count].copy()
            for _ in range(n):
                w = angular_velocity(short_axis(current), jeffery_lambda, g)
                current = advance_orientation(current, w, h)
            q[k, :rbc_count] = current
            spin = .5*np.array([s.vorticity_s_inv for s in samples[rbc_count:]])
            q[k, rbc_count:] = advance_orientation(q0[rbc_count:], spin, t)
            step_count += n
        omega[k, :rbc_count] = angular_velocity(short_axis(q[k, :rbc_count]), jeffery_lambda, g)
        omega[k, rbc_count:] = .5*np.array([s.vorticity_s_inv for s in samples[rbc_count:]])
    return q, omega, step_count
