"""Reference-only variable-radius no-slip translation; never imported by production.

Garbin 2009 Eq (1), Cattaneo & Supponen arXiv:2506.16413v2 Eqs (11),(14).
Zero-Re, no gas inertia. Conservative added momentum M(t)v(t), M=rho*V/2.
The Abel kernel is integrated exactly over piecewise-linear q=R*v in theta.
BDF2 momentum marching (backward Euler first step); no memory truncation.
"""

import numpy as np
from scipy.integrate import cumulative_trapezoid


def history_integral(theta, q):
    """Integral dq/dtheta / sqrt(theta_n-theta) dtheta at each node."""
    result = np.zeros(len(q))
    dq = np.diff(q)
    for n in range(1, len(q)):
        weights = 2 / (np.sqrt(theta[n] - theta[:n]) + np.sqrt(theta[n] - theta[1 : n + 1]))
        result[n] = weights @ dq[:n]
    return result


def solve(time_s, radius_m, force_N, *, density=1000.0, viscosity=0.001, history=False):
    t, r, force = map(lambda x: np.asarray(x, dtype=float), (time_s, radius_m, force_N))
    if t.ndim != 1 or t.shape != r.shape or t.shape != force.shape or len(t) < 3:
        raise ValueError("Matched one-dimensional time, radius, force required")
    dt = t[1] - t[0]
    if (
        not np.isfinite([*t, *r, *force, density, viscosity]).all()
        or np.any(r <= 0)
        or density <= 0
        or viscosity <= 0
        or dt <= 0
        or not np.allclose(np.diff(t), dt, rtol=1e-7, atol=0)
    ):
        raise ValueError("Finite, positive radius/medium and uniform increasing clock required")
    mass = 0.5 * density * 4 * np.pi * r**3 / 3
    drag = 6 * np.pi * viscosity * r
    theta = cumulative_trapezoid(1 / r**2, t, initial=0)
    coeff = 6 * np.pi * viscosity / np.sqrt(np.pi * viscosity / density)
    v = np.zeros(len(t))
    q = np.zeros(len(t))
    dq = np.zeros(len(t) - 1)
    fh = np.zeros(len(t))
    fam = np.zeros(len(t))
    for n in range(1, len(t)):
        alpha = 1.0 if n == 1 else 1.5
        rhs_m = (
            mass[n - 1] * v[n - 1] / dt
            if n == 1
            else (2 * mass[n - 1] * v[n - 1] - 0.5 * mass[n - 2] * v[n - 2]) / dt
        )
        past = last = 0.0
        if history:
            weights = 2 / (np.sqrt(theta[n] - theta[:n]) + np.sqrt(theta[n] - theta[1 : n + 1]))
            past = weights[:-1] @ dq[: n - 1]
            last = weights[-1]
        v[n] = (force[n] + rhs_m - coeff * past + coeff * last * q[n - 1]) / (
            alpha * mass[n] / dt + drag[n] + coeff * last * r[n]
        )
        q[n] = r[n] * v[n]
        dq[n - 1] = q[n] - q[n - 1]
        fh[n] = -coeff * (past + last * dq[n - 1])
        fam[n] = rhs_m - alpha * mass[n] * v[n] / dt
    return dict(
        time_s=t,
        radius_m=r,
        force_N=force,
        velocity_m_s=v,
        displacement_m=cumulative_trapezoid(v, t, initial=0),
        quasisteady_drag_N=-drag * v,
        added_mass_N=fam,
        history_force_N=fh,
    )
