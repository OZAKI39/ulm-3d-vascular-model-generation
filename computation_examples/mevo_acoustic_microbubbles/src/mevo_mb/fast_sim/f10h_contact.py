"""Bounded nonlinear geometric overdamped penalty, fixed initial-step resistance.

Only the spring distance/normal/active set is refreshed. k_c is the F10G 2kc
candidate. This is a CPU diagnostic and does not change any production default.
"""

from time import perf_counter
import numpy as np
from scipy import linalg, optimize
from .f10g_contact import contact_geometry, normals, overlaps


def geometric_contact(x, free, radius, zeta, dt, kc, A, max_iterations=24, atol_m=1e-13):
    start = perf_counter()
    if min(radius, zeta, dt, kc) <= 0:
        raise ValueError("Positive physical parameters required")
    x = np.asarray(x, float)
    free = np.asarray(free, float)
    A = np.asarray(A, float)
    if (
        not np.all(np.isfinite(x))
        or not np.all(np.isfinite(free))
        or not np.allclose(A, A.T, rtol=0, atol=1e-12)
    ):
        raise ValueError("Finite data and symmetric resistance required")
    factor = linalg.cho_factor(A)
    u = linalg.cho_solve(factor, free.ravel()).reshape(-1, 3)
    origin = x.mean(0)
    prediction = x + dt * u
    yp = (x - origin + dt * u) / radius
    alpha = dt * kc / zeta
    q = np.zeros(x.size)
    records = []
    initialization = perf_counter() - start

    def evaluate(q):
        y = yp + q.reshape(-1, 3)
        B, pairs, gaps = contact_geometry(y, 1.0)
        delta = np.maximum(0, -gaps)
        grad = A @ q - alpha * (B.T @ delta)
        energy = 0.5 * (q @ A @ q + alpha * (delta @ delta))
        return energy, grad, B, delta, pairs

    begin = perf_counter()
    converged = False
    for iteration in range(max_iterations + 1):
        energy, grad, B, delta, pairs = evaluate(q)
        residual = linalg.cho_solve(factor, grad) * radius
        res = float(np.linalg.norm(residual.reshape(-1, 3), axis=1).max())
        if res <= atol_m:
            converged = True
            break
        if iteration == max_iterations:
            break
        # SPD Gauss-Newton metric, NOT an assertion that the nonlinear energy is convex.
        active = delta > 0
        H = A + alpha * B[active].T @ B[active]
        direction = linalg.solve(H, -grad, assume_a="pos")
        rate = 1.0
        for _ in range(20):
            if evaluate(q + rate * direction)[0] <= energy + 1e-4 * rate * np.dot(grad, direction):
                break
            rate *= 0.5
        else:
            break
        dx = rate * direction * radius
        q += rate * direction
        records.append(
            dict(
                iteration=iteration + 1,
                residual_before_m=res,
                line_search_rate=rate,
                max_correction_m=float(np.linalg.norm(dx.reshape(-1, 3), axis=1).max()),
                correction_m=dx.reshape(-1, 3),
                active_contacts=int(active.sum()),
                objective=evaluate(q)[0],
            )
        )
    contact_dx = q.reshape(-1, 3) * radius
    positions = prediction + contact_dx
    # Re-evaluate on returned *actual* global coordinates, not a clipped linearized gap.
    B, pairs, gap = contact_geometry(positions, radius)
    delta = np.maximum(-gap, 0)
    pair = kc * delta
    force = (B.T @ pair).reshape(-1, 3)
    velocity = u + linalg.cho_solve(factor, force.ravel() / zeta).reshape(-1, 3)
    consistent = positions - x - dt * velocity
    actual_res = float(np.linalg.norm(consistent, axis=1).max())
    B0, _, g0 = contact_geometry(x, radius)
    return dict(
        positions_m=positions,
        prediction_m=prediction,
        free_m_s=free,
        mobility_free_m_s=u,
        corrected_velocity_m_s=velocity,
        contact_displacement_m=contact_dx,
        contact_force_N=force,
        pair_force_N=pair,
        pairs=pairs,
        geometric_gap_m=gap,
        geometric_overlap_m=overlaps(positions, radius),
        frozen_normal_overlap_m=np.maximum(-g0 - B0 @ (positions - x).ravel(), 0),
        current_normals=normals(positions)[3],
        initial_normals=normals(x)[3],
        relative_normal_velocity_m_s=B @ velocity.ravel(),
        mobility_equation_residual_m=actual_res,
        force_position_consistency_m=float(np.max(abs(consistent))),
        net_contact_force_N=force.sum(0),
        contact_torque_N_m=np.sum(np.cross(positions - positions.mean(0), force), axis=0),
        centroid_contact_shift_m=contact_dx.mean(0),
        converged=converged and actual_res <= atol_m * 1.01,
        iterations=records,
        max_single_iteration_m=max((v["max_correction_m"] for v in records), default=0),
        max_total_contact_correction_m=float(np.linalg.norm(contact_dx, axis=1).max()),
        initialization_s=initialization,
        contact_s=perf_counter() - begin,
        total_s=perf_counter() - start,
        kc_N_m=kc,
        dt_s=dt,
        model="Actual final-position geometric unilateral spring; initial-step resistance fixed; bounded nonlinear Gauss-Newton",
    )


def geometric_reference(x, free, radius, zeta, dt, kc, A):
    """Independent force-loop residual and finite-difference Powell hybrid solve.

    Does not use candidate B matrices, gradients, active set, or Newton metric.
    A is still the same archived resistance model; no new physical reference.
    """
    start = perf_counter()
    x = np.asarray(x, float)
    D = np.linalg.inv(A)
    u = (D @ np.asarray(free).ravel()).reshape(-1, 3)
    yp = (x - x.mean(0) + dt * u) / radius
    alpha = dt * kc / zeta

    def force(y):
        f = np.zeros_like(y)
        for i in range(len(y)):
            for j in range(i + 1, len(y)):
                dx = y[j] - y[i]
                d = np.linalg.norm(dx)
                if d < 2:
                    v = (2 - d) * dx / d
                    f[i] -= v
                    f[j] += v
        return f.ravel()

    def residual(q):
        return q - alpha * (D @ force(yp + q.reshape(-1, 3)))

    result = optimize.root(residual, np.zeros(x.size), method="hybr", options=dict(xtol=1e-10, maxfev=2000))
    res = float(np.linalg.norm(residual(result.x).reshape(-1, 3) * radius, axis=1).max())
    pos = x + dt * u + result.x.reshape(-1, 3) * radius
    return dict(
        positions_m=pos,
        contact_force_N=force(yp + result.x.reshape(-1, 3)).reshape(-1, 3) * kc * radius,
        residual_m=res,
        solver_success=bool(result.success),
        solver_message=str(result.message),
        function_evaluations=int(result.nfev),
        reference_pass=res < 1e-14,
        total_s=perf_counter() - start,
    )
