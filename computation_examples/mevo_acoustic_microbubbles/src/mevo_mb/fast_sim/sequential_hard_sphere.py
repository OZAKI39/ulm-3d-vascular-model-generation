"""Opt-in bounded true-sphere SQP, fixed original F/A and physical time.

s=(x_new-x_n)/R. At candidate q: B(q)s+[g(q)+B(q)(x_n-q)]/R>=0.
Only contact geometry is implicit. The existing F10P active-set QP is reused.
No spring, position projection, force refresh, or accumulation of physical steps.
"""

from time import perf_counter
import numpy as np
from scipy import linalg
from .f10g_contact import contact_geometry
from .hard_sphere import active_set_qp


def sequential_step(x, force, radius, zeta, dt, A, policy):
    begin = perf_counter()
    x, force, A = [np.asarray(a, dtype=float) for a in (x, force, A)]
    if (
        min(radius, zeta, dt) <= 0
        or x.shape != force.shape
        or not all(np.isfinite(a).all() for a in (x, force, A))
    ):
        raise ValueError("Finite matching inputs and positive scales required")
    if abs(A - A.T).max() > 1e-12:
        raise ValueError("Nonsymmetric resistance")
    factor = linalg.cho_factor(A, lower=True)
    _, pairs, g0 = contact_geometry(x, radius)
    if g0.min(initial=0) < -policy["hard_geometry_m"]:
        raise ValueError("Illegal initial hard-sphere state")
    limit = int(policy.get("max_sequential_qp", 3))
    if not 1 <= limit <= 5:
        raise ValueError("Sequential QP hard cap is five")
    b = dt * force.ravel() / (zeta * radius)
    free = linalg.cho_solve(factor, b)
    q = x.copy()
    history = []
    for iteration in range(limit):
        B, _, g = contact_geometry(q, radius)
        offset = (g + B @ (x - q).ravel()) / radius
        if np.all(B @ free + offset >= 0):
            s, lam, status, its = free, np.zeros(len(pairs)), "unconstrained_verified", 0
        else:
            s, lam, status, its = active_set_qp(A, factor, b, B, offset)
        y = x + radius * s.reshape(x.shape)
        By, py, gy = contact_geometry(y, radius)
        assert np.array_equal(pairs, py)
        slack = radius * (B @ s + offset)
        reaction = zeta * radius / dt * lam
        local_stat = A @ s - b - B.T @ lam
        true_stat = A @ s - b - By.T @ lam

        def resid(stat):
            return (
                float(
                    np.linalg.norm((radius * linalg.cho_solve(factor, stat)).reshape(x.shape), axis=1).max()
                ),
                float(np.linalg.norm((zeta * radius / dt * stat).reshape(x.shape), axis=1).max()),
            )

        local_m, local_N = resid(local_stat)
        true_m, true_N = resid(true_stat)
        scale = max(float(np.linalg.norm(force, axis=1).max()), float(abs(reaction).max(initial=0)), 1e-18)
        force_tol = policy["force_residual_atol_N"] + policy["force_residual_rtol"] * scale
        primal = float(np.maximum(0, -slack).max(initial=0))
        dual = float(np.maximum(0, -reaction).max(initial=0))
        comp_Nm = float(abs(reaction * slack).max(initial=0))
        comp = comp_Nm / scale
        truecomp = float(abs(reaction * gy).max(initial=0)) / scale
        overlap = np.maximum(0, -gy)
        dx = radius * (s - free).reshape(x.shape)
        correction = float(np.linalg.norm(dx, axis=1).max())
        gates = dict(
            LOCAL_QP_PASS=bool(
                status in ["active_set_solved", "unconstrained_verified"]
                and local_m <= policy["mobility_residual_m"]
                and local_N <= force_tol
                and primal <= policy["primal_violation_m"]
                and dual <= policy["dual_negative_N"]
                and comp <= policy["complementarity_equivalent_gap_m"]
            ),
            TRUE_GEOMETRY_PASS=bool(
                overlap.max(initial=0) <= min(policy["hard_geometry_m"], policy["contact_overlap_m"])
            ),
            TRUE_GAP_COMPLEMENTARITY_PASS=bool(truecomp <= policy["actual_geometry_complementarity_m"]),
            CONTACT_REACTION_BALANCE_PASS=bool(
                true_m <= policy["mobility_residual_m"] and true_N <= force_tol
            ),
            CONTACT_CORRECTION=bool(correction <= policy["contact_correction_m"]),
        )
        history.append(
            dict(
                index=iteration + 1,
                linearization_m=q.copy(),
                positions_m=y.copy(),
                pair_force_N=reaction.copy(),
                linear_gap_m=slack.copy(),
                actual_gap_m=gy.copy(),
                active_count=int(np.sum(reaction > policy["dual_negative_N"])),
                local_residual_m=local_m,
                local_force_residual_N=local_N,
                true_residual_m=true_m,
                true_force_residual_N=true_N,
                true_complementarity_m=truecomp,
                overlap_m=float(overlap.max(initial=0)),
                working_set_iterations=its,
                gates=gates.copy(),
            )
        )
        if all(gates.values()) or not gates["LOCAL_QP_PASS"]:
            break
        q = y.copy()
    contact = (By.T @ reaction).reshape(x.shape)
    return dict(
        positions_m=y,
        prediction_m=x + radius * free.reshape(x.shape),
        velocity_m_s=radius / dt * s.reshape(x.shape),
        contact_displacement_m=dx,
        contact_force_N=contact,
        local_contact_force_N=(B.T @ reaction).reshape(x.shape),
        pair_force_N=reaction,
        pairs=pairs,
        initial_gap_m=g0,
        linear_gap_m=slack,
        actual_gap_m=gy,
        geometric_overlap_m=overlap,
        mobility_equation_residual_m=true_m,
        force_equation_residual_N=true_N,
        local_mobility_residual_m=local_m,
        local_force_residual_N=local_N,
        primal_violation_m=primal,
        dual_negative_N=dual,
        complementarity_N_m=comp_Nm,
        complementarity_equivalent_gap_m=comp,
        actual_complementarity_equivalent_gap_m=truecomp,
        net_contact_force_N=contact.sum(0),
        symmetry_error=float(abs(A - A.T).max()),
        SPD_verified=True,
        max_single_iteration_m=correction,
        max_total_contact_correction_m=correction,
        iterations=history,
        qp_iterations=len(history),
        working_set_iterations=sum(h["working_set_iterations"] for h in history),
        active_count=history[-1]["active_count"],
        solver_status=status,
        total_s=perf_counter() - begin,
        gates=gates,
        converged=all(gates.values()),
        dt_s=dt,
        force_scale_N=scale,
        model="F10Q sequential hard spheres; F10G resistance; explicit frozen forces/resistance; true-normal contact balance",
    )
