"""Opt-in local hard-sphere velocity QP; no spring and no position projection.

s=dt*v/R, b=dt*F/(zeta*R). Minimize .5*s.T*A*s-b.T*s subject
to B*s+g/R >= 0. A is the unchanged F10G resistance. Local normals only.
The SPD local problem has a unique velocity, not a unique nonlinear trajectory.
"""

from time import perf_counter
from itertools import product
import numpy as np
from scipy import linalg, sparse
from scipy.spatial.distance import pdist
import osqp
from .f10g_contact import contact_geometry


def active_set_qp(A, factor, b, B, gap_over_R, max_iterations=96):
    """Small exact equality-QP working sets; all inequalities rechecked.

    This solves the same resistance-metric velocity QP, never projects positions.
    Linear dependent working sets use a least-squares dual; inconsistent or cycling
    sets fail closed. No uniqueness of redundant reaction multipliers is claimed.
    """
    free = linalg.cho_solve(factor, b)
    inverse_BT = linalg.cho_solve(factor, B.T)
    active = set(np.flatnonzero(B @ free + gap_over_R < 0).tolist())
    visited = set()
    lam = np.zeros(len(B))
    for iteration in range(max_iterations):
        key = tuple(sorted(active))
        if key in visited:
            return free, lam, "active_set_cycle", iteration
        visited.add(key)
        ids = np.array(key, dtype=int)
        lam = np.zeros(len(B))
        if len(ids):
            S = B[ids] @ inverse_BT[:, ids]
            rhs = -gap_over_R[ids] - B[ids] @ free
            lam[ids] = linalg.lstsq(S, rhs, cond=1e-12, lapack_driver="gelsy")[0]
        s = free + inverse_BT @ lam
        slack = B @ s + gap_over_R
        if lam.min(initial=0) < -1e-12:
            active.remove(int(np.argmin(lam)))
            continue
        inactive = np.ones(len(B), dtype=bool)
        inactive[ids] = False
        violations = np.where(inactive, slack, np.inf)
        if violations.min(initial=0) < -1e-12:
            active.add(int(np.argmin(violations)))
            continue
        if abs(slack[ids]).max(initial=0) > 1e-10:
            return s, lam, "active_set_inconsistent_equalities", iteration + 1
        return s, lam, "active_set_solved", iteration + 1
    return s, lam, "active_set_iteration_limit", max_iterations


def hard_step(x, force, radius, zeta, dt, A, policy):
    start = perf_counter()
    x, force, A = [np.asarray(k, dtype=float) for k in (x, force, A)]
    if min(radius, zeta, dt) <= 0 or x.shape != force.shape:
        raise ValueError("Positive scales and matching positions/forces required")
    if not all(np.isfinite(k).all() for k in (x, force, A)):
        raise ValueError("Finite inputs required")
    symmetry = float(abs(A - A.T).max())
    if symmetry > 1e-12:
        raise ValueError("Nonsymmetric resistance")
    factor = linalg.cho_factor(A, lower=True)  # independently verifies SPD each step
    B, pairs, gap = contact_geometry(x, radius)
    if gap.min(initial=0) < -policy["hard_geometry_m"]:
        raise ValueError("Initial state violates hard-sphere diagnostic")
    b = dt * force.ravel() / (zeta * radius)
    free = linalg.cho_solve(factor, b)
    if np.all(B @ free + gap / radius >= 0):
        s, lam = free, np.zeros(len(gap))
        status, iterations, solve_s, polish = "unconstrained_verified", 0, 0.0, None
    elif policy.get("hard_solver") == "active_set":
        s, lam, status, iterations = active_set_qp(A, factor, b, B, gap / radius)
        solve_s, polish = 0.0, None
    else:
        qp = osqp.OSQP()
        qp.setup(
            P=sparse.csc_matrix(np.triu(A)),
            q=-b,
            A=sparse.csc_matrix(B),
            l=-gap / radius,
            u=np.full(len(gap), np.inf),
            **policy["osqp"],
        )
        result = qp.solve(raise_error=False)
        s, lam = result.x, -result.y  # lower-bound OSQP dual has negative sign
        status, iterations = result.info.status, result.info.iter
        solve_s, polish = result.info.run_time, result.info.status_polish
    if s is None or not np.isfinite(s).all() or not np.isfinite(lam).all():
        raise RuntimeError(f"QP returned no finite proposal: {status}")
    slack = B @ s + gap / radius
    stationarity = A @ s - b - B.T @ lam
    reaction = zeta * radius / dt * lam
    contact_force = (B.T @ reaction).reshape(x.shape)
    velocity = radius / dt * s.reshape(x.shape)
    new = x + dt * velocity
    actual_gap = pdist(new) - 2 * radius
    contact_dx = radius * (s - free).reshape(x.shape)
    force_residual = zeta * radius / dt * stationarity
    force_scale = max(float(np.linalg.norm(force, axis=1).max()), float(abs(reaction).max(initial=0)), 1e-18)
    residual_m = float(
        np.linalg.norm((radius * linalg.cho_solve(factor, stationarity)).reshape(x.shape), axis=1).max()
    )
    primal = float(np.maximum(0, -slack).max(initial=0) * radius)
    dual = float(np.maximum(0, -reaction).max(initial=0))
    comp_Nm = float(abs(reaction * slack * radius).max(initial=0))
    comp = comp_Nm / force_scale
    actual_comp = float(abs(reaction * actual_gap).max(initial=0)) / force_scale
    overlap = np.maximum(0, -actual_gap)
    correction = float(np.linalg.norm(contact_dx, axis=1).max())
    force_error = float(np.linalg.norm(force_residual.reshape(x.shape), axis=1).max())
    gates = dict(
        HARD_SPHERE_GEOMETRY=bool(
            overlap.max(initial=0) <= policy["hard_geometry_m"]
            and overlap.max(initial=0) <= policy["contact_overlap_m"]
            and primal <= policy["primal_violation_m"]
        ),
        SOLVER_RESIDUAL=bool(
            status in ["solved", "unconstrained_verified", "active_set_solved"]
            and residual_m <= policy["mobility_residual_m"]
            and force_error <= policy["force_residual_atol_N"] + policy["force_residual_rtol"] * force_scale
        ),
        CONTACT_COMPLEMENTARITY=bool(
            dual <= policy["dual_negative_N"]
            and comp <= policy["complementarity_equivalent_gap_m"]
            and actual_comp <= policy["actual_geometry_complementarity_m"]
        ),
        CONTACT_CORRECTION=bool(correction <= policy["contact_correction_m"]),
    )
    return dict(
        positions_m=new,
        prediction_m=x + radius * free.reshape(x.shape),
        velocity_m_s=velocity,
        contact_displacement_m=contact_dx,
        contact_force_N=contact_force,
        pair_force_N=reaction,
        pairs=pairs,
        initial_gap_m=gap,
        linear_gap_m=slack * radius,
        actual_gap_m=actual_gap,
        geometric_overlap_m=overlap,
        mobility_equation_residual_m=residual_m,
        force_equation_residual_N=force_error,
        primal_violation_m=primal,
        dual_negative_N=dual,
        complementarity_N_m=comp_Nm,
        complementarity_equivalent_gap_m=comp,
        actual_complementarity_equivalent_gap_m=actual_comp,
        net_contact_force_N=contact_force.sum(0),
        symmetry_error=symmetry,
        SPD_verified=True,
        max_single_iteration_m=correction,
        max_total_contact_correction_m=correction,
        iterations=[dict(count=int(iterations))] if iterations else [],
        qp_iterations=int(iterations),
        solver_status=status,
        polish_status=polish,
        qp_internal_s=solve_s,
        total_s=perf_counter() - start,
        gates=gates,
        converged=all(gates.values()),
        dt_s=dt,
        model="F10P local frictionless nonadhesive hard spheres; F10G resistance; no spring",
    )


def enumerated_hard_reference(x, force, radius, zeta, dt, A):
    """Independent <=4-sphere dual active-set enumeration, all <=64 subsets."""
    if len(x) > 4:
        raise ValueError("Reference limited to four particles")
    B, pairs, gap = contact_geometry(x, radius)
    inverse = np.linalg.inv(A)
    free = inverse @ (dt * force.ravel() / (zeta * radius))
    S = B @ inverse @ B.T
    w = B @ free + gap / radius
    for flags in product([False, True], repeat=len(pairs)):
        ids = np.flatnonzero(flags)
        lam = np.zeros(len(pairs))
        if len(ids):
            sub = S[np.ix_(ids, ids)]
            # Singular subsets are skipped; redundant zero multipliers need not be active.
            try:
                lam[ids] = np.linalg.solve(sub, -w[ids])
            except np.linalg.LinAlgError:
                continue
        slack = w + S @ lam
        if (
            lam.min(initial=0) >= -1e-10
            and slack.min(initial=0) >= -1e-10
            and abs(lam * slack).max(initial=0) < 1e-10
        ):
            s = free + inverse @ B.T @ lam
            return dict(positions_m=x + radius * s.reshape(x.shape), reaction_N=lam * zeta * radius / dt)
    raise RuntimeError("No independent feasible complementary active set")
