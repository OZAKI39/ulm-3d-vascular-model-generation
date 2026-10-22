"""R3: original nearest-wall resistance, multiple finite WALL contact features.

The R2 S1 equations/iteration are retained below as a versioned implementation.
Only contact geometry and the continuous certificate are changed in V1. Additional
wall rows are finite-triangle distance supporting tangents, queried at the actual
start state within the activation/reach band; no wall lubrication terms are added.
"""

from dataclasses import replace
from time import perf_counter
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu
from .motion_r2 import MotionR2, SemiImplicit, anchored_constraints, endpoint_wall_constraints
from .constrained_motion import Constraints, sparse_rows, SubdivideInterval, kkt_residual
from .motion_sweep_r3 import closest_point


class MotionR3(MotionR2):
    name = "constrained_lubrication_r3"


class SemiImplicitR3(SemiImplicit):
    def feature_constraints(self, x, r, wall, dt, ids, velocity):
        base = anchored_constraints(x, x, r, wall, wall, dt, ids, self.tolerance, velocity)
        extra_rows = []
        extra_cols = []
        extra_values = []
        lower = []
        keys = []
        for i in np.flatnonzero(wall[0] <= 0.1 * r + np.linalg.norm(velocity, axis=1) * dt):
            reach = 0.1 * r[i] + np.linalg.norm(velocity[i]) * dt
            cells = self.sweep.candidates(x[i], x[i], np.nextafter(r[i] + reach, np.inf))
            for cell in cells:
                point = closest_point(x[i], self.sweep.triangles[cell])
                delta = x[i] - point
                distance = np.linalg.norm(delta)
                gap = distance - r[i]
                if gap > reach or distance == 0:
                    continue
                inward = delta / distance
                row = len(lower)
                extra_rows.extend([row] * 3)
                extra_cols.extend((3 * i + np.arange(3)).tolist())
                extra_values.extend(inward)
                allowed = gap if gap <= self.tolerance else 0.0
                lower.append((allowed - gap) / dt)
                keys.append(("wall_feature", int(ids[i]), int(cell)))
        if not lower:
            return base
        extra = sparse.coo_matrix(
            (extra_values, (extra_rows, extra_cols)), shape=(len(lower), x.size)
        ).tocsc()
        return Constraints(
            sparse.vstack([base.matrix, extra], format="csc"),
            np.r_[base.lower, lower],
            base.ids + tuple(keys),
            base.pairs,
        )

    def advance(self, x, r, ids, t, dt, stage):
        self.last_candidate = None
        start = perf_counter()
        f = self.forcing(x, r, t)
        self.timings["forcing"] += perf_counter() - start
        force = f["primary_force_N"] + f["secondary_force_N"]
        flow = f["flow_velocity_m_s"]
        w0 = self.wall(x, r)
        xk = x.copy()
        uk = np.zeros_like(x)
        diagnostics = []
        last_delta = None
        history_x = []
        history_f = []
        fixed_c = None
        scratch = None
        prior_wall_normals = None
        prior_wall_gaps = None
        wall_normal_jump = 0.0
        jump_particle_id = None
        for iteration in range(1 if self.method == "S0" else self.max_iterations):
            wk = self.wall(xk, r)
            if prior_wall_normals is not None:
                changes = np.linalg.norm(wk[1] - prior_wall_normals, axis=1)
                relevant = np.minimum(wk[0], prior_wall_gaps) < self.backend.activation * r
                changes = np.where(relevant, changes, 0.0)
                j = int(np.argmax(changes))
                if changes[j] > wall_normal_jump:
                    wall_normal_jump = float(changes[j])
                    jump_particle_id = int(ids[j])
            prior_wall_normals = wk[1].copy()
            prior_wall_gaps = wk[0].copy()
            start = perf_counter()
            a = self.backend.solve(xk, r, force, flow, wk)
            self.timings["assembly"] += perf_counter() - start
            if fixed_c is None:
                fixed_c = self.feature_constraints(x, r, w0, dt, ids, a.velocity)
            c = fixed_c
            if self.wall_contact_mode == "endpoint":
                c = endpoint_wall_constraints(fixed_c, x, xk, r, w0, wk, dt, ids, self.tolerance, a.velocity)
            # Normal gap derivative in the nonlinear SOLVER Jacobian only.
            # Positive correction preserves an SPD QP metric; the fixed point is
            # still audited against the unmodified A(x1) equations below.
            numerical = a
            if self.method == "S1" and self.nonlinear_strategy == "quasi_newton":
                begin = perf_counter()

                def derivative(gap, kind):
                    e = np.asarray(gap) / r[0]
                    v = np.clip(e, 0.001, 0.1)
                    raw = (
                        3 * np.pi * self.backend.mu * (-0.5 / v**2 - 0.45 / v - (3 / 56) * (np.log(v) + 1))
                        if kind == "pair"
                        else 6 * np.pi * self.backend.mu * (-1 / v**2 - 0.2 / v - (1 / 21) * (np.log(v) + 1))
                    )
                    return np.where((e > 0.001) & (e < 0.1), raw, 0.0)

                wi = np.flatnonzero(a.wall_coefficients_N_s_m > 0)
                rows = sparse_rows(a.pairs, a.normals, wi, wk[1][wi], len(r))
                rate = rows @ uk.ravel()
                correction = np.maximum(
                    dt * np.r_[derivative(a.gaps, "pair"), derivative(wk[0][wi], "wall")] * rate, 0.0
                )
                extra = rows.T @ sparse.diags(correction) @ rows
                jac = (a.matrix + extra).tocsc()
                rhs = a.matrix @ a.velocity.ravel() + extra @ uk.ravel()
                factor = splu(jac / a.scale)
                target = factor.solve(rhs / a.scale).reshape(x.shape)
                numerical = replace(a, matrix=jac, velocity=target, factor=factor)
                self.timings["assembly"] += perf_counter() - begin
            start = perf_counter()
            try:
                u, cache, diag = self.backend.constrained_stage(
                    numerical, c, ids, stage, feasible=uk.ravel(), scratch=scratch
                )
                scratch = cache
            finally:
                self.timings["qp"] += perf_counter() - start
                self.calls += 1
            y = x + dt * u
            self.last_candidate = dict(
                start_m=x.copy(),
                end_m=y.copy(),
                velocity_m_s=u.copy(),
                radii_m=r.copy(),
                particle_id=ids.copy(),
                time_s=t,
                dt_s=dt,
                stage=stage,
                method=self.method,
            )
            delta = np.max(np.linalg.norm(u - uk, axis=1), initial=0)
            scale = max(np.max(np.linalg.norm(u, axis=1), initial=0), 1e-12)
            diag.update(linear_residual=a.diagnostics["relative_residual"], nonlinear_delta=delta / scale)
            diagnostics.append(diag)
            if self.method == "S0" or delta / scale <= 1e-7:
                # Re-evaluate the unlagged endpoint equations with the returned U,
                # multipliers and the SAME fixed Stage5A normalization.
                if self.method == "S1":
                    we = self.wall(y, r)
                    beg = perf_counter()
                    ae = self.backend.solve(y, r, force, flow, we)
                    self.timings["assembly"] += perf_counter() - beg
                    ce = c
                    if self.wall_contact_mode == "endpoint":
                        ce = endpoint_wall_constraints(
                            fixed_c, x, y, r, w0, we, dt, ids, self.tolerance, ae.velocity
                        )
                    vs = max(float(np.max(abs(ae.velocity))), 1e-30)
                    an = ae.matrix / ae.scale
                    bn = an @ (ae.velocity.ravel() / vs)
                    lam = np.array([cache.dual.get(key, 0.0) / (ae.scale * vs) for key in ce.ids])
                    k, terms = kkt_residual(an, bn, ce.matrix, ce.lower / vs, u.ravel() / vs, lam)
                    diag["nonlinear_kkt"] = k
                    diag["nonlinear_kkt_terms"] = terms
                    if k > self.nonlinear_rtol:
                        uk = u
                        xk = y
                        continue
                    diag["endpoint_linear_residual"] = ae.diagnostics["relative_residual"]
                    physical_lambda = np.array([cache.dual.get(key, 0.0) for key in ce.ids])
                    physical_rhs = ae.matrix @ ae.velocity.ravel()
                    reaction = ce.matrix.T @ physical_lambda
                    force_residual = ae.matrix @ u.ravel() - physical_rhs - reaction
                    diag["physical_force_residual_l2_N"] = float(np.linalg.norm(force_residual))
                    diag["physical_force_residual_relative"] = float(
                        np.linalg.norm(force_residual)
                        / max(np.linalg.norm(physical_rhs), np.linalg.norm(reaction), 1e-300)
                    )
                    diag["physical_constraint_velocity_error_m_s"] = float(
                        max(0, -np.min(ce.matrix @ u.ravel() - ce.lower, initial=0))
                    )
                else:
                    diag["nonlinear_kkt"] = None
                begin = perf_counter()
                safety = self.sweep(x, y, r, self.tolerance)
                self.timings["geometry"] += perf_counter() - begin
                if not safety["pass_gate"]:
                    raise SubdivideInterval(
                        "SWEEP_NOT_CERTIFIED",
                        dict(
                            dt_s=dt,
                            safety=safety,
                            qp=diag,
                            substep_method=self.method,
                            nonlinear_iterations=iteration + 1,
                        ),
                    )
                return (
                    y,
                    u,
                    cache,
                    dict(
                        substep_method=self.method,
                        resistance_evaluation="endpoint" if self.method == "S1" else "start",
                        qp=diag,
                        iterations=iteration + 1,
                        kkt=max(z["kkt"] for z in diagnostics),
                        linear_residual=max(z["linear_residual"] for z in diagnostics),
                        qp_iterations=sum(z["iterations"] for z in diagnostics),
                        fallbacks=sum(s.get("fallbacks", 0) for z in diagnostics for s in z["solves"]),
                        safety=safety,
                        qp_timings={
                            key: sum(s.get(key, 0.0) for z in diagnostics for s in z["solves"])
                            for key in ["setup_s", "update_s", "solve_s", "polish_s"]
                        },
                        workspace_hits=sum(
                            bool(s.get("workspace_reused", False)) for z in diagnostics for s in z["solves"]
                        ),
                        component_solves=sum(len(z["solves"]) for z in diagnostics),
                    ),
                )
            # Relaxation is ONLY nonlinear iteration, not a force or motion term.
            alpha = self.relaxation
            if last_delta is not None and delta > last_delta * 1.2:
                alpha = min(alpha, 0.5)
            residual = (u - uk).ravel()
            history_x.append(uk.ravel().copy())
            history_f.append(residual)
            if len(history_x) > 7:
                history_x.pop(0)
                history_f.pop(0)
            candidate = (1 - alpha) * uk + alpha * u
            if len(history_x) > 1:
                dx = np.diff(np.array(history_x), axis=0).T
                df = np.diff(np.array(history_f), axis=0).T
                gamma = np.linalg.lstsq(df, residual, rcond=1e-10)[0]
                accelerated = u.ravel() - (dx + df) @ gamma
                # Safeguard ONLY trial nonlinear iterates. Final physical U is
                # always the original-equation certified QP result above.
                if np.isfinite(accelerated).all() and np.linalg.norm(accelerated - uk.ravel()) <= 4 * max(
                    np.linalg.norm(residual), 1e-30
                ):
                    candidate = accelerated.reshape(uk.shape)
            uk = candidate
            xk = x + dt * uk
            last_delta = delta
        end_wall = self.wall(y, r)
        begin = perf_counter()
        end_result = self.backend.solve(y, r, force, flow, end_wall)
        self.timings["assembly"] += perf_counter() - begin
        end_c = (
            fixed_c
            if self.wall_contact_mode == "start"
            else endpoint_wall_constraints(
                fixed_c, x, y, r, w0, end_wall, dt, ids, self.tolerance, end_result.velocity
            )
        )
        vs = max(float(np.max(abs(end_result.velocity))), 1e-30)
        an = end_result.matrix / end_result.scale
        bn = an @ (end_result.velocity.ravel() / vs)
        lam = np.array([cache.dual.get(key, 0.0) / (end_result.scale * vs) for key in end_c.ids])
        terminal_kkt, terminal_terms = kkt_residual(
            an, bn, end_c.matrix, end_c.lower / vs, u.ravel() / vs, lam
        )
        terminal_force_residual = (
            end_result.matrix @ u.ravel()
            - end_result.matrix @ end_result.velocity.ravel()
            - end_c.matrix.T @ (lam * end_result.scale * vs)
        )
        raise SubdivideInterval(
            "NONLINEAR_NOT_CONVERGED",
            dict(
                dt_s=dt,
                iterations=len(diagnostics),
                last_delta=diagnostics[-1]["nonlinear_delta"],
                wall_normal_jump=wall_normal_jump,
                jump_particle_id=jump_particle_id,
                terminal_original_nonlinear_kkt=terminal_kkt,
                terminal_kkt_terms=terminal_terms,
                terminal_physical_force_residual_l2_N=float(np.linalg.norm(terminal_force_residual)),
                last_qp_kkt=diagnostics[-1]["kkt"],
            ),
        )
