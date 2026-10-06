"""Versioned R2 motion: S0 diagnostic and converged S1 semi-implicit Euler.

S1 V2 freezes F_acoustic, u_flow and contact normals at each substep start.
It solves A(x0+dt U) U = F0 + D u_flow0 + J(x0)' lambda.
The single start tangent set guarantees pair clearance along the adopted straight
path, including initially negative tolerance-band gaps. End geometry and real WALL
sweeps are separately certified. Historical tangent rows are not stacked.
The resistance is never differentiated as a position energy.
"""

from collections import OrderedDict
from dataclasses import replace
from time import perf_counter
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu
from scipy.spatial.distance import pdist
from .constrained_motion import (
    ConstrainedLubrication,
    Constraints,
    sparse_rows,
    MotionSolverFailure,
    SubdivideInterval,
    kkt_residual,
)
from .translation_backends import _pairs, _gap_vectors
from .motion_qp_r2 import NativeQP

BACKEND = "constrained_lubrication_r2"


class MotionR2(ConstrainedLubrication):
    name = BACKEND

    def __init__(self, viscosity_Pa_s, *, solver="prox", **kwargs):
        super().__init__(viscosity_Pa_s, **kwargs)
        self.qp_factory = lambda: NativeQP(solver)
        self.qp = self.qp_factory()
        self.qp_solver = solver
        self.method = "STAGE5A_R2_S1_SEMI_IMPLICIT_EULER_V3"
        self.stage_workspaces = OrderedDict()

    def constrained_stage(self, result, constraints, ids, stage, feasible=None, scratch=None):
        if stage not in self.stage_workspaces:
            self.stage_workspaces[stage] = OrderedDict()
        self.stage_workspaces.move_to_end(stage)
        self.component_workspaces = self.stage_workspaces[stage]
        accepted = self.cache
        try:
            if scratch is not None:
                self.cache = scratch
            answer = super().constrained(result, constraints, ids, feasible=feasible)
        finally:
            self.cache = accepted
        u, cache, d = answer
        reaction = constraints.matrix.T @ np.array([cache.dual[k] for k in constraints.ids])
        rhs = result.matrix @ result.velocity.ravel()
        residual = result.matrix @ u.ravel() - rhs - reaction
        d["physical_force_residual_l2_N"] = float(np.linalg.norm(residual))
        d["physical_force_residual_relative"] = float(
            np.linalg.norm(residual) / max(np.linalg.norm(rhs), np.linalg.norm(reaction), 1e-300)
        )
        d["physical_constraint_velocity_error_m_s"] = float(
            max(0, -np.min(constraints.matrix @ u.ravel() - constraints.lower, initial=0))
        )
        while len(self.component_workspaces) > 128:
            self.component_workspaces.popitem(last=False)
        while len(self.stage_workspaces) > 12:
            self.stage_workspaces.popitem(last=False)
        return answer


def anchored_constraints(x0, xk, r, wall0, wallk, dt, ids, tolerance, reach_velocity):
    """Every pair potentially reached by the known candidate plus activation band.

    Broad phase is augmented by both endpoints; final all-pair path check rejects
    any omission. Initial tolerance-band gaps are retained, never clamped in data.
    """
    radius = r[0]
    reach = 2 * np.max(np.linalg.norm(reach_velocity - reach_velocity.mean(0), axis=1), initial=0) * dt
    pairs = np.unique(np.vstack([_pairs(x0, 2.1 * radius + reach), _pairs(xk, 2.1 * radius)]), axis=0)
    hk, nk = _gap_vectors(xk, pairs, radius)
    h0, _ = _gap_vectors(x0, pairs, radius)
    wg0, _, valid0 = wall0
    wg, nw, valid = wallk
    if not valid.all() or not valid0.all():
        raise MotionSolverFailure(dict(reason="GEOMETRY_DATA_CORRUPTED"))
    wi = np.flatnonzero(
        np.minimum(wg, wg0) <= np.maximum(0.1 * radius, np.linalg.norm(reach_velocity, axis=1) * dt)
    )
    g = sparse_rows(pairs, nk, wi, nw[wi], len(r))
    initial = np.r_[h0, wg0[wi]]
    current = np.r_[hk, wg[wi]]
    allowed = np.where(initial <= tolerance, initial, 0.0)
    lower = (allowed - current + g @ (xk - x0).ravel()) / dt
    keys = tuple(("pair", int(ids[i]), int(ids[j])) for i, j in pairs)
    keys += tuple(("wall", int(ids[i])) for i in wi)
    return Constraints(g, lower, keys, pairs)


class SemiImplicit:
    def __init__(
        self,
        backend,
        geometry,
        forcing,
        sweep,
        *,
        method="S1",
        tolerance=1e-11,
        max_iterations=32,
        nonlinear_rtol=1e-9,
        atol_R=0.002,
        rtol=0.01,
        gap_rtol=0.02,
        relaxation=1.0,
        nonlinear_strategy="picard",
    ):
        self.backend, self.geometry, self.forcing, self.sweep = backend, geometry, forcing, sweep
        self.method, self.tolerance = method, tolerance
        self.max_iterations, self.nonlinear_rtol = max_iterations, nonlinear_rtol
        self.atol_R, self.rtol, self.gap_rtol = atol_R, rtol, gap_rtol
        self.relaxation = relaxation
        self.nonlinear_strategy = nonlinear_strategy
        self.timings = {k: 0.0 for k in ["forcing", "wall", "assembly", "qp", "geometry", "error"]}
        self.calls = 0

    def wall(self, x, r):
        start = perf_counter()
        w = self.geometry.wall_gap_and_normal(x, r)
        self.timings["wall"] += perf_counter() - start
        return w

    def advance(self, x, r, ids, t, dt, stage):
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
        for iteration in range(1 if self.method == "S0" else self.max_iterations):
            wk = self.wall(xk, r)
            start = perf_counter()
            a = self.backend.solve(xk, r, force, flow, wk)
            self.timings["assembly"] += perf_counter() - start
            if fixed_c is None:
                fixed_c = anchored_constraints(x, x, r, w0, w0, dt, ids, self.tolerance, a.velocity)
            c = fixed_c
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
                else:
                    diag["nonlinear_kkt"] = None
                begin = perf_counter()
                safety = self.sweep(x, y, r, self.tolerance)
                self.timings["geometry"] += perf_counter() - begin
                if not safety["pass_gate"]:
                    raise SubdivideInterval("SWEEP_NOT_CERTIFIED", dict(dt_s=dt, safety=safety))
                return (
                    y,
                    u,
                    cache,
                    dict(
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
        raise SubdivideInterval(
            "NONLINEAR_NOT_CONVERGED",
            dict(dt_s=dt, iterations=len(diagnostics), last_delta=diagnostics[-1]["nonlinear_delta"]),
        )

    def trial(self, x, r, ids, t, dt):
        full, _, _, df = self.advance(x, r, ids, t, dt, "full")
        half, uh, ch, dh = self.advance(x, r, ids, t, dt / 2, "half1")
        fine, u, c, d = self.advance(half, r, ids, t + dt / 2, dt / 2, "half2")
        begin = perf_counter()
        diff = fine - full
        position = float(np.max(np.linalg.norm(diff, axis=1), initial=0))
        relative = float(np.max(np.linalg.norm(diff - diff.mean(0), axis=1), initial=0))
        displacement = np.linalg.norm(fine - x, axis=1).max(initial=0)
        posbudget = self.atol_R * r[0] + self.rtol * displacement
        h0 = pdist(x) - 2 * r[0]
        hf = pdist(full) - 2 * r[0]
        hh = pdist(fine) - 2 * r[0]
        pairratio = abs(hh - hf) / (self.atol_R * r[0] + self.gap_rtol * np.maximum(h0, 0.001 * r[0]))
        w0 = self.wall(x, r)[0]
        wf = self.wall(full, r)[0]
        wh = self.wall(fine, r)[0]
        valid = np.isfinite(w0)
        wallratio = abs(wf[valid] - wh[valid]) / (
            self.atol_R * r[0] + self.gap_rtol * np.maximum(w0[valid], 0.001 * r[0])
        )
        ratios = dict(
            position=position / posbudget,
            relative=relative / posbudget,
            pair=float(pairratio.max(initial=0)),
            wall=float(wallratio.max(initial=0)),
        )
        error = max(ratios.values())
        self.timings["error"] += perf_counter() - begin
        diag = dict(
            error_ratio=error,
            error_terms=ratios,
            position_error_m=position,
            kkt=max(q["kkt"] for q in [df, dh, d]),
            linear_residual=max(q["linear_residual"] for q in [df, dh, d]),
            nonlinear_iterations=sum(q["iterations"] for q in [df, dh, d]),
            qp_iterations=sum(q["qp_iterations"] for q in [df, dh, d]),
            fallbacks=sum(q["fallbacks"] for q in [df, dh, d]),
        )
        if error > 1:
            raise SubdivideInterval("LOCAL_ERROR", dict(diag, dt_s=dt))
        # Actual adopted path consists of TWO certified segments, both recorded.
        return [(half, uh, ch, dh), (fine, u, c, d)], diag
