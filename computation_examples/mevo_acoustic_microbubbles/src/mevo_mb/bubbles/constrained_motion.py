"""Stage5A opt-in sparse, frictionless motion; frozen Stage4L resistance law.

No acoustic/flow/response physics is defined in this module. SI throughout.
The active-set QP uses sparse saddle systems, never a dense contact Gram matrix.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
from time import perf_counter
import warnings
from concurrent.futures import ThreadPoolExecutor

import numpy as np
from scipy import sparse
from scipy.optimize import minimize
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import splu, lsmr, MatrixRankWarning
from scipy.spatial.distance import pdist

from .translation_backends import (
    LubricationCoupledTranslation,
    LegacyStokesTranslation,
    TranslationResult,
    _pairs,
    _gap_vectors,
)

METHOD = "STAGE5A_SPARSE_ACTIVE_SET_PC_V1"
BACKEND = "constrained_lubrication_v1"


class MotionSolverFailure(RuntimeError):
    def __init__(self, diagnostics):
        self.diagnostics = diagnostics
        super().__init__(str(diagnostics))


class SubdivideInterval(RuntimeError):
    def __init__(self, reason, diagnostics):
        self.reason, self.diagnostics = reason, diagnostics
        super().__init__(reason)


@dataclass
class AcceptedCache:
    particle_ids: tuple = ()
    primal: dict = field(default_factory=dict)
    dual: dict = field(default_factory=dict)  # dimensional N, canonical constraint IDs
    active: tuple = ()
    sparsity_pattern: str = ""


@dataclass
class Constraints:
    matrix: sparse.csc_matrix
    lower: np.ndarray
    ids: tuple
    pairs: np.ndarray


def sparse_rows(pairs, normals, wall_ids, wall_normals, n):
    """<=6 pair entries and <=3 wall entries per row; zero entries eliminated."""
    pairs = np.asarray(pairs, int).reshape(-1, 2)
    m = len(pairs)
    k = len(wall_ids)
    rows = np.r_[np.repeat(np.arange(m), 6), np.repeat(m + np.arange(k), 3)]
    columns = np.r_[
        (pairs[:, :, None] * 3 + np.arange(3)).ravel(),
        (np.asarray(wall_ids)[:, None] * 3 + np.arange(3)).ravel(),
    ]
    values = np.r_[np.stack([-normals, normals], axis=1).ravel(), -wall_normals.ravel()]
    out = sparse.coo_matrix((values, (rows, columns)), shape=(m + k, 3 * n)).tocsc()
    out.eliminate_zeros()
    return out


def horizon_constraints(x, r, wall, velocity, dt, tolerance, ids, stage="start"):
    radius = float(r[0])
    centered = velocity - velocity.mean(axis=0)
    reach = 2 * np.linalg.norm(centered, axis=1).max(initial=0) * dt
    pairs = _pairs(x, (2 + 0.1) * radius + reach)
    gap, normal = _gap_vectors(x, pairs, radius)
    wg, wn, valid = wall
    if not np.asarray(valid).all():
        raise MotionSolverFailure({"reason": "INVALID_WALL"})
    wi = np.flatnonzero(wg <= np.maximum(0.1 * radius, np.linalg.norm(velocity, axis=1) * dt))
    gaps = np.r_[gap, wg[wi]]
    lower = -np.maximum(gaps, 0) / dt
    lower[gaps <= tolerance] = 0
    keys = tuple((stage, "pair", int(ids[i]), int(ids[j])) for i, j in pairs)
    keys += tuple((stage, "wall", int(ids[i])) for i in wi)
    return Constraints(sparse_rows(pairs, normal, wi, wn[wi], len(r)), lower, keys, pairs)


def kkt_residual(a, b, g, lower, u, lam):
    slack = g @ u - lower
    station = float(
        np.linalg.norm(a @ u - b - g.T @ lam) / max(1.0, np.linalg.norm(b), np.linalg.norm(g.T @ lam))
    )
    feasibility = float(max(0.0, -np.min(slack, initial=0.0)))
    dual = float(max(0.0, -np.min(lam, initial=0.0)) / max(1.0, np.max(np.abs(lam), initial=0.0)))
    positive = lam > 1e-10 * max(1.0, np.max(np.abs(lam), initial=0.0))
    complementary = float(np.max(np.abs(slack[positive]), initial=0.0))
    product = float(np.max(np.abs(slack * lam), initial=0.0) / max(1.0, np.max(np.abs(lam), initial=0.0)))
    return max(station, feasibility, dual, complementary, product), dict(
        stationarity=station,
        feasibility=feasibility,
        dual_feasibility=dual,
        active_slack=complementary,
        complementarity=product,
    )


class SparseActiveSetQP:
    """Matrix-free bounded dual seed + deterministic sparse KKT active set.

    Reuses the factorized resistance throughout a solve. Accepted primal/dual
    caches are external: rejected trials never change their contents. Sparse LU
    of the active saddle system is polished, with sparse least squares only for
    dependent equalities. Every original inequality remains in the KKT audit.
    """

    def __init__(self):
        self.workspace_key = None
        self.workspace_factor = None

    def solve(self, a, b, g, lower, *, warm_dual=None, warm_primal=None, feasible=None):
        started = perf_counter()
        a = a.tocsc()
        g = g.tocsc()
        m = g.shape[0]
        pattern = hashlib.sha256(a.indptr.tobytes() + a.indices.tobytes() + a.data.tobytes()).hexdigest()
        reused = pattern == self.workspace_key
        if not reused:
            self.workspace_key = pattern
            self.workspace_factor = splu(a)
        lu = self.workspace_factor
        u0 = lu.solve(b)
        setup = perf_counter() - started
        diag = dict(
            solver="SCIPY_MATRIX_FREE_DUAL_SPARSE_ACTIVE_SET",
            setup_s=setup,
            solve_s=0.0,
            polish_s=0.0,
            iterations=0,
            cold_start=warm_dual is None,
            primal_warm_start=warm_primal is not None,
            workspace_reuse="resistance LU across dual iterations and identical matrices",
            workspace_reused=reused,
            matrix_dimension=len(b),
            constraint_count=m,
            A_nnz=a.nnz,
            G_nnz=g.nnz,
        )
        if warm_primal is not None and warm_dual is not None:
            res, terms = kkt_residual(a, b, g, lower, warm_primal, warm_dual)
            if res <= 1e-11:
                return (
                    warm_primal.copy(),
                    warm_dual.copy(),
                    dict(
                        diag,
                        kkt=res,
                        kkt_terms=terms,
                        warm_solution_certified=True,
                        objective=float(0.5 * warm_primal @ (a @ warm_primal) - b @ warm_primal),
                    ),
                )
        if not m or np.min(g @ u0 - lower, initial=0) >= -1e-12:
            lam = np.zeros(m)
            res, terms = kkt_residual(a, b, g, lower, u0, lam)
            return (
                u0,
                lam,
                dict(diag, kkt=res, kkt_terms=terms, objective=float(0.5 * u0 @ (a @ u0) - b @ u0)),
            )
        # Diagonal Jacobi scaling only, no dense inverse or contact matrix.
        d = np.sqrt(np.asarray(g.power(2) @ (1 / a.diagonal())).ravel())
        c = (sparse.diags(1 / d) @ g).tocsc()
        bound = lower / d
        q = c @ u0 - bound

        def value(z):
            change = lu.solve(c.T @ z)
            grad = c @ change + q
            return 0.5 * z @ (grad + q), grad

        seed = np.zeros(m) if warm_dual is None else np.maximum(warm_dual * d, 0)
        solve_start = perf_counter()
        fit = minimize(
            value,
            seed,
            method="L-BFGS-B",
            jac=True,
            bounds=[(0, None)] * m,
            options=dict(ftol=1e-15, gtol=1e-12, maxiter=1500, maxls=50),
        )
        z = np.maximum(fit.x, 0)
        active = z > 1e-10 * max(1.0, z.max(initial=0))
        diag.update(
            iterations=int(fit.nit), dual_seed_message=str(fit.message), solve_s=perf_counter() - solve_start
        )
        polish_start = perf_counter()
        seen = set()
        sparse_lstsq = 0
        for iteration in range(4 * m + 20):
            key = active.tobytes()
            if key in seen:
                break
            seen.add(key)
            ix = np.flatnonzero(active)
            ca = c[ix]
            saddle = sparse.bmat([[a, -ca.T], [-ca, None]], format="csc")
            rhs = np.r_[b, -bound[ix]]
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("error", MatrixRankWarning)
                    solution = splu(saddle).solve(rhs)
                if not np.isfinite(solution).all() or np.linalg.norm(saddle @ solution - rhs) > 1e-10 * max(
                    1.0, np.linalg.norm(rhs)
                ):
                    raise RuntimeError("ill-conditioned saddle solve")
            except (RuntimeError, MatrixRankWarning):
                solution = lsmr(saddle, rhs, atol=2e-15, btol=2e-15, conlim=1e14, maxiter=20 * len(rhs))[0]
                sparse_lstsq += 1
            u = solution[: len(b)]
            z = np.zeros(m)
            z[ix] = solution[len(b) :]
            slack = c @ u - bound
            # Dependent equalities may be inconsistent; release a strictly slack
            # equality, but keep the inequality in all subsequent audits.
            positive = np.where(active, slack, -np.inf)
            if np.max(positive, initial=-np.inf) > 1e-10:
                active[np.argmax(positive)] = False
                continue
            if np.min(z, initial=0) < -1e-11 * max(1.0, np.max(np.abs(z))):
                active[np.argmin(z)] = False
                continue
            violation = np.where(~active, slack, np.inf)
            if np.min(violation, initial=np.inf) < -1e-11:
                active[np.argmin(violation)] = True
                continue
            lam = np.maximum(z, 0) / d
            residual, terms = kkt_residual(a, b, g, lower, u, lam)
            if residual <= 1e-9 and np.isfinite(u).all():
                diag.update(
                    polish_s=perf_counter() - polish_start,
                    active_iterations=iteration + 1,
                    sparse_lstsq_calls=sparse_lstsq,
                    kkt=residual,
                    kkt_terms=terms,
                    objective=float(0.5 * u @ (a @ u) - b @ u),
                )
                return u, lam, diag
            break
        # Deterministic feasible primal active set. A feasible prior predictor
        # handles positive corrector bounds; zero works for all start bounds.
        v = np.zeros(len(b)) if feasible is None else np.array(feasible, copy=True)
        if warm_primal is not None and np.min(g @ warm_primal - lower, initial=0) >= -1e-12:
            v = warm_primal.copy()
        if np.min(g @ v - lower, initial=0) < -1e-10:
            raise MotionSolverFailure(dict(diag, reason="NO_FEASIBLE_PRIMAL_SEED"))
        working = []
        for count in range(20 * m + 100):
            cw = c[working]
            gradient = a @ v - b
            saddle = sparse.bmat([[a, -cw.T], [-cw, None]], format="csc")
            rhs = np.r_[-gradient, np.zeros(len(working))]
            try:
                answer = splu(saddle).solve(rhs)
            except RuntimeError:
                answer = lsmr(saddle, rhs, atol=1e-15, btol=1e-15, maxiter=20 * len(rhs))[0]
            direction = answer[: len(b)]
            mult = answer[len(b) :]
            if np.linalg.norm(direction) <= 1e-11 * max(1.0, np.linalg.norm(v)):
                if len(mult) and np.min(mult) < -1e-11 * max(1.0, np.max(abs(mult))):
                    working.pop(int(np.argmin(mult)))
                    continue
                z = np.zeros(m)
                z[working] = np.maximum(mult, 0)
                lam = z / d
                residual, terms = kkt_residual(a, b, g, lower, v, lam)
                if residual <= 1e-9 and np.isfinite(v).all():
                    diag.update(
                        polish_s=perf_counter() - polish_start,
                        active_iterations=count + 1,
                        primal_fallback=True,
                        kkt=residual,
                        kkt_terms=terms,
                        objective=float(0.5 * v @ (a @ v) - b @ v),
                    )
                    return v, lam, diag
                break
            slack = c @ v - bound
            closing = c @ direction
            eligible = closing < -1e-12 * max(1.0, np.linalg.norm(direction))
            eligible[working] = False
            ratios = np.full(m, np.inf)
            ratios[eligible] = np.maximum(slack[eligible], 0) / (-closing[eligible])
            hit = int(np.argmin(ratios))
            alpha = min(1.0, float(ratios[hit]))
            v += alpha * direction
            if alpha < 1.0:
                working.append(hit)
        residual, terms = kkt_residual(a, b, g, lower, v, np.zeros(m))
        raise MotionSolverFailure(
            dict(diag, reason="SPARSE_CONTACT_KKT_FAILURE", kkt=residual, kkt_terms=terms)
        )


class ConstrainedLubrication(LubricationCoupledTranslation):
    name = BACKEND

    def __init__(
        self,
        viscosity_Pa_s,
        *,
        gap_limiter_mode="constrained_adaptive",
        components=True,
        warm_start=True,
        component_workers=1,
        qp_solver="sparse_active_set",
        **kw,
    ):
        super().__init__(viscosity_Pa_s, **kw)
        if gap_limiter_mode not in ["legacy_half_gap", "constrained_adaptive"]:
            raise ValueError("Unknown gap limiter")
        self.gap_limiter_mode = gap_limiter_mode
        self.components = components
        self.warm_start = warm_start
        self.cache = AcceptedCache()
        if qp_solver == "sparse_active_set":
            self.qp_factory = SparseActiveSetQP
            self.method = METHOD
        elif qp_solver == "osqp":
            from .osqp_motion import OSQPMotionQP, METHOD as OSQP_METHOD

            self.qp_factory = OSQPMotionQP
            self.method = OSQP_METHOD
        else:
            raise ValueError(f"Unknown constrained QP solver: {qp_solver}")
        self.qp_solver = qp_solver
        self.qp = self.qp_factory()
        self.component_workers = int(component_workers)
        self.component_workspaces = {}
        self.executor = None

    def solve(self, x, r, force, flow, wall=None):
        """Vectorized CSC assembly; coefficients inherited byte-for-byte."""
        start = perf_counter()
        x = np.asarray(x)
        r = np.asarray(r)
        radius = self.mono_radius(r)
        base = LegacyStokesTranslation.solve(self, x, r, force, flow, wall)
        n = len(r)
        pairs = _pairs(x, (2 + self.activation) * radius)
        gaps, normals = _gap_vectors(x, pairs, radius)
        coeff = self.coefficient(gaps, radius)
        take = coeff > 0
        pairs, gaps, normals, coeff = pairs[take], gaps[take], normals[take], coeff[take]
        if wall is None:
            wall = (np.full(n, np.inf), np.zeros((n, 3)), np.ones(n, bool))
        wg, wn, valid = wall
        if not np.asarray(valid).all():
            raise MotionSolverFailure({"reason": "INVALID_WALL"})
        wc = self.coefficient(wg, radius, "wall")
        wi = np.flatnonzero(wc > 0)
        rows = sparse_rows(pairs, normals, wi, wn[wi], n)
        a = (base.matrix + rows.T @ sparse.diags(np.r_[coeff, wc[wi]]) @ rows).tocsc()
        a.eliminate_zeros()
        asym = a - a.T
        sym = float(np.max(np.abs(asym.data), initial=0))
        if sym > 1e-13 * max(a.diagonal()) or np.any(coeff < 0) or np.any(wc < 0):
            raise MotionSolverFailure({"reason": "SPD_CERTIFICATE_FAILED"})
        rhs = (force + base.drag[:, None] * flow).ravel()
        scale = float(np.median(base.drag))
        beg = perf_counter()
        factor = splu(a / scale) if rows.shape[0] else None
        u = factor.solve(rhs / scale).reshape(n, 3) if factor else base.velocity
        residual = float(np.linalg.norm(a @ u.ravel() - rhs) / max(np.linalg.norm(rhs), np.finfo(float).tiny))
        if not np.isfinite(u).all() or not np.isfinite(residual) or residual > self.residual_tolerance:
            raise MotionSolverFailure(
                {"reason": "LINEAR_RESIDUAL", "residual": residual if np.isfinite(residual) else None}
            )
        diag = dict(
            backend=self.name,
            relative_residual=residual,
            solve_time_s=perf_counter() - beg,
            assembly_and_solve_time_s=perf_counter() - start,
            matrix_size=3 * n,
            nnz=a.nnz,
            symmetry_error=sym,
            spd_certificate="D>0 and nonnegative weighted normal outer products",
            pair_active=len(pairs),
            wall_active=len(wi),
        )
        return TranslationResult(
            u, base.stokes_velocity, a, base.drag, pairs, gaps, normals, wall, diag, factor, scale, coeff, wc
        )

    def constrained(self, result, constraints, ids, *, target=None, feasible=None):
        start = perf_counter()
        ids = np.asarray(ids)
        n = len(ids)
        g = constraints.matrix
        lower = constraints.lower
        free = result.velocity if target is None else np.asarray(target)
        vs = max(float(np.max(np.abs(free))), 1e-30)
        a = result.matrix / result.scale
        b = a @ (free.ravel() / vs)
        # Components include every A and G coupling, not merely hard contacts.
        ap = result.pairs
        gp = constraints.pairs
        pairs = np.vstack([ap, gp])
        graph = sparse.coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(n, n))
        count, labels = connected_components(graph, directed=False)
        sizes = np.bincount(labels)
        if not self.components:
            labels = np.zeros(n, int)
        u = np.empty(3 * n)
        dual = np.zeros(g.shape[0])
        diagnostics = []
        csr = g.tocsr()
        owners = np.full(g.shape[0], -1, int)
        for row in range(g.shape[0]):
            involved = csr.indices[csr.indptr[row] : csr.indptr[row + 1]] // 3
            labs = np.unique(labels[involved])
            if len(labs) != 1:
                raise MotionSolverFailure({"reason": "CROSS_COMPONENT_CONSTRAINT"})
            owners[row] = labs[0]

        def component(label):
            ix = np.flatnonzero(labels == label)
            cols = (3 * ix[:, None] + np.arange(3)).ravel()
            rr = np.flatnonzero(owners == label)
            aa = a[cols][:, cols].tocsc()
            gg = g[rr][:, cols].tocsc()
            keys = [constraints.ids[k] for k in rr]
            have = self.warm_start and bool(self.cache.particle_ids)
            warm = (
                np.array([self.cache.dual.get(key, 0.0) / (result.scale * vs) for key in keys])
                if have
                else None
            )
            wp = (
                np.concatenate([self.cache.primal.get(int(ids[i]), np.zeros(3)) for i in ix]) / vs
                if have
                else None
            )
            vv = None if feasible is None else feasible.ravel()[cols] / vs
            key = tuple(map(int, ids[ix]))
            if key not in self.component_workspaces:
                self.component_workspaces[key] = self.qp_factory()
            solver = self.component_workspaces[key]
            answer, lam, d = solver.solve(
                aa, b[cols], gg, lower[rr] / vs, warm_dual=warm, warm_primal=wp, feasible=vv
            )
            return cols, rr, answer, lam, d

        groups = np.unique(labels)
        if self.component_workers > 1 and len(groups) >= 4:
            if self.executor is None:
                self.executor = ThreadPoolExecutor(max_workers=self.component_workers)
            outputs = list(self.executor.map(component, groups))
        else:
            outputs = [component(k) for k in groups]
        for cols, rr, answer, lam, d in outputs:
            u[cols] = answer
            dual[rr] = lam
            diagnostics.append(d)
        residual, terms = kkt_residual(a, b, g, lower / vs, u, dual)
        if residual > 1e-9 or not np.isfinite(u).all():
            raise MotionSolverFailure({"reason": "GLOBAL_KKT", "kkt": residual})
        pattern = hashlib.sha256(
            a.indptr.tobytes() + a.indices.tobytes() + g.indptr.tobytes() + g.indices.tobytes()
        ).hexdigest()
        candidate = AcceptedCache(
            tuple(map(int, ids)),
            dict(zip(map(int, ids), u.reshape(-1, 3) * vs)),
            dict(zip(constraints.ids, dual * result.scale * vs)),
            tuple(k for k, v in zip(constraints.ids, dual) if v > 1e-10),
            pattern,
        )
        diag = dict(
            kkt=residual,
            kkt_terms=terms,
            qp_seconds=perf_counter() - start,
            component_count=int(count),
            component_sizes=sizes.tolist(),
            largest_component=int(sizes.max()),
            median_component=float(np.median(sizes)),
            isolated_fraction=float(np.sum(sizes == 1) / n),
            largest_component_fraction=float(sizes.max() / n),
            constraint_count=g.shape[0],
            G_nnz=g.nnz,
            global_dimension=3 * n,
            components_used=self.components,
            solves=diagnostics,
            setup_s=sum(d["setup_s"] for d in diagnostics),
            solve_s=sum(d["solve_s"] for d in diagnostics),
            polish_s=sum(d["polish_s"] for d in diagnostics),
            iterations=sum(d["iterations"] for d in diagnostics),
        )
        return u.reshape(-1, 3) * vs, candidate, diag

    def timestep(self, x, r, velocity, wall, remaining):
        if self.gap_limiter_mode == "legacy_half_gap":
            return super().timestep(x, r, velocity, wall, remaining)
        return remaining


class PredictorCorrector:
    """Trial has no accepted-cache mutation. The driver commits accepted children."""

    def __init__(
        self,
        backend,
        geometry,
        forcing,
        *,
        contact_tolerance_m=1e-11,
        position_atol_R=0.002,
        position_rtol=0.01,
        gap_atol_R=0.002,
        gap_rtol=0.02,
    ):
        self.backend, self.geometry, self.forcing = backend, geometry, forcing
        self.tolerance = contact_tolerance_m
        self.position_atol_R, self.position_rtol = position_atol_R, position_rtol
        self.gap_atol_R, self.gap_rtol = gap_atol_R, gap_rtol

    def trial(self, x, r, ids, t, dt):
        start = perf_counter()
        wall_s = 0.0
        force_s = 0.0

        def wall(p):
            nonlocal wall_s
            begin = perf_counter()
            v = self.geometry.wall_gap_and_normal(p, r)
            wall_s += perf_counter() - begin
            return v

        def force(p, at):
            nonlocal force_s
            begin = perf_counter()
            v = self.forcing(p, r, at)
            force_s += perf_counter() - begin
            return v

        f0 = force(x, t)
        w0 = wall(x)
        b = self.backend
        a0 = b.solve(x, r, f0["primary_force_N"] + f0["secondary_force_N"], f0["flow_velocity_m_s"], w0)
        c0 = horizon_constraints(x, r, w0, a0.velocity, dt, self.tolerance, ids)
        up, _, dp = b.constrained(a0, c0, ids)
        xp = x + dt * up
        wp = wall(xp)
        if not np.isfinite(xp).all() or np.min(wp[0]) < -self.tolerance:
            raise SubdivideInterval("PREDICTOR_WALL_GEOMETRY", {"dt_s": dt})
        fp = force(xp, t + dt)
        ap = b.solve(xp, r, fp["primary_force_N"] + fp["secondary_force_N"], fp["flow_velocity_m_s"], wp)
        cp = horizon_constraints(xp, r, wp, ap.velocity, dt, self.tolerance, ids, stage="predictor")
        cp.lower = cp.lower + cp.matrix @ up.ravel()
        # Constraints at both anchors are expressed in the same x_n displacement.
        cc = Constraints(
            sparse.vstack([c0.matrix, cp.matrix], format="csc"),
            np.r_[c0.lower, cp.lower],
            c0.ids + cp.ids,
            np.vstack([c0.pairs, cp.pairs]),
        )
        target = 0.5 * (up + ap.velocity)
        uc, cache, dc = b.constrained(ap, cc, ids, target=target, feasible=up)
        xc = x + dt * uc
        wc = wall(xc)
        pairs = np.unique(np.vstack([c0.pairs, cp.pairs]), axis=0)
        hp, _ = _gap_vectors(xp, pairs, r[0])
        hc, _ = _gap_vectors(xc, pairs, r[0])
        h0, _ = _gap_vectors(x, pairs, r[0])
        finite_wall = np.isfinite(w0[0])
        gap_error = np.r_[abs(hc - hp), abs(wc[0][finite_wall] - wp[0][finite_wall])]
        original_gap = np.r_[h0, w0[0][finite_wall]]
        gap_budget = self.gap_atol_R * r[0] + self.gap_rtol * np.maximum(original_gap, 0.001 * r[0])
        ex = float(np.max(np.linalg.norm(xc - xp, axis=1)))
        xb = self.position_atol_R * r[0] + self.position_rtol * np.max(np.linalg.norm(xc - x, axis=1))
        error = max(ex / xb, float(np.max(gap_error / gap_budget, initial=0)))
        diag = dict(
            method=b.method,
            time_s=t + dt,
            dt_s=dt,
            error_ratio=error,
            position_error_m=ex,
            gap_error_max_m=float(np.max(gap_error, initial=0)),
            predictor=dp,
            corrector=dc,
            linear_residual=max(a0.diagnostics["relative_residual"], ap.diagnostics["relative_residual"]),
            assembly_s=a0.diagnostics["assembly_and_solve_time_s"]
            + ap.diagnostics["assembly_and_solve_time_s"],
            wall_s=wall_s,
            forcing_s=force_s,
            kkt=max(dp["kkt"], dc["kkt"]),
        )
        if error > 1:
            raise SubdivideInterval("LOCAL_ERROR", diag)
        # Independent all-pair straight-segment clearance: endpoint checks alone
        # cannot detect tunnelling when broad-phase velocities were insufficient.
        ii, jj = np.triu_indices(len(x), 1)
        delta = x[jj] - x[ii]
        change = (xc - x)[jj] - (xc - x)[ii]
        denom = np.einsum("ij,ij->i", change, change)
        fraction = np.clip(
            -np.einsum("ij,ij->i", delta, change) / np.maximum(denom, np.finfo(float).tiny), 0, 1
        )
        swept_gap = np.linalg.norm(delta + fraction[:, None] * change, axis=1) - r[ii] - r[jj]
        minimum_pair = float(np.min(pdist(xc) - 2 * r[0], initial=np.inf))
        minimum_wall = float(np.min(wc[0], initial=np.inf))
        if (
            not np.isfinite(xc).all()
            or min(minimum_pair, minimum_wall, float(np.min(swept_gap, initial=np.inf))) < -self.tolerance
        ):
            raise SubdivideInterval("FINAL_OR_SWEPT_GEOMETRY", diag)
        # No projection is needed for a tolerance-feasible state. Retain explicit
        # zero diagnostics so speedup cannot hide geometric corrections.
        diag.update(
            minimum_pair_gap_m=minimum_pair if np.isfinite(minimum_pair) else None,
            minimum_wall_gap_m=minimum_wall if np.isfinite(minimum_wall) else None,
            projection_calls=0,
            projection_corrected_calls=0,
            projection_max_m=0.0,
            projection_rms_m=0.0,
            multi_radius_corrections=0,
            wall_s=wall_s,
            total_s=perf_counter() - start,
        )
        return xc, uc, cache, diag

    def interval(self, x, r, ids, t, end, *, depth=0, on_accept=None, on_reject=None):
        try:
            y, u, cache, diag = self.trial(x, r, ids, t, end - t)
        except SubdivideInterval as exc:
            if on_reject:
                on_reject(
                    dict(time_s=t, end_s=end, depth=depth, reason=exc.reason, diagnostics=exc.diagnostics)
                )
            mid = t + (end - t) / 2
            if depth >= 24 or not t < mid < end:
                raise MotionSolverFailure(
                    {"reason": "PHYSICAL_TIME_SUBDIVISION_LIMIT", "time_s": t, "end_s": end}
                ) from exc
            y = self.interval(x, r, ids, t, mid, depth=depth + 1, on_accept=on_accept, on_reject=on_reject)
            return self.interval(
                y, r, ids, mid, end, depth=depth + 1, on_accept=on_accept, on_reject=on_reject
            )
        self.backend.cache = cache
        if on_accept:
            on_accept(y, u, diag)
        return y
