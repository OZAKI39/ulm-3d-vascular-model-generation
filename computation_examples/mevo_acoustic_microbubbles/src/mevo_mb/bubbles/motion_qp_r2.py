"""Native, bounded QP adapters; certificate always uses the unmodified input.

min .5 u'A u-b'u, G u>=l, project convention lambda>=0 and Au-b-G'lambda=0.
OSQP and ProxQP lower-bound multipliers are negated; Clarabel (-G) uses +z.
No original inequality is removed and no one-shot guessed active set is used.
"""

from time import perf_counter
import hashlib
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu
from .constrained_motion import kkt_residual, MotionSolverFailure


def certificate(a, b, g, lower, u, lam, *, velocity_scale=None, force_scale=None, dt=None):
    if not all(np.isfinite(v).all() for v in [a.data, b, g.data, lower, u, lam]):
        return dict(pass_gate=False, reason="NONFINITE_INPUT", kkt=None)
    k, terms = kkt_residual(a, b, g, lower, u, lam)
    residual = a @ u - b - g.T @ lam
    slack = g @ u - lower
    objective = float(0.5 * u @ (a @ u) - b @ u)
    # Exact dual objective, using SPD resistance solve, no dense inverse.
    dual_rhs = b + g.T @ lam
    dual_u = splu(a.tocsc()).solve(dual_rhs)
    dual_objective = float(-0.5 * dual_rhs @ dual_u + lower @ lam)
    return dict(
        pass_gate=bool(k <= 1e-9),
        kkt=k,
        kkt_terms=terms,
        objective=objective,
        dual_objective=dual_objective,
        primal_dual_gap=objective - dual_objective,
        normalized_force_residual_l2=float(np.linalg.norm(residual)),
        force_residual_l2_N=float(np.linalg.norm(residual) * force_scale)
        if force_scale is not None
        else None,
        velocity_feasibility_m_s=float(max(0, -slack.min(initial=0)) * velocity_scale)
        if velocity_scale is not None
        else None,
        dt_geometry_error_bound_m=float(max(0, -slack.min(initial=0)) * velocity_scale * dt)
        if velocity_scale is not None and dt is not None
        else None,
    )


def pattern(a, g):
    return hashlib.sha256(
        np.array([*a.shape, *g.shape]).tobytes()
        + a.indptr.tobytes()
        + a.indices.tobytes()
        + g.indptr.tobytes()
        + g.indices.tobytes()
    ).hexdigest()


class NativeQP:
    """One solver/stage/component/pattern workspace; numeric values always updated.

    ProxQP dense defaults to <=300 unknowns; the separately benchmarked bounded
    dense-global option allows <=600. Both cap rows at 1200 and input at 8 MiB.
    Larger systems use its sparse primal API. Calls have bounded iteration counts;
    the external watchdog additionally bounds native-call wall time.
    """

    def __init__(self, name="prox", eps=1e-12):
        self.name, self.eps = name, eps
        self.workspace = None
        self.key = None

    def candidate(self, a, b, g, lower, warm_primal=None, warm_dual=None):
        start = perf_counter()
        a, g = a.tocsc(), g.tocsc()
        n = len(b)
        m = len(lower)
        key = pattern(a, g)
        reused = key == self.key and self.workspace is not None
        setup = update = 0.0
        if self.name.startswith("prox"):
            import proxsuite as ps

            cap = 600 if self.name == "prox_dense_bounded" else 300
            dense = self.name != "prox_sparse" and n <= cap and m <= 1200 and 8 * (n * n + n * m) <= 8 * 2**20
            aa, gg = (a.toarray(), g.toarray()) if dense else (a, g)
            if not reused:
                q = (ps.proxqp.dense if dense else ps.proxqp.sparse).QP(n, 0, m)
                q.settings.eps_abs = self.eps
                q.settings.eps_rel = 0
                q.settings.max_iter = 1000
                q.settings.max_iter_in = 1500
                q.settings.nb_iterative_refinement = 10
                q.settings.compute_timings = True
                q.init(aa, -b, None, None, gg, lower, np.full(m, np.inf))
                self.workspace = q
                self.key = key
                setup = perf_counter() - start
            else:
                q = self.workspace
                q.update(H=aa, g=-b, C=gg, l=lower, u=np.full(m, np.inf), update_preconditioner=True)
                update = perf_counter() - start
            if warm_primal is None:
                q.settings.initial_guess = ps.proxqp.InitialGuess.EQUALITY_CONSTRAINED_INITIAL_GUESS
                q.solve()
            else:
                q.settings.initial_guess = ps.proxqp.InitialGuess.WARM_START
                q.solve(warm_primal, np.empty(0), -warm_dual if warm_dual is not None else np.zeros(m))
            u, lam = q.results.x.copy(), -q.results.z.copy()
            status = str(q.results.info.status)
            iterations = int(q.results.info.iter)
            version = ps.__version__
            backend = "dense" if dense else "sparse"
            native_time = q.results.info.solve_time / 1e6
        elif self.name == "osqp":
            import osqp

            p = sparse.triu(a, format="csc")
            if not reused:
                q = osqp.OSQP()
                q.setup(
                    P=p,
                    q=-b,
                    A=g,
                    l=lower,
                    u=np.full(m, np.inf),
                    verbose=False,
                    eps_abs=self.eps,
                    eps_rel=0,
                    max_iter=300000,
                    polishing=True,
                    polish_refine_iter=10,
                    time_limit=5.0,
                    adaptive_rho_interval=50,
                )
                self.workspace = q
                self.key = key
                setup = perf_counter() - start
            else:
                q = self.workspace
                q.update(Px=p.data, q=-b, Ax=g.data, l=lower, u=np.full(m, np.inf))
                update = perf_counter() - start
            q.warm_start(
                x=np.zeros(n) if warm_primal is None else warm_primal,
                y=np.zeros(m) if warm_dual is None else -warm_dual,
            )
            result = q.solve(raise_error=False)
            u, lam = result.x, -result.y
            status = result.info.status
            iterations = int(result.info.iter)
            version = osqp.__version__
            backend = "sparse"
            native_time = result.info.solve_time
        elif self.name == "clarabel":
            import clarabel

            settings = clarabel.DefaultSettings()
            settings.verbose = False
            settings.max_iter = 300
            settings.tol_gap_abs = self.eps
            settings.tol_gap_rel = self.eps
            settings.tol_feas = self.eps
            settings.time_limit = 5.0
            # This installed API has no user-provided primal/dual warm start.
            q = clarabel.DefaultSolver(
                sparse.triu(a, format="csc"), -b, -g, -lower, [clarabel.NonnegativeConeT(m)], settings
            )
            setup = perf_counter() - start
            result = q.solve()
            u, lam = np.array(result.x), np.array(result.z)
            status = str(result.status)
            iterations = int(result.iterations)
            version = clarabel.__version__
            backend = "sparse"
            native_time = result.solve_time
            reused = False
        else:
            raise ValueError(self.name)
        elapsed = perf_counter() - start
        audit = certificate(a, b, g, lower, u, lam)
        d = dict(
            audit,
            solver=self.name,
            version=version,
            backend=backend,
            status=status,
            setup_s=setup,
            update_s=update,
            solve_s=elapsed - setup - update,
            polish_s=0.0,
            native_solve_s=native_time,
            iterations=iterations,
            workspace_reused=reused,
            warm_start=warm_primal is not None and self.name != "clarabel",
            wall_s=elapsed,
            factorization_s=None,
            refinement_s=None,
            timing_note="factorization/refinement included in native solve; API does not expose separate costs",
        )
        return u, lam, d

    def solve(self, a, b, g, lower, *, warm_primal=None, warm_dual=None, feasible=None):
        a, g = a.tocsc(), g.tocsc()
        if not all(np.isfinite(v).all() for v in [a.data, b, g.data, lower]):
            raise MotionSolverFailure(dict(reason="NONFINITE_INPUT"))
        start = perf_counter()
        try:
            free = splu(a).solve(b)
        except RuntimeError as exc:
            raise MotionSolverFailure(dict(reason="INVALID_MODEL", native_error=str(exc))) from exc
        if np.min(g @ free - lower, initial=0) >= 0:
            lam = np.zeros(len(lower))
            d = certificate(a, b, g, lower, free, lam)
            if d["pass_gate"]:
                return (
                    free,
                    lam,
                    dict(
                        d,
                        setup_s=perf_counter() - start,
                        update_s=0.0,
                        solve_s=0.0,
                        polish_s=0.0,
                        iterations=0,
                        solver="certified_free",
                        attempts=[],
                        workspace_reused=False,
                    ),
                )
        attempts = []
        wp = (
            feasible
            if feasible is not None and np.min(g @ feasible - lower, initial=0) >= -1e-12
            else warm_primal
        )
        candidates = []
        for solver in [self, NativeQP("osqp", 1e-12), NativeQP("prox_sparse", 1e-13)]:
            u, lam, d = solver.candidate(a, b, g, lower, wp, warm_dual)
            attempts.append(d)
            candidates.append((u, lam))
            if "nonconvex" in d["status"].lower().replace(" ", ""):
                raise MotionSolverFailure(dict(reason="INVALID_MODEL", attempts=attempts))
            solved = d["status"] in ["QPSolverOutput.PROXQP_SOLVED", "solved", "Solved"]
            if solved and d["pass_gate"]:
                return (
                    u,
                    lam,
                    dict(
                        d,
                        attempts=attempts,
                        fallbacks=len(attempts) - 1,
                        setup_s=sum(z["setup_s"] for z in attempts),
                        update_s=sum(z["update_s"] for z in attempts),
                        solve_s=sum(z["solve_s"] for z in attempts),
                        iterations=sum(z["iterations"] for z in attempts),
                        total_s=perf_counter() - start,
                    ),
                )
            # Inaccurate candidates are retained, never promoted by status alone.
            # The alternative mature solver supplies the final certified candidate.
        infeasibility = None
        reason = "QP_CERTIFICATE_REJECTED"
        if any("infeasible" in d["status"].lower() for d in attempts):
            known = [("zero", np.zeros(len(b))), ("free", free), ("trial", feasible), ("warm", warm_primal)]
            witnesses = [
                {"source": name, "minimum_slack": float(np.min(g @ v - lower, initial=0))}
                for name, v in known
                if v is not None and np.isfinite(v).all()
            ]
            feasible_witness = next((w for w in witnesses if w["minimum_slack"] >= -1e-12), None)
            infeasibility = dict(known_candidates=witnesses, feasible_witness=feasible_witness)
            reason = "QP_INFEASIBILITY_CONTRADICTED" if feasible_witness else "QP_INFEASIBILITY_UNRESOLVED"
        fingerprint = hashlib.sha256(
            a.data.tobytes() + g.data.tobytes() + b.tobytes() + lower.tobytes() + pattern(a, g).encode()
        ).hexdigest()
        exc = MotionSolverFailure(
            dict(
                reason=reason,
                attempts=attempts,
                wall_s=perf_counter() - start,
                infeasibility_audit=infeasibility,
                input_sha256=fingerprint,
            )
        )
        exc.qp_state = dict(
            A_data=a.data,
            A_indices=a.indices,
            A_indptr=a.indptr,
            A_shape=a.shape,
            G_data=g.data,
            G_indices=g.indices,
            G_indptr=g.indptr,
            G_shape=g.shape,
            b=b,
            lower=lower,
            candidates_U=np.stack([v[0] for v in candidates]),
            candidates_lambda=np.stack([v[1] for v in candidates]),
        )
        raise exc
