"""Optional CPU OSQP adapter. Sparse original QP, unchanged acceptance audit."""

from time import perf_counter
import hashlib
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu

from .constrained_motion import kkt_residual, MotionSolverFailure

METHOD = "STAGE5A_OSQP_PC_V2"
SETTINGS = dict(
    verbose=False,
    eps_abs=1e-10,
    eps_rel=0.0,
    max_iter=100000,
    polishing=True,
    polish_refine_iter=10,
    adaptive_rho_interval=50,
    check_termination=25,
    warm_starting=True,
)


class OSQPMotionQP:
    """Reuse CSC workspaces; explicitly reset starts to the accepted cache.

    A failed warm solve may be retried once cold. Both attempts and their cost
    are retained. A sparse KKT polish may refine OSQP's candidate; the independent
    1e-9 audit of every original inequality must pass.
    Numerical ADMM regularization does not change the QP objective.
    """

    def __init__(self):
        import osqp

        self.osqp = osqp
        self.workspace = None
        self.pattern = None

    def solve(self, a, b, g, lower, *, warm_dual=None, warm_primal=None, feasible=None):
        begin = perf_counter()
        a, g = a.tocsc(), g.tocsc()
        p = sparse.triu(a, format="csc")
        p.sort_indices()
        g.sort_indices()
        m, n = g.shape
        diag = dict(
            solver="OSQP_BUILTIN_CPU",
            version=self.osqp.__version__,
            matrix_dimension=n,
            constraint_count=m,
            A_nnz=a.nnz,
            G_nnz=g.nnz,
            cold_start=warm_dual is None,
            primal_warm_start=warm_primal is not None,
            setup_s=0.0,
            solve_s=0.0,
            polish_s=0.0,
            iterations=0,
            workspace_reused=False,
            settings=SETTINGS,
            attempts=[],
        )
        if not m:
            # Isolated, unrestrained Stokes coordinates avoid OSQP entirely.
            u = b / a.diagonal() if a.nnz == n else splu(a).solve(b)
            lam = np.empty(0)
            kkt, terms = kkt_residual(a, b, g, lower, u, lam)
            if not np.isfinite(u).all() or not np.isfinite(kkt) or kkt > 1e-9:
                raise MotionSolverFailure(dict(diag, reason="UNCONSTRAINED_KKT", kkt=kkt))
            diag.update(
                setup_s=perf_counter() - begin,
                kkt=kkt,
                kkt_terms=terms,
                objective=float(0.5 * u @ (a @ u) - b @ u),
                isolated_bypass=True,
            )
            return u, lam, diag
        free = splu(a).solve(b)
        if np.min(g @ free - lower, initial=0) >= -1e-12:
            lam = np.zeros(m)
            kkt, terms = kkt_residual(a, b, g, lower, free, lam)
            if np.isfinite(free).all() and np.isfinite(kkt) and kkt <= 1e-9:
                diag.update(
                    setup_s=perf_counter() - begin,
                    kkt=kkt,
                    kkt_terms=terms,
                    objective=float(0.5 * free @ (a @ free) - b @ free),
                    free_feasible_bypass=True,
                )
                return free, lam, diag
        pattern = hashlib.sha256(
            np.array([m, n]).tobytes()
            + p.indptr.tobytes()
            + p.indices.tobytes()
            + g.indptr.tobytes()
            + g.indices.tobytes()
        ).hexdigest()
        reused = self.workspace is not None and self.pattern == pattern
        if reused:
            self.workspace.update(Px=p.data, Ax=g.data, q=-b, l=lower, u=np.full(m, np.inf))
        else:
            self.workspace = self.osqp.OSQP(algebra="builtin")
            self.workspace.setup(P=p, q=-b, A=g, l=lower, u=np.full(m, np.inf), **SETTINGS)
            self.pattern = pattern
        diag.update(setup_s=perf_counter() - begin, workspace_reused=reused)
        warm = warm_primal is not None and warm_dual is not None
        for mode in ["accepted_warm", "cold_retry"] if warm else ["cold"]:
            if mode == "accepted_warm":
                self.workspace.warm_start(x=warm_primal, y=-warm_dual)
            else:
                # Also reset rho after a failed warm solve: identical cold policy.
                self.workspace.update_settings(rho=0.1)
                self.workspace.warm_start(x=np.zeros(n), y=np.zeros(m))
            result = self.workspace.solve(raise_error=False)
            u, lam = result.x, -result.y
            finite = bool(np.isfinite(u).all() and np.isfinite(lam).all())
            kkt, terms = kkt_residual(a, b, g, lower, u, lam) if finite else (None, {})
            good = bool(finite and np.isfinite(kkt) and kkt <= 1e-9 and result.info.status_val == 1)
            # One sparse saddle solve can polish an ADMM candidate even if its
            # internal stopping criterion stalled. All inequalities remain in
            # the final KKT audit; no dense contact Gram matrix is constructed.
            polish_begin = perf_counter()
            external = dict(attempted=False, accepted=False)
            if not good and finite:
                ix = np.flatnonzero(lam > 1e-10 * max(1.0, np.max(np.abs(lam), initial=0)))
                if len(ix):
                    external["attempted"] = True
                    ca = g[ix]
                    saddle = sparse.bmat([[a, -ca.T], [-ca, None]], format="csc")
                    rhs = np.r_[b, -lower[ix]]
                    try:
                        answer = splu(saddle).solve(rhs)
                        linear = float(np.linalg.norm(saddle @ answer - rhs) / max(1.0, np.linalg.norm(rhs)))
                        v, mult = answer[:n], np.zeros(m)
                        mult[ix] = answer[n:]
                        err, part = kkt_residual(a, b, g, lower, v, mult)
                        external.update(
                            linear_residual=linear if np.isfinite(linear) else None,
                            kkt=err if np.isfinite(err) else None,
                            accepted=bool(
                                np.isfinite(answer).all()
                                and np.isfinite(err)
                                and linear <= 1e-11
                                and err <= 1e-9
                            ),
                        )
                        if external["accepted"]:
                            u, lam, kkt, terms, good = v, mult, err, part, True
                    except RuntimeError as exc:
                        external["error"] = str(exc)
            external["seconds"] = perf_counter() - polish_begin
            row = dict(
                mode=mode,
                status=result.info.status,
                status_val=int(result.info.status_val),
                iterations=int(result.info.iter),
                kkt=kkt,
                kkt_terms=terms,
                solve_s=float(result.info.solve_time),
                polish_s=float(result.info.polish_time) + external["seconds"],
                polish_status=int(result.info.status_polish),
                external_sparse_polish=external,
                accepted=good,
            )
            diag["attempts"].append(row)
            for key in ["solve_s", "polish_s", "iterations"]:
                diag[key] += row[key]
            if good:
                diag.update(
                    kkt=kkt,
                    kkt_terms=terms,
                    objective=float(0.5 * u @ (a @ u) - b @ u),
                    cold_retry_used=mode == "cold_retry",
                    wall_s=perf_counter() - begin,
                )
                return u, lam, diag
        raise MotionSolverFailure(
            dict(diag, reason="OSQP_CONTACT_KKT_FAILURE", wall_s=perf_counter() - begin)
        )
