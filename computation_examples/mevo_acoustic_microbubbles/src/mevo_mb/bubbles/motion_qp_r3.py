"""R3 QP adapter: same R2 native solvers and independent final certificates.

An already computed SAME-matrix free solution avoids redundant factorization.
Its original-equation KKT is checked before the free fast path. No matrix factor
or numerical values persist across changed inputs. Fallback behavior stays R2.
"""

from time import perf_counter
import hashlib
import numpy as np
from scipy.sparse.linalg import splu
from .motion_qp_r2 import NativeQP, certificate, pattern
from .constrained_motion import MotionSolverFailure, kkt_residual


class NativeQPR3(NativeQP):
    def solve(self, a, b, g, lower, *, warm_primal=None, warm_dual=None, feasible=None, free_hint=None):
        a, g = a.tocsc(), g.tocsc()
        if not all(np.isfinite(v).all() for v in [a.data, b, g.data, lower]):
            raise MotionSolverFailure(dict(reason="NONFINITE_INPUT"))
        start = perf_counter()
        try:
            free = splu(a).solve(b) if free_hint is None else np.asarray(free_hint).copy()
        except RuntimeError as exc:
            raise MotionSolverFailure(dict(reason="INVALID_MODEL", native_error=str(exc))) from exc
        free_kkt, _ = kkt_residual(a, b, g, lower, free, np.zeros(len(lower)))
        if free_kkt <= 1e-12:
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
