"""Opt-in F10J correction; inherits the complete F10I fixed-state algorithm.

Same A, H, forces, tolerances, 24 outer / 20 Armijo limits. Only PCG's bounded
budget (256), explicit true-residual admission, and translated stable merit differ.
No production import changes. CPU independent validation remains mandatory.
"""

import numpy as np
import warp as wp
from warp.optim.linear import cg
from .f10i_gpu import FixedStateB2, vec, apply, product


@wp.kernel
def energy_increment(
    pred: wp.array(dtype=wp.float64),
    q: wp.array(dtype=wp.float64),
    direction: wp.array(dtype=wp.float64),
    c: wp.array2d(dtype=wp.float64),
    n: wp.array2d(dtype=wp.vec3d),
    rate: wp.float64,
    alpha: wp.float64,
    work: wp.array(dtype=wp.float64),
):
    i = wp.tid()
    s = rate * vec(direction, i)
    As = rate * (vec(direction, i) + apply(i, direction, c, n))
    value = wp.dot(vec(q, i), As) + wp.float64(0.5) * wp.dot(s, As)
    for j in range(c.shape[1]):
        if i != j:
            dx = (vec(pred, i) - vec(pred, j)) + (vec(q, i) - vec(q, j))
            ds = s - rate * vec(direction, j)
            d = wp.length(dx)
            dn = wp.length(dx + ds)
            dd = (wp.float64(2.0) * wp.dot(dx, ds) + wp.dot(ds, ds)) / (dn + d)
            signed = wp.float64(2.0) - d
            old = wp.max(signed, wp.float64(0.0))
            new = wp.max(signed - dd, wp.float64(0.0))
            change = new - old
            if signed > wp.float64(0.0) and signed - dd > wp.float64(0.0):
                change = -dd
            value += wp.float64(0.25) * alpha * change * (new + old)
    work[i] = value


@wp.kernel
def sum_increment(work: wp.array(dtype=wp.float64), stats: wp.array(dtype=wp.float64)):
    total = wp.float64(0.0)
    for i in range(work.shape[0]):
        total += work[i]
    stats[0] = total


@wp.kernel
def check_linear(
    b: wp.array(dtype=wp.float64),
    rec: wp.array(dtype=wp.float64),
    true: wp.array(dtype=wp.float64),
    inv: wp.array(dtype=wp.mat33d),
    counter: wp.array(dtype=wp.int32),
    rows: wp.array2d(dtype=wp.float64),
):
    rr = wp.float64(0.0)
    tt = wp.float64(0.0)
    bb = wp.float64(0.0)
    pre = wp.float64(0.0)
    for i in range(inv.shape[0]):
        rr += wp.dot(vec(rec, i), vec(rec, i))
        tt += wp.dot(vec(true, i), vec(true, i))
        bb += wp.dot(vec(b, i), vec(b, i))
        pre += wp.dot(vec(true, i), inv[i] * vec(true, i))
    k = counter[0]
    rows[k, 0] = wp.sqrt(rr)
    rows[k, 1] = wp.sqrt(tt)
    rows[k, 2] = wp.max(wp.float64(1.0e-14), wp.float64(1.0e-12) * wp.sqrt(bb))
    rows[k, 3] = wp.sqrt(wp.max(pre, wp.float64(0.0)))
    rows[k, 4] = wp.sqrt(bb)
    counter[0] = k + 1


class CorrectedB2(FixedStateB2):
    def __init__(self, *args, **kwargs):
        self.conditional_graph = kwargs.pop("conditional_graph", False)
        super().__init__(*args, **kwargs)
        from time import perf_counter

        t = perf_counter()
        self.increment = wp.zeros(self.n, dtype=wp.float64, device=self.device)
        self.linear_records = wp.zeros((25, 5), dtype=wp.float64, device=self.device)
        self.linear_count = wp.zeros(1, dtype=wp.int32, device=self.device)
        self.true_residual = wp.zeros(3 * self.n, dtype=wp.float64, device=self.device)
        self.candidate_trace = []
        self.rate = 1.0
        self.outer_index = 0
        for kind, contact, inv in [("A", 0, self.invA), ("H", 1, self.invH)]:
            old = getattr(self, "cg" + kind)
            solver = cg(
                old._A,
                old._b,
                old._x,
                M=old._M,
                tol=1e-12,
                atol=1e-14,
                maxiter=256,
                check_every=0,
                use_cuda_graph=self.conditional_graph,
                run=False,
            )
            setattr(self, "cg" + kind, solver)
            with wp.ScopedCapture(device=self.device) as cap:
                result = solver()
                self.launch(
                    product,
                    self.n,
                    [
                        solver._x,
                        solver._b,
                        self.true_residual,
                        self.c,
                        self.normal,
                        self.a,
                        self.current,
                        contact,
                        -1.0,
                        1.0,
                    ],
                )
                self.launch(
                    check_linear,
                    1,
                    [
                        solver._b,
                        solver._r_and_z_buf[0],
                        self.true_residual,
                        inv,
                        self.linear_count,
                        self.linear_records,
                    ],
                )
            setattr(self, "linear" + kind, result)
            setattr(self, "graph" + kind, cap.graph)
        self.correction_graph_build_s = perf_counter() - t

    def launch(self, k, dim, args):
        from .f10i_gpu import add

        if k is add and args[2] is self.trial:
            self.rate = float(args[3])
        super().launch(k, dim, args)

    def evaluate(self, q):
        # Translate the merit origin to E(q_current)=0. The inherited Armijo test
        # then compares Delta E <= c * rate * g.d, without subtracting total energies.
        from .f10i_gpu import evaluate, reduce
        from time import perf_counter

        self.launch(
            evaluate,
            self.n,
            [self.pred, q, self.direction, self.c, self.normal, self.alpha, self.rhs, self.work, self.force],
        )
        self.launch(reduce, 1, [self.work, self.stats])
        if q is self.trial:
            self.launch(
                energy_increment,
                self.n,
                [
                    self.pred,
                    self.q,
                    self.direction,
                    self.c,
                    self.normal,
                    self.rate,
                    self.alpha,
                    self.increment,
                ],
            )
            self.launch(sum_increment, 1, [self.increment, self.stats])
        t = perf_counter()
        stats = self.stats.numpy()
        self.scalar_read_s += perf_counter() - t
        self.scalar_reads += 1
        if q is self.q:
            stats[0] = 0.0
            self.old_derivative = stats[2]
        else:
            self.candidate_trace.append(
                dict(
                    rate=self.rate,
                    energy_change=float(stats[0]),
                    armijo_rhs=1e-4 * self.rate * self.old_derivative,
                    gradient_bound_m=float(self.r * stats[1]),
                    overlap_m=float(self.r * stats[4]),
                    accepted=bool(stats[0] <= 1e-4 * self.rate * self.old_derivative),
                )
            )
        return stats

    def solve(self, x, free):
        from time import perf_counter

        start = perf_counter()
        x, free = np.asarray(x, float), np.asarray(free, float)
        if x.shape != (self.n, 3) or free.shape != x.shape or not np.all(np.isfinite(x + free)):
            raise ValueError("Finite N x 3 inputs required")
        distances = np.linalg.norm(x[:, None, :] - x[None, :, :], axis=2)
        if np.any(distances[np.triu_indices(self.n, 1)] <= 0):
            raise ValueError("Coincident centers have undefined contact normals")
        self.linear_count.zero_()
        self.candidate_trace = []
        result = super().solve(x, free)
        t = perf_counter()
        count = int(self.linear_count.numpy()[0])
        rows = self.linear_records.numpy()[:count]
        result["linear_audit"] = rows
        result["candidate_trace"] = self.candidate_trace
        result["mathematical_stop_before_quality"] = result["converged"]
        if not np.all(np.isfinite(rows)) or not np.all(np.isfinite(result["positions_m"])):
            result.update(converged=False, reason="NONFINITE")
        elif not np.all(rows[:, :2] <= rows[:, 2, None]):
            result.update(converged=False, reason="TRUE_OR_RECURSIVE_LINEAR_RESIDUAL")
        elif result["max_overlap_m"] > 17.8e-9:
            result.update(converged=False, reason="OVERLAP_BUDGET")
        elif result["max_single_iteration_m"] > 178e-9:
            result.update(converged=False, reason="CORRECTION_BUDGET")
        result["timing"]["linear_audit_read_wall_s"] = perf_counter() - t
        result["timing"]["whole_step_wall_s"] = perf_counter() - start
        result["solver_version"] = "F10J_STABLE_MERIT_PCG256_V1"
        return result


class OptimizedB2(CorrectedB2):
    """Reuse accepted evaluation at the next outer check; PCG uses Warp's bounded
    captured conditional loop. Both optimizations preserve the accepted q exactly.
    Unaccepted trials are never cached. No post-convergence q update is introduced.
    """

    def __init__(self, *args, **kwargs):
        kwargs["conditional_graph"] = True
        super().__init__(*args, **kwargs)
        self.accepted_stats = None
        self.reused_evaluations = 0

    def evaluate(self, q):
        if q is self.q and self.accepted_stats is not None:
            stats = self.accepted_stats
            self.accepted_stats = None
            self.reused_evaluations += 1
            return stats.copy()
        stats = super().evaluate(q)
        if q is self.trial and self.candidate_trace[-1]["accepted"]:
            self.accepted_stats = stats.copy()
            self.accepted_stats[0] = 0.0
        return stats

    def solve(self, x, free):
        self.accepted_stats = None
        self.reused_evaluations = 0
        result = super().solve(x, free)
        result["reused_evaluations"] = self.reused_evaluations
        result["solver_version"] = "F10J_STABLE_MERIT_CONDITIONAL_CG_CACHE_V2"
        return result
