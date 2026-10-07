"""R3 V2 overhead reduction; resistance/contact/geometry equations unchanged.

Versioned from R2's component solve: vectorized row ownership and certified
same-input free-solution reuse. Native ProxQP, alternatives, KKT, trial isolation
and persistent per-stage/component workspaces remain intact.
"""

from collections import OrderedDict
from concurrent.futures import ThreadPoolExecutor
from time import perf_counter, thread_time
from threading import get_native_id
import hashlib
import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from .motion_r3 import MotionR3
from .motion_qp_r3 import NativeQPR3
from .constrained_motion import MotionSolverFailure, AcceptedCache, kkt_residual


class MotionR3Optimized(MotionR3):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.qp_factory = lambda: NativeQPR3(self.qp_solver)
        self.qp = self.qp_factory()
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
            answer = self.constrained(result, constraints, ids, feasible=feasible)
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
        if np.any(np.diff(csr.indptr) == 0):
            raise MotionSolverFailure({"reason": "EMPTY_CONSTRAINT_ROW"})
        if g.shape[0]:
            labs = labels[csr.indices // 3]
            owners = np.minimum.reduceat(labs, csr.indptr[:-1])
            maximum = np.maximum.reduceat(labs, csr.indptr[:-1])
            if np.any(owners != maximum):
                raise MotionSolverFailure({"reason": "CROSS_COMPONENT_CONSTRAINT"})
        else:
            owners = np.empty(0, int)

        def component(label):
            component_start = perf_counter()
            cpu_start = thread_time()
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
                aa,
                b[cols],
                gg,
                lower[rr] / vs,
                warm_dual=warm,
                warm_primal=wp,
                feasible=vv,
                free_hint=free.ravel()[cols] / vs,
            )
            d.update(
                component_wall_s=perf_counter() - component_start,
                component_thread_cpu_s=thread_time() - cpu_start,
                native_thread_id=get_native_id(),
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
