"""Opt-in full-particle execution; the inherited F6 M0 arithmetic is unchanged.

One control family and clock per batch. Global device errors are intentionally
fail-stop for the whole batch, not silently reinterpreted as per-member success.
"""

from time import perf_counter
import warp as wp
from .engine import FastSimulator
from .f6_engine import F6Simulator
from . import kernels as k


@wp.kernel
def decompose(
    x: wp.array(dtype=wp.vec3d),
    velocity: wp.array(dtype=wp.vec3d),
    parent: wp.array(dtype=int),
    correction: wp.array(dtype=wp.vec3d),
    count: wp.array(dtype=int),
    n: int,
    center: wp.array(dtype=wp.vec3d),
    mean: wp.array(dtype=wp.vec3d),
    relative_x: wp.array(dtype=wp.vec3d),
    relative_u: wp.array(dtype=wp.vec3d),
    residual: wp.array2d(dtype=wp.float64),
):
    i = wp.tid()
    start = i // n * n
    xc = wp.vec3d(wp.float64(0.0))
    um = xc
    for j in range(start, start + n):
        root = k.root_of(parent, j)
        xc = xc + x[j] / wp.float64(n)
        um = um + (velocity[j] + correction[root] / wp.float64(count[root])) / wp.float64(n)
    root = k.root_of(parent, i)
    u = velocity[i] + correction[root] / wp.float64(count[root])
    r = x[i] - xc
    ur = u - um
    center[i] = xc
    mean[i] = um
    relative_x[i] = r
    relative_u[i] = ur
    residual[i, 0] = wp.max(residual[i, 0], wp.length((xc + r) - x[i]))
    residual[i, 1] = wp.max(residual[i, 1], wp.length((um + ur) - u))


class FullParticleBatch(F6Simulator):
    def __init__(self, *args, sample_ids=None, decomposition=False, **kwargs):
        self.decomposition_enabled = False
        super().__init__(*args, **kwargs)
        if self.mobility_mode != "f5_local":
            raise ValueError("F10C admits only full-particle M0")
        self.sample_ids = list(sample_ids or [f"sample_{b:04d}" for b in range(self.batch)])
        if len(self.sample_ids) != self.batch or len(set(self.sample_ids)) != self.batch:
            raise ValueError("Unique sample_id required for every batch member")
        # Capacity is bounded and explicit. Download only the used prefix.
        self.event_ids = wp.empty((max(200000, 20000 * self.batch), 4), dtype=int, device=self.device)
        self.event_time = wp.empty(self.event_ids.shape[0], dtype=wp.float64, device=self.device)
        self.decomposition_enabled = bool(decomposition)
        if decomposition:
            self.decomp_center = wp.zeros_like(self.x)
            self.decomp_mean = wp.zeros_like(self.x)
            self.relative_x = wp.zeros_like(self.x)
            self.relative_u = wp.zeros_like(self.x)
            self.decomp_residual = wp.zeros((self.total, 2), dtype=wp.float64, device=self.device)
            begin = perf_counter()
            wp.load_module(module=decompose.module, device=self.device)
            wp.synchronize_device(self.device)
            self.compile_s += perf_counter() - begin

    def launch(self, kernel, dim, args):
        if self.decomposition_enabled and kernel is k.predict:
            FastSimulator.launch(
                self,
                decompose,
                self.total,
                [
                    self.x,
                    args[1],
                    self.component_parent,
                    self.component_delta,
                    self.component_count,
                    self.n,
                    self.decomp_center,
                    self.decomp_mean,
                    self.relative_x,
                    self.relative_u,
                    self.decomp_residual,
                ],
            )
        return super().launch(kernel, dim, args)

    def failure_members(self):
        """Local evidence plus explicit global abort; never assign a false PASS."""
        diag = self.diag.numpy().reshape(self.batch, self.n, 11)
        flow = self.flow_query.stats.numpy().reshape(self.batch, self.n, 4)
        wall = self.wall_stats.numpy().reshape(self.batch, self.n, 12)
        rows = []
        for b, sample in enumerate(self.sample_ids):
            causes = []
            if diag[b, :, 3].max() > 0.01 * self.p.radius_m:
                causes.append("CONTACT_RESIDUAL")
            if diag[b, :, 6].max() > 0.1 * self.p.radius_m:
                causes.append("CONTACT_PASS")
            if flow[b, :, 3].sum():
                causes.append("CFD_DOMAIN")
            if wall[b, :, 1].max() > 0.1 * self.p.radius_m:
                causes.append("WALL_PROPOSAL")
            rows.append(
                dict(
                    sample_id=sample,
                    status="FAILED_MEMBER" if causes else "ABORTED_BATCH",
                    local_causes=causes,
                    batch_error_bits=int(self.error.numpy()[0]),
                    training_ready=False,
                )
            )
        return rows
