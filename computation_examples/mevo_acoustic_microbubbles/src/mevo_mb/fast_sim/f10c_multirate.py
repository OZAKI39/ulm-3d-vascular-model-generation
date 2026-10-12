"""Guarded experimental external-field caching. All particle mechanics stay F6.

CFD proposal validity is still queried EVERY physical step. Cached flow is used
only for the free-velocity RHS. Acoustic caching includes pressure, eight volume
harmonics, consistent gradient and primary force; F2 geometry is never cached.
"""

from math import lcm
from time import perf_counter
import warp as wp
from .f10c_engine import FullParticleBatch
from .engine import FastSimulator
from .f5_engine import flow_free_velocity
from .flow_gpu import contains, barycentric
from . import kernels as k


@wp.kernel
def field_guard(
    x: wp.array(dtype=wp.vec3d),
    old: wp.array(dtype=wp.vec3d),
    active: wp.array(dtype=int),
    mask: wp.array(dtype=int),
    grid: wp.array2d(dtype=wp.vec2d),
    origin: wp.vec3d,
    spacing: wp.vec3d,
    nx: int,
    ny: int,
    nz: int,
    knots: wp.array(dtype=wp.float64),
    valid: wp.array(dtype=int),
    kvalid: wp.array(dtype=int),
    cache_x: wp.array(dtype=wp.vec3d),
    cache_d: wp.array(dtype=wp.float64),
    radius: wp.float64,
    events: wp.array(dtype=int),
    stamp: wp.array(dtype=int),
    error: wp.array(dtype=int),
    stats: wp.array2d(dtype=int),
    diag: wp.array2d(dtype=wp.float64),
):
    i = wp.tid()
    mask[i] = active[i]
    if error[0] != 0:
        return
    if active[i] == 0:
        wp.atomic_or(error, 0, 65536)
        return
    q = wp.cw_div(x[i] - origin, spacing)
    q0 = wp.cw_div(old[i] - origin, spacing)
    for axis in range(3):
        limit = nx - 1
        if axis == 1:
            limit = ny - 1
        if axis == 2:
            limit = nz - 1
        if not wp.isfinite(q[axis]) or q[axis] < wp.float64(0.0) or q[axis] > wp.float64(limit):
            wp.atomic_or(error, 0, 1)
            stats[i, 3] = stats[i, 3] + 1
            return
    # Reject close-wall cache use; physical per-step projection/checks remain on.
    if cache_d[i] - wp.length(x[i] - cache_x[i]) - radius < wp.float64(1e-4):
        wp.atomic_or(error, 0, 65536)
        stats[i, 3] = stats[i, 3] + 1
        return
    ix = wp.min(int(wp.floor(q[0])), nx - 2)
    iy = wp.min(int(wp.floor(q[1])), ny - 2)
    iz = wp.min(int(wp.floor(q[2])), nz - 2)
    w = q - wp.vec3d(wp.float64(ix), wp.float64(iy), wp.float64(iz))
    p = wp.vec2d(wp.float64(0.0))
    for a in range(2):
        for b in range(2):
            for c in range(2):
                wa = wp.float64(1.0) - w[0]
                wb = wp.float64(1.0) - w[1]
                wc = wp.float64(1.0) - w[2]
                if a == 1:
                    wa = w[0]
                if b == 1:
                    wb = w[1]
                if c == 1:
                    wc = w[2]
                p = p + (wa * wb * wc) * grid[((ix + a) * ny + iy + b) * nz + iz + c, 0]
    amp = wp.length(p)
    ok = int(amp == wp.float64(0.0))
    diag[i, 0] = wp.min(diag[i, 0], amp)
    diag[i, 1] = wp.max(diag[i, 1], amp)
    for j in range(knots.shape[0]):
        if amp == knots[j]:
            ok = kvalid[j]
        if j < knots.shape[0] - 1 and amp > knots[j] and amp < knots[j + 1]:
            ok = valid[j]
    if ok == 0:
        wp.atomic_or(error, 0, 2)
        stats[i, 3] = stats[i, 3] + 1
        return
    # Cell-face derivatives can jump. Any adjacency change invalidates response.
    changed = (
        int(wp.floor(q0[0])) != ix
        or int(wp.floor(q0[1])) != iy
        or int(wp.floor(q0[2])) != iz
        or events[0] != stamp[i]
    )
    if not changed:
        mask[i] = 0
        stats[i, 0] = stats[i, 0] + 1
    else:
        stats[i, 1] = stats[i, 1] + 1


@wp.kernel
def field_cache(
    mask: wp.array(dtype=int),
    x: wp.array(dtype=wp.vec3d),
    force: wp.array(dtype=wp.vec3d),
    pressure: wp.array(dtype=wp.vec2d),
    harmonics: wp.array2d(dtype=wp.vec2d),
    grad: wp.array2d(dtype=wp.vec2d),
    old: wp.array(dtype=wp.vec3d),
    cf: wp.array(dtype=wp.vec3d),
    cp: wp.array(dtype=wp.vec2d),
    ch: wp.array2d(dtype=wp.vec2d),
    cg: wp.array2d(dtype=wp.vec2d),
    stamp: wp.array(dtype=int),
    events: wp.array(dtype=int),
):
    i = wp.tid()
    if mask[i] == 1:
        old[i] = x[i]
        cf[i] = force[i]
        cp[i] = pressure[i]
        stamp[i] = events[0]
        for j in range(8):
            ch[i, j] = harmonics[i, j]
        for j in range(3):
            cg[i, j] = grad[i, j]
    else:
        force[i] = cf[i]
        pressure[i] = cp[i]
        for j in range(8):
            harmonics[i, j] = ch[i, j]
        for j in range(3):
            grad[i, j] = cg[i, j]


@wp.kernel
def save_flow(
    velocity: wp.array(dtype=wp.vec3d),
    cache: wp.array(dtype=int),
    saved: wp.array(dtype=wp.vec3d),
    saved_cell: wp.array(dtype=int),
    stamp: wp.array(dtype=int),
    events: wp.array(dtype=int),
):
    i = wp.tid()
    saved[i] = velocity[i]
    saved_cell[i] = cache[i]
    stamp[i] = events[0]


@wp.kernel
def select_flow(
    x: wp.array(dtype=wp.vec3d),
    active: wp.array(dtype=int),
    velocity: wp.array(dtype=wp.vec3d),
    valid: wp.array(dtype=int),
    cache: wp.array(dtype=int),
    origins: wp.array(dtype=wp.vec3d),
    inverses: wp.array(dtype=wp.mat33d),
    gradient: wp.array(dtype=wp.mat33d),
    saved: wp.array(dtype=wp.vec3d),
    saved_cell: wp.array(dtype=int),
    stamp: wp.array(dtype=int),
    events: wp.array(dtype=int),
    out: wp.array(dtype=wp.vec3d),
    cache_x: wp.array(dtype=wp.vec3d),
    cache_d: wp.array(dtype=wp.float64),
    radius: wp.float64,
    error: wp.array(dtype=int),
    stats: wp.array2d(dtype=int),
):
    i = wp.tid()
    out[i] = velocity[i]
    if error[0] != 0:
        return
    if active[i] == 0:
        wp.atomic_or(error, 0, 65536)
        return
    cell = cache[i]
    if valid[i] == 0 or cell < 0:
        wp.atomic_or(error, 0, 16384)
        return
    if not contains(barycentric(x[i], origins[cell], inverses[cell])):
        wp.atomic_or(error, 0, 65536)
        return
    shear = wp.float64(0.0)
    for a in range(3):
        for b in range(3):
            d = wp.float64(0.5) * (gradient[i][a, b] + gradient[i][b, a])
            shear = shear + wp.float64(2.0) * d * d
    if wp.sqrt(shear) > wp.float64(2.0) or cache_d[i] - wp.length(x[i] - cache_x[i]) - radius < wp.float64(
        1e-4
    ):
        wp.atomic_or(error, 0, 65536)
        return
    if cell == saved_cell[i] and stamp[i] == events[0]:
        out[i] = saved[i]
        stats[i, 2] = stats[i, 2] + 1
    else:
        # Previous proposal lookup is at EXACT current accepted x; safe refresh.
        saved[i] = velocity[i]
        saved_cell[i] = cell
        stamp[i] = events[0]
        stats[i, 4] = stats[i, 4] + 1


class MultirateBatch(FullParticleBatch):
    def __init__(self, *args, field_stride=1, flow_stride=1, **kwargs):
        if field_stride not in (1, 2, 4) or flow_stride not in (1, 2, 4):
            raise ValueError("Unregistered stride")
        self.field_stride = field_stride
        self.flow_stride = flow_stride
        self.cache_ready = False
        self.phase = 0
        self.build_phase = 0
        self.last_dt = None
        super().__init__(*args, **kwargs)
        self.refresh_mask = wp.ones(self.total, dtype=int, device=self.device)
        self.saved_x = wp.empty_like(self.x)
        self.saved_primary = wp.zeros_like(self.force)
        self.saved_pressure = wp.zeros_like(self.pressure)
        self.saved_harmonics = wp.zeros_like(self.harmonics)
        self.saved_grad = wp.zeros_like(self.consistent_gradients)
        self.field_stamp = wp.zeros(self.total, dtype=int, device=self.device)
        self.flow_stamp = wp.zeros_like(self.field_stamp)
        self.saved_flow = wp.zeros_like(self.free)
        self.saved_cell = wp.full(self.total, -1, dtype=int, device=self.device)
        self.flow_for_free = wp.zeros_like(self.free)
        self.cache_stats = wp.zeros((self.total, 5), dtype=int, device=self.device)
        self.cache_ready = True
        begin = perf_counter()
        wp.load_module(module=field_guard.module, device=self.device)
        wp.synchronize_device(self.device)
        self.compile_s += perf_counter() - begin

    def reset(self):
        super().reset()
        self.phase = 0
        self.build_phase = 0
        self.last_dt = None
        if self.cache_ready:
            self.cache_stats.zero_()

    def control(self, event):
        record = super().control(event)
        self.phase = 0
        self.build_phase = 0
        self.graphs.clear()
        return dict(record, multirate_cache_invalidated=True)

    def fields(self):
        if self.field_stride == 1 or not self.acoustic_enabled:
            return super().fields()
        mask = self.active
        if self.build_phase % self.field_stride:
            self.launch(
                field_guard,
                self.total,
                [
                    self.x,
                    self.saved_x,
                    self.active,
                    self.refresh_mask,
                    self.grid,
                    self.origin,
                    self.spacing,
                    *self.shape,
                    self.knots,
                    self.valid,
                    self.kval,
                    self.cache_x,
                    self.cache_d,
                    self.p.radius_m,
                    self.event_count,
                    self.field_stamp,
                    self.error,
                    self.cache_stats,
                    self.diag,
                ],
            )
            mask = self.refresh_mask
        actual_active = self.active
        self.active = mask
        try:
            super().fields()
        finally:
            self.active = actual_active
        self.launch(
            field_cache,
            self.total,
            [
                mask,
                self.x,
                self.force,
                self.pressure,
                self.harmonics,
                self.consistent_gradients,
                self.saved_x,
                self.saved_primary,
                self.saved_pressure,
                self.saved_harmonics,
                self.saved_grad,
                self.field_stamp,
                self.event_count,
            ],
        )

    def launch(self, kernel, dim, args):
        if self.cache_ready and kernel is k.free_velocity and self.flow_enabled and self.flow_stride > 1:
            if self.build_phase % self.flow_stride == 0:
                result = super().launch(kernel, dim, args)
                FastSimulator.launch(
                    self,
                    save_flow,
                    self.total,
                    [
                        self.flow_query.velocity,
                        self.flow_query.cache,
                        self.saved_flow,
                        self.saved_cell,
                        self.flow_stamp,
                        self.event_count,
                    ],
                )
                return result
            f = self.flow
            q = self.flow_query
            FastSimulator.launch(
                self,
                select_flow,
                self.total,
                [
                    self.x,
                    self.active,
                    q.velocity,
                    q.valid,
                    q.cache,
                    f.origins,
                    f.inverses,
                    q.gradient,
                    self.saved_flow,
                    self.saved_cell,
                    self.flow_stamp,
                    self.event_count,
                    self.flow_for_free,
                    self.cache_x,
                    self.cache_d,
                    self.p.radius_m,
                    self.error,
                    self.cache_stats,
                ],
            )
            FastSimulator.launch(
                self,
                flow_free_velocity,
                self.total,
                [
                    self.force,
                    self.flow_for_free,
                    self.free,
                    self.active,
                    args[2],
                    1,
                    self._current_dt,
                    self.free_path,
                ],
            )
            return
        return super().launch(kernel, dim, args)

    def step(self, dt):
        super().step(dt)
        self.build_phase += 1

    def graph(self, dt, count):
        phase = self.phase % lcm(self.field_stride, self.flow_stride)
        key = (float(dt), int(count), phase)
        if key not in self.graphs:
            begin = perf_counter()
            saved = self.build_phase
            self.build_phase = phase
            with wp.ScopedCapture(device=self.device) as capture:
                for _ in range(count):
                    self.step(dt)
            self.build_phase = saved
            self.graphs[key] = capture.graph
            self.capture_s += perf_counter() - begin
        return self.graphs[key]

    def advance(self, dt, count, graphs=True):
        if self.last_dt is not None and dt != self.last_dt:
            self.phase = 0
            self.graphs.clear()
        self.last_dt = dt
        remaining = count
        while remaining:
            size = min(remaining, self.graph_chunk)
            if graphs:
                wp.capture_launch(self.graph(dt, size))
            else:
                self.build_phase = self.phase
                for _ in range(size):
                    self.step(dt)
            self.phase += size
            remaining -= size
        q = self.quality()
        if (
            q["max_accepted_penetration_m"] > 0.01 * self.p.radius_m
            or q["max_contact_pass_correction_m"] > 0.1 * self.p.radius_m
        ):
            raise RuntimeError("F6_REGISTERED_CONTACT_BUDGET")
        if q["error_bits"]:
            raise RuntimeError("MULTIRATE_NOT_ADMITTED_ERROR_" + str(q["error_bits"]))
