"""B=1 fixed-state Warp prototype of the unchanged geometric B2 equation.

Matrix-free pair-normal resistance and SPD Gauss-Newton metric. PCG is used
only on A=I+sum c B'B and H=A+alpha B_active'B_active, never on the full
nonlinear Hessian. No production engine import or default is changed.
"""

from time import perf_counter
import numpy as np
import warp as wp
from warp.optim.linear import LinearOperator, cg


@wp.func
def vec(a: wp.array(dtype=wp.float64), i: int):
    return wp.vec3d(a[3 * i], a[3 * i + 1], a[3 * i + 2])


@wp.func
def phi(e: wp.float64):
    return wp.float64(0.5) / e - wp.float64(0.45) * wp.log(e) - wp.float64(3.0 / 56.0) * e * wp.log(e)


@wp.kernel
def resistance_pairs(
    x: wp.array(dtype=wp.float64), c: wp.array2d(dtype=wp.float64), normal: wp.array2d(dtype=wp.vec3d)
):
    i, j = wp.tid()
    c[i, j] = wp.float64(0.0)
    normal[i, j] = wp.vec3d(0.0)
    if i != j:
        dx = vec(x, j) - vec(x, i)
        d = wp.length(dx)
        e = d - wp.float64(2.0)
        normal[i, j] = dx / d
        if e < wp.float64(0.1):
            c[i, j] = wp.float64(0.5) * (phi(wp.max(e, wp.float64(0.001))) - phi(wp.float64(0.1)))


@wp.kernel
def current_geometry(
    pred: wp.array(dtype=wp.float64),
    q: wp.array(dtype=wp.float64),
    active: wp.array2d(dtype=wp.float64),
    normal: wp.array2d(dtype=wp.vec3d),
    alpha: wp.float64,
):
    i, j = wp.tid()
    active[i, j] = wp.float64(0.0)
    normal[i, j] = wp.vec3d(0.0)
    if i != j:
        dx = vec(pred, j) + vec(q, j) - vec(pred, i) - vec(q, i)
        d = wp.length(dx)
        normal[i, j] = dx / d
        if d < wp.float64(2.0):
            active[i, j] = alpha


@wp.func
def apply(
    i: int, v: wp.array(dtype=wp.float64), c: wp.array2d(dtype=wp.float64), n: wp.array2d(dtype=wp.vec3d)
):
    out = wp.vec3d(0.0)
    for j in range(c.shape[1]):
        diff = vec(v, i) - vec(v, j)
        out += c[i, j] * wp.dot(n[i, j], diff) * n[i, j]
    return out


@wp.kernel
def product(
    x: wp.array(dtype=wp.float64),
    y: wp.array(dtype=wp.float64),
    z: wp.array(dtype=wp.float64),
    c: wp.array2d(dtype=wp.float64),
    n: wp.array2d(dtype=wp.vec3d),
    a: wp.array2d(dtype=wp.float64),
    nn: wp.array2d(dtype=wp.vec3d),
    contact: int,
    alpha: wp.float64,
    beta: wp.float64,
):
    i = wp.tid()
    v = vec(x, i) + apply(i, x, c, n)
    if contact != 0:
        v += apply(i, x, a, nn)
    for k in range(3):
        previous = wp.float64(0.0)
        if beta != wp.float64(0.0):
            previous = y[3 * i + k]
        z[3 * i + k] = alpha * v[k] + beta * previous


@wp.kernel
def blocks(
    c: wp.array2d(dtype=wp.float64),
    n: wp.array2d(dtype=wp.vec3d),
    a: wp.array2d(dtype=wp.float64),
    nn: wp.array2d(dtype=wp.vec3d),
    contact: int,
    inv: wp.array(dtype=wp.mat33d),
):
    i = wp.tid()
    m = wp.identity(n=3, dtype=wp.float64)
    for j in range(c.shape[1]):
        m += c[i, j] * wp.outer(n[i, j], n[i, j])
        if contact != 0:
            m += a[i, j] * wp.outer(nn[i, j], nn[i, j])
    inv[i] = wp.inverse(m)


@wp.kernel
def precondition(
    x: wp.array(dtype=wp.float64),
    y: wp.array(dtype=wp.float64),
    z: wp.array(dtype=wp.float64),
    inv: wp.array(dtype=wp.mat33d),
    alpha: wp.float64,
    beta: wp.float64,
):
    i = wp.tid()
    v = inv[i] * vec(x, i)
    for k in range(3):
        old = wp.float64(0.0)
        if beta != wp.float64(0.0):
            old = y[3 * i + k]
        z[3 * i + k] = alpha * v[k] + beta * old


@wp.kernel
def add(
    a: wp.array(dtype=wp.float64),
    b: wp.array(dtype=wp.float64),
    out: wp.array(dtype=wp.float64),
    rate: wp.float64,
):
    i = wp.tid()
    out[i] = a[i] + rate * b[i]


@wp.kernel
def evaluate(
    pred: wp.array(dtype=wp.float64),
    q: wp.array(dtype=wp.float64),
    direction: wp.array(dtype=wp.float64),
    c: wp.array2d(dtype=wp.float64),
    n: wp.array2d(dtype=wp.vec3d),
    alpha: wp.float64,
    rhs: wp.array(dtype=wp.float64),
    work: wp.array2d(dtype=wp.float64),
    forces: wp.array(dtype=wp.float64),
):
    i = wp.tid()
    y = vec(pred, i) + vec(q, i)
    aq = vec(q, i) + apply(i, q, c, n)
    f = wp.vec3d(0.0)
    pen = wp.float64(0.0)
    ov = wp.float64(0.0)
    for j in range(c.shape[1]):
        if j != i:
            dx = y - vec(pred, j) - vec(q, j)
            d = wp.length(dx)
            delta = wp.max(wp.float64(2.0) - d, wp.float64(0.0))
            f += delta * dx / d
            pen += delta * delta
            ov = wp.max(ov, delta)
    g = aq - alpha * f
    for k in range(3):
        rhs[3 * i + k] = -g[k]
        forces[3 * i + k] = f[k]
    work[i, 0] = wp.float64(0.5) * wp.dot(vec(q, i), aq) + wp.float64(0.25) * alpha * pen
    work[i, 1] = wp.dot(g, g)
    work[i, 2] = wp.dot(g, vec(direction, i))
    work[i, 3] = wp.length(vec(direction, i))
    work[i, 4] = ov


@wp.kernel
def reduce(work: wp.array2d(dtype=wp.float64), out: wp.array(dtype=wp.float64)):
    energy = wp.float64(0.0)
    norm = wp.float64(0.0)
    gd = wp.float64(0.0)
    dx = wp.float64(0.0)
    ov = wp.float64(0.0)
    for i in range(work.shape[0]):
        energy += work[i, 0]
        norm += work[i, 1]
        gd += work[i, 2]
        dx = wp.max(dx, work[i, 3])
        ov = wp.max(ov, work[i, 4])
    out[0] = energy
    out[1] = wp.sqrt(norm)
    out[2] = gd
    out[3] = dx
    out[4] = ov


class FixedStateB2:
    """Reusable buffers and two captured bounded PCG solves; host Armijo loop.

    Full vectors are read only once at the end for validation. Each outer/line
    iteration reads five diagnostic scalars. No CPU resistance assembly/solve.
    """

    # Original public research limit stays32. F10K's explicitly opted-in subclass
    # raises only allocation scope; all equations, kernels and solver targets stay unchanged.
    max_fixed_state_bubbles = 32

    def __init__(self, n, radius=1.78e-6, dt=1e-5, kc=0.04641319518378658, device="cuda:0"):
        if not 2 <= n <= self.max_fixed_state_bubbles:
            raise ValueError("Only 2..32 bubbles, B=1 fixed-state research")
        self.n = n
        self.r = radius
        self.dt = dt
        self.kc = kc
        self.zeta = 6 * np.pi * 0.001 * radius
        self.alpha = dt * kc / self.zeta
        self.device = wp.get_device(device)
        t = perf_counter()
        for name in ["x", "b", "u", "pred", "q", "trial", "rhs", "direction", "force"]:
            setattr(self, name, wp.zeros(3 * n, dtype=wp.float64, device=self.device))
        self.c = wp.zeros((n, n), dtype=wp.float64, device=self.device)
        self.a = wp.zeros_like(self.c)
        self.normal = wp.zeros((n, n), dtype=wp.vec3d, device=self.device)
        self.current = wp.zeros_like(self.normal)
        self.invA = wp.zeros(n, dtype=wp.mat33d, device=self.device)
        self.invH = wp.zeros_like(self.invA)
        self.work = wp.zeros((n, 5), dtype=wp.float64, device=self.device)
        self.stats = wp.zeros(5, dtype=wp.float64, device=self.device)

        def operator(contact):
            def mv(x, y, z, alpha, beta):
                self.launch(
                    product,
                    self.n,
                    [x, y, z, self.c, self.normal, self.a, self.current, contact, alpha, beta],
                )

            return LinearOperator((3 * n, 3 * n), wp.float64, self.device, mv)

        def pre(inv):
            def mv(x, y, z, alpha, beta):
                self.launch(precondition, self.n, [x, y, z, inv, alpha, beta])

            return LinearOperator((3 * n, 3 * n), wp.float64, self.device, mv)

        self.A = operator(0)
        self.H = operator(1)
        self.cgA = cg(
            self.A,
            self.b,
            self.u,
            tol=1e-12,
            atol=1e-14,
            maxiter=128,
            M=pre(self.invA),
            check_every=0,
            use_cuda_graph=False,
            run=False,
        )
        self.cgH = cg(
            self.H,
            self.rhs,
            self.direction,
            tol=1e-12,
            atol=1e-14,
            maxiter=128,
            M=pre(self.invH),
            check_every=0,
            use_cuda_graph=False,
            run=False,
        )
        wp.synchronize_device(self.device)
        self.allocation_s = perf_counter() - t
        t = perf_counter()
        with wp.ScopedCapture(device=self.device) as cap:
            self.linearA = self.cgA()
        self.graphA = cap.graph
        with wp.ScopedCapture(device=self.device) as cap:
            self.linearH = self.cgH()
        self.graphH = cap.graph
        self.graph_build_s = perf_counter() - t
        self.buffer_bytes = sum(v.capacity for v in self.__dict__.values() if isinstance(v, wp.array))

    def launch(self, k, dim, args):
        wp.launch(k, dim=dim, inputs=args, device=self.device)

    def evaluate(self, q):
        self.launch(
            evaluate,
            self.n,
            [self.pred, q, self.direction, self.c, self.normal, self.alpha, self.rhs, self.work, self.force],
        )
        self.launch(reduce, 1, [self.work, self.stats])
        t = perf_counter()
        s = self.stats.numpy()
        self.scalar_read_s += perf_counter() - t
        self.scalar_reads += 1
        return s

    def solve(self, x, free):
        start = perf_counter()
        self.scalar_read_s = 0.0
        self.scalar_reads = 0
        x = np.asarray(x, float)
        free = np.asarray(free, float)
        self.origin = x.mean(0)
        if x.shape != (self.n, 3) or free.shape != x.shape or not np.all(np.isfinite(x + free)):
            raise ValueError("Finite N x 3 inputs required")
        self.x.assign(((x - self.origin) / self.r).ravel())
        self.b.assign((free * self.dt / self.r).ravel())
        for a in [self.u, self.q, self.trial, self.direction]:
            a.zero_()
        self.launch(resistance_pairs, (self.n, self.n), [self.x, self.c, self.normal])
        self.launch(blocks, self.n, [self.c, self.normal, self.a, self.current, 0, self.invA])
        wp.synchronize_device(self.device)
        prep = perf_counter() - start
        ev0 = wp.Event(device=self.device, enable_timing=True)
        ev1 = wp.Event(device=self.device, enable_timing=True)
        t = perf_counter()
        wp.record_event(ev0)
        wp.capture_launch(self.graphA)
        self.launch(add, 3 * self.n, [self.x, self.u, self.pred, 1.0])
        rows = []
        converged = False
        reason = "OUTER_LIMIT"
        for outer in range(25):
            stats = self.evaluate(self.q)
            if self.r * stats[1] <= 1e-14:
                converged = True
                reason = "CONSERVATIVE_GRADIENT_BOUND"
                break
            if outer == 24:
                break
            self.launch(
                current_geometry, (self.n, self.n), [self.pred, self.q, self.a, self.current, self.alpha]
            )
            self.launch(blocks, self.n, [self.c, self.normal, self.a, self.current, 1, self.invH])
            self.direction.zero_()
            wp.capture_launch(self.graphH)
            old = self.evaluate(self.q)
            rate = 1.0
            if old[2] >= 0:
                reason = "NON_DESCENT"
                break
            for line in range(20):
                self.launch(add, 3 * self.n, [self.q, self.direction, self.trial, rate])
                trial = self.evaluate(self.trial)
                if trial[0] <= old[0] + 1e-4 * rate * old[2]:
                    break
                rate *= 0.5
            else:
                reason = "LINE_SEARCH_LIMIT"
                break
            wp.copy(self.q, self.trial)
            correction = self.r * rate * old[3]
            rows.append(
                dict(
                    outer=outer + 1,
                    gradient_bound_before_m=self.r * stats[1],
                    line_search_rate=rate,
                    max_correction_m=correction,
                )
            )
            if correction > 178e-9:
                reason = "CORRECTION_BUDGET"
                break
        final = self.evaluate(self.q)
        wp.record_event(ev1)
        wp.synchronize_device(self.device)
        solve_s = perf_counter() - t
        t = perf_counter()
        pos = (self.pred.numpy() + self.q.numpy()).reshape(-1, 3) * self.r + self.origin
        force = self.force.numpy().reshape(-1, 3) * self.kc * self.r
        c = self.c.numpy()
        n = self.normal.numpy()
        a = self.a.numpy()
        linear = [
            dict(
                iterations=int(v[0].numpy()[0]),
                recursive_residual_squared=float(v[1].numpy()[0]),
                tolerance_squared=float(v[2].numpy()[0]),
            )
            for v in [self.linearA, self.linearH]
        ]
        read_s = perf_counter() - t
        return dict(
            positions_m=pos,
            contact_force_N=force,
            converged=bool(converged),
            reason=reason,
            conservative_equation_bound_m=float(self.r * final[1]),
            max_overlap_m=float(self.r * final[4]),
            max_single_iteration_m=max((z["max_correction_m"] for z in rows), default=0),
            iterations=rows,
            linear_solves=linear,
            resistance_coefficients=c,
            resistance_normals=n,
            active_metric_coefficients=a,
            timing=dict(
                preparation_wall_s=prep,
                solve_wall_s=solve_s,
                scalar_sync_read_wall_s=self.scalar_read_s,
                scalar_reads=self.scalar_reads,
                final_sync_read_wall_s=read_s,
                whole_step_wall_s=perf_counter() - start,
                cuda_activity_span_s=wp.get_event_elapsed_time(ev0, ev1) / 1000,
                interpretation="Solve includes host Armijo decisions and scalar waits; CUDA span contains those gaps, not pure kernel time. Final read includes operator evidence for audit.",
            ),
        )
