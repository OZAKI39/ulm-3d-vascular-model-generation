"""Independent M0 centroid proxies: no pair-force, PBD, or mobility solver.

Only the centroid evolves. Virtual samples are never returned as particle states.
Known-state probes and independent rollouts are deliberately separate APIs.
"""

from time import perf_counter
import numpy as np
import warp as wp
from . import kernels as k
from . import consistent_field as cf
from .engine import complex_vectors
from .flow_gpu import QueryBuffers
from .metrics import frame_metrics

VERSION = "f10b_centroid_v1"
LABELS = {"G0": "CENTROID_SURROGATE_G0", "G1": "FROZEN_TEMPLATE_SURROGATE_G1"}


def instantaneous(x, event, objects, mode):
    """Read-only local prediction at a supplied state; never used in rollouts."""
    p, b, r, g, f = objects
    samples = x.mean(0)[None] if mode == "G0" else np.asarray(x)
    field, _, _ = b.combine(event["amplitude_Pa"], event["phase_rad"], event["enabled"])
    pressure, grad = cf.query_numpy(field.pressure_Pa, b.origin, b.spacing, samples)
    force = -0.5 * np.real(r.query(pressure)[:, 0, None] * grad.conj())
    flow = f.sample(samples)
    if flow["error"]:
        raise ValueError("CFD domain invalid")
    return (flow["velocity_m_s"] + force / (6 * np.pi * p.viscosity_Pa_s * p.radius_m)).mean(0)


def admit(
    x,
    event,
    objects,
    limits,
    *,
    mobility,
    duration_s,
    maturity,
    validated_initial,
    validated_control,
    all_active=True,
    morphology_requested=False,
    known_contact_rejection=False,
):
    """Conservative static admission. Flags carry audited history, not predictions."""
    p, b, r, g, f = objects
    x = np.asarray(x)
    reasons = []
    if x.shape != (32, 3) or not np.isfinite(x).all():
        return dict(status="SURROGATE_NOT_ADMITTED", reasons=["INVALID_32_STATE"])
    m = frame_metrics(x, p.radius_m)
    wall = g.query(x)
    gap = float((wall["distance_m"] - p.radius_m).min())
    flow = f.sample(x)
    flow_enabled = bool(event.get("flow_enabled", True))
    if not flow_enabled:
        # Static scope reports the actual scenario, not an inactive CFD field.
        flow["velocity_m_s"] = np.zeros_like(x)
        flow["gradient_s_inv"] = np.zeros((len(x), 3, 3))
    tensor = flow["gradient_s_inv"]
    strain = (tensor + tensor.transpose(0, 2, 1)) / 2
    shear = float(np.sqrt(2 * np.sum(strain**2, axis=(1, 2))).max())
    pressure_max = 0.0
    table_valid = True
    try:
        field, _, _ = b.combine(event["amplitude_Pa"], event["phase_rad"], event["enabled"])
        pressure, _ = cf.query_numpy(field.pressure_Pa, b.origin, b.spacing, x)
        r.query(pressure)
        pressure_max = float(abs(pressure).max())
    except ValueError:
        table_valid = False
        reasons.append("FIELD_OR_TABLE_INVALID")
    speed = float(np.linalg.norm(flow["velocity_m_s"], axis=1).max())
    # A conservative forecast plus a strict accumulated-path runtime bound.
    forecast_gap = gap - 1.5 * speed * duration_s
    checks = {
        "MOBILITY_NOT_M0": mobility == "f5_local",
        "FLOW_CONFIGURATION_UNVALIDATED": flow_enabled,
        "RADIUS_NOT_EQUAL_PRODUCTION": p.radius_m == limits["radius_m"],
        "HORIZON_UNVALIDATED": 0 < duration_s <= limits["max_horizon_s"],
        "NO_MATURITY_HISTORY": bool(maturity),
        "INITIAL_STATE_UNVALIDATED": bool(validated_initial),
        "CONTROL_UNVALIDATED": bool(validated_control),
        "INACTIVE_OR_EXIT": bool(all_active),
        "MORPHOLOGY_REQUIRED": not morphology_requested,
        "KNOWN_CONTACT_REJECTION": not known_contact_rejection,
        "DISPERSED": m["largest_size"] == 32,
        "SIZE_OUTSIDE_SCOPE": limits["rg_range_m"][0] <= m["global_rg_m"] <= limits["rg_range_m"][1],
        "WALL_NEAR": gap >= limits["minimum_wall_gap_m"] and wall["inside"].all(),
        "WALL_FORECAST_UNSAFE": forecast_gap >= limits["minimum_wall_gap_m"],
        "HIGH_SHEAR": shear <= limits["maximum_shear_s_inv"],
        "HIGH_PRESSURE": pressure_max <= limits["maximum_pressure_Pa"],
        "CFD_INVALID": not flow["error"] and flow["valid"].all(),
    }
    reasons += [key for key, ok in checks.items() if not ok]
    return dict(
        status="ADMITTED_FOR_SHADOW_TEST" if not reasons else "SURROGATE_NOT_ADMITTED",
        reasons=reasons,
        minimum_wall_gap_m=gap,
        forecast_gap_m=forecast_gap,
        maximum_shear_s_inv=shear,
        maximum_pressure_Pa=pressure_max,
        initial_rg_m=m["global_rg_m"],
        largest_fraction=m["largest_size"] / 32,
        table_valid=table_valid,
        checks=checks,
        scope="STATIC_ADMISSION; no trajectory generated",
        training_ready=False,
        flow_enabled=flow_enabled,
    )


@wp.kernel
def virtual_samples(
    center: wp.array(dtype=wp.vec3d), offsets: wp.array(dtype=wp.vec3d), x: wp.array(dtype=wp.vec3d)
):
    i = wp.tid()
    x[i] = center[0] + offsets[i]


@wp.kernel
def centroid_step(
    center: wp.array(dtype=wp.vec3d),
    force: wp.array(dtype=wp.vec3d),
    flow: wp.array(dtype=wp.vec3d),
    gradients: wp.array(dtype=wp.mat33d),
    pressure: wp.array(dtype=wp.vec2d),
    n: int,
    drag: wp.float64,
    dt: wp.float64,
    error: wp.array(dtype=int),
    path: wp.array(dtype=wp.float64),
    path_limit: wp.float64,
    pmax: wp.float64,
    shearmax: wp.float64,
    lower: wp.vec3d,
    upper: wp.vec3d,
    diagnostics: wp.array(dtype=wp.float64),
):
    if error[0] != 0:
        return
    velocity = wp.vec3d(wp.float64(0.0))
    for i in range(n):
        amp = wp.length(pressure[i])
        diagnostics[0] = wp.max(diagnostics[0], amp)
        strain = (gradients[i] + wp.transpose(gradients[i])) * wp.float64(0.5)
        shear = wp.sqrt(wp.float64(2.0) * wp.ddot(strain, strain))
        diagnostics[1] = wp.max(diagnostics[1], shear)
        if amp > pmax or shear > shearmax:
            wp.atomic_or(error, 0, 65536)
        velocity = velocity + (flow[i] + force[i] / drag) / wp.float64(n)
    q = center[0] + dt * velocity
    travel = path[0] + wp.length(dt * velocity)
    for axis in range(3):
        if not wp.isfinite(q[axis]) or q[axis] < lower[axis] or q[axis] > upper[axis]:
            wp.atomic_or(error, 0, 131072)
    if travel > path_limit:
        wp.atomic_or(error, 0, 262144)
    if error[0] == 0:
        center[0] = q
        path[0] = travel


class CentroidSurrogate:
    """Independent Euler integration; constructor accepts initial state only."""

    def __init__(self, x, physics, response, basis, flow, event, mode, limits, admission, device="cuda:0"):
        if mode not in LABELS or admission["status"] != "ADMITTED_FOR_SHADOW_TEST":
            raise ValueError("SURROGATE_NOT_ADMITTED")
        self.device = device
        self.p = physics
        self.r = response
        self.b = basis
        self.flow = flow
        self.mode = mode
        self.limits = limits
        x = np.asarray(x)
        self.initial_center = x.mean(0)
        self.template = x - self.initial_center
        offsets = np.zeros((1, 3)) if mode == "G0" else self.template
        self.n = len(offsets)

        def arr(a, dtype):
            return wp.array(np.ascontiguousarray(a), dtype=dtype, device=device)

        self.center = arr(self.initial_center[None], wp.vec3d)
        self.offsets = arr(offsets, wp.vec3d)
        self.x = wp.empty(self.n, dtype=wp.vec3d, device=device)
        self.force = wp.zeros_like(self.x)
        self.harmonics = wp.zeros((self.n, 8), dtype=wp.vec2d, device=device)
        self.pressure = wp.zeros(self.n, dtype=wp.vec2d, device=device)
        self.gradient = wp.zeros((self.n, 3), dtype=wp.vec2d, device=device)
        self.diag = wp.zeros((self.n, 11), dtype=wp.float64, device=device)
        self.error = wp.zeros(1, dtype=int, device=device)
        self.active = wp.ones(self.n, dtype=int, device=device)
        self.pending = wp.zeros(self.n, dtype=int, device=device)
        self.path = wp.zeros(1, dtype=wp.float64, device=device)
        self.diagnostics = wp.zeros(2, dtype=wp.float64, device=device)
        self.flow_query = QueryBuffers(self.n, device)
        self.uniform = wp.zeros(1, dtype=wp.vec2d, device=device)
        self.knots = arr(response.knots, wp.float64)
        self.values = arr(complex_vectors(response.values), wp.vec2d)
        self.valid = arr(response.valid, wp.int32)
        self.kv = arr(complex_vectors(response.knot_values), wp.vec2d)
        self.kval = arr(response.knot_valid, wp.int32)
        _, grid, _ = basis.combine(event["amplitude_Pa"], event["phase_rad"], event["enabled"])
        self.grid = arr(grid.reshape(-1, 4, 2), wp.vec2d)
        self.origin = wp.vec3d(*basis.origin)
        self.spacing = wp.vec3d(*basis.spacing)
        self.upper = wp.vec3d(*(basis.origin + basis.spacing * (np.array(basis.shape) - 1)))
        extent = np.abs(self.template).max(0) + physics.radius_m + 2e-8
        self.center_lower = wp.vec3d(*(basis.origin + extent))
        self.center_upper = wp.vec3d(*(np.array(self.upper) - extent))
        self.path_limit = min(
            limits["maximum_centroid_travel_m"],
            admission["minimum_wall_gap_m"] - limits["minimum_wall_gap_m"],
        )
        for module in [k, cf, virtual_samples.module]:
            wp.load_module(module=module, device=device)
        wp.synchronize_device(device)

    def launch(self, kernel, dim, args):
        wp.launch(kernel, dim=dim, inputs=args, device=self.device, block_dim=32)

    def sample(self):
        self.launch(virtual_samples, self.n, [self.center, self.offsets, self.x])
        self.launch(
            k.field_response,
            self.n,
            [
                self.x,
                self.grid,
                self.origin,
                self.spacing,
                *self.b.shape,
                0,
                self.uniform,
                self.knots,
                self.values,
                self.valid,
                self.kv,
                self.kval,
                self.harmonics,
                self.force,
                self.pressure,
                self.diag,
                self.error,
            ],
        )
        self.launch(
            cf.consistent_primary,
            self.n,
            [
                self.x,
                self.grid,
                self.origin,
                self.spacing,
                self.upper,
                *self.b.shape,
                self.active,
                self.harmonics,
                self.force,
                self.pressure,
                self.gradient,
                self.error,
            ],
        )
        self.flow.launch(self.x, self.active, self.pending, self.flow_query, self.error)

    def step(self, dt):
        self.sample()
        self.launch(
            centroid_step,
            1,
            [
                self.center,
                self.force,
                self.flow_query.velocity,
                self.flow_query.gradient,
                self.pressure,
                self.n,
                6 * np.pi * self.p.viscosity_Pa_s * self.p.radius_m,
                dt,
                self.error,
                self.path,
                self.path_limit,
                self.limits["maximum_pressure_Pa"],
                self.limits["maximum_shear_s_inv"],
                self.center_lower,
                self.center_upper,
                self.diagnostics,
            ],
        )

    def prepare(self, dt, steps):
        with wp.ScopedCapture(device=self.device) as cap:
            for _ in range(steps):
                self.step(dt)
        self.graph = cap.graph

    def read_center(self):
        code = int(self.error.numpy()[0])
        if code:
            raise RuntimeError("SURROGATE_NOT_ADMITTED runtime error=" + str(code))
        return self.center.numpy()[0].copy()

    def rollout(self, duration, dt, sample_interval):
        steps = round(sample_interval / dt)
        count = round(duration / sample_interval)
        if (
            not np.isclose(steps * dt, sample_interval, rtol=0, atol=1e-15)
            or duration > self.limits["max_horizon_s"]
        ):
            raise ValueError("Unsupported time protocol")
        started = perf_counter()
        self.prepare(dt, steps)
        prepare = perf_counter() - started
        values = [self.read_center()]
        sync = 0.0
        integration = 0.0
        for _ in range(count):
            start = perf_counter()
            wp.capture_launch(self.graph)
            wp.synchronize_device(self.device)
            integration += perf_counter() - start
            start = perf_counter()
            values.append(self.read_center())
            sync += perf_counter() - start
        self.sample()
        self.read_center()  # final endpoint domain/table validation
        return dict(
            time_s=np.arange(count + 1) * sample_interval,
            centroid_m=np.asarray(values),
            label=LABELS[self.mode],
            training_ready=False,
            template_label="FROZEN_TEMPLATE_ASSUMPTION",
            timing=dict(graph_prepare_s=prepare, integration_s=integration, sync_readback_s=sync),
            quality=dict(
                error_bits=int(self.error.numpy()[0]),
                centroid_path_m=float(self.path.numpy()[0]),
                maximum_pressure_Pa=float(self.diagnostics.numpy()[0]),
                maximum_shear_s_inv=float(self.diagnostics.numpy()[1]),
            ),
        )
