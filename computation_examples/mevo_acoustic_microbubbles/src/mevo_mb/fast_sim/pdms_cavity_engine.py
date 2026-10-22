"""Opt-in four-source M0 adapter with exposed-union walls, no open-X semantics.

No changes to bubble forces, lubrication or PBD. Intersection curves use bounded
1D minimization, checked against the independent dense CPU oracle. Top scope
failures retain the valid prefix; the numerical ceiling is never projected.
"""

import numpy as np
import warp as wp
from .f10c_engine import FullParticleBatch
from .engine import FastSimulator
from .straight_tube_engine import GPUZeroFlow
from . import kernels as k


@wp.func
def intersection(theta: wp.float64, c: wp.float64, R: wp.float64, a: wp.float64):
    y = a * wp.sin(theta)
    return wp.vec3d(c + a * wp.cos(theta), y, wp.sqrt(R * R - y * y))


@wp.func
def curve_d2(theta: wp.float64, c: wp.float64, R: wp.float64, a: wp.float64, q: wp.vec3d):
    v = intersection(theta, c, R, a) - q
    return wp.dot(v, v)


@wp.func
def curve_nearest(c: wp.float64, R: wp.float64, a: wp.float64, q: wp.vec3d):
    step = wp.float64(6.283185307179586 / 128.0)
    best = intersection(wp.float64(0.0), c, R, a)
    delta = best - q
    best_d = wp.dot(delta, delta)
    for j in range(128):
        mid = wp.float64(j) * step
        val = curve_d2(mid, c, R, a, q)
        if val <= curve_d2(mid - step, c, R, a, q) and val <= curve_d2(mid + step, c, R, a, q):
            lo = mid - step
            hi = mid + step
            for it in range(40):
                t1 = hi - wp.float64(0.6180339887498949) * (hi - lo)
                t2 = lo + wp.float64(0.6180339887498949) * (hi - lo)
                if curve_d2(t1, c, R, a, q) <= curve_d2(t2, c, R, a, q):
                    hi = t2
                else:
                    lo = t1
            v = intersection((lo + hi) * wp.float64(0.5), c, R, a)
            d = v - q
            dd = wp.dot(d, d)
            if dd < best_d:
                best_d = dd
                best = v
    return best


@wp.func
def wall_query(q: wp.vec3d, R: wp.float64, a: wp.float64, L: wp.float64, H: wp.float64, c_abs: wp.float64):
    # Return signed distance, inward normal, wall point. All exposed patches
    # compete; a maximum of component SDFs would underestimate union clearance.
    rad = wp.sqrt(q[1] * q[1] + q[2] * q[2])
    best = wp.vec3d(-L, wp.float64(0.0), wp.float64(0.0))
    if rad > wp.float64(0.0):
        scale = wp.min(wp.float64(1.0), R / rad)
        best[1] = q[1] * scale
        best[2] = q[2] * scale
    best_d = wp.length(q - best)
    kind = int(1)
    v = best
    v[0] = L
    if wp.length(q - v) < best_d:
        best = v
        best_d = wp.length(q - v)
        kind = 2
    v = wp.vec3d(wp.clamp(q[0], -L, L), wp.float64(0.0), -R)
    if rad > wp.float64(0.0):
        v[1] = q[1] * R / rad
        v[2] = q[2] * R / rad
    hidden = bool(False)
    for i in range(2):
        c = wp.float64(2 * i - 1) * c_abs
        if v[2] >= wp.float64(0.0) and wp.sqrt((v[0] - c) * (v[0] - c) + v[1] * v[1]) < a - wp.float64(1e-14):
            hidden = True
    if not hidden and wp.length(q - v) < best_d:
        best = v
        best_d = wp.length(q - v)
        kind = 0
    inside = rad <= R
    for i in range(2):
        c = wp.float64(2 * i - 1) * c_abs
        rr = wp.sqrt((q[0] - c) * (q[0] - c) + q[1] * q[1])
        if rr <= a and q[2] >= wp.float64(0.0) and q[2] <= H:
            inside = True
        v = wp.vec3d(c + a, wp.float64(0.0), wp.clamp(q[2], wp.float64(0.0), H))
        if rr > wp.float64(0.0):
            v[0] = c + (q[0] - c) * a / rr
            v[1] = q[1] * a / rr
        if wp.sqrt(v[1] * v[1] + v[2] * v[2]) >= R - wp.float64(1e-14):
            dd = wp.length(q - v)
            if dd < best_d:
                best = v
                best_d = dd
                kind = 3 + i
        if wp.max(wp.abs(q[0] - c) - a, wp.float64(0.0)) <= best_d:
            v = curve_nearest(c, R, a, q)
            dd = wp.length(q - v)
            if dd < best_d:
                best = v
                best_d = dd
                kind = 5 + i
    inside = inside and wp.abs(q[0]) <= L
    sign = wp.float64(1.0)
    if not inside:
        sign = wp.float64(-1.0)
    normal = wp.vec3d(wp.float64(0.0), wp.float64(0.0), wp.float64(-1.0))
    if best_d > wp.float64(1e-16):
        normal = sign * (q - best) / best_d
    elif kind == 1:
        normal = wp.vec3d(wp.float64(1.0), wp.float64(0.0), wp.float64(0.0))
    elif kind == 2:
        normal = wp.vec3d(wp.float64(-1.0), wp.float64(0.0), wp.float64(0.0))
    elif kind == 0:
        normal = wp.vec3d(wp.float64(0.0), -best[1] / R, -best[2] / R)
    else:
        c = -c_abs
        if kind == 4 or kind == 6:
            c = c_abs
        normal = wp.vec3d((c - best[0]) / a, -best[1] / a, wp.float64(0.0))
        if kind >= 5:
            normal = wp.normalize(normal + wp.vec3d(wp.float64(0.0), -best[1] / R, -best[2] / R))
    return sign * best_d, normal, best


@wp.kernel
def query_walls(
    x: wp.array(dtype=wp.vec3d),
    R: wp.float64,
    a: wp.float64,
    L: wp.float64,
    H: wp.float64,
    c_abs: wp.float64,
    distances: wp.array(dtype=wp.float64),
    normals: wp.array(dtype=wp.vec3d),
    closest: wp.array(dtype=wp.vec3d),
):
    i = wp.tid()
    d, n, v = wall_query(x[i], R, a, L, H, c_abs)
    distances[i] = d
    normals[i] = n
    closest[i] = v


@wp.kernel
def cavity_project(
    x: wp.array(dtype=wp.vec3d),
    center: wp.vec3d,
    rotation: wp.mat33d,
    R: wp.float64,
    a: wp.float64,
    L: wp.float64,
    H: wp.float64,
    c_abs: wp.float64,
    radius: wp.float64,
    ceiling: wp.float64,
    stats: wp.array2d(dtype=wp.float64),
    error: wp.array(dtype=int),
):
    i = wp.tid()
    if error[0] != 0:
        return
    original = wp.transpose(rotation) * (x[i] - center)
    q = original
    if q[2] + radius >= ceiling:
        wp.atomic_or(error, 0, 131072)
        return
    converged = bool(False)
    for it in range(8):
        d, n, v = wall_query(q, R, a, L, H, c_abs)
        gap = d - radius
        stats[i, 0] = wp.min(stats[i, 0], gap)
        if gap >= wp.float64(-1e-12):
            converged = True
            break
        q = q + (-gap + wp.float64(2e-13)) * n
        if q[2] + radius >= ceiling:
            wp.atomic_or(error, 0, 131072)
            return
    correction = wp.length(q - original)
    stats[i, 1] = wp.max(stats[i, 1], correction)
    if not converged:
        wp.atomic_or(error, 0, 262144)
        return
    if correction > wp.float64(0.1) * radius:
        wp.atomic_or(error, 0, 1024)
        return
    if correction > wp.float64(0.0):
        x[i] = center + rotation * q
        stats[i, 2] = wp.max(stats[i, 2], correction)
        stats[i, 3] = stats[i, 3] + correction
        stats[i, 4] = stats[i, 4] + wp.float64(1.0)


class PDMSCavitySimulator(FullParticleBatch):
    def __init__(self, positions_m, physics, response, basis, *, cavity, **kwargs):
        # Formal ON is intentionally gated: unknown end seals and free surfaces.
        if basis is not None or kwargs.get("acoustic_enabled", True):
            raise ValueError("R3 apparatus ON not admitted: unknown interfaces and no four-field basis")
        if not all(
            cavity.query(q, physics.radius_m)["valid"] for q in np.asarray(positions_m).reshape(-1, 3)
        ):
            raise ValueError("Initial finite spheres violate union wall or interior validation scope")
        self.cavity = cavity
        self._cavity_ready = False
        self._cavity_frames = None
        super().__init__(
            positions_m,
            physics,
            response,
            basis,
            vessel_enabled=False,
            flow=GPUZeroFlow(),
            flow_enabled=False,
            gradient_mode="consistent_trilinear",
            mobility_mode="f5_local",
            parameterized_m0=True,
            **kwargs,
        )
        self.cavity_stats = wp.zeros((self.total, 5), dtype=wp.float64, device=self.device)
        self._cavity_ready = True
        self.reset()
        wp.load_module(module=cavity_project.module, device=self.device)

    def reset(self):
        super().reset()
        if self._cavity_ready:
            values = np.zeros((self.total, 5))
            values[:, 0] = [
                self.cavity.query(q, self.p.radius_m)["finite_radius_gap_m"]
                for q in self.initial.reshape(-1, 3)
            ]
            self.cavity_stats.assign(values)

    def control(self, event):
        if event.get("flow_enabled", False) or event.get("pressure_Pa", 0) != 0:
            raise ValueError("R3 validation admits strictly OFF and zero background flow")
        for key in ["amplitude_Pa", "phase_rad", "enabled"]:
            if key in event and len(event[key]) != 4:
                raise ValueError("Four independent source controls required")
        if any(event.get("amplitude_Pa", [0] * 4)):
            raise ValueError("Apparatus ON not admitted")
        # Uniform OFF parent accepts a scalar phase; array controls are validated
        # above and stored as four zero controls without passing a list to float().
        record = super().control(
            dict(time_s=event["time_s"], pressure_Pa=0.0, phase_rad=0.0, flow_enabled=False)
        )
        record.update(control_real_Pa=[0.0] * 4, control_imag_Pa=[0.0] * 4)
        return record

    def launch(self, kernel, dim, args):
        result = super().launch(kernel, dim, args)
        if self._cavity_ready and kernel in (k.predict, k.contact):
            buffer = args[2] if kernel is k.predict else args[1]
            c = self.cavity
            FastSimulator.launch(
                self,
                cavity_project,
                self.total,
                [
                    buffer,
                    wp.vec3d(*c.lumen.center),
                    wp.mat33d(*c.lumen.rotation.ravel()),
                    c.R,
                    c.a,
                    c.L,
                    c.H,
                    c.ports[1],
                    self.p.radius_m,
                    c.ceiling,
                    self.cavity_stats,
                    self.error,
                ],
            )
        return result

    def quality(self, sample=None):
        q = super().quality(sample)
        if self._cavity_ready:
            d = self.cavity_stats.numpy()
            if sample is not None:
                d = d[sample * self.n : (sample + 1) * self.n]
            q.update(
                minimum_proposed_union_gap_m=float(d[:, 0].min()),
                maximum_wall_proposal_m=float(d[:, 1].max()),
                maximum_wall_correction_m=float(d[:, 2].max()),
                cumulative_wall_correction_m=float(d[:, 3].sum()),
            )
        return q

    def state(self):
        x = super().state()
        if self._cavity_ready and self._cavity_frames is not None:
            self._cavity_frames.append(self.cavity_stats.numpy().copy())
        return x

    def run(self, *args, **kwargs):
        self._cavity_frames = []
        result = super().run(*args, **kwargs)
        if len(self._cavity_frames) != len(result["time_s"]):
            raise RuntimeError("Union diagnostic/frame mismatch")
        result["cavity_stats"] = np.asarray(self._cavity_frames)
        self._cavity_frames = None
        result["summary"].update(
            solver_version="F10M_R3_exposed_union_M0",
            end_conditions=["closed", "closed"],
            flow_physics="STRICT_ZERO; not residual flow, streaming or evaporation",
            no_CFD_loaded=True,
            wall_physics="EXPOSED_CYLINDER_PATCHES_AND_INTERSECTION_CURVES",
            scope="DEMO_INTERIOR_WETTING_ONLY",
        )
        result["summary"]["model"].update(
            wall="union with two sealed X ends; unknown top surface is fail-stop, never a wall", flow="zero"
        )
        return result
