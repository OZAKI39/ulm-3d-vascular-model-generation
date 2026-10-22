"""Opt-in analytic cylinder adapter around the unchanged full-particle M0 path."""

import numpy as np
import warp as wp
from .f10c_engine import FullParticleBatch
from .engine import FastSimulator
from . import kernels as k


@wp.kernel
def zero_flow(
    velocity: wp.array(dtype=wp.vec3d),
    gradient: wp.array(dtype=wp.mat33d),
    valid: wp.array(dtype=int),
    cache: wp.array(dtype=int),
):
    i = wp.tid()
    velocity[i] = wp.vec3d(wp.float64(0.0))
    gradient[i] = wp.mat33d(wp.float64(0.0))
    valid[i] = 1
    cache[i] = -1


class GPUZeroFlow:
    metadata = dict(mode="zero", cfd_files=[], exact=True, gradient="zero")

    def launch(self, x, active, pending, buffers, error, allow_terminal=False):
        wp.launch(
            zero_flow,
            dim=len(x),
            inputs=[buffers.velocity, buffers.gradient, buffers.valid, buffers.cache],
            device=x.device,
        )


@wp.kernel
def cylinder_project(
    x: wp.array(dtype=wp.vec3d),
    center: wp.vec3d,
    rotation: wp.mat33d,
    lumen_radius: wp.float64,
    half_length: wp.float64,
    bubble_radius: wp.float64,
    closed_minus: int,
    closed_plus: int,
    stats: wp.array2d(dtype=wp.float64),
    error: wp.array(dtype=int),
):
    i = wp.tid()
    if error[0] != 0:
        return
    q = wp.transpose(rotation) * (x[i] - center)
    r = wp.sqrt(q[1] * q[1] + q[2] * q[2])
    gap = lumen_radius - bubble_radius - r
    stats[i, 0] = wp.min(stats[i, 0], gap)
    next_q = q
    if gap < wp.float64(0.0):
        next_q[1] = q[1] * (lumen_radius - bubble_radius) / r
        next_q[2] = q[2] * (lumen_radius - bubble_radius) / r
    if closed_minus == 1:
        next_q[0] = wp.max(next_q[0], -half_length + bubble_radius)
    elif q[0] < -half_length:
        # Open endpoint: terminate at the last valid saved frame; never reflect.
        wp.atomic_or(error, 0, 32768)
        return
    if closed_plus == 1:
        next_q[0] = wp.min(next_q[0], half_length - bubble_radius)
    elif q[0] > half_length:
        wp.atomic_or(error, 0, 32768)
        return
    correction = wp.length(next_q - q)
    stats[i, 1] = wp.max(stats[i, 1], correction)
    if correction > wp.float64(0.1) * bubble_radius:
        wp.atomic_or(error, 0, 1024)
        return
    if correction > wp.float64(0.0):
        x[i] = center + rotation * next_q
        stats[i, 2] = wp.max(stats[i, 2], correction)
        stats[i, 3] = stats[i, 3] + correction
        stats[i, 4] = stats[i, 4] + wp.float64(1.0)


class StraightTubeSimulator(FullParticleBatch):
    def __init__(self, positions_m, physics, response, basis, *, tube, **kwargs):
        if basis is not None and len(basis.db.pressure_basis) != 3:
            raise ValueError("Exactly three newly computed PZT fields required")
        if not tube.contains(np.asarray(positions_m), physics.radius_m).all():
            raise ValueError("Initial finite-radius bubble violates tube")
        self.tube = tube
        self._straight_ready = False
        self._cylinder_frames = None
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
        self.cylinder_stats = wp.zeros((self.total, 5), dtype=wp.float64, device=self.device)
        self._straight_ready = True
        self.reset()
        wp.load_module(module=cylinder_project.module, device=self.device)

    def reset(self):
        super().reset()
        if self._straight_ready:
            values = np.zeros((self.total, 5))
            values[:, 0] = self.tube.gaps(self.initial.reshape(-1, 3), self.p.radius_m)[:, 0]
            self.cylinder_stats.assign(values)

    def control(self, event):
        if event.get("flow_enabled", False):
            raise ValueError("Straight-tube protocol forbids enabling background flow")
        for name in ("amplitude_Pa", "phase_rad", "enabled"):
            if name in event and len(event[name]) != 3:
                raise ValueError("Control dimension must be three")
        record = super().control(dict(event, flow_enabled=False))
        if not self.acoustic_enabled:
            record.update(control_real_Pa=[0.0] * 3, control_imag_Pa=[0.0] * 3)
        return record

    def launch(self, kernel, dim, args):
        result = super().launch(kernel, dim, args)
        if self._straight_ready and kernel in (k.predict, k.contact):
            buffer = args[2] if kernel is k.predict else args[1]
            t = self.tube
            FastSimulator.launch(
                self,
                cylinder_project,
                self.total,
                [
                    buffer,
                    wp.vec3d(*t.center),
                    wp.mat33d(*t.rotation.ravel()),
                    t.radius_m,
                    t.length_m / 2,
                    self.p.radius_m,
                    int(t.ends[0] == "closed"),
                    int(t.ends[1] == "closed"),
                    self.cylinder_stats,
                    self.error,
                ],
            )
        return result

    def quality(self, sample=None):
        q = super().quality(sample)
        if self._straight_ready:
            d = self.cylinder_stats.numpy()
            if sample is not None:
                d = d[sample * self.n : (sample + 1) * self.n]
            q.update(
                minimum_proposed_radial_gap_m=float(d[:, 0].min()),
                maximum_wall_proposal_m=float(d[:, 1].max()),
                maximum_wall_correction_m=float(d[:, 2].max()),
                cumulative_wall_correction_m=float(d[:, 3].sum()),
            )
        return q

    def run(self, *args, **kwargs):
        self._cylinder_frames = []
        result = super().run(*args, **kwargs)
        if len(self._cylinder_frames) != len(result["time_s"]):
            raise RuntimeError("Cylinder diagnostic/frame misalignment")
        result["cylinder_stats"] = np.asarray(self._cylinder_frames)
        self._cylinder_frames = None
        result["summary"].update(
            solver_version="F10L_analytic_tube_M0",
            flow_physics="STRICT_ZERO",
            wall_physics="ANALYTIC_FINITE_RADIUS_CYLINDER",
            end_conditions=list(self.tube.ends),
            no_CFD_loaded=True,
        )
        result["summary"]["model"].update(
            flow="zero", wall="analytic cylinder; open endpoints fail-stop, no reflection"
        )
        return result

    def state(self):
        x = super().state()
        if self._straight_ready and self._cylinder_frames is not None:
            self._cylinder_frames.append(self.cylinder_stats.numpy().copy())
        return x
