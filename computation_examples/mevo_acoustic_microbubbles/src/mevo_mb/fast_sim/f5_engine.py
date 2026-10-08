"""F4 confinement plus prescribed flow before local mobility and contact.

The inherited local Jacobi mobility approximates (D + L_pp) U = D u_f + F_ac.
Its component mean balance is retained. No rigid body constraint, artificial
repulsion, inertial motion, or post-contact advection is introduced.
"""

import numpy as np
import warp as wp
from .f4_engine import F4VesselSimulator
from .f3_engine import F3Simulator
from .engine import FastSimulator
from . import kernels as k
from . import vessel_kernels as vk
from . import vessel_active_kernels as masked
from .flow_gpu import QueryBuffers


@wp.kernel
def acoustic_zero(
    force: wp.array(dtype=wp.vec3d),
    pressure: wp.array(dtype=wp.vec2d),
    harmonics: wp.array2d(dtype=wp.vec2d),
    diag: wp.array2d(dtype=wp.float64),
):
    i = wp.tid()
    force[i] = wp.vec3d(wp.float64(0.0))
    pressure[i] = wp.vec2d(wp.float64(0.0))
    for h in range(8):
        harmonics[i, h] = wp.vec2d(wp.float64(0.0))
    diag[i, 0] = wp.float64(0.0)


@wp.kernel
def flow_free_velocity(
    force: wp.array(dtype=wp.vec3d),
    flow: wp.array(dtype=wp.vec3d),
    out: wp.array(dtype=wp.vec3d),
    active: wp.array(dtype=int),
    drag: wp.float64,
    enabled: int,
    dt: wp.float64,
    free_path: wp.array(dtype=wp.float64),
):
    i = wp.tid()
    u = force[i] / drag
    if enabled == 1:
        u = u + flow[i]
    if active[i] == 0:
        u = wp.vec3d(wp.float64(0.0))
    out[i] = u
    free_path[i] = free_path[i] + wp.length(u) * dt


class F5FlowSimulator(F4VesselSimulator):
    def __init__(
        self,
        positions_m,
        physics,
        response,
        basis=None,
        *,
        flow=None,
        acoustic_enabled=True,
        flow_enabled=True,
        **kwargs,
    ):
        self.flow = flow
        self.flow_enabled = bool(flow_enabled and flow is not None)
        self.acoustic_enabled = bool(acoustic_enabled)
        self.secondary_enabled = bool(acoustic_enabled)
        self._flow_frames = None
        self._f5_ready = False
        self.switch_records = []
        # OFF-only continuations may start outside the acoustic grid. Build F4
        # without an active field check, then attach the static basis storage.
        super().__init__(positions_m, physics, response, basis if acoustic_enabled else None, **kwargs)
        if not acoustic_enabled and basis is not None:
            self.basis = basis
            self.shape = basis.shape
            self.origin = wp.vec3d(*basis.origin)
            self.spacing = wp.vec3d(*basis.spacing)
            self.grid = wp.zeros((int(np.prod(self.shape)), 4), dtype=wp.vec2d, device=self.device)
        if not self.vessel_enabled:
            self.active = wp.ones(self.total, dtype=int, device=self.device)
            self.pending = wp.zeros(self.total, dtype=int, device=self.device)
        self.field_enabled = int(self.acoustic_enabled and basis is not None)
        if basis is not None:
            margin = self.p.radius_m + 2e-8
            self.field_lower = wp.vec3d(*(basis.origin + margin))
            self.field_upper = wp.vec3d(
                *(basis.origin + basis.spacing * (np.array(basis.shape) - 1) - margin)
            )
        self.flow_query = QueryBuffers(self.total, self.device)
        self.flow_observation = QueryBuffers(self.total, self.device)
        self._f5_ready = True

    def reset(self):
        super().reset()
        if self._f5_ready:
            self.flow_query.reset()
            self.flow_observation.reset()

    def control(self, event):
        if self.basis is not None:
            amplitudes = np.asarray(event["amplitude_Pa"])
            enabled = np.asarray(event.get("enabled", [True] * len(amplitudes)))
            acoustic = bool(np.any(amplitudes * enabled != 0))
        else:
            acoustic = event.get("pressure_Pa", 0) != 0
        secondary = bool(acoustic and event.get("secondary_enabled", True))
        flow = bool(event.get("flow_enabled", self.flow_enabled))
        if flow and self.flow is None:
            raise ValueError("Flow ON requires a resident provider")
        mode = (
            "ACOUSTIC_ON_FLOW_" + ("ON" if flow else "OFF")
            if secondary
            else "SECONDARY_OFF_PZT_ON"
            if acoustic
            else "ULTRASOUND_OFF_FLOW_" + ("ON" if flow else "OFF")
        )
        if "mode" in event and event["mode"] != mode:
            raise ValueError("Control label disagrees with actual enabled forces")
        previous = (self.acoustic_enabled, self.secondary_enabled, self.flow_enabled)
        current = (acoustic, secondary, flow)
        before = self.x.numpy().copy()
        if current != previous:
            self.graphs.clear()
            self.switch_records.append(
                dict(
                    time_s=event["time_s"],
                    previous=list(previous),
                    current=list(current),
                    mode=mode,
                    CUDA_graphs_invalidated=True,
                )
            )
        self.acoustic_enabled, self.secondary_enabled, self.flow_enabled = current
        self.field_enabled = int(acoustic and self.basis is not None)
        if acoustic:
            if self.basis is not None:
                active = self.active.numpy().astype(bool)
                x = before[active]
                if (x < np.asarray(self.field_lower)).any() or (x > np.asarray(self.field_upper)).any():
                    raise RuntimeError("Cannot enable acoustics outside stored field; no extrapolation")
            record = super().control(event)
        else:
            record = dict(event, control_real_Pa=[0.0] * 4, control_imag_Pa=[0.0] * 4)
        np.testing.assert_array_equal(self.x.numpy(), before)
        return dict(
            record,
            mode=mode,
            acoustic_enabled=acoustic,
            secondary_enabled=secondary,
            flow_enabled=flow,
            acoustic_domain_required=bool(self.field_enabled),
        )

    def fields(self):
        if self.acoustic_enabled:
            return super().fields()
        FastSimulator.launch(
            self, acoustic_zero, self.total, [self.force, self.pressure, self.harmonics, self.diag]
        )

    def pairs(self):
        if self.secondary_enabled:
            return super().pairs()

    def launch(self, kernel, dim, args):
        if not self._f5_ready:
            return super().launch(kernel, dim, args)
        if kernel is k.free_velocity:
            if self.flow_enabled:
                self.flow.launch(self.x, self.active, self.pending, self.flow_query, self.error)
            FastSimulator.launch(
                self,
                flow_free_velocity,
                self.total,
                [
                    self.force,
                    self.flow_query.velocity,
                    self.free,
                    self.active,
                    args[2],
                    int(self.flow_enabled),
                    self._current_dt,
                    self.free_path,
                ],
            )
            return
        if kernel is k.accept:
            x = args[0]
            if self.vessel_enabled:
                g = self.geometry
                # F4's exact open-terminal and field-exit sequence is preserved.
                F3Simulator.launch(
                    self,
                    vk.clip_boundaries,
                    self.total,
                    [
                        g.closed.id,
                        g.tags,
                        g.origin,
                        g.padding,
                        self.x,
                        x,
                        self.active,
                        self.error,
                        g.cap_origins,
                        g.cap_normals,
                        self.field_lower,
                        self.field_upper,
                        self.field_enabled,
                        self.pending,
                        self.fractions,
                    ],
                )
                F3Simulator.launch(
                    self,
                    vk.final_geometry_check,
                    self.total,
                    [
                        *g.args(),
                        self.x,
                        x,
                        self.n,
                        self.p.radius_m,
                        self.active,
                        self.cache_x,
                        self.cache_d,
                        self.wall_stats,
                        self.error,
                        g.cap_origins,
                        g.cap_normals,
                        self.field_lower,
                        self.field_upper,
                        self.field_enabled,
                        int(self.wall_contact_enabled),
                        self.pending,
                        self.fractions,
                        self.wall_gap,
                    ],
                )
            # Validate proposed accepted positions too: a final-step CFD domain
            # failure must not be hidden by ending the run. Only already verified
            # physical terminal crossings may bypass this lookup.
            if self.flow_enabled:
                self.flow.launch(
                    x,
                    self.active,
                    self.pending,
                    self.flow_query,
                    self.error,
                    allow_terminal=self.vessel_enabled,
                )
            if self.vessel_enabled:
                F3Simulator.launch(
                    self,
                    vk.commit_events,
                    self.total,
                    [
                        self.active,
                        self.pending,
                        self.fractions,
                        self.exit_label,
                        self.exit_time,
                        self.wall_stats,
                        self.error,
                        self.step_count,
                        self.time_origin,
                        self._current_dt,
                    ],
                )
            F3Simulator.launch(self, masked.accept, dim, [*args, self.active])
            return
        return super().launch(kernel, dim, args)

    def state(self):
        try:
            x = super().state()
        except RuntimeError as exc:
            raise RuntimeError(str(exc) + "; F5 bit 16384=FLOW_DOMAIN_INVALID (no extrapolation)") from exc
        if self._f5_ready and self._flow_frames is not None:
            if self.flow is not None:
                self.flow.launch(self.x, self.active, self.pending, self.flow_observation, self.error)
                code = int(self.error.numpy()[0])
                if code:
                    raise RuntimeError("Invalid saved flow query, bit=" + str(code))
            b = self.flow_observation
            u = b.velocity.numpy()
            valid = b.valid.numpy().astype(bool)
            u[~valid] = np.nan
            self._flow_frames.append(
                dict(
                    velocity_m_s=u.reshape(self.batch, self.n, 3),
                    gradient_s_inv=b.gradient.numpy().reshape(self.batch, self.n, 3, 3),
                    valid=valid.reshape(self.batch, self.n),
                    tetra_id=b.cache.numpy().reshape(self.batch, self.n),
                    stats=self.flow_query.stats.numpy().reshape(self.batch, self.n, 4),
                    flow_enabled=self.flow_enabled,
                    acoustic_enabled=self.acoustic_enabled,
                    secondary_enabled=self.secondary_enabled,
                )
            )
        return x

    def run(self, *args, **kwargs):
        self._flow_frames = []
        result = super().run(*args, **kwargs)
        # Invalid final diagnostics can fail before an observation is appended;
        # keep only matched successful frames, never invent missing samples.
        if len(self._flow_frames) != len(result["time_s"]):
            raise RuntimeError("Flow/position sampling misalignment")
        result["flow"] = {k: np.asarray([r[k] for r in self._flow_frames]) for k in self._flow_frames[0]}
        # Event-aligned labels at an exact switch refer to the subsequent state;
        # the positions there are continuous and still the same accepted sample.
        controls = result["summary"]["controls"]
        for key in ["flow_enabled", "acoustic_enabled", "secondary_enabled"]:
            result["flow"][key] = np.array([controls[i][key] for i in result["control_index"]])
        result["summary"].update(
            flow_physics="PRESCRIBED_STEADY_FLOW",
            physics_scope="EXPLORATORY",
            flow_provider=None if self.flow is None else self.flow.metadata,
            switch_records=self.switch_records,
            flow_query_counts=self.flow_query.stats.numpy().sum(0).tolist(),
            free_path_definition="Integral of |u_flow + F_ac/(6*pi*mu*R)| before local mobility/contact",
            mobility_scope="Two inherited local Jacobi iterations plus component mean balance; pair lubrication from ambient differences is approximate, not a global Stokes solve",
            wall_physics="GEOMETRIC_WALL_ONLY",
        )
        result["summary"]["model"].update(
            flow="prescribed, event controlled",
            secondary="event controlled",
            acoustic_domain="required only with nonzero enabled PZT drive",
        )
        if int(self.error.numpy()[0]) & 16384:
            result["flow_domain_failure"] = dict(
                attempted_positions_m=self.flow_query.failed_x.numpy(),
                particle_failure_count=self.flow_query.stats.numpy()[:, 3],
                time_bracket_s=[
                    float(result["time_s"][-1]),
                    float(result["time_s"][-1] + kwargs.get("sample_interval_s", 0.0005)),
                ],
            )
        self._flow_frames = None
        return result
