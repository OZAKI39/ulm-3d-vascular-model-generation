"""Explicit F6 switches. Legacy default path remains byte-for-byte F5 execution."""

import numpy as np
import warp as wp
from .f5_engine import F5FlowSimulator
from .engine import FastSimulator
from . import kernels as k
from . import consistent_field as cf
from . import rpy


class F6Simulator(F5FlowSimulator):
    def __init__(
        self,
        *args,
        gradient_mode="consistent_trilinear",
        mobility_mode="f5_local",
        correction_mode="none",
        **kwargs,
    ):
        if gradient_mode not in ("legacy_nodal", "consistent_trilinear"):
            raise ValueError("Unknown gradient mode")
        if mobility_mode not in ("f5_local", "rpy_free_space", "diagonal_stokes"):
            raise ValueError("Unknown mobility mode")
        if correction_mode != "none":
            raise ValueError("UNCALIBRATED: no validated translational correction supplied")
        self.gradient_mode, self.mobility_mode = gradient_mode, mobility_mode
        self.correction_mode = correction_mode
        self._adaptive = None
        super().__init__(*args, **kwargs)
        self.consistent_gradients = wp.zeros((self.total, 3), dtype=wp.vec2d, device=self.device)
        wp.load_module(module=cf, device=self.device)
        wp.load_module(module=rpy, device=self.device)

    def fields(self):
        super().fields()
        if self.gradient_mode == "consistent_trilinear" and self.acoustic_enabled and self.basis is not None:
            FastSimulator.launch(
                self,
                cf.consistent_primary,
                self.total,
                [
                    self.x,
                    self.grid,
                    self.origin,
                    self.spacing,
                    wp.vec3d(*(self.basis.origin + self.basis.spacing * (np.array(self.shape) - 1))),
                    *self.shape,
                    self.active,
                    self.harmonics,
                    self.force,
                    self.pressure,
                    self.consistent_gradients,
                    self.error,
                ],
            )

    def launch(self, kernel, dim, args):
        if getattr(self, "_f5_ready", False):
            if kernel is k.free_velocity and self.mobility_mode == "rpy_free_space":
                if self.flow_enabled:
                    self.flow.launch(self.x, self.active, self.pending, self.flow_query, self.error)
                FastSimulator.launch(
                    self,
                    rpy.velocity,
                    self.total,
                    [
                        self.x,
                        self.force,
                        self.flow_query.velocity,
                        self.active,
                        self.free,
                        self.n,
                        self.p.radius_m,
                        self.p.viscosity_Pa_s,
                        int(self.flow_enabled),
                        self._current_dt,
                        self.free_path,
                        self.error,
                    ],
                )
                return
            if self.mobility_mode != "f5_local":
                if kernel is k.local_mobility:
                    FastSimulator.launch(self, rpy.copy_velocity, self.total, [self.free, args[3]])
                    return
                if kernel is k.component_links:
                    return  # independent roots; balance shift exactly zero
        return super().launch(kernel, dim, args)

    def advance(self, dt, count, graphs=True):
        if self._adaptive is not None:
            return self._adaptive.advance(dt, count, graphs)
        result = super().advance(dt, count, graphs)
        q = self.quality()
        if (
            q["max_accepted_penetration_m"] > 0.01 * self.p.radius_m
            or q["max_contact_pass_correction_m"] > 0.1 * self.p.radius_m
        ):
            self.failed_proposal = dict(positions_m=self.x.numpy(), attempted_steps=count, dt_s=dt)
            raise RuntimeError(
                "F6_REGISTERED_CONTACT_BUDGET: rejected sampled interval; last valid frame retained"
            )
        return result

    def run(self, *args, adaptive=False, adaptive_options=None, **kwargs):
        if adaptive:
            from .adaptive import StepDoubling

            self._adaptive = StepDoubling(self, **(adaptive_options or {}))
        result = super().run(*args, **kwargs)
        if self._adaptive is not None:
            result["physical_dt_s"] = np.array(self._adaptive.steps)
            result["adaptive_attempts"] = np.array(self._adaptive.attempts)
            result["summary"].update(
                actual_substep_count=len(self._adaptive.steps),
                macro_step_count=None,
                actual_trial_step_count=self._adaptive.trial_steps,
                rejected_attempts=sum(1 for v in self._adaptive.attempts if not v[4]),
                timestep_mode="adaptive_step_doubling",
                adaptive_error="maximum same-ID coarse/fine position difference; no error estimator claims at unresolved acoustic cells",
            )
            self._adaptive = None
        else:
            # Infer exact segment remainder steps, rather than duplicating PBD work.
            duration = args[0] if args else kwargs["duration_s"]
            dt = args[1] if len(args) > 1 else kwargs["dt_s"]
            events = args[2] if len(args) > 2 else kwargs["events"]
            times = []
            for i, e in enumerate(events):
                end = events[i + 1]["time_s"] if i + 1 < len(events) else duration
                span = end - e["time_s"]
                n = int(np.floor(span / dt + 1e-10))
                times.extend([dt] * n)
                if span - n * dt > max(1e-15, dt * 1e-10):
                    times.append(span - n * dt)
            result["physical_dt_s"] = np.array(times[: result["summary"]["actual_substep_count"]])
            result["summary"].update(timestep_mode="fixed_event_aligned", rejected_attempts=0)
        if hasattr(self, "failed_proposal"):
            result["failure_proposal"] = self.failed_proposal
        s = result["summary"]
        s.update(
            solver_version="fast_pzt_f6_v1",
            gradient_mode=self.gradient_mode,
            mobility_mode=self.mobility_mode,
            correction_mode=self.correction_mode,
            evidence_level="EXPLORATORY",
            translation_calibration="UNCALIBRATED",
        )
        s["model"].update(
            gradient=self.gradient_mode,
            mobility=self.mobility_mode,
            local_lubrication_applied=self.mobility_mode == "f5_local",
            contact_metric="inherited diagonal central PBD; not full RPY constrained dynamics",
        )
        if self.mobility_mode == "rpy_free_space":
            s["mobility_scope"] = (
                "Free-space equal-radius no-slip RPY; old lubrication bypassed; no wall images or stresslets; no validated dense bubble hydrodynamics"
            )
            s["free_path_definition"] = "Integral of |u_flow + M_RPY F_ac| before contact"
        return result
