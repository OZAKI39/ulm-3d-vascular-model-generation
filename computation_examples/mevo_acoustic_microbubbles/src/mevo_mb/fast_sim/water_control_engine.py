"""Explicit F10N water-proxy admission, preserving R3 walls and full M0 dynamics.

R3 remains OFF-gated. This separately named adapter admits only the documented
water proxy and bounded interior domain; it does not modify historical physics.
"""

import numpy as np
import warp as wp
from .pdms_cavity_engine import PDMSCavitySimulator, cavity_project
from .f10c_engine import FullParticleBatch
from .engine import FastSimulator
from .straight_tube_engine import GPUZeroFlow
from .f5_engine import F5FlowSimulator
from .cubic_field import cubic_response


class WaterControlSimulator(PDMSCavitySimulator):
    def __init__(self, positions_m, physics, response, basis, *, cavity, **kwargs):
        if basis.metadata["model"] != "HOMOGENEOUS_WATER_CONTROL_PROXY":
            raise ValueError("F10N admits only the explicitly labelled water proxy")
        if basis.metadata["frequency_Hz"] != physics.frequency_Hz:
            raise ValueError("Frequency mismatch")
        if not all(
            cavity.query(q, physics.radius_m)["valid"] for q in np.asarray(positions_m).reshape(-1, 3)
        ):
            raise ValueError("Initial finite spheres outside cavity")
        self.cavity = cavity
        self._cavity_ready = False
        self._cavity_frames = None
        FullParticleBatch.__init__(
            self,
            positions_m,
            physics,
            response,
            basis,
            vessel_enabled=False,
            flow=GPUZeroFlow(),
            flow_enabled=False,
            acoustic_enabled=True,
            gradient_mode="consistent_trilinear",
            mobility_mode="f5_local",
            parameterized_m0=True,
            **kwargs,
        )
        self.cavity_stats = wp.zeros((self.total, 5), dtype=wp.float64, device=self.device)
        self._cavity_ready = True
        self.reset()
        wp.load_module(module=cavity_project.module, device=self.device)
        wp.load_module(module=cubic_response.module, device=self.device)

    def control(self, event):
        if event.get("flow_enabled", False):
            raise ValueError("F10N strictly zero background flow")
        self.basis.vector(event["amplitude_Pa"], event["phase_rad"], event.get("enabled"))
        return FullParticleBatch.control(self, event)

    def fields(self):
        if not self.acoustic_enabled:
            self.consistent_gradients.zero_()
            return F5FlowSimulator.fields(self)
        FastSimulator.launch(
            self,
            cubic_response,
            self.total,
            [
                self.x,
                self.grid,
                self.origin,
                self.spacing,
                *self.shape,
                self.knots,
                self.values,
                self.valid,
                self.kv,
                self.kval,
                self.harmonics,
                self.force,
                self.pressure,
                self.consistent_gradients,
                self.diag,
                self.error,
            ],
        )

    def run(self, duration_s, dt_s, events, **kwargs):
        for i, e in enumerate(events):
            end = events[i + 1]["time_s"] if i + 1 < len(events) else duration_s
            if "duration_s" not in e or not np.isclose(
                e["duration_s"], end - e["time_s"], rtol=0, atol=1e-14
            ):
                raise ValueError("Every control segment needs an exact positive duration")
            if e["duration_s"] <= 0:
                raise ValueError("Nonpositive duration")
        result = super().run(duration_s, dt_s, events, **kwargs)
        s = result["summary"]
        s.update(
            solver_version="F10N_water_proxy_R3_walls_full_M0",
            scope="HOMOGENEOUS_WATER_CONTROL_PROXY; all motion cavity assumed water; bounded interior",
            gradient_mode="C1_CATMULL_ROM_PRESSURE_EXACT_POLYNOMIAL_DERIVATIVE",
            training_ready=False,
            experiment_matched=False,
            statuses=["EXPLORATORY", "PHYSICS_UNCALIBRATED"],
        )
        s["model"].update(
            gradient=s["gradient_mode"],
            wall="R3 exposed union walls; sealed X; numerical scope exit, not a top wall",
            flow="strict zero",
        )
        return result
