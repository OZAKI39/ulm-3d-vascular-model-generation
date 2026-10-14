"""F10O interior adapter: original F10H CPU B2, F10N forces, R3 GPU walls.

No new contact or resistance implementation. A wall correction rejects the step,
exactly as in the historical interior-only B2 reference. No PBD after B2.
"""

import time
import numpy as np
import warp as wp
from .water_control_engine import WaterControlSimulator
from .pdms_cavity_engine import cavity_project
from .engine import FastSimulator
from .f10g_contact import resistance
from .f10h_contact import geometric_contact

KC = 0.04641319518378658
ATOL = 1e-13  # F10K trajectory reference, not the tighter fixed-state scaling test


class WaterB2:
    def __init__(self, initial, physics, response, basis, cavity):
        self.sim = WaterControlSimulator(initial, physics, response, basis, cavity=cavity, graph_chunk=1)
        self.radius = physics.radius_m
        self.zeta = 6 * np.pi * physics.viscosity_Pa_s * self.radius
        self.basis = basis
        self.cavity = cavity
        self.x = np.array(initial, copy=True)

    def control(self, event):
        if event.get("flow_enabled", False):
            raise ValueError("Strict zero flow required")
        if event["duration_s"] <= 0:
            raise ValueError("Positive control duration required")
        self.sim.control(event)

    def fields(self, x):
        # F10N CPU query validates exactly the same strict interpolation domain.
        coord = (x - self.basis.origin) / self.basis.spacing
        if np.any(coord < 1) or np.any(coord >= np.array(self.basis.shape) - 2):
            raise ValueError("Outside stored cubic field interior; no clipping")
        self.sim.x.assign(x)
        self.sim.fields()
        p = self.sim.pressure.numpy()
        pressure = p[:, 0] + 1j * p[:, 1]
        primary = self.sim.force.numpy().copy()
        self.sim.pairs()
        total = self.sim.force.numpy().copy()
        bits = int(self.sim.error.numpy()[0])
        if bits:
            raise ValueError(f"Original acoustic/APECSS fail-stop bits={bits}")
        if self.sim.flow_enabled:
            raise ValueError("Nonzero-flow configuration")
        return dict(
            pressure_Pa=pressure,
            primary_force_N=primary,
            secondary_force_N=total - primary,
            total_force_N=total,
        )

    def wall_check(self, x):
        s, c = self.sim, self.cavity
        s.tmp1.assign(x)
        FastSimulator.launch(
            s,
            cavity_project,
            s.total,
            [
                s.tmp1,
                wp.vec3d(*c.lumen.center),
                wp.mat33d(*c.lumen.rotation.ravel()),
                c.R,
                c.a,
                c.L,
                c.H,
                c.ports[1],
                self.radius,
                c.ceiling,
                s.cavity_stats,
                s.error,
            ],
        )
        correction = s.tmp1.numpy() - x
        return dict(
            correction_m=float(np.linalg.norm(correction, axis=1).max()), error_bits=int(s.error.numpy()[0])
        )

    def step(self, dt):
        start = time.perf_counter()
        f = self.fields(self.x)
        t = time.perf_counter()
        op = resistance(self.x, self.radius)
        assembly = time.perf_counter() - t
        result = geometric_contact(
            self.x, f["total_force_N"] / self.zeta, self.radius, self.zeta, dt, KC, op["A"], atol_m=ATOL
        )
        walls = [self.wall_check(result[k]) for k in ["prediction_m", "positions_m"]]
        reasons = []
        if not result["converged"]:
            reasons.append("B2_MATHEMATICAL_RESIDUAL")
        if result["geometric_overlap_m"].max(initial=0) > 0.01 * self.radius:
            reasons.append("RESIDUAL_OVERLAP_BUDGET")
        if result["max_single_iteration_m"] > 0.1 * self.radius:
            reasons.append("CONTACT_ITERATION_BUDGET")
        if any(w["error_bits"] or w["correction_m"] != 0 for w in walls):
            reasons.append("INTERIOR_B2_WALL_CORRECTION_OR_SCOPE")
        if not reasons:
            self.x = result["positions_m"].copy()
        return dict(
            accepted=not reasons,
            reasons=reasons,
            result=result,
            fields=f,
            wall_correction_m=max(w["correction_m"] for w in walls),
            wall_error_bits=max(w["error_bits"] for w in walls),
            resistance_s=assembly,
            full_step_s=time.perf_counter() - start,
        )
