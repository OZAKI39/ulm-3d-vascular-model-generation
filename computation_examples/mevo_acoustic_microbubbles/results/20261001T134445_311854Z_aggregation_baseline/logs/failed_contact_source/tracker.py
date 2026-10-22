"""Slow-time RK integration with explicit geometric feasibility guards and audit logs."""
from __future__ import annotations

import json
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from mevo_mb.bubbles.contact import constrain_pair_velocities, project_contacts
from mevo_mb.bubbles.forces import NumpyPairBackend, primary_bjerknes
from mevo_mb.bubbles.metrics import aggregation_metrics
from mevo_mb.bubbles.population import BubblePopulation
from mevo_mb.bubbles.response import LinearResonantBubbleModel
from mevo_mb.coupling.geometry import Cylinder
from mevo_mb.flow.fields import FlowField
from mevo_mb.io.artifacts import write_json


class ParticleTracker:
    """Fixed nominal RK2/RK4 clock, optional displacement-limited internal substeps.

    Every accepted substep is saved. Displacement guards establish geometry safety,
    not integration-error convergence. Contact stages use geometric projections.
    """
    def __init__(self, config: dict, acoustic_field, flow: FlowField):
        self.config, self.field, self.flow = config, acoustic_field, flow
        self.p = config["particles"]
        self.cylinder = Cylinder(**config["vessel"])
        response = config["bubble_response"]
        if response["model"] != "linear_surrogate":
            raise NotImplementedError("Only the explicitly declared linear surrogate is enabled")
        self.response = LinearResonantBubbleModel(response["resonance_frequency_Hz"], response["damping_ratio"])
        if self.p["pair_backend"] != "numpy":
            raise ValueError("This version enables the tested NumPy pair backend for N < 500")
        self.pair = NumpyPairBackend(self.p["minimum_separation_m"], self.p["collision_regularization_m"], self.p["secondary_force_cutoff_m"])

    def evaluate(self, x, radii, t):
        pressure = self.field.interpolate_complex_pressure(x)
        gradient = self.field.interpolate_complex_gradient(x)
        volume = self.response.response(radii, self.field.frequency_Hz, pressure, self.config["medium"])
        primary = primary_bjerknes(volume, gradient) if self.p["primary_enabled"] else np.zeros_like(x)
        secondary = self.pair.forces(x, radii, volume, self.field.frequency_Hz, self.config["medium"]["density_kg_m3"]) if self.p["secondary_enabled"] else np.zeros_like(x)
        flow = self.flow.velocity(x, t)
        free = flow + (primary+secondary)/(6*np.pi*self.config["medium"]["viscosity_Pa_s"]*radii[:, None])
        velocity = self.cylinder.constrain_velocity(x, radii, free)
        velocity = constrain_pair_velocities(x, radii, velocity, self.p["contact_tolerance_m"], self.p["contact_max_iterations"])
        velocity = self.cylinder.constrain_velocity(x, radii, velocity)
        if not np.isfinite(velocity).all():
            raise RuntimeError("Nonfinite translation velocity")
        return {"primary_force_N": primary, "secondary_force_N": secondary, "flow_velocity_m_s": flow,
                "free_total_velocity_m_s": free, "total_velocity_m_s": velocity,
                "local_pressure_complex_Pa": pressure, "volume_response_complex_m3": volume}

    def _project(self, x, radii):
        return project_contacts(x, radii, self.cylinder, self.p["contact_tolerance_m"], self.p["contact_max_iterations"])

    def step(self, x, radii, t, dt, first):
        def velocity(position, at):
            feasible, _, _, _ = self._project(position, radii)
            return self.evaluate(feasible, radii, at)["total_velocity_m_s"]
        k1 = first["total_velocity_m_s"]
        k2 = velocity(x+dt*k1/2, t+dt/2)
        if self.p["integrator"] == "RK2":
            proposed = x+dt*k2
        elif self.p["integrator"] == "RK4":
            k3 = velocity(x+dt*k2/2, t+dt/2)
            k4 = velocity(x+dt*k3, t+dt)
            proposed = x+dt*(k1+2*k2+2*k3+k4)/6
        else:
            raise ValueError("Integrator must be RK2 or RK4")
        corrected, wall, contacts, iterations = self._project(proposed, radii)
        return corrected, wall, contacts, float(np.max(np.linalg.norm(corrected-proposed, axis=1))), iterations

    def run(self, population: BubblePopulation, run: Path) -> dict:
        start = time.perf_counter()
        x, radii = population.positions_m.copy(), population.radii_m
        if len(radii) < 2:
            raise ValueError("Population metrics require at least two bubbles")
        duration, nominal = self.p["duration_s"], self.p["dt_s"]
        time_s, accepted, subdivisions, contact_events, wall_events = 0., 0, 0, 0, 0
        max_projection, max_radial_fraction = 0., 0.
        minimum_dt = nominal
        metric_rows = []
        file = run/"particles/trajectories.h5"
        with h5py.File(file, "x") as h5, (run/"particles/events.jsonl").open("x") as events:
            h5.attrs["velocity_semantics"] = "Instantaneous force/drag velocity at saved state; not the preceding displacement divided by dt. All accepted substeps saved."
            h5.attrs["model"] = "Linear radial surrogate, quasisteady Stokes translation, far-field secondary forces, geometric contact only"
            h5.create_dataset("radii_m", data=radii)
            h5.create_dataset("particle_id", data=population.ids)
            datasets = {}

            def record(at, position, data, wall, actual_dt):
                nonlocal max_radial_fraction
                metric, labels, hist, edges = aggregation_metrics(position, radii, self.config["metrics"])
                metric.update(time_s=at, mean_displacement_m=float(np.linalg.norm(position-population.positions_m, axis=1).mean()),
                              mean_lateral_displacement_m=float(np.linalg.norm(position[:, 1:]-population.positions_m[:, 1:], axis=1).mean()))
                metric_rows.append(metric)
                max_radial_fraction = max(max_radial_fraction, float(np.max(np.abs(data["volume_response_complex_m3"])/(4*np.pi*radii**3))))
                items = {"time_s": np.array(at), "accepted_dt_s": np.array(actual_dt), "positions_m": position,
                         "cluster_labels": labels, "wall_hit": wall, "pair_distance_histogram": hist, **data}
                for name, value in items.items():
                    value = np.asarray(value)
                    if name not in datasets:
                        datasets[name] = h5.create_dataset(name, shape=(0,)+value.shape, maxshape=(None,)+value.shape,
                                                          chunks=(1,)+value.shape, dtype=value.dtype, compression="gzip", compression_opts=1)
                    ds = datasets[name]
                    ds.resize(len(metric_rows), axis=0)
                    ds[-1] = value
                if "pair_distance_histogram_edges_m" not in h5:
                    h5.create_dataset("pair_distance_histogram_edges_m", data=edges)

            state = self.evaluate(x, radii, time_s)
            record(0., x, state, np.zeros(len(radii), bool), 0.)
            macro_count = int(np.ceil(duration/nominal))
            for macro in range(1, macro_count+1):
                target = min(macro*nominal, duration)
                local_count = 0
                while time_s < target-1e-15:
                    speed = np.linalg.norm(state["total_velocity_m_s"], axis=1)
                    geometric_dt = np.min(self.p["max_displacement_per_radius"]*radii/np.maximum(speed, 1e-30))
                    dt = min(target-time_s, geometric_dt)
                    if dt <= 1e-15:
                        raise RuntimeError("Translation step collapsed; preserve failure instead of fabricating motion")
                    x, wall, contacts, correction, _ = self.step(x, radii, time_s, dt, state)
                    time_s = min(time_s+dt, target)
                    state = self.evaluate(x, radii, time_s)
                    accepted += 1
                    local_count += 1
                    minimum_dt = min(minimum_dt, dt)
                    max_projection = max(max_projection, correction)
                    wall_events += int(wall.sum())
                    contact_events += len(contacts)
                    if contacts or wall.any():
                        events.write(json.dumps({"time_s": time_s, "contact_pairs": contacts, "wall_particle_ids": population.ids[wall].tolist(), "maximum_projection_m": correction})+"\n")
                    record(time_s, x, state, wall, dt)
                    if local_count > self.p["max_substeps"]:
                        raise RuntimeError("Substep budget exceeded")
                subdivisions += max(0, local_count-1)
                if macro % 40 == 0:
                    h5.flush()
                    write_json(run/"PROGRESS.json", {"time_s": time_s, "duration_s": duration, "accepted_steps": accepted, "elapsed_s": time.perf_counter()-start})
            h5.attrs["complete"] = True
        pd.DataFrame(metric_rows).to_csv(run/"metrics.csv", index=False)
        summary = {"status": "COMPLETE", "particle_count": len(radii), "duration_s": time_s,
                   "nominal_dt_s": nominal, "minimum_accepted_dt_s": minimum_dt, "accepted_substeps": accepted,
                   "additional_geometric_substeps": subdivisions, "contact_projection_pair_events": contact_events,
                   "wall_hit_events": wall_events, "maximum_projection_m": max_projection,
                   "maximum_linear_radial_amplitude_over_R0": max_radial_fraction,
                   "radial_surrogate_validated_for_SonoVue": False, "initial": metric_rows[0], "final": metric_rows[-1],
                   "elapsed_s": time.perf_counter()-start, "trajectory_file": str(file),
                   "aggregation_definition": "spatial clustering, no coalescence", "all_accepted_substeps_recorded": True}
        write_json(run/"particles/SUMMARY.json", summary)
        return summary
