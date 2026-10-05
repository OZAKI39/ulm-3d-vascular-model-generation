"""Slow-time RK integration with explicit geometric feasibility guards and audit logs."""

from __future__ import annotations

import json
from collections import deque
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from mevo_mb.bubbles.translation_backends import (
    make_translation_backend, LUBRICATION_PHYSICS, INTERFACE_SCOPE, NOT_INCLUDED, LEGACY_LABEL,
)
from mevo_mb.bubbles.contact import constrain_pair_velocities, project_contacts
from mevo_mb.bubbles.forces import NumpyPairBackend, ClassicalFarFieldBjerknes, primary_bjerknes
from mevo_mb.bubbles.metrics import aggregation_metrics
from mevo_mb.bubbles.population import BubblePopulation
from mevo_mb.bubbles.response import LinearResonantBubbleModel
from mevo_mb.bubbles.response_table import ResponseTableBubbleModel
from mevo_mb.bubbles.adaptive_table import AdaptiveResponseTableBubbleModel
from mevo_mb.bubbles.validity import PhysicsValidityMonitor, PairInteractionValidity, PhysicsValidityError
from mevo_mb.coupling.geometry import geometry_from_config
from mevo_mb.flow.fields import FlowField
from mevo_mb.io.artifacts import write_json


class ParticleTracker:
    """RK2/RK4 with an explicit translation backend and geometric guards.

    Legacy uses the original displacement guard and diagonal contact metric.
    Lubrication uses closing-gap guards and A-metric normal contact constraints.
    These guards do not establish integration-error convergence.
    """

    def __init__(self, config: dict, acoustic_field, flow: FlowField):
        self.config, self.field, self.flow = config, acoustic_field, flow
        self.p = config["particles"]
        if self.p.get("history_force") or config.get("history_force"):
            raise ValueError("History force is validation-only and cannot be enabled in ParticleTracker")
        self.timings = {k: 0.0 for k in ("response_queries_s", "pair_force_s", "contact_search_s", "particle_integration_s", "flow_interpolation_s", "vtk_interpolation_s")}
        self.timings["contact_search_calls"] = 0
        self.contact_backend = self.p.get("contact_search_backend", "brute")
        self.cylinder = geometry_from_config(config)
        self.translation = make_translation_backend(config)
        self.lubrication_enabled = self.translation.name == "lubrication_coupled"
        self.translation_diagnostics = deque(maxlen=256)
        self.translation_log = None
        self.translation_summary = dict(linear_solves=0, contact_solves=0,
            max_linear_relative_residual=0., max_contact_kkt_residual=0.,
            total_linear_solve_time_s=0., total_contact_solve_time_s=0.)
        self.projection_diagnostics = dict(calls=0, corrected_calls=0, max_correction_m=0.,
                                           squared_correction_m2=0., particle_queries=0,
                                           multi_radius_corrections=0)

        response = config["bubble_response"]
        # Old explicitly-surrogate configurations retain their baseline behavior,
        # always labelled DEMO ONLY. New response-table runs default to strict.
        self.policy = self.p.get(
            "near_field_policy", "legacy_demo" if response["model"] == "linear_surrogate" else "strict"
        )
        self.monitor = PhysicsValidityMonitor(
            response["model"],
            policy=self.policy,
            regime=config.get("physics_validity", {}).get("regime"),
            pair=PairInteractionValidity(**config.get("physics_validity", {}).get("pair", {})),
            min_wall_distance_over_radius=config.get("physics_validity", {}).get(
                "min_wall_distance_over_radius", 5.0
            ),
            regime_map_path=config.get('physics_validity', {}).get('regime_map_path'),
            max_invalid_pair_fraction=config.get('physics_validity', {}).get('max_invalid_pair_fraction', 0.0),
            max_transition_pair_fraction=config.get('physics_validity', {}).get('max_transition_pair_fraction', 0.05),
        )
        if self.policy == "pair_reference":
            # Full coupled reference has its own radial/translation integrator.
            # It cannot be inserted as an N-body force sum or combined with flow.
            self.response = self.pair = None
            return
        if response["model"] in ("response_table", "adaptive_response_table"):
            table_class = AdaptiveResponseTableBubbleModel if response['model']=='adaptive_response_table' else ResponseTableBubbleModel
            self.response = table_class(
                response["path"], backend=self.p.get("response_backend", 'numpy' if self.p.get('pair_backend')=='auto' else self.p.get("pair_backend", "numpy"))
            )
            self.monitor.response_context = {
                "table": str(self.response.path),
                "metadata": self.response.metadata,
                "axis_bounds_SI": {
                    k: [float(a[0]), float(a[-1])]
                    for k, a in zip(("radius_m", "frequency_Hz", "pressure_Pa"), self.response.axes)
                },
                "extrapolation_allowed": False,
                "validated_cells_required": True,
                "shell_parameters": "Fixed table parameter manifold; other shell parameters require another table",
                "experimental_shell_calibration_validated": False,
            }
        elif response["model"] == "linear_surrogate" and self.policy == "legacy_demo":
            self.response = LinearResonantBubbleModel(
                response["resonance_frequency_Hz"], response["damping_ratio"]
            )
        else:
            raise ValueError(
                "Strict trajectories require a validated response table; linear surrogate is demo only"
            )
        if self.policy == "legacy_demo":
            self.pair = NumpyPairBackend(
                self.p["minimum_separation_m"],
                self.p["collision_regularization_m"],
                self.p["secondary_force_cutoff_m"],
            )
        else:
            if self.p.get('pair_backend') == 'auto':
                from mevo_mb.bubbles.auto_backend import AutoPairBackend
                self.pair = AutoPairBackend(benchmark_path=self.p.get('pair_benchmark_path', 'results/stage2/performance/STATE.json'),
                    crossover=self.p.get('pair_crossover'),min_distance_over_sum=self.monitor.pair.min_distance_over_sum, enforce_far_field=self.policy != "simulator")
            else:
                self.pair = ClassicalFarFieldBjerknes(
                    self.p.get("pair_backend", "numpy"),
                    min_distance_over_sum=self.monitor.pair.min_distance_over_sum,
                    enforce_far_field=self.policy != "simulator",
                )

    def evaluate(self, x, radii, t):
        self._last_query = dict(
            positions_m=np.array(x, copy=True), radii_m=np.array(radii, copy=True), time_s=float(t)
        )
        pressure = self.field.interpolate_complex_pressure(x)
        gradient = self.field.interpolate_complex_gradient(x)

        def cpu(value):
            return value.get() if hasattr(value, "get") else np.asarray(value)

        pair_volume = None
        resolved = True
        query_start = time.perf_counter()
        try:
            if not self.p["primary_enabled"] and not self.p["secondary_enabled"]:
                # Exactly zero acoustic force: no radial table extrapolation at P=0.
                volume = np.zeros_like(radii, dtype=complex)
                radial_fraction = np.zeros_like(radii)
            elif isinstance(self.response, ResponseTableBubbleModel):
                self.response.validate_medium(self.config["medium"])
                q = self.response.harmonic_response(radii, self.field.frequency_Hz, pressure)
                pair_volume = q["harmonics_m3"]
                volume = cpu(pair_volume[..., 0])
                radial_fraction = np.maximum(
                    cpu(q["corner_envelope_max_R_over_R0"]) - 1, 1 - cpu(q["corner_envelope_min_R_over_R0"])
                )
                resolved = bool(cpu(q["harmonic_resolved"]).all())
            else:
                volume = cpu(
                    self.response.response(radii, self.field.frequency_Hz, pressure, self.config["medium"])
                )
                radial_fraction = np.abs(volume) / (4 * np.pi * radii**3)
        except ValueError as error:
            self.monitor.observe(
                x,
                radii,
                t,
                radial_fraction=np.zeros_like(radii),
                frequency=self.field.frequency_Hz,
                sound_speed=self.config["medium"]["sound_speed_m_s"],
                response_valid=False,
                secondary_enabled=self.p["secondary_enabled"],
            )
            raise PhysicsValidityError(str(error)) from error
        self.timings["response_queries_s"] += time.perf_counter()-query_start
        wall_distance = self.cylinder.wall_distance(x)
        self.monitor.observe(
            x,
            radii,
            t,
            radial_fraction=radial_fraction,
            frequency=self.field.frequency_Hz,
            sound_speed=self.config["medium"]["sound_speed_m_s"],
            wall_distance_m=wall_distance,
            harmonic_resolved=resolved,
            secondary_enabled=self.p["secondary_enabled"],
            pressure_amplitude=np.abs(pressure),
            pressure_phase_rad=np.angle(pressure),
            shell_model=self.config['bubble_response'].get('shell_model', self.config['bubble_response']['model']),
        )
        primary = primary_bjerknes(volume, gradient) if self.p["primary_enabled"] else np.zeros_like(x)
        pair_start = time.perf_counter()
        secondary = (
            cpu(
                self.pair.forces(
                    x,
                    radii,
                    pair_volume if pair_volume is not None and self.policy != "legacy_demo" else volume,
                    self.field.frequency_Hz,
                    self.config["medium"]["density_kg_m3"],
                )
            )
            if self.p["secondary_enabled"]
            else np.zeros_like(x)
        )
        self.timings["pair_force_s"] += time.perf_counter()-pair_start
        flow_start = time.perf_counter()
        flow = self.flow.velocity(x, t)
        flow_elapsed = time.perf_counter()-flow_start
        self.timings["flow_interpolation_s"] += flow_elapsed
        from mevo_mb.flow.fields import VTKFlowField
        if isinstance(self.flow, VTKFlowField):
            self.timings["vtk_interpolation_s"] += flow_elapsed
        wall = self.cylinder.wall_gap_and_normal(x, radii) if self.lubrication_enabled else None
        translation = self.translation.solve(x, radii, primary + secondary, flow, wall)
        self.last_translation = translation
        free = translation.velocity
        if self.lubrication_enabled:
            velocity, contact_diag = self.translation.contact(
                translation, x, radii, free, self.p["contact_tolerance_m"])
            self._record_translation(dict(translation.diagnostics, **contact_diag, time_s=float(t)))
        else:
            velocity = self.cylinder.constrain_velocity(x, radii, free)
            velocity = constrain_pair_velocities(
                x, radii, velocity, self.p["contact_tolerance_m"], self.p["contact_max_iterations"],
                self.contact_backend, self.timings)
            velocity = self.cylinder.constrain_velocity(x, radii, velocity)
        if not np.isfinite(velocity).all():
            raise RuntimeError("Nonfinite translation velocity")
        return {
            "primary_force_N": primary,
            "secondary_force_N": secondary,
            "flow_velocity_m_s": flow,
            "free_total_velocity_m_s": free,
            "total_velocity_m_s": velocity,
            "local_pressure_complex_Pa": pressure,
            "volume_response_complex_m3": volume,
        }

    def _project(self, x, radii):
        result = project_contacts(
            x, radii, self.cylinder, self.p["contact_tolerance_m"], self.p["contact_max_iterations"], self.contact_backend, self.timings
        )
        correction = np.linalg.norm(result[0] - x, axis=1)
        diag = self.projection_diagnostics
        diag["calls"] += 1
        diag["corrected_calls"] += int(np.any(correction > 0))
        diag["max_correction_m"] = max(diag["max_correction_m"], float(np.max(correction)))
        diag["squared_correction_m2"] += float(correction @ correction)
        diag["particle_queries"] += len(x)
        diag["multi_radius_corrections"] += int(np.count_nonzero(correction > 2*radii))
        return result

    def _record_translation(self, diagnostic):
        self.translation_diagnostics.append(diagnostic)
        summary = self.translation_summary
        if "solver" in diagnostic:
            summary["linear_solves"] += 1
        summary["contact_solves"] += int(diagnostic.get("contact_solve_time_s", 0) > 0)
        summary["max_linear_relative_residual"] = max(summary["max_linear_relative_residual"],
            diagnostic.get("relative_residual", 0.))
        summary["max_contact_kkt_residual"] = max(summary["max_contact_kkt_residual"],
            diagnostic.get("contact_kkt_residual", 0.))
        summary["total_linear_solve_time_s"] += diagnostic.get("solve_time_s", 0.)
        summary["total_contact_solve_time_s"] += diagnostic.get("contact_solve_time_s", 0.)
        if self.translation_log is not None:
            self.translation_log.write(json.dumps(diagnostic)+"\n")

    def select_timestep(self, x, radii, state, remaining):
        if self.lubrication_enabled:
            return self.translation.timestep(x, radii, state["total_velocity_m_s"],
                                            self.last_translation.wall, remaining)
        speed = np.linalg.norm(state["total_velocity_m_s"], axis=1)
        return min(remaining, float(np.min(self.p["max_displacement_per_radius"] * radii /
                                          np.maximum(speed, 1e-30))))

    def step(self, x, radii, t, dt, first):
        integration_start = time.perf_counter()
        initial_translation = self.last_translation if self.lubrication_enabled else None

        def feasible_step(candidate, horizon):
            if not self.lubrication_enabled:
                return candidate
            average_velocity = (candidate - x) / horizon
            safe, diag = self.translation.contact(initial_translation, x, radii, average_velocity,
                                                   self.p["contact_tolerance_m"], horizon=horizon)
            self._record_translation(dict(diag, time_s=float(t), stage="finite_horizon_contact"))
            return x + horizon * safe

        def velocity(position, at):
            position = feasible_step(position, at-t)
            feasible, _, _, _ = self._project(position, radii)
            return self.evaluate(feasible, radii, at)["total_velocity_m_s"]

        k1 = first["total_velocity_m_s"]
        k2 = velocity(x + dt * k1 / 2, t + dt / 2)
        if self.p["integrator"] == "RK2":
            proposed = x + dt * k2
        elif self.p["integrator"] == "RK4":
            k3 = velocity(x + dt * k2 / 2, t + dt / 2)
            k4 = velocity(x + dt * k3, t + dt)
            proposed = x + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6
        else:
            raise ValueError("Integrator must be RK2 or RK4")
        proposed = feasible_step(proposed, dt)
        corrected, wall, contacts, iterations = self._project(proposed, radii)
        self.timings["particle_integration_s"] += time.perf_counter()-integration_start
        return (
            corrected,
            wall,
            contacts,
            float(np.max(np.linalg.norm(corrected - proposed, axis=1))),
            iterations,
        )

    def run(self, population: BubblePopulation, run: Path) -> dict:
        if self.policy == "pair_reference":
            from mevo_mb.validation.pair_policy import run_pair_reference

            return run_pair_reference(self.config, population, run, self.field, self.flow)
        self.banner = {
            "production_physics": ["linear external acoustics", self.config["bubble_response"]["model"], "cycle-averaged primary Bjerknes", "classical secondary Bjerknes" if self.policy != "legacy_demo" else "deprecated regularized secondary DEMO", "quasi-steady Stokes translation", "one-way flow coupling", "geometric no-penetration"],
            "not_included": ["hydrodynamic history force", "high-fidelity near-contact dynamics", "coalescence", "adhesion", "fusion", "nonspherical bubble deformation", "acoustic streaming", "bubble feedback on flow", "bubble feedback on acoustic field"],
            "experimental_calibration": False,
            "primary_enabled": self.p["primary_enabled"],
            "secondary_enabled": self.p["secondary_enabled"],
            "contact_search_backend": self.contact_backend,
            "contact_semantics": "Geometric nonpenetration only; contact is not coalescence, adhesion, fusion or validated near-field dynamics",
        }
        self.banner["translation_backend"] = self.translation.name
        if self.lubrication_enabled:
            self.banner["production_physics"].remove("quasi-steady Stokes translation")
            self.banner["production_physics"].extend([LUBRICATION_PHYSICS, INTERFACE_SCOPE])
            self.banner["not_included"].extend(NOT_INCLUDED)
            self.banner["contact_semantics"] = "A-metric frictionless normal nonpenetration; geometric residual projection"
        else:
            self.banner["baseline_identity"] = LEGACY_LABEL
        write_json(run / "PHYSICS_BANNER.json", self.banner)
        print("PRODUCTION_PHYSICS " + json.dumps(self.banner), flush=True)
        if self.lubrication_enabled:
            self.translation_log = (run / "TRANSLATION_SOLVES.jsonl").open("x")
        try:
            return self._run(population, run)
        except PhysicsValidityError as error:
            # Preserve both the last accepted trajectory and rejected RK-stage
            # coordinates; do not mark an incomplete trajectory complete.
            report = {
                "status": "OUTSIDE_VALIDITY",
                "error": str(error),
                "physics_validity": self.monitor.summary(),
            }
            write_json(run / "particles/SUMMARY.json", report)
            write_json(run / "physics_validity.json", self.monitor.summary())
            if hasattr(self, "_last_query"):
                np.savez(run / "particles/rejected_state.npz", **self._last_query)
            raise
        finally:
            if self.translation_log is not None:
                self.translation_log.close()
                self.translation_log = None
                write_json(run / "TRANSLATION_DIAGNOSTICS.json", dict(
                    **self.translation_summary, projection=self.projection_diagnostics,
                    failure=getattr(self.translation, "last_failure", None)))

    def _run(self, population: BubblePopulation, run: Path) -> dict:
        start = time.perf_counter()
        last_progress = start
        x, radii = population.positions_m.copy(), population.radii_m
        if len(radii) < 2:
            raise ValueError("Population metrics require at least two bubbles")
        duration, nominal = self.p["duration_s"], self.p["dt_s"]
        time_s, accepted, subdivisions, contact_events, wall_events = 0.0, 0, 0, 0, 0
        max_projection, max_radial_fraction = 0.0, 0.0
        minimum_dt = nominal
        metric_rows = []
        file = run / "particles/trajectories.h5"
        with h5py.File(file, "x") as h5, (run / "particles/events.jsonl").open("x") as events:
            h5.attrs["translation_backend"] = self.translation.name
            h5.attrs["velocity_semantics"] = (
                "Instantaneous force/drag velocity at saved state; not the preceding displacement divided by dt. All accepted substeps saved."
            )
            h5.attrs["model"] = (
                f"{self.config['bubble_response']['model']}; {self.policy}; " +
                (LUBRICATION_PHYSICS + "; A-metric normal contact" if self.lubrication_enabled else
                 "quasisteady Stokes; geometric contact only")
            )
            h5.attrs["physics_status"] = (
                "NOT PHYSICALLY VALIDATED" if self.policy == "legacy_demo" else "partially validated"
            )
            h5.create_dataset("radii_m", data=radii)
            h5.create_dataset("particle_id", data=population.ids)
            datasets = {}

            def record(at, position, data, wall, actual_dt):
                nonlocal max_radial_fraction
                self.monitor.record_accepted_interval(actual_dt)
                metric, labels, hist, edges = aggregation_metrics(position, radii, self.config["metrics"])
                metric.update(
                    time_s=at,
                    mean_displacement_m=float(
                        np.linalg.norm(position - population.positions_m, axis=1).mean()
                    ),
                    mean_lateral_displacement_m=float(
                        np.linalg.norm(position[:, 1:] - population.positions_m[:, 1:], axis=1).mean()
                    ),
                )
                metric.update({'pair_regime_'+k:v for k,v in self.monitor.last['pair_regime'].items() if k.endswith('_fraction')})
                metric_rows.append(metric)
                max_radial_fraction = max(
                    max_radial_fraction,
                    float(np.max(np.abs(data["volume_response_complex_m3"]) / (4 * np.pi * radii**3))),
                )
                items = {
                    "time_s": np.array(at),
                    "accepted_dt_s": np.array(actual_dt),
                    "positions_m": position,
                    "cluster_labels": labels,
                    "wall_hit": wall,
                    "pair_distance_histogram": hist,
                    **data,
                }
                for name, value in items.items():
                    value = np.asarray(value)
                    if name not in datasets:
                        datasets[name] = h5.create_dataset(
                            name,
                            shape=(0,) + value.shape,
                            maxshape=(None,) + value.shape,
                            chunks=(1,) + value.shape,
                            dtype=value.dtype,
                            compression="gzip",
                            compression_opts=1,
                        )
                    ds = datasets[name]
                    ds.resize(len(metric_rows), axis=0)
                    ds[-1] = value
                if "pair_distance_histogram_edges_m" not in h5:
                    h5.create_dataset("pair_distance_histogram_edges_m", data=edges)

            state = self.evaluate(x, radii, time_s)
            record(0.0, x, state, np.zeros(len(radii), bool), 0.0)
            macro_count = int(np.ceil(duration / nominal))
            for macro in range(1, macro_count + 1):
                target = min(macro * nominal, duration)
                local_count = 0
                while time_s < target - 1e-15:
                    dt = self.select_timestep(x, radii, state, target - time_s)
                    if dt <= 1e-15:
                        raise RuntimeError(
                            "Translation step collapsed; preserve failure instead of fabricating motion"
                        )
                    x, wall, contacts, correction, _ = self.step(x, radii, time_s, dt, state)
                    time_s = min(time_s + dt, target)
                    state = self.evaluate(x, radii, time_s)
                    accepted += 1
                    local_count += 1
                    minimum_dt = min(minimum_dt, dt)
                    max_projection = max(max_projection, correction)
                    wall_events += int(wall.sum())
                    contact_events += len(contacts)
                    if contacts or wall.any():
                        events.write(
                            json.dumps(
                                {
                                    "time_s": time_s,
                                    "contact_pairs": contacts,
                                    "contact_event": bool(contacts),
                                    "wall_particle_ids": population.ids[wall].tolist(),
                                    "maximum_projection_m": correction,
                                }
                            )
                            + "\n"
                        )
                    record(time_s, x, state, wall, dt)
                    if time.perf_counter() - last_progress >= 5:
                        h5.flush()
                        write_json(
                            run / "PROGRESS.json",
                            {
                                "time_s": time_s,
                                "duration_s": duration,
                                "accepted_steps": accepted,
                                "elapsed_s": time.perf_counter() - start,
                            },
                        )
                        last_progress = time.perf_counter()
                    if local_count > self.p["max_substeps"]:
                        raise RuntimeError("Substep budget exceeded")
                subdivisions += max(0, local_count - 1)
                if macro % 40 == 0:
                    h5.flush()
                    write_json(
                        run / "PROGRESS.json",
                        {
                            "time_s": time_s,
                            "duration_s": duration,
                            "accepted_steps": accepted,
                            "elapsed_s": time.perf_counter() - start,
                        },
                    )
            h5.attrs["complete"] = True
        pd.DataFrame(metric_rows).to_csv(run / "metrics.csv", index=False)
        summary = {
            "status": "COMPLETE",
            "particle_count": len(radii),
            "duration_s": time_s,
            "nominal_dt_s": nominal,
            "minimum_accepted_dt_s": minimum_dt,
            "accepted_substeps": accepted,
            "additional_geometric_substeps": subdivisions,
            "contact_projection_pair_events": contact_events,
            "wall_hit_events": wall_events,
            "maximum_projection_m": max_projection,
            "maximum_linear_radial_amplitude_over_R0": max_radial_fraction,
            "radial_surrogate_validated_for_SonoVue": False,
            "initial": metric_rows[0],
            "final": metric_rows[-1],
            "elapsed_s": time.perf_counter() - start,
            "trajectory_file": str(file),
            "aggregation_definition": "spatial clustering, no coalescence",
            "all_accepted_substeps_recorded": True,
        }
        summary["physics_banner"] = self.banner
        summary["timings"] = dict(self.timings, note="Wall seconds; integration includes its RK-stage response/force/contact/flow work; component times must not be summed with integration")
        summary["physics_validity"] = self.monitor.summary()
        summary["NEAR_CONTACT_PHYSICS_APPROXIMATED"] = summary["physics_validity"]["NEAR_CONTACT_PHYSICS_APPROXIMATED"]
        summary['pair_regime'] = summary['physics_validity']['pair_regime']
        summary['validity_status'] = summary['physics_validity']['validity_status']
        if hasattr(self.pair, 'decisions'):
            summary['pair_backend_decisions'] = self.pair.decisions
        if isinstance(self.response, ResponseTableBubbleModel):
            summary["maximum_fundamental_volume_equivalent_radial_fraction"] = summary.pop(
                "maximum_linear_radial_amplitude_over_R0"
            )
            summary["maximum_radial_excursion_over_R0"] = max(
                item["max_radial_fraction"] for item in self.monitor.history
            )
            summary["SonoVue_experimental_calibration_validated"] = False
        write_json(run / "physics_validity.json", self.monitor.summary())
        write_json(run / "particles/SUMMARY.json", summary)
        return summary
