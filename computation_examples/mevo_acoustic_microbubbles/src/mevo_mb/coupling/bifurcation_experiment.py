"""Open-boundary experiment recorder using the unchanged Stage 3 force/RK methods."""
import json
import time
from pathlib import Path

import h5py
import numpy as np

from mevo_mb.bubbles.metrics import aggregation_metrics
from mevo_mb.bubbles.tracker import ParticleTracker
from mevo_mb.coupling.bifurcation import BranchClassifier
from mevo_mb.io.artifacts import write_json


class BifurcationExperimentTracker(ParticleTracker):
    """Same constructor, force law and stepper; adds first-exit recording/removal.

    Exited particles are no longer force/flow query targets. Stored positions freeze
    at the measured crossing and carry active=False, never pretending to advect.
    All crossing/contact events are examined at every accepted substep. Visualization
    states are saved at a separate configured cadence to bound output size.
    """
    def run(self, population, run):
        try:
            return super().run(population, run)
        except Exception as error:
            run = Path(run)
            if hasattr(self, "_last_accepted"):
                np.savez_compressed(run / "particles/last_accepted_state.npz", **self._last_accepted)
            if hasattr(self, "classifier"):
                write_json(run / "BRANCH_METRICS_INCOMPLETE.json", self.classifier.summary())
            write_json(run / "EXECUTION_FAILURE.json", dict(status="INCOMPLETE",error=repr(error),
                       completed=False, saved_visualization_states_are_partial=True))
            raise

    def _run(self, population, run):
        start = time.perf_counter()
        output_seconds = 0.0
        run = Path(run)
        ids = population.ids.copy()
        if not np.array_equal(ids, np.arange(len(ids))):
            raise ValueError("Experiment recorder requires stable contiguous population IDs")
        n = len(ids)
        active = ids.copy()
        x = population.positions_m.copy()
        radii = population.radii_m.copy()
        frozen = x.copy()
        classifier = BranchClassifier(self.cylinder, ids, self.config.get("classification_tolerance_m"))
        self.classifier = classifier
        nominal = float(self.p["dt_s"])
        duration = float(self.p["duration_s"])
        cadence = float(self.config.get("output", {}).get("save_interval_s", .01))
        t, next_save, accepted, minimum_dt = 0., cadence, 0, nominal
        wall_hits = np.zeros(n, np.int64)
        contacts_count = 0
        maximum_projection = 0.0
        metrics = []
        saved_count, last_saved_time = 0, None
        trajectory_path = run / "particles/trajectories.h5"
        with h5py.File(trajectory_path,"x") as h:
            h.attrs["complete"] = False
            h.attrs["model"] = "Stage3 forces/RK; open classification planes; geometric wall/contact"
            h.attrs["all_accepted_substeps_recorded"] = False
            h.attrs["event_semantics"] = "Every accepted substep tested; inactive coordinates freeze at first crossing"
            h.attrs["inactive_force_velocity_semantics"] = "NaN: inactive particles are not evaluated"
            h.create_dataset("particle_id", data=ids)
            h.create_dataset("radii_m", data=radii)
        roi_volume = self.cylinder.concentration_roi_volume()
        last_progress = time.perf_counter()

        def record(state, actual_dt):
            nonlocal output_seconds, saved_count, last_saved_time
            begin = time.perf_counter()
            frozen[active] = x
            full = {key: np.full((n,)+np.asarray(value).shape[1:], np.nan,
                                dtype=complex if np.iscomplexobj(value) else float)
                    for key, value in state.items()}
            for key, value in state.items():
                full[key][active] = value
            mask = np.zeros(n, bool)
            mask[active] = True
            snapshot = dict(time_s=t, positions_m=frozen.copy(), active=mask,
                            accepted_dt_s=actual_dt, **full)
            with h5py.File(trajectory_path,"a") as h:
                for key,value in snapshot.items():
                    value = np.asarray(value)
                    if key not in h:
                        h.create_dataset(key,shape=(0,)+value.shape,maxshape=(None,)+value.shape,
                                         dtype=value.dtype,chunks=(1,)+value.shape,compression="gzip",compression_opts=1)
                    ds = h[key]
                    ds.resize(saved_count+1,axis=0)
                    ds[saved_count] = value
            saved_count += 1
            last_saved_time = t
            if len(active) > 1:
                row, _, _, _ = aggregation_metrics(x, radii[active], self.config["metrics"])
            else:
                row = dict(mean_nearest_neighbor_m=None, median_nearest_neighbor_m=None,
                           cluster_count=len(active), largest_cluster_fraction=1. if len(active) else 0.,
                           near_contact_fraction=0., minimum_surface_gap_m=None)
            row.update(time_s=t, active_particles=len(active),
                       pair_regime=self.monitor.last.get("pair_regime", {}))
            for name, direction in self.cylinder.directions.items():
                along = x @ direction
                transverse = x - along[:, None]*direction
                lower = self.cylinder.config["concentration_roi_start_m"]
                upper = self.cylinder.config["concentration_roi_end_m"]
                selected = ((along >= lower) & (along <= upper)
                            & (np.linalg.norm(transverse, axis=1) <= self.cylinder.daughter_radius)
                            & (x[:,1] * direction[1] > 0))
                row[name+"_roi_count"] = int(selected.sum())
                row[name+"_normalized_occupancy"] = float(selected.sum()/n)
                row[name+"_particles_per_nominal_roi_m3"] = float(selected.sum()/roi_volume)
                row[name+"_normalized_number_density_m_inverse3"] = float(selected.sum()/n/roi_volume)
                branch_mask = (along >= lower) & (x[:,1]*direction[1] > 0)
                if branch_mask.sum() > 1:
                    branch_metrics, _, _, _ = aggregation_metrics(x[branch_mask], radii[active][branch_mask], self.config["metrics"])
                    row[name+"_spatial_aggregation"] = branch_metrics
            a_count, b_count = row["BRANCH_A_roi_count"], row["BRANCH_B_roi_count"]
            row["A_to_B_concentration_ratio"] = a_count/b_count if b_count else None
            row["B_to_A_concentration_ratio"] = b_count/a_count if a_count else None
            metrics.append(row)
            output_seconds += time.perf_counter()-begin

        state = self.evaluate(x, radii, 0.)
        self._last_accepted = dict(time_s=t,positions_m=x.copy(),particle_id=active.copy(),radii_m=radii.copy())
        record(state, 0.)
        with (run / "particles/events.jsonl").open("x") as events:
            for macro in range(1, int(np.ceil(duration/nominal))+1):
                target = min(macro*nominal, duration)
                local_steps = 0
                while len(active) and t < target-1e-15:
                    speed = np.linalg.norm(state["total_velocity_m_s"], axis=1)
                    geometric_dt = np.min(self.p["max_displacement_per_radius"]*radii[active]
                                          /np.maximum(speed, 1e-30))
                    dt = min(target-t, float(geometric_dt))
                    if dt <= 1e-15:
                        raise RuntimeError("Translation step collapsed")
                    old = x.copy()
                    old_state = state
                    x, wall, contacts, correction, _ = self.step(x, radii[active], t, dt, state)
                    end_time = min(t+dt, target)
                    state = self.evaluate(x, radii[active], end_time)
                    crossings = classifier.observe(active, old, x, old_state["total_velocity_m_s"],
                                                   state["total_velocity_m_s"], t, end_time)
                    self.monitor.record_accepted_interval(dt)
                    for event in crossings:
                        events.write(json.dumps(event)+"\n")
                    if wall.any() or contacts:
                        events.write(json.dumps(dict(time_s=end_time, wall_particle_ids=active[wall].tolist(),
                                                     contact_particle_ids=[[int(active[i]), int(active[j])]
                                                                           for i,j in contacts],
                                                     maximum_projection_m=correction))+"\n")
                    wall_hits[active[wall]] += 1
                    contacts_count += len(contacts)
                    maximum_projection = max(maximum_projection, correction)
                    frozen[active] = x
                    leaving = {event["particle_id"] for event in crossings}
                    for event in crossings:
                        frozen[event["particle_id"]] = event["position_m"]
                    keep = np.array([int(pid) not in leaving for pid in active], bool)
                    active = active[keep]
                    x = x[keep]
                    state = {key:value[keep] for key,value in state.items()}
                    t = end_time
                    self._last_accepted = dict(time_s=t,positions_m=x.copy(),particle_id=active.copy(),radii_m=radii[active].copy())
                    if leaving and len(active):
                        state = self.evaluate(x, radii[active], t)
                    accepted += 1
                    local_steps += 1
                    minimum_dt = min(minimum_dt, dt)
                    if t >= next_save-1e-15 or not len(active):
                        record(state, dt)
                        next_save = t+cadence
                    if local_steps > self.p["max_substeps"]:
                        raise RuntimeError("Substep budget exceeded")
                    if time.perf_counter()-last_progress > 10:
                        events.flush()
                        write_json(run / "PROGRESS.json", dict(time_s=t, duration_s=duration,
                                   active_particles=len(active), accepted_steps=accepted,
                                   elapsed_s=time.perf_counter()-start))
                        last_progress = time.perf_counter()
                if not len(active):
                    break
        if last_saved_time != t:
            record(state, minimum_dt)
        output_start = time.perf_counter()
        with h5py.File(trajectory_path, "a") as h:
            h.attrs["complete"] = True
        write_json(run / "metrics.json", dict(rows=metrics,nominal_roi_volume_m3=roi_volume,
                   concentration_scope="Finite initial ensemble, no injection concentration; normalized occupancy/density, not clinical particles/mL. Nominal cylindrical ROI intersected with branch y half-space."))
        branch = classifier.summary()
        assert sum(branch[k] for k in ("N_A", "N_B", "N_unassigned", "N_wall_loss")) == n
        write_json(run / "BRANCH_METRICS.json", branch)
        validity = self.monitor.summary()
        write_json(run / "physics_validity.json", validity)
        output_seconds += time.perf_counter()-output_start
        self.timings["output_s"] = output_seconds
        summary = dict(status="COMPLETE", particle_count=n, duration_s=t, observation_limit_s=duration,
                       accepted_steps=accepted, nominal_dt_s=nominal, minimum_dt_s=minimum_dt,
                       wall_hit_events=int(wall_hits.sum()), wall_hit_particle_count=int((wall_hits>0).sum()),
                       contact_pair_events=contacts_count, maximum_projection_m=maximum_projection,
                       branch_metrics=branch, physics_banner=self.banner, physics_validity=validity,
                       timings=self.timings, elapsed_s=time.perf_counter()-start,
                       geometry_cache=dict(hits=self.cylinder.distance_cache_hits,misses=self.cylinder.distance_cache_misses,
                                           method="Exact repeated coordinates only; no quantization or field approximation"),
                       all_crossings_tested_at_accepted_substeps=True, saved_visualization_states=saved_count,
                       scope=self.config.get("experiment_scope", "STAGE4_SCOPED_CONTROL_STUDY"))
        write_json(run / "particles/SUMMARY.json", summary)
        return summary
