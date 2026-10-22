"""Bounded, sequential GPU control trials, fresh replay and dt/cache comparisons."""

import os
import time
import traceback
import numpy as np
import h5py
from scipy.optimize import minimize
from mevo_mb.acoustics.water_disc_proxy import CachedWaterBasis
from mevo_mb.acoustics.kwave_backend import GPUTelemetry
from mevo_mb.coupling.pdms_cavity import apparatus
from mevo_mb.fast_sim.water_control_engine import WaterControlSimulator
from mevo_mb.fast_sim.models import ResponseSection
from mevo_mb.fast_sim.nearfield import F2Physics, pair_reference
from mevo_mb.fast_sim.metrics import trajectory_metrics
from .common import BASE, configs, read, write, sha, identity


def shape_arrays(x):
    rel = x - x.mean(axis=-2, keepdims=True)
    Q = np.einsum("tni,tnj->tij", rel, rel) / x.shape[1]
    ev, axes = np.linalg.eigh(Q)
    principal = axes[:, :, -1]
    identified = ev[:, -1] - ev[:, -2] > 1e-3 * ev[:, -1]
    return dict(
        centroid_m=x.mean(1),
        tensor_m2=Q,
        rg_m=np.sqrt(np.trace(Q, axis1=1, axis2=2)),
        aspect=np.sqrt(ev[:, -1] / ev[:, 0]),
        directional_xy=np.sqrt(Q[:, 0, 0] / Q[:, 1, 1]),
        principal_axis=principal,
        axis_identifiable=identified,
        angle_to_x_deg=np.where(
            identified, np.degrees(np.arccos(np.clip(abs(principal[:, 0]), 0, 1))), np.nan
        ),
    )


class Experiment:
    def __init__(self, out):
        self.out = out
        self.started = time.perf_counter()
        self.c = configs()
        self.p = self.c["protocol"]
        assert self.c == read(BASE / "CONFIG_SNAPSHOT.json")
        self.reg = read(BASE / "PREREGISTRATION.json")
        assert self.reg["protocol"] == self.p
        assert sha(BASE / "initial_states.npz") == self.reg["initial_sha256"]
        self.initial = np.load(BASE / "initial_states.npz")
        self.cavity, _ = apparatus(self.c["geometry"])
        self.physics = F2Physics(radius_m=self.p["radius_m"], frequency_Hz=self.p["frequency_Hz"])
        self.response = ResponseSection(BASE / "inputs/response.h5", self.physics)
        assert sha(BASE / "inputs/response.h5") == self.p["response_sha256"]
        self.bases = {n: CachedWaterBasis(BASE / "cache" / f"{n}.h5") for n in ["coarse", "fine"]}
        for n in self.bases:
            assert sha(BASE / "cache" / f"{n}.h5") == self.reg["cache_sha256"][n]
        self.records = []
        self.static_log = []

    def event(self, v, dt=None):
        v = np.asarray(v, float)
        if v.shape != (7,) or np.any(v[:4] < 0) or np.any(v[:4] > 1):
            raise ValueError("Outside source bounds")
        return dict(
            time_s=0.0,
            duration_s=self.p["duration_s"],
            amplitude_Pa=(v[:4] * self.p["maximum_source_amplitude_Pa"]).tolist(),
            phase_rad=np.r_[0.0, v[4:]].tolist(),
            enabled=[True] * 4,
            flow_enabled=False,
        )

    def baseline(self):
        return np.r_[
            np.array(self.p["baseline_amplitude_Pa"]) / self.p["maximum_source_amplitude_Pa"],
            self.p["baseline_phase_rad"][1:],
        ]

    def sim(self, x, cache="coarse"):
        return WaterControlSimulator(
            x, self.physics, self.response, self.bases[cache], cavity=self.cavity, graph_chunk=32
        )

    def checks(self):
        start = time.perf_counter()
        x = self.initial["aggregate_m"]
        sim = self.sim(x)
        event = self.event(self.baseline())
        sim.control(event)
        sim.fields()
        basis = self.bases["coarse"]
        pb, gb = basis.query_basis(x)
        u = basis.vector(event["amplitude_Pa"], event["phase_rad"])
        p = u @ pb
        g = np.einsum("i,inj->nj", u, gb)
        v = self.response.query(p)
        f = -0.5 * np.real(v[:, 0, None] * np.conj(g))
        fp = sim.force.numpy()
        pg = sim.pressure.numpy()
        pg = pg[:, 0] + 1j * pg[:, 1]
        gg = sim.consistent_gradients.numpy()
        gg = gg[:, :, 0] + 1j * gg[:, :, 1]
        checks = dict(
            GPU_cubic_pressure=np.allclose(pg, p, rtol=1e-12, atol=1e-10),
            GPU_cubic_gradient=np.allclose(gg, g, rtol=1e-10, atol=1e-7),
            GPU_primary_matches_CPU=np.max(abs(fp - f)) < max(np.max(abs(f)) * 1e-10, 1e-26),
            strict_zero_flow=not sim.flow_enabled,
        )
        sim.pairs()
        secondary = sim.force.numpy() - fp
        f2 = pair_reference(x, v, p, self.physics, sim.chi_zero)
        checks["GPU_secondary_unchanged"] = np.max(abs(f2 - secondary)) < max(np.max(abs(f2)) * 1e-10, 1e-26)
        for label, amplitude in [("above_table", 1e6), ("below_table", 1e-9)]:
            sim.reset()
            bad = dict(event, amplitude_Pa=[amplitude] * 4)
            sim.control(bad)
            sim.fields()
            checks[label + "_strictly_rejected"] = int(sim.error.numpy()[0]) & 2 != 0
        sim.reset()
        sim.control(self.event(np.zeros(7)))
        sim.fields()
        sim.pairs()
        checks["OFF_exact_zero"] = (
            np.count_nonzero(sim.force.numpy()) == np.count_nonzero(sim.pressure.numpy()) == 0
        )
        checks["contact_and_wall_budgets_unchanged"] = (
            self.p["contact_residual_over_R"] == 0.01 and self.p["wall_correction_over_R"] == 0.1
        )
        # Control duration and no-flow checks need not launch a trajectory.
        for name, ev in [
            ("flow", dict(event, flow_enabled=True)),
            ("three_sources", dict(event, amplitude_Pa=[1] * 3)),
        ]:
            try:
                sim.control(ev)
                checks["reject_" + name] = False
            except ValueError:
                checks["reject_" + name] = True
        write(
            self.out / "GPU_CHECKS.json",
            dict(
                status="PASS" if all(checks.values()) else "FAIL",
                checks=checks,
                pressure_max_error_Pa=float(np.max(abs(pg - p))),
                gradient_max_error_Pa_m=float(np.max(abs(gg - g))),
                wall_s=time.perf_counter() - start,
                identity=identity(),
            ),
        )
        assert all(checks.values()), checks
        # Initial force response, not a claim of globally independent controllability.
        v0 = self.baseline()
        J = []
        for j in range(7):
            plus = v0.copy()
            minus = v0.copy()
            plus[j] += 1e-4
            minus[j] -= 1e-4
            J.append((self.static_velocity(plus)[0].mean(0) - self.static_velocity(minus)[0].mean(0)) / 2e-4)
        J = np.array(J).T
        write(
            self.out / "CONTROLLABILITY.json",
            dict(
                velocity_J_m_s=J,
                singular_values_m_s=np.linalg.svd(J, compute_uv=False),
                z_reflection_symmetry=True,
                scope="local COM force-response at registered symmetric initial state; not global controllability",
            ),
        )

    def static_velocity(self, v):
        e = self.event(v)
        basis = self.bases["coarse"]
        u = basis.vector(e["amplitude_Pa"], e["phase_rad"])
        pb, gb = self.static_fields
        p = u @ pb
        g = np.einsum("i,inj->nj", u, gb)
        if (
            np.any((abs(p) > 0) & (abs(p) < 0.1))
            or np.max(abs(p)) > self.p["maximum_registered_particle_pressure_Pa"]
        ):
            raise ValueError("Actual particle pressure outside registered admission")
        vv = self.response.query(p)
        f = -0.5 * np.real(vv[:, 0, None] * np.conj(g))
        return f / (6 * np.pi * self.physics.viscosity_Pa_s * self.physics.radius_m), p

    def optimize_static(self, start, mode):
        count = 0
        best = (np.inf, start.copy())

        def fun(v):
            nonlocal count, best
            count += 1
            if count > self.p["max_static_optimizer_evaluations"]:
                raise StopIteration
            try:
                vel, pressure = self.static_velocity(v)
                x = self.initial["aggregate_m"] + self.p["duration_s"] * vel
                if mode == "navigation":
                    score = float(np.linalg.norm(x.mean(0) - self.p["target_position_m"]) / 1e-6)
                else:
                    m = shape_arrays(x[None])
                    score = float(
                        abs(m["directional_xy"][0] - self.p["target_shape"]["directional_aspect_x_over_y"])
                        + np.linalg.norm(x.mean(0)) / self.p["target_shape"]["centroid_penalty_scale_m"]
                    )
                score += 1e-5 * np.sum(v[:4] ** 2)
                if score < best[0]:
                    best = (score, v.copy())
                why = None
            except ValueError as ex:
                score = 1e6
                pressure = []
                why = str(ex)
            self.static_log.append(
                dict(
                    mode=mode, parameters=v.copy(), objective=score, pressure_Pa=np.abs(pressure), invalid=why
                )
            )
            return score

        try:
            res = minimize(
                fun,
                start,
                method="SLSQP",
                bounds=[(0, 1)] * 4 + [(-np.pi, np.pi)] * 3,
                options=dict(maxiter=30, ftol=1e-10),
            )
            status = str(res.message)
        except StopIteration:
            status = "REGISTERED_STATIC_EVALUATION_LIMIT"
        return best[1], dict(
            mode=mode,
            calls=count,
            best_score=best[0],
            stop=status,
            scope="uncoupled initial primary-force predictor only; every selected control validated by full-particle rollout",
        )

    def trial(self, name, v, state="aggregate", cache="coarse", dt=None):
        if (
            len(self.records) >= self.p["max_total_trajectories"]
            or time.perf_counter() - self.started > self.p["campaign_wall_budget_s"]
        ):
            raise TimeoutError("Registered campaign budget")
        folder = self.out / "runs" / name
        folder.mkdir(parents=True, exist_ok=False)
        start = time.perf_counter()
        dt = self.p["dt_s"] if dt is None else dt
        ev = self.event(v)
        write(
            folder / "INPUT.json",
            dict(
                state=state,
                events=[ev],
                dt_s=dt,
                cache=cache,
                cache_sha256=self.reg["cache_sha256"][cache],
                preregistration_sha256=sha(BASE / "PREREGISTRATION.json"),
                identity=identity(),
            ),
        )
        x0 = self.initial[state + "_m"]
        sim = self.sim(x0, cache)
        result = sim.run(
            self.p["duration_s"],
            dt,
            [ev],
            sample_interval_s=self.p["sample_interval_s"],
            wall_budget_s=self.p["trajectory_wall_budget_s"],
        )
        summary = result.pop("summary")
        metrics, edges = trajectory_metrics(result["positions_m"], self.physics.radius_m, result["time_s"])
        shapes = shape_arrays(result["positions_m"][:, 0])
        pressures = []
        gradients = []
        primary = []
        secondary = []
        # Independent diagnostic reevaluation on saved states; no position edits.
        sim.error.zero_()
        for x in result["positions_m"][:, 0]:
            sim.x.assign(x)
            sim.control(ev)
            sim.fields()
            f = sim.force.numpy().copy()
            sim.pairs()
            pp = sim.pressure.numpy()
            gg = sim.consistent_gradients.numpy()
            pressures.append(pp[:, 0] + 1j * pp[:, 1])
            gradients.append(gg[:, :, 0] + 1j * gg[:, :, 1])
            primary.append(f)
            secondary.append(sim.force.numpy() - f)
        with h5py.File(folder / "trajectory.h5", "x") as h:
            for key, val in result.items():
                if isinstance(val, np.ndarray):
                    h.create_dataset(key, data=val, compression="gzip")
            for key, val in result.get("flow", {}).items():
                if isinstance(val, np.ndarray):
                    h.create_dataset("flow/" + key, data=val, compression="gzip")
            if "failure_proposal" in result:
                for key, val in result["failure_proposal"].items():
                    h.create_dataset("failure_proposal/" + key, data=val)
            for key, val in dict(
                pressure_Pa=pressures,
                gradient_Pa_m=gradients,
                primary_force_N=primary,
                secondary_force_N=secondary,
            ).items():
                h.create_dataset(key, data=np.asarray(val), compression="gzip")
            h.create_dataset("particle_id", data=np.arange(32))
            h.attrs["radius_m"] = self.physics.radius_m
            h.attrs["execution_commit"] = identity()["execution_commit"]
            h.attrs["model"] = self.p["model"]
        np.savez_compressed(folder / "metrics.npz", **metrics, **{"shape_" + k: v for k, v in shapes.items()})
        write(folder / "EDGE_EVENTS.json", edges)
        end = shapes["centroid_m"][-1]
        target = np.array(self.p["target_position_m"])
        quality = summary["quality"]
        eligible = (
            summary["status"] == "COMPLETE"
            and quality["error_bits"] == 0
            and quality["pressure_max_Pa"] <= self.p["maximum_registered_particle_pressure_Pa"]
        )
        score = (
            float(np.linalg.norm(end - target) / 1e-6)
            + 1000 * (not eligible)
            + 10 * (metrics["cluster_count"][-1, 0] - 1)
        )
        score += (
            quality["max_accepted_penetration_m"] / self.physics.radius_m
            + quality["maximum_wall_correction_m"] / self.physics.radius_m
        )
        record = dict(
            name=name,
            state=state,
            cache=cache,
            dt_s=dt,
            parameters=np.asarray(v),
            eligible=bool(eligible),
            status=summary["status"],
            objective=score,
            initial_error_m=float(np.linalg.norm(x0.mean(0) - target)),
            endpoint_error_m=float(np.linalg.norm(end - target)),
            displacement_m=end - x0.mean(0),
            cluster_count_final=int(metrics["cluster_count"][-1, 0]),
            cluster_count_max=int(metrics["cluster_count"][:, 0].max()),
            directional_xy_initial=float(shapes["directional_xy"][0]),
            directional_xy_final=float(shapes["directional_xy"][-1]),
            rg_initial_m=float(shapes["rg_m"][0]),
            rg_final_m=float(shapes["rg_m"][-1]),
            source_cost_Pa2_s=float(self.p["duration_s"] * np.sum(np.array(ev["amplitude_Pa"]) ** 2)),
            quality=quality,
            integration_wall_s=summary["wall_s"],
            end_to_end_wall_s=time.perf_counter() - start,
            trajectory_sha256=sha(folder / "trajectory.h5"),
        )
        self.records.append(record)
        write(folder / "SUMMARY.json", summary)
        write(folder / "RESULT.json", record)
        write(self.out / "TRIALS.json", self.records)
        print(
            name,
            summary["status"],
            "target_um",
            record["endpoint_error_m"] * 1e6,
            "d_um",
            np.round(record["displacement_m"] * 1e6, 6),
            "clusters",
            record["cluster_count_final"],
            "wall_s",
            record["end_to_end_wall_s"],
            flush=True,
        )
        return record

    def run(self):
        self.static_fields = self.bases["coarse"].query_basis(self.initial["aggregate_m"])
        with GPUTelemetry() as telemetry:
            self.checks()
            self.trial("OFF", np.zeros(7))
            baseline = self.trial("BASELINE", self.baseline())
            candidates = [baseline]
            optimizers = []
            for index, start in enumerate([self.baseline(), np.r_[[0.7] * 4, [-np.pi / 2, 0, np.pi / 2]]]):
                v, info = self.optimize_static(start, "navigation")
                optimizers.append(info)
                candidates.append(self.trial("NAV_STATIC_" + str(index), v))
            # At most four sequential actual-trajectory local perturbations.
            for i, (axis, delta) in enumerate([(4, 0.3), (5, 0.3), (6, -0.3), (0, -0.15)]):
                best = min(candidates, key=lambda r: r["objective"])
                v = best["parameters"].copy()
                lo, hi = (0, 1) if axis < 4 else (-np.pi, np.pi)
                v[axis] = np.clip(v[axis] + delta, lo, hi)
                candidates.append(self.trial("NAV_LOCAL_" + str(i), v))
            best = min(candidates, key=lambda r: r["objective"])
            write(
                self.out / "NAV_SELECTION.json",
                dict(
                    selected=best["name"],
                    objective=best["objective"],
                    optimizers=optimizers,
                    target_unchanged=self.p["target_position_m"],
                    candidate_count=len(candidates),
                ),
            )
            self.trial("OPTIMIZED_REPLAY", best["parameters"])
            self.trial("OPTIMIZED_HALF_DT", best["parameters"], dt=self.p["dt_s"] / 2)
            self.trial("OPTIMIZED_FINE_CACHE", best["parameters"], cache="fine")
            self.trial("BASELINE_HALF_DT", self.baseline(), dt=self.p["dt_s"] / 2)
            self.trial("DISPERSED_OFF", np.zeros(7), state="dispersed")
            self.trial("DISPERSED_OPTIMIZED", best["parameters"], state="dispersed")
            # Shape trials follow completed position comparisons, not parallel claims.
            shape = []
            for i, start in enumerate([self.baseline(), best["parameters"]]):
                v, info = self.optimize_static(start, "shape")
                optimizers.append(info)
                shape.append(self.trial("SHAPE_STATIC_" + str(i), v))
            for i, phase in enumerate([-np.pi / 2, np.pi / 2]):
                v = self.baseline()
                v[4:] = [0, phase, phase]
                shape.append(self.trial("SHAPE_DIRECTION_" + str(i), v))

            def shape_objective(r):
                return (
                    abs(r["directional_xy_final"] - self.p["target_shape"]["directional_aspect_x_over_y"])
                    + np.linalg.norm(r["displacement_m"]) / self.p["target_shape"]["centroid_penalty_scale_m"]
                    + 1000 * (not r["eligible"])
                    + 10 * (r["cluster_count_final"] - 1)
                )

            sb = min(shape, key=shape_objective)
            write(
                self.out / "SHAPE_SELECTION.json",
                dict(
                    selected=sb["name"],
                    objective=shape_objective(sb),
                    optimizers=optimizers,
                    candidate_count=len(shape),
                ),
            )
            self.trial("SHAPE_REPLAY", sb["parameters"])
            self.trial("SHAPE_HALF_DT", sb["parameters"], dt=self.p["dt_s"] / 2)
            self.trial("SHAPE_FINE_CACHE", sb["parameters"], cache="fine")
        write(self.out / "STATIC_SEARCH.json", self.static_log)
        write(self.out / "GPU_TELEMETRY.json", dict(samples=telemetry.samples, errors=telemetry.errors))
        record = dict(
            status="COMPLETE",
            trajectory_count=len(self.records),
            wall_s=time.perf_counter() - self.started,
            integration_wall_s=sum(r["integration_wall_s"] for r in self.records),
            trial_end_to_end_wall_s=sum(r["end_to_end_wall_s"] for r in self.records),
            sampled_device_peak_MiB=max(s["used_MiB"] for s in telemetry.samples),
            basis_combinations={n: b.combine_count for n, b in self.bases.items()},
            actual_wave_solves=0,
            identity=identity(),
            training_ready=False,
        )
        write(self.out / "CAMPAIGN.json", record)
        write(
            BASE / "ACCEPTED_CAMPAIGN.json",
            dict(path=str(self.out.relative_to(BASE)), sha256=sha(self.out / "CAMPAIGN.json")),
        )
        print(record, flush=True)


if __name__ == "__main__":
    label = os.environ.get("F10N_ATTEMPT", "campaign_v1")
    if not label.replace("_", "").isalnum():
        raise ValueError("Simple label required")
    out = BASE / label
    out.mkdir(exist_ok=False)
    try:
        Experiment(out).run()
    except Exception:
        write(out / "FAILURE.json", dict(traceback=traceback.format_exc(), identity=identity()))
        raise
