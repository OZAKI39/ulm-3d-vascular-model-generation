"""Bounded hard-sphere experiments. Actual forces and cavity checks use existing GPU kernels."""

import os
import time
import traceback
import resource
from dataclasses import replace
import numpy as np
import h5py
from mevo_mb.acoustics.water_disc_proxy import CachedWaterBasis
from mevo_mb.coupling.pdms_cavity import apparatus
from mevo_mb.fast_sim.models import ResponseSection
from mevo_mb.fast_sim.nearfield import F2Physics, pair_reference
from mevo_mb.fast_sim.water_hard_sphere import WaterHardSphere
from .common import BASE, OLD, CONFIG, read, write, sha, identity, event


class Campaign:
    def __init__(self, out):
        self.started = time.perf_counter()
        self.prior_s = sum(read(p)["wall_s"] for p in BASE.glob("*/CAMPAIGN.json"))
        self.out = out
        self.p = read(CONFIG)
        self.p.update(read(CONFIG.parent / "algorithm_revision_2.json"))
        self.reg = read(BASE / "PREREGISTRATION.json")
        for rel, digest in self.reg["input_sha256"].items():
            assert sha(OLD / rel) == digest, rel
        self.x0 = np.load(OLD / "initial_states.npz")["aggregate_m"]
        self.physics = F2Physics(radius_m=self.p["radius_m"], frequency_Hz=self.p["frequency_Hz"])
        self.response = ResponseSection(OLD / "inputs/response.h5", self.physics)
        self.cavity, _ = apparatus(read(OLD / "CONFIG_SNAPSHOT.json")["geometry"])
        self.bases = {k: CachedWaterBasis(OLD / "cache" / f"{k}.h5") for k in ["coarse", "fine"]}
        self.records = []
        self.ident = identity()

    def engine(self, cache):
        return WaterHardSphere(
            self.x0, self.physics, self.response, self.bases[cache], self.cavity, policy=self.p
        )

    def checks(self):
        import warp as wp
        from mevo_mb.fast_sim import vessel_kernels

        start = time.perf_counter()
        sim = self.engine("coarse")
        # Compile actual cleaned module; absence of f64 cannot break launched kernels.
        wp.load_module(module=vessel_kernels, device=sim.sim.device)
        checks = {}
        evidence = {}
        for name in ["OFF", "SHAPE_REPLAY", "OPTIMIZED_REPLAY"]:
            e = event(name, 0.002)
            sim.control(e)
            f = sim.fields(self.x0)
            pb, gb = self.bases["coarse"].query_basis(self.x0)
            u = self.bases["coarse"].vector(e["amplitude_Pa"], e["phase_rad"], e["enabled"])
            p = u @ pb
            g = np.einsum("i,inj->nj", u, gb)
            v = self.response.query(p)
            primary = -0.5 * np.real(v[:, 0, None] * g.conj())
            secondary = pair_reference(self.x0, v, p, self.physics, sim.sim.chi_zero)
            far = pair_reference(
                self.x0, v, p, replace(self.physics, force_model="f1_smooth_near_field"), sim.sim.chi_zero
            )
            errors = {}
            for key, ref in [
                ("pressure_Pa", p),
                ("primary_force_N", primary),
                ("secondary_force_N", secondary),
            ]:
                err = float(np.max(abs(f[key] - ref)))
                errors[key] = err
                checks[name + "_" + key] = err < max(float(np.max(abs(ref))) * 1e-9, 1e-24)
            evidence[name] = dict(
                pressure_Pa=p,
                primary_force_N=primary,
                secondary_force_N=secondary,
                no_feedback_regularized_secondary_N=far,
                nearfield_feedback_increment_N=secondary - far,
                errors=errors,
                controls=e,
            )
        checks["valid_center_wall_no_correction"] = sim.wall_check(self.x0) == dict(
            correction_m=0.0, error_bits=0
        )
        # Artificial wall-only input is an engineering rejection fixture, never a trajectory.
        outside = self.x0.copy()
        outside[0, 0] = self.cavity.L + 2 * self.p["radius_m"]
        w = sim.wall_check(outside)
        checks["sealed_X_rejects"] = bool(w["error_bits"] or w["correction_m"])
        sim.sim.error.zero_()
        top = self.x0.copy()
        top[0] = [self.cavity.ports[1], 0, self.cavity.ceiling]
        w = sim.wall_check(top)
        checks["numerical_top_scope_rejects"] = bool(w["error_bits"])
        sim.sim.error.zero_()
        for amp in [1e-9, 1e6]:
            try:
                sim.control(dict(event("BASELINE", 0.002), amplitude_Pa=[amp] * 4))
                sim.fields(self.x0)
                ok = False
            except ValueError:
                ok = True
            checks[f"pressure_{amp}_rejects"] = ok
            sim.sim.error.zero_()
        write(
            self.out / "CHECKS.json",
            dict(
                status="PASS" if all(checks.values()) else "FAIL",
                checks=checks,
                evidence=evidence,
                wall_s=time.perf_counter() - start,
                identity=self.ident,
            ),
        )
        assert all(checks.values()), checks

    def run(self, name, source, duration=0.002, dt=1e-6, cache="coarse", control=None):
        folder = self.out / "runs" / name
        folder.mkdir(parents=True, exist_ok=False)
        e = event(source, duration) if control is None else dict(control, duration_s=duration)
        write(
            folder / "INPUT.json",
            dict(
                control=e,
                source_F10N=source,
                dt_s=dt,
                duration_s=duration,
                cache=cache,
                preregistration_sha256=sha(BASE / "PREREGISTRATION.json"),
                identity=self.ident,
            ),
        )
        begin = time.perf_counter()
        sim = self.engine(cache)
        sim.control(e)
        x = [self.x0.copy()]
        times = [0.0]
        step_rows = []
        fields = []
        reason = []
        last = None
        steps = round(duration / dt)
        integration = time.perf_counter()
        try:
            for j in range(steps):
                if time.perf_counter() - begin > self.p["per_run_budget_s"]:
                    reason = ["PER_RUN_WALL_BUDGET"]
                    break
                if self.prior_s + time.perf_counter() - self.started > self.p["total_numerical_budget_s"]:
                    reason = ["CAMPAIGN_WALL_BUDGET"]
                    break
                last = sim.step(dt)
                r = last["result"]
                if j % max(1, round(1e-4 / dt)) == 0 or not last["accepted"]:
                    fields.append(dict(time_s=j * dt, **last["fields"], contact_force_N=r["contact_force_N"]))
                step_rows.append(
                    [
                        last["accepted"],
                        r["mobility_equation_residual_m"],
                        r["geometric_overlap_m"].max(initial=0),
                        r["max_single_iteration_m"],
                        r["max_total_contact_correction_m"],
                        last["wall_correction_m"],
                        last["full_step_s"],
                        r["total_s"],
                        last["resistance_s"],
                        r["qp_iterations"],
                        np.linalg.norm(r["prediction_m"] - x[-1], axis=1).sum(),
                        np.linalg.norm(r["contact_displacement_m"], axis=1).sum(),
                        r["primal_violation_m"],
                        r["dual_negative_N"],
                        r["force_equation_residual_N"],
                        r["complementarity_equivalent_gap_m"],
                        r["complementarity_N_m"],
                        r["actual_complementarity_equivalent_gap_m"],
                    ]
                )
                if not last["accepted"]:
                    reason = last["reasons"]
                    break
                x.append(sim.x.copy())
                times.append((j + 1) * dt)
            if last is not None:
                write(folder / "LAST_STEP.json", last)
        except Exception as ex:
            reason = ["EXCEPTION: " + repr(ex)]
            (folder / "exception.txt").write_text(traceback.format_exc())
        integration_s = time.perf_counter() - integration
        rows = np.array(step_rows).reshape(-1, 18)
        x = np.array(x)
        with h5py.File(folder / "trajectory.h5", "w") as h:
            h.create_dataset("positions_m", data=x, compression="gzip")
            h.create_dataset("time_s", data=times)
            h.create_dataset("particle_id", data=np.arange(32))
            h.create_dataset("step_metrics", data=rows)
            h["step_metrics"].attrs["columns"] = (
                "accepted,residual_m,overlap_m,max_iteration_m,total_contact_m,wall_m,full_step_s,solver_s,resistance_s,iterations,physical_path_sum_m,contact_displacement_sum_m,primal_m,dual_N,force_residual_N,complementarity_m,complementarity_N_m,actual_complementarity_m"
            )
            h.attrs["flow"] = "STRICT_ZERO"
            h.attrs["execution_commit"] = self.ident["execution_commit"]
            if reason and last is not None:
                h.create_dataset("failure_proposal_m", data=last["result"]["positions_m"])
            ff = h.create_group("force_samples")
            if fields:
                for key in fields[0]:
                    ff.create_dataset(key, data=np.array([f[key] for f in fields]))
        disp = x[-1].mean(0) - x[0].mean(0)
        result = dict(
            name=name,
            source=source,
            status="PASS" if not reason and len(x) == steps + 1 else "FAIL",
            reasons=reason,
            dt_s=dt,
            cache=cache,
            requested_s=duration,
            completed_s=times[-1],
            accepted_steps=len(x) - 1,
            attempted_steps=len(rows),
            displacement_m=disp,
            integration_s=integration_s,
            wall_s=time.perf_counter() - begin,
            mean_step_s=float(rows[:, 6].mean()) if len(rows) else None,
            max_residual_m=float(rows[:, 1].max(initial=0)),
            max_overlap_m=float(rows[:, 2].max(initial=0)),
            max_contact_iteration_m=float(rows[:, 3].max(initial=0)),
            max_contact_total_m=float(rows[:, 4].max(initial=0)),
            max_wall_correction_m=float(rows[:, 5].max(initial=0)),
            max_primal_m=float(rows[:, 12].max(initial=0)),
            max_dual_N=float(rows[:, 13].max(initial=0)),
            max_force_residual_N=float(rows[:, 14].max(initial=0)),
            max_complementarity_m=float(rows[:, 15].max(initial=0)),
            max_complementarity_N_m=float(rows[:, 16].max(initial=0)),
            max_actual_complementarity_m=float(rows[:, 17].max(initial=0)),
            rss_peak_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024,
            trajectory_sha256=sha(folder / "trajectory.h5"),
            training_ready=False,
        )
        write(folder / "RESULT.json", result)
        self.records.append(result)
        write(self.out / "RUNS.json", self.records)
        print(
            name,
            result["status"],
            "t_ms=",
            times[-1] * 1e3,
            "step_ms=",
            (result["mean_step_s"] or 0) * 1e3,
            reason,
            flush=True,
        )
        return result


def main():
    mode = os.environ.get("F10P_MODE", "pilot")
    label = os.environ["F10P_ATTEMPT"]
    out = BASE / label
    out.mkdir(parents=True, exist_ok=False)
    start = time.perf_counter()
    try:
        c = Campaign(out)
        if mode == "pilot":
            from .foundations import run as foundations

            checks = foundations()
            write(out / "FOUNDATIONS.json", checks)
            assert checks["status"] == "PASS"
            c.checks()
            pilot = c.run("STRONG_150", "OPTIMIZED_REPLAY", duration=150e-6)
            if pilot["status"] == "PASS" and pilot["mean_step_s"] * 8000 < 350:
                c.run("STRONG_SHORT", "OPTIMIZED_REPLAY", duration=0.00025)
        elif mode == "validation":
            assert read(BASE / "pilot_v2/runs/STRONG_SHORT/RESULT.json")["status"] == "PASS"
            r = c.run("STRONG_DT1", "OPTIMIZED_REPLAY")
            if r["status"] == "PASS":
                r = c.run("STRONG_DT05", "OPTIMIZED_REPLAY", dt=0.5e-6)
                if r["status"] == "PASS":
                    c.run("STRONG_FINE", "OPTIMIZED_REPLAY", cache="fine")
            for source, tag in [("OFF", "OFF"), ("SHAPE_REPLAY", "WEAK")]:
                r = c.run(tag + "_DT1", source)
                if r["status"] == "PASS":
                    c.run(tag + "_DT05", source, dt=0.5e-6)
                    c.run(tag + "_FINE", source, cache="fine")
        elif mode == "diagnostic":
            c.run("STRONG_PREFIX_DT05", "OPTIMIZED_REPLAY", duration=0.0005, dt=0.5e-6)
            c.run("STRONG_PREFIX_FINE", "OPTIMIZED_REPLAY", duration=0.0005, cache="fine")
        elif mode == "navigation":
            from .navigation import run_navigation

            run_navigation(c)
        else:
            raise ValueError(mode)
        write(
            out / "CAMPAIGN.json",
            dict(
                status="COMPLETED",
                mode=mode,
                wall_s=time.perf_counter() - start,
                runs=c.records,
                identity=c.ident,
                peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024,
            ),
        )
    except Exception:
        (out / "EXCEPTION.txt").write_text(traceback.format_exc())
        write(
            out / "CAMPAIGN.json",
            dict(status="EXCEPTION", wall_s=time.perf_counter() - start, identity=identity()),
        )
        raise


if __name__ == "__main__":
    main()
