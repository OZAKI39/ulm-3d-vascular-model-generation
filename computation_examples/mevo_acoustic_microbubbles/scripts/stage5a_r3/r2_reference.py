"""Preserved R2 numerical driver, new R3 output root and read-only original inputs.

Only paths/provenance differ; R2 kernel source hashes are checked before execution.
"""

import argparse
from collections import Counter
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import signal
import time
import traceback
import zipfile
import h5py
import numpy as np
from scripts.stage4l.dense_probe import load_case, sha, write
from scripts.stage5a_r2.resources import ResourceMonitor
from scripts.stage5a_r3.common import ROOT, BASE
from mevo_mb.bubbles.tracker import ParticleTracker
from mevo_mb.bubbles.motion_r2 import SemiImplicit
from mevo_mb.bubbles.motion_recovery_r2 import (
    selected_method,
    recovery_ceiling,
    begin_recovery,
    accepted_recovery,
)
from mevo_mb.bubbles.motion_sweep_r2 import WallSweep
from mevo_mb.bubbles.constrained_motion import SubdivideInterval, MotionSolverFailure, AcceptedCache


def checkpoint(path, x, u, r, ids, mask, t, controller, cache, identity):
    value = dict(
        time_s=t,
        positions_m=x,
        velocity_m_s=u,
        radii_m=r,
        particle_id=ids,
        active=mask,
        metadata_json=json.dumps(
            dict(
                controller=controller,
                identity=identity,
                cache=dict(
                    particle_ids=cache.particle_ids,
                    primal=[[k, np.asarray(v).tolist()] for k, v in cache.primal.items()],
                    dual=[[list(k), v] for k, v in cache.dual.items()],
                    active=cache.active,
                    sparsity_pattern=cache.sparsity_pattern,
                ),
            ),
            allow_nan=False,
        ),
    )
    tmp = path.with_suffix(".tmp")
    with tmp.open("wb") as f:
        np.savez_compressed(f, **value)
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, path)


def restore(path):
    with np.load(path) as f:
        values = {k: f[k] for k in f.files if k != "metadata_json"}
        meta = json.loads(str(f["metadata_json"]))
    c = meta["cache"]
    cache = AcceptedCache(
        tuple(c["particle_ids"]),
        {int(k): np.array(v) for k, v in c["primal"]},
        {tuple(k): v for k, v in c["dual"]},
        tuple(tuple(k) for k in c["active"]),
        c["sparsity_pattern"],
    )
    return values, meta, cache


def load(case, method, pair="numpy", components=True):
    config, field, flow, gold, inputs = load_case(Path("/workspace/mevo_acoustic_microbubbles"), case)
    selected = json.loads(
        (
            Path("/workspace/mevo_acoustic_microbubbles_stage5a_r2")
            / "results/stage5/motion_solver/audit_v2/BASELINE_SELECTION.json"
        ).read_text()
    )
    source = Path(next(r["run"] for r in selected if r["key"] == case + "_1ms"))
    snapshot = source / "INPUT_SNAPSHOT.npz"
    inputs[str(snapshot)] = sha(snapshot)
    with np.load(snapshot) as f:
        state = {k: f[k] for k in f.files}
    config["particles"].update(
        translation_backend="constrained_lubrication_r2",
        pair_backend=pair,
        response_backend="numpy",
        lubrication=dict(activation=0.1, floor=0.001, residual_tolerance=1e-11),
        constrained_solver=dict(solver="prox", components=components, component_workers=1),
    )
    tracker = ParticleTracker(config, field, flow)
    wall = Path(config["geometry"]["mesh_surface"]).parent / "mesh-surfaces/WALL.vtp"
    if not wall.exists():
        wall = Path(
            "/workspace/mevo_acoustic_microbubbles/results/stage4/simvascular/mesh/20261002T141441884393Z/SV_MESH/mesh-surfaces/WALL.vtp"
        )
    inputs[str(wall)] = sha(wall)
    sweep = WallSweep(tracker.cylinder, wall)
    return tracker, sweep, config, state, inputs


def run(
    case,
    duration,
    method="S1",
    max_dt=20e-6,
    tol_factor=1.0,
    pair="numpy",
    budget=1800.0,
    label="development",
    initial_dt=1e-8,
    components=True,
    resume=None,
    nonlinear_strategy="picard",
    nonsmooth_recovery=False,
):
    for filename in [
        "motion_r2.py",
        "motion_qp_r2.py",
        "motion_sweep_r2.py",
        "motion_recovery_r2.py",
        "constrained_motion.py",
    ]:
        local = ROOT / "src/mevo_mb/bubbles" / filename
        original = Path("/workspace/mevo_acoustic_microbubbles_stage5a_r2/src/mevo_mb/bubbles") / filename
        assert sha(local) == sha(original), "FROZEN_R2_KERNEL_CHANGED"
    if label.startswith("formal"):
        from scripts.stage5a_r2.freeze import verify

        protocol = verify()
        assert method == "S1" and pair == "cupy" and components
        assert nonlinear_strategy == protocol["integrator"]["nonlinear_strategy"]
        assert nonsmooth_recovery == protocol["integrator"].get("nonsmooth_S0_recovery", False)
        assert (max_dt, tol_factor) in [(1e-5, 1.0), (5e-6, 0.25)]
        assert budget <= protocol["resources"]["per_run_s"]
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    out = BASE / "benchmarks/r2_protected" / f"{stamp}_{case}_{method}_{duration:g}s_{label}"
    out.mkdir(parents=True)
    start = time.perf_counter()
    cpu = time.process_time()
    status = "INCOMPLETE"
    error = None
    failure = None
    tracker, sweep, config, state, inputs = load(case, method, pair, components)
    initialization = time.perf_counter() - start
    full, rfull, ids, mask = (
        state["positions_m"].copy(),
        state["radii_m"],
        state["particle_id"],
        state["active"],
    )
    t0 = float(state["time_s"])
    t = t0
    elapsed = 0.0
    ufull = np.zeros_like(full)
    controller = dict(
        next_dt=initial_dt, last_error=None, last_reason=None, consecutive_rejections=0, elapsed_s=0.0
    )
    if resume:
        data, meta, cache = restore(Path(resume))
        full = data["positions_m"]
        ufull = data["velocity_m_s"]
        if not all(
            np.array_equal(data[k], state[k], equal_nan=True) for k in ["radii_m", "particle_id", "active"]
        ):
            raise ValueError("CONFIG_ERROR: restart particle identity differs")
        previous_identity = Path(resume).parent / "IDENTITY.json"
        if sha(previous_identity) != meta["identity"]:
            raise ValueError("CONFIG_ERROR: checkpoint identity hash mismatch")
        old_identity = json.loads(previous_identity.read_text())
        for name, digest in old_identity["source_sha256"].items():
            if name.startswith("src/") and sha(ROOT / name) != digest:
                raise ValueError("CONFIG_ERROR: restart solver source differs")
        if old_identity["environment_sha256"] != sha(BASE / "environment/IDENTITY.json"):
            raise ValueError("CONFIG_ERROR: restart environment differs")
        if (
            old_identity["config"] != config
            or old_identity["method"] != method
            or old_identity["tolerance_factor"] != tol_factor
            or old_identity.get("nonlinear_strategy", "picard") != nonlinear_strategy
            or old_identity.get("nonsmooth_recovery", False) != nonsmooth_recovery
        ):
            raise ValueError("CONFIG_ERROR: exact continuation configuration differs")
        t = float(data["time_s"])
        controller = meta["controller"]
        elapsed = controller.get("elapsed_s", t - t0)
        tracker.translation.cache = cache
    x = full[mask].copy()
    r = rfull[mask]
    active = ids[mask]
    step = SemiImplicit(
        tracker.translation,
        tracker.cylinder,
        tracker.evaluate,
        sweep,
        method=method,
        atol_R=0.002 * tol_factor,
        rtol=0.01 * tol_factor,
        gap_rtol=0.02 * tol_factor,
        nonlinear_strategy=nonlinear_strategy,
    )
    sources = {
        str(p.relative_to(ROOT)): sha(p)
        for folder in ["src", "scripts/stage5a_r2", "scripts/stage5a_r3"]
        for p in (ROOT / folder).rglob("*.py")
    }
    with zipfile.ZipFile(out / "SOURCE.zip", "x", zipfile.ZIP_DEFLATED) as z:
        for name in sources:
            z.write(ROOT / name, name)
    identity = dict(
        case=case,
        method=method,
        duration_s=duration,
        integrator_variant="S1_V6_WITH_BOUNDED_S0_WALL_RECOVERY" if nonsmooth_recovery else "S1_V4_OR_S0",
        config=config,
        max_dt_s=max_dt,
        tolerance_factor=tol_factor,
        initial_dt_s=initial_dt,
        active_count=len(r),
        source_sha256=sources,
        input_sha256=inputs,
        environment_sha256=sha(BASE / "environment/IDENTITY.json"),
        plan_sha256=sha(ROOT / "docs/STAGE5A_R2_MOTION_PLAN.md"),
        source_commit=(ROOT / "SOURCE_COMMIT_R3.txt").read_text().strip()
        if (ROOT / "SOURCE_COMMIT_R3.txt").exists()
        else "UNCOMMITTED_DEVELOPMENT",
        resume=str(resume) if resume else None,
        wall_budget_s=budget,
        pair_backend=pair,
        components=components,
        nonlinear_strategy=nonlinear_strategy,
        nonsmooth_recovery=nonsmooth_recovery,
    )
    write(out / "IDENTITY.json", identity)
    np.savez_compressed(out / "INPUT_SNAPSHOT.npz", **state)
    write(out / "CONFIG.json", config)
    statepath = out / "CHECKPOINT.npz"
    checkpoint(
        statepath,
        full,
        ufull,
        rfull,
        ids,
        mask,
        t,
        controller,
        tracker.translation.cache,
        sha(out / "IDENTITY.json"),
    )
    times = [t]
    positions = [full.copy()]
    accepted = []
    rejected = []
    trials = []
    io_seconds = 0.0
    interrupted = False
    repeated_failures = Counter()
    history = h5py.File(out / "trajectories.h5", "x")
    history.attrs.update(
        complete=False,
        method=method,
        all_accepted_substeps_recorded=True,
        actual_path="piecewise linear accepted half steps; no coordinate interpolation",
    )
    history.create_dataset("time_s", data=[t], maxshape=(None,), chunks=(128,))
    history.create_dataset(
        "positions_m", data=[full], maxshape=(None, *full.shape), chunks=(1, *full.shape), compression="gzip"
    )
    history.create_dataset(
        "active", data=[mask], maxshape=(None, len(mask)), chunks=(1, len(mask)), compression="gzip"
    )
    history["particle_id"] = ids
    history["radii_m"] = rfull
    history.flush()
    last_report = time.perf_counter()

    def check_budget():
        if time.perf_counter() - start > budget or len(accepted) >= 100000 or len(rejected) >= 10000:
            raise TimeoutError("PREREGISTERED_RUNTIME_BUDGET")

    def alarm(signum, frame):
        raise TimeoutError("WATCHDOG_WALL_BUDGET")

    prior = signal.signal(signal.SIGALRM, alarm)
    signal.setitimer(signal.ITIMER_REAL, budget + 10)
    with (
        ResourceMonitor() as monitor,
        (out / "ACCEPTED.jsonl").open("x") as al,
        (out / "REJECTED.jsonl").open("x") as rl,
    ):
        try:
            while elapsed < duration - 1e-15:
                check_budget()
                next_clock = (np.floor((elapsed + 1e-14) / 1e-5) + 1) * 1e-5
                dt = min(controller["next_dt"], max_dt, duration - elapsed, next_clock - elapsed)
                step.method = selected_method(controller, method) if nonsmooth_recovery else method
                if step.method == "S0" and nonsmooth_recovery:
                    dt = min(dt, recovery_ceiling(tol_factor))
                if dt < 1e-12 or controller["consecutive_rejections"] >= 30:
                    raise RuntimeError("MINIMUM_PHYSICAL_PROGRESS")
                trial_start = time.perf_counter()
                try:
                    states, d = step.trial(x, r, active, t, dt)
                except (SubdivideInterval, MotionSolverFailure) as exc:
                    reason = getattr(exc, "reason", exc.diagnostics.get("reason", "UNKNOWN"))
                    row = dict(
                        trial_method=step.method,
                        time_s=t,
                        elapsed_s=elapsed,
                        dt_s=dt,
                        reason=reason,
                        diagnostics=exc.diagnostics,
                        wall_s=time.perf_counter() - start,
                        trial_wall_s=time.perf_counter() - trial_start,
                    )
                    rejected.append(row)
                    qp_state = getattr(exc, "qp_state", None)
                    if qp_state is not None:
                        np.savez_compressed(out / f"FAILED_QP_{len(rejected):05d}.npz", **qp_state)
                    fingerprint = exc.diagnostics.get("input_sha256")
                    if fingerprint:
                        repeated_failures[fingerprint] += 1
                        if repeated_failures[fingerprint] >= 3:
                            raise RuntimeError(
                                "REPEATED_IDENTICAL_QP_FAILURE: stop ineffective interval retries"
                            ) from exc
                    if reason == "QP_INFEASIBILITY_UNRESOLVED":
                        raise RuntimeError(
                            "INFEASIBILITY_REPORTED: independent diagnosis required; no mechanical bisection"
                        ) from exc
                    rl.write(json.dumps(row, allow_nan=False) + "\n")
                    rl.flush()
                    if reason in [
                        "NONFINITE_INPUT",
                        "INVALID_MODEL",
                        "CONFIG_ERROR",
                        "GEOMETRY_DATA_CORRUPTED",
                    ]:
                        raise
                    controller.update(
                        next_dt=dt * 0.5,
                        last_reason=reason,
                        consecutive_rejections=controller["consecutive_rejections"] + 1,
                    )
                    if nonsmooth_recovery and step.method == "S1":
                        begin_recovery(
                            controller,
                            exc,
                            dt,
                            min(max_dt, duration - elapsed, next_clock - elapsed),
                            tol_factor,
                        )
                    check_budget()
                    continue
                trials.append(d)
                io_start = time.perf_counter()
                for k, (y, u, cache, sub) in enumerate(states):
                    time_sub = t + dt * (k + 1) / 2
                    full[mask] = y
                    ufull[mask] = u
                    row = dict(
                        substep_method=sub["substep_method"],
                        resistance_evaluation=sub["resistance_evaluation"],
                        time_s=time_sub,
                        elapsed_s=elapsed + dt * (k + 1) / 2,
                        dt_s=dt / 2,
                        wall_s=time.perf_counter() - start,
                        error_ratio=d["error_ratio"],
                        error_terms=d["error_terms"],
                        kkt=sub["kkt"],
                        linear_residual=sub["linear_residual"],
                        nonlinear_kkt=sub["qp"]["nonlinear_kkt"],
                        nonlinear_iterations=sub["iterations"],
                        qp_iterations=sub["qp_iterations"],
                        fallbacks=sub["fallbacks"],
                        force_residual_l2_N=sub["qp"]["physical_force_residual_l2_N"],
                        force_residual_relative=sub["qp"]["physical_force_residual_relative"],
                        dt_geometry_error_bound_m=dt
                        / 2
                        * sub["qp"]["physical_constraint_velocity_error_m_s"],
                        qp_timings=sub["qp_timings"],
                        workspace_hits=sub["workspace_hits"],
                        component_solves=sub["component_solves"],
                        component_sizes=sub["qp"]["component_sizes"],
                        safety=sub["safety"],
                    )
                    accepted.append(row)
                    al.write(json.dumps(row, allow_nan=False) + "\n")
                    times.append(time_sub)
                    positions.append(full.copy())
                    count = len(history["time_s"])
                    for key, value in [("time_s", time_sub), ("positions_m", full), ("active", mask)]:
                        dataset = history[key]
                        dataset.resize(count + 1, axis=0)
                        dataset[count] = value
                x = states[-1][0]
                repeated_failures.clear()
                tracker.translation.cache = states[-1][2]
                elapsed += dt
                t = t0 + elapsed
                factor = min(2.0, max(0.5, 0.85 / max(d["error_ratio"], 1e-8) ** 0.5))
                controller.update(
                    next_dt=accepted_recovery(controller, step.method, dt, dt * factor, max_dt),
                    last_error=d["error_ratio"],
                    last_reason="ACCEPTED",
                    consecutive_rejections=0,
                    elapsed_s=elapsed,
                )
                controller["total_accepted_steps"] = controller.get("total_accepted_steps", 0) + 2
                checkpoint(
                    statepath,
                    full,
                    ufull,
                    rfull,
                    ids,
                    mask,
                    t,
                    controller,
                    tracker.translation.cache,
                    sha(out / "IDENTITY.json"),
                )
                history.flush()
                al.flush()
                rl.flush()
                io_seconds += time.perf_counter() - io_start
                if time.perf_counter() - last_report > 5:
                    checkpoint(
                        statepath,
                        full,
                        ufull,
                        rfull,
                        ids,
                        mask,
                        t,
                        controller,
                        tracker.translation.cache,
                        sha(out / "IDENTITY.json"),
                    )
                    al.flush()
                    rl.flush()
                    progress = dict(
                        run=out.name,
                        case=case,
                        elapsed_s=elapsed,
                        target_s=duration,
                        accepted=len(accepted),
                        rejected=len(rejected),
                        last_dt_s=dt,
                        wall_s=time.perf_counter() - start,
                    )
                    write(out / "PROGRESS.json", progress)
                    print(json.dumps(progress), flush=True)
                    last_report = time.perf_counter()
            status = "COMPLETE"
        except (Exception, KeyboardInterrupt) as exc:
            error = repr(exc)
            failure = getattr(exc, "diagnostics", None)
            status = "INCOMPLETE" if isinstance(exc, (TimeoutError, KeyboardInterrupt)) else "FAIL"
            interrupted = isinstance(exc, KeyboardInterrupt)
            (out / "ERROR.txt").write_text(traceback.format_exc())
        finally:
            signal.setitimer(signal.ITIMER_REAL, 0)
            signal.signal(signal.SIGALRM, prior)
            checkpoint(
                statepath,
                full,
                ufull,
                rfull,
                ids,
                mask,
                t,
                controller,
                tracker.translation.cache,
                sha(out / "IDENTITY.json"),
            )
    history.attrs["complete"] = status == "COMPLETE"
    history.flush()
    history.close()
    write(out / "RESOURCE_SAMPLES.json", monitor.rows)
    write(out / "COMPLETED_TRIALS.json", trials)
    dtvalues = [a["dt_s"] for a in accepted]
    wall = time.perf_counter() - start
    summary = dict(
        status=status,
        error=error,
        failure=failure,
        run=str(out),
        case=case,
        method=method,
        target_s=duration,
        integrator_variant="S1_V6_WITH_BOUNDED_S0_WALL_RECOVERY" if nonsmooth_recovery else "S1_V4_OR_S0",
        completed_duration_s=elapsed,
        accepted_steps=len(accepted),
        rejected_trials=len(rejected),
        rejection_reasons=dict(Counter(z["reason"] for z in rejected)),
        actual_dt_quantiles_s=np.quantile(dtvalues, [0, 0.5, 1]).tolist() if dtvalues else None,
        accepted_kkt_max=max((a["kkt"] for a in accepted), default=None),
        accepted_linear_residual_max=max((a["linear_residual"] for a in accepted), default=None),
        accepted_force_residual_l2_max_N=max((a["force_residual_l2_N"] for a in accepted), default=None),
        accepted_force_residual_relative_max=max(
            (a["force_residual_relative"] for a in accepted), default=None
        ),
        accepted_nonlinear_kkt_max=max(
            (a["nonlinear_kkt"] for a in accepted if a["nonlinear_kkt"] is not None), default=None
        ),
        accepted_pair_swept_gap_min_m=min(
            (a["safety"]["minimum_pair_swept_gap_m"] for a in accepted), default=None
        ),
        accepted_wall_swept_gap_bound_min_m=min(
            (a["safety"]["minimum_wall_swept_gap_bound_m"] for a in accepted), default=None
        ),
        projection_max_m=0.0 if accepted else None,
        wall_s=wall,
        initialization_s=initialization,
        process_cpu_s=time.process_time() - cpu,
        timings_exclusive_s=step.timings,
        tracker_timings_inclusive_s=tracker.timings,
        qp_calls=step.calls,
        pair_backend=pair,
        peak_rss_bytes=max((z["rss_bytes"] for z in monitor.rows), default=None),
        io_seconds=io_seconds,
        monitor_inclusive_seconds=monitor.overhead_s,
        segment_initial_time_s=times[0],
        segment_completed_duration_s=t - times[0],
        cumulative_accepted_steps=controller.get("total_accepted_steps", len(accepted)),
        accepted_method_counts=dict(Counter(a["substep_method"] for a in accepted)),
        explicit_recovery_duration_s=controller.get("recovery_duration_s", 0.0),
        explicit_recovery_bursts=controller.get("recovery_bursts", 0),
        explicit_recovery_max_dt_s=max(
            (a["dt_s"] * 2 for a in accepted if a["substep_method"] == "S0"), default=None
        ),
        recovery_policy="S0 start-state QP is audited separately; no claim of endpoint implicit equilibrium for S0 steps"
        if nonsmooth_recovery
        else None,
        source_inputs_unchanged=all(sha(p) == s for p, s in inputs.items()),
    )
    write(out / "SUMMARY.json", summary)
    print(json.dumps(summary), flush=True)
    if interrupted:
        raise KeyboardInterrupt("Safe checkpoint and evidence saved; stop campaign chain")
    return out, summary


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--case", choices=["B2", "B3"], required=True)
    p.add_argument("--duration", type=float, required=True)
    p.add_argument("--method", choices=["S0", "S1"], default="S1")
    p.add_argument("--max-dt", type=float, default=20e-6)
    p.add_argument("--tol-factor", type=float, default=1.0)
    p.add_argument("--pair", default="numpy")
    p.add_argument("--budget", type=float, default=1800)
    p.add_argument("--label", default="development")
    p.add_argument("--initial-dt", type=float, default=1e-8)
    p.add_argument("--global-qp", action="store_true")
    p.add_argument("--resume")
    p.add_argument("--nonsmooth-recovery", action="store_true")
    p.add_argument("--nonlinear-strategy", choices=["picard", "quasi_newton"], default="picard")
    a = p.parse_args()
    run(
        a.case,
        a.duration,
        a.method,
        a.max_dt,
        a.tol_factor,
        a.pair,
        a.budget,
        a.label,
        a.initial_dt,
        not a.global_qp,
        a.resume,
        a.nonlinear_strategy,
        a.nonsmooth_recovery,
    )
