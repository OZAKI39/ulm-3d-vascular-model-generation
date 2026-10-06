"""Exclusive finite Stage5A probes. Archived inputs/results are read-only."""

from __future__ import annotations
import argparse
import json
import time
import traceback
import threading
import subprocess
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
import h5py
import psutil
import xml.etree.ElementTree as ET
from scripts.stage4l.dense_probe import load_case, sha, write
from mevo_mb.bubbles.tracker import ParticleTracker
from mevo_mb.bubbles.constrained_motion import PredictorCorrector, METHOD, BACKEND
from mevo_mb.bubbles.translation_backends import LUBRICATION_PHYSICS, INTERFACE_SCOPE, NOT_INCLUDED

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "results/stage5/motion_solver"


class ResourceMonitor:
    def __init__(self):
        self.rows = []
        self.stop = threading.Event()
        self.process = psutil.Process()
        self.thread = threading.Thread(target=self.collect, daemon=True)

    def collect(self):
        self.process.cpu_percent()
        while not self.stop.is_set():
            row = dict(
                monotonic_s=time.perf_counter(),
                rss_bytes=self.process.memory_info().rss,
                cpu_percent=self.process.cpu_percent(),
            )
            try:
                values = (
                    subprocess.check_output(
                        [
                            "nvidia-smi",
                            "--query-gpu=utilization.gpu,memory.used",
                            "--format=csv,noheader,nounits",
                        ],
                        text=True,
                        timeout=3,
                    )
                    .strip()
                    .split(",")
                )
                row.update(gpu_percent=float(values[0]), gpu_memory_MiB=float(values[1]))
            except Exception as exc:
                row["gpu_observation_error"] = type(exc).__name__
            self.rows.append(row)
            self.stop.wait(0.5)

    def __enter__(self):
        self.thread.start()
        return self

    def __exit__(self, *args):
        self.stop.set()
        self.thread.join(timeout=5)


def run(case, duration, dt, qp_solver="sparse_active_set"):
    if duration not in [0.001, 0.005] or dt not in [20e-6, 10e-6]:
        raise ValueError("Only preregistered probes allowed")
    test_path = BASE / (
        "environment/osqp_20261006/combined_tests.xml" if qp_solver == "osqp" else "tests/attempt02.xml"
    )
    tests = ET.parse(test_path).getroot()
    if any(int(s.get(k, 0)) for s in tests.iter("testsuite") for k in ["failures", "errors"]):
        raise ValueError("Unit-test gate failed")
    bench_path = (
        BASE
        / "benchmarks"
        / ("OSQP_BENCHMARK_SUMMARY.json" if qp_solver == "osqp" else "BENCHMARK_SUMMARY.json")
    )
    bench = json.load(open(bench_path))
    if bench["status"] != "PASS":
        raise ValueError("Fixed-state equivalence gate required")
    if duration == 0.005:
        gate = json.load(
            open(BASE / "comparisons" / ("OSQP_GATE_1ms.json" if qp_solver == "osqp" else "GATE_1ms.json"))
        )
        if not gate["pass_gate"]:
            raise ValueError("Joint B2/B3 1 ms gate required")
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    suffix = "_osqp" if qp_solver == "osqp" else ""
    out = BASE / "runs" / f"{stamp}_{case}_{duration * 1000:g}ms_dt{dt * 1e6:g}us{suffix}"
    out.mkdir(parents=True, exist_ok=False)
    start = time.perf_counter()
    process_start = time.process_time()
    selection = json.load(open(BASE / "audit_v2/BASELINE_SELECTION.json"))
    baseline = Path(next(row["run"] for row in selection if row["key"] == f"{case}_{duration * 1000:g}ms"))
    baseline_identity = json.load(open(baseline / "IDENTITY.json"))
    config, field, flow, _, inputs = load_case(Path("/workspace/mevo_acoustic_microbubbles"), case)
    for path, digest in inputs.items():
        if baseline_identity["input_sha256"].get(path) != digest:
            raise ValueError("FROZEN_INPUT_HASH_MISMATCH " + path)
    snapshot = baseline / "INPUT_SNAPSHOT.npz"
    if sha(snapshot) != baseline_identity["snapshot_sha256"]:
        raise ValueError("FROZEN_SNAPSHOT_HASH_MISMATCH")
    inputs[str(snapshot)] = sha(snapshot)
    with np.load(snapshot) as h:
        full = h["positions_m"].copy()
        rfull = h["radii_m"].copy()
        ids = h["particle_id"].copy()
        mask = h["active"].copy()
        t0 = float(h["time_s"])
    active = ids[mask]
    x = full[mask].copy()
    r = rfull[mask]
    config["particles"].update(
        translation_backend=BACKEND,
        dt_s=dt,
        duration_s=duration,
        pair_backend=bench["selected_pair_backend"],
        lubrication=dict(activation=0.1, floor=0.001, residual_tolerance=1e-11),
        constrained_solver=dict(
            qp_solver=qp_solver,
            gap_limiter_mode="constrained_adaptive",
            components=True,
            warm_start=True,
            component_workers=bench["selected_component_workers"],
        ),
    )
    sources = {
        str(p.relative_to(ROOT)): sha(p)
        for folder in ["src", "scripts/stage5a"]
        for p in (ROOT / folder).rglob("*.py")
    }
    plan = ROOT / "docs/STAGE5A_MOTION_SOLVER_PLAN.md"
    method = METHOD
    if qp_solver == "osqp":
        from mevo_mb.bubbles.osqp_motion import METHOD as method
    supplement = ROOT / "docs/STAGE5A_OSQP_SUPPLEMENT.md"
    identity = dict(
        run=out.name,
        case=case,
        backend=BACKEND,
        method=method,
        qp_solver=qp_solver,
        numerical_supplement_sha256=sha(supplement) if qp_solver == "osqp" else None,
        config=config,
        dt_s=dt,
        duration_s=duration,
        input_sha256=inputs,
        source_sha256=sources,
        plan_sha256=sha(plan),
        baseline_run=str(baseline),
        initial_time_s=t0,
        branch="dev/stage5a-motion-constrained-20261006",
        source_commit=(ROOT / ("SOURCE_COMMIT_OSQP.txt" if qp_solver == "osqp" else "SOURCE_COMMIT.txt"))
        .read_text()
        .strip(),
        wall_budget_s=1800,
        step_budget=100000,
    )
    write(out / "IDENTITY.json", identity)
    write(out / "CONFIG.json", config)
    write(
        out / "PHYSICS_BANNER.json",
        dict(
            backend=BACKEND,
            production_physics=[
                LUBRICATION_PHYSICS,
                INTERFACE_SCOPE,
                "unchanged classical secondary Bjerknes A0",
            ],
            not_included=NOT_INCLUDED,
            experimental_numerics=True,
            full_near_contact_bubble_validation=False,
        ),
    )
    write(
        out / "SOLVER_BANNER.json",
        dict(
            method=method,
            device="CPU OSQP sparse direct" if qp_solver == "osqp" else "CPU sparse active set",
            pair_backend=bench["selected_pair_backend"],
            predictor_corrector_solves_per_trial=2,
            cache_commit="accepted only",
            gap_limiter_mode="constrained_adaptive",
            coordinate_interpolation=False,
        ),
    )
    tracker = ParticleTracker(config, field, flow)
    wall_distance = tracker.cylinder.wall_distance
    tracker.timings["wall_distance_s"] = 0.0

    def timed_wall_distance(x):
        begin = time.perf_counter()
        try:
            return wall_distance(x)
        finally:
            tracker.timings["wall_distance_s"] += time.perf_counter() - begin

    tracker.cylinder.wall_distance = timed_wall_distance
    stepper = PredictorCorrector(
        tracker.translation,
        tracker.cylinder,
        tracker.evaluate,
        contact_tolerance_m=config["particles"]["contact_tolerance_m"],
    )
    times = [t0]
    positions = [full.copy()]
    steps = []
    rejects = []
    last_time = t0
    last_report = time.perf_counter()
    status = "RUNNING"
    error = None
    failure = None
    # Source state and fixed identities are saved even if the first trial fails.
    np.savez_compressed(
        out / "INPUT_SNAPSHOT.npz", time_s=t0, positions_m=full, radii_m=rfull, particle_id=ids, active=mask
    )

    def checkpoint():
        np.savez_compressed(
            out / "LAST_ACCEPTED_STATE.npz",
            time_s=last_time,
            positions_m=positions[-1],
            radii_m=rfull,
            particle_id=ids,
            active=mask,
        )

    checkpoint()
    with (
        (out / "ACCEPTED.jsonl").open("x") as accepted_log,
        (out / "REJECTED.jsonl").open("x") as rejected_log,
        ResourceMonitor() as monitor,
    ):

        def accepted(y, u, diag):
            nonlocal last_time, last_report
            last_time = diag["time_s"]
            full[mask] = y
            times.append(last_time)
            positions.append(full.copy())
            steps.append(diag)
            accepted_log.write(json.dumps(diag, allow_nan=False) + "\n")
            if time.perf_counter() - last_report > 5:
                accepted_log.flush()
                rejected_log.flush()
                checkpoint()
                progress = dict(
                    case=case,
                    run=out.name,
                    time_s=last_time,
                    elapsed_physical_s=last_time - t0,
                    target_s=duration,
                    accepted=len(steps),
                    rejected=len(rejects),
                    wall_s=time.perf_counter() - start,
                    last_dt_s=diag["dt_s"],
                )
                write(out / "PROGRESS.json", progress)
                print(json.dumps(progress), flush=True)
                last_report = time.perf_counter()
            if time.perf_counter() - start > 1800 or len(steps) > 100000:
                raise RuntimeError("PREREGISTERED_RUNTIME_BUDGET")

        def rejected(diag):
            rejects.append(diag)
            rejected_log.write(json.dumps(diag, allow_nan=False) + "\n")

        try:
            f = tracker.evaluate(x, r, t0)
            np.savez_compressed(out / "INITIAL_FORCES.npz", **f)
            for k in range(1, int(round(duration / dt)) + 1):
                end = min(t0 + k * dt, t0 + duration)
                x = stepper.interval(x, r, active, last_time, end, on_accept=accepted, on_reject=rejected)
            status = "COMPLETE"
        except (Exception, KeyboardInterrupt) as exc:
            status = (
                "STOPPED_BY_USER"
                if isinstance(exc, KeyboardInterrupt)
                else ("INCOMPLETE_RUNTIME" if "PREREGISTERED_RUNTIME_BUDGET" in str(exc) else "FAILED")
            )
            error = repr(exc)
            failure = getattr(exc, "diagnostics", None)
            (out / "ERROR.txt").write_text(traceback.format_exc())
        finally:
            checkpoint()
    with h5py.File(out / "trajectories.h5", "x") as h:
        h.attrs.update(
            complete=status == "COMPLETE",
            method=method,
            translation_backend=BACKEND,
            all_accepted_substeps_recorded=True,
            velocity_semantics="constrained interval-average motion; raw coordinates retained",
        )
        h.create_dataset("time_s", data=times)
        h.create_dataset("positions_m", data=positions, compression="gzip")
        h.create_dataset("active", data=np.tile(mask, (len(times), 1)), compression="gzip")
        h.create_dataset("particle_id", data=ids)
        h.create_dataset("radii_m", data=rfull)
    write(out / "RESOURCE_SAMPLES.json", monitor.rows)
    dts = np.array([s["dt_s"] for s in steps])
    all_trials = steps + [r["diagnostics"] for r in rejects if "predictor" in r["diagnostics"]]
    # Timings cover all successful QP attempts, including local-error rejections.
    timing = {
        key: sum(row.get(key, 0) for row in all_trials) for key in ["assembly_s", "wall_s", "forcing_s"]
    }
    for key in ["setup_s", "solve_s", "polish_s", "qp_seconds"]:
        timing["qp_" + key] = sum(
            row[s].get(key, 0) for row in all_trials for s in ["predictor", "corrector"]
        )
    wall_s = time.perf_counter() - start
    summary = dict(
        status=status,
        error=error,
        solver_failure=failure,
        case=case,
        backend=BACKEND,
        qp_solver=qp_solver,
        method=method,
        run=str(out),
        baseline_run=str(baseline),
        duration_s=duration,
        completed_duration_s=last_time - t0,
        macro_dt_s=dt,
        accepted_steps=len(steps),
        rejected_intervals=len(rejects),
        wall_s=wall_s,
        process_cpu_s=time.process_time() - process_start,
        cpu_utilization_mean_percent=100 * (time.process_time() - process_start) / wall_s,
        dt_quantiles_s=np.quantile(dts, [0, 0.05, 0.5, 0.95, 1]).tolist() if len(dts) else None,
        kkt_max=max((r["kkt"] for r in all_trials), default=None),
        linear_residual_max=max((r["linear_residual"] for r in all_trials), default=None),
        minimum_pair_gap_m=min(
            (r["minimum_pair_gap_m"] for r in steps if r["minimum_pair_gap_m"] is not None), default=None
        ),
        minimum_wall_gap_m=min(
            (r["minimum_wall_gap_m"] for r in steps if r["minimum_wall_gap_m"] is not None), default=None
        ),
        projection_calls=0,
        projection_corrected_calls=0,
        projection_max_m=0.0,
        projection_rms_m=0.0,
        multi_radius_corrections=0,
        timings=timing,
        tracker_timings=tracker.timings,
        peak_ram_bytes=max((s["rss_bytes"] for s in monitor.rows), default=0),
        peak_gpu_memory_MiB=max((s.get("gpu_memory_MiB", 0) for s in monitor.rows), default=0),
        gpu_utilization_mean_percent=float(
            np.mean([s["gpu_percent"] for s in monitor.rows if "gpu_percent" in s])
        )
        if monitor.rows
        else None,
        source_inputs_unchanged=all(sha(p) == s for p, s in inputs.items()),
        pair_backend=bench["selected_pair_backend"],
        physical_scope="Fixed active dense checkpoint cohort; no long-transport promotion",
        plan_sha256=sha(plan),
    )
    write(out / "SUMMARY.json", summary)
    if status != "COMPLETE":
        write(
            out / "FAILURE.json",
            dict(
                status=status,
                error=error,
                diagnostics=failure,
                last_accepted_time_s=last_time,
                identity_sha256=sha(out / "IDENTITY.json"),
            ),
        )
    print(json.dumps(summary), flush=True)
    return out, summary


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--case", required=True, choices=["B2", "B3"])
    p.add_argument("--duration", type=float, required=True, choices=[0.001, 0.005])
    p.add_argument("--dt", type=float, default=20e-6, choices=[20e-6, 10e-6])
    p.add_argument("--qp-solver", choices=["sparse_active_set", "osqp"], default="sparse_active_set")
    a = p.parse_args()
    run(a.case, a.duration, a.dt, a.qp_solver)


if __name__ == "__main__":
    main()
