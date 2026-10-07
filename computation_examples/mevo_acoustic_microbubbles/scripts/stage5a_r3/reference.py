"""Rerun actual frozen Stage4L source in a separate process and new output root.

Only execution resources, pair device and exact 10 us output-clock truncation
are selected here. Old solver/evaluator/contact/projection methods are imported
from the frozen server tree, not copied or altered.
"""
# ruff: noqa: E402
# Imports after the explicit frozen-source path switch are intentional.

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
FROZEN = Path("/workspace/mevo_acoustic_microbubbles_stage4l")
sys.path[:0] = [str(FROZEN / "src"), str(FROZEN)]
import argparse
from datetime import datetime, timezone
import json
import time
import traceback
import h5py
import numpy as np
from mevo_mb.coupling.integrators import GlobalSynchronousIntegrator
from mevo_mb.bubbles.ordered_contact_sweep import OrderedContactSweep
import mevo_mb.bubbles.tracker as tracking
from scripts.stage4l.dense_probe import load_case, sha, write
from scripts.stage5a_r3.resources import ResourceMonitor
import signal


def run(case, duration, pair="cupy", repeat=0):
    assert str(Path(tracking.__file__).resolve()).startswith(str(FROZEN / "src"))
    base = ROOT / "results/stage5/motion_solver_r3"
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    out = base / "benchmarks/frozen_reference" / f"{stamp}_{case}_{duration:g}s_{pair}_repeat{repeat}"
    out.mkdir(parents=True)
    start = time.perf_counter()
    cpu = time.process_time()
    cfg, field, flow, gold, inputs = load_case(Path("/workspace/mevo_acoustic_microbubbles"), case)
    cfg["particles"].update(
        translation_backend="lubrication_coupled",
        pair_backend=pair,
        response_backend="numpy",
        lubrication=dict(activation=0.1, floor=0.001, residual_tolerance=1e-11),
    )
    selected = json.loads(
        (
            Path("/workspace/mevo_acoustic_microbubbles_stage5a_r2")
            / "results/stage5/motion_solver/audit_v2/BASELINE_SELECTION.json"
        ).read_text()
    )
    src = Path(next(r["run"] for r in selected if r["key"] == case + "_1ms")) / "INPUT_SNAPSHOT.npz"
    inputs[str(src)] = sha(src)
    with np.load(src) as f:
        state = {k: f[k] for k in f.files}
    full = state["positions_m"].copy()
    mask = state["active"]
    r = state["radii_m"]
    ids = state["particle_id"]
    t0 = float(state["time_s"])
    x = full[mask].copy()
    rr = r[mask]
    t = t0
    elapsed = 0.0
    tracker = GlobalSynchronousIntegrator(cfg, field, flow)
    ordered = OrderedContactSweep(out / "compiled_ordered_contact")
    tracking.project_contacts = ordered.project
    write(out / "PROJECTION_RUNTIME.json", ordered.identity)
    identity = dict(
        case=case,
        duration_s=duration,
        pair_backend=pair,
        repeat=repeat,
        config=cfg,
        input_sha256=inputs,
        frozen_source_root=str(FROZEN),
        source_sha256={str(p.relative_to(FROZEN)): sha(p) for p in (FROZEN / "src").rglob("*.py")},
        driver_sha256=sha(Path(__file__)),
        environment_sha256=sha(base / "environment/IDENTITY.json"),
        output_clock_s=1e-5,
        clock_alignment_cost_included=True,
        reference_lineage="fresh frozen Stage4L source, as separately audited; not claimed identical to archived reference",
        wall_budget_s=4500,
        threads=1,
    )
    write(out / "IDENTITY.json", identity)
    initialization = time.perf_counter() - start
    h = h5py.File(out / "trajectories.h5", "x")
    for key, value in [
        ("time_s", np.array([t])),
        ("positions_m", np.array([full])),
        ("active", np.array([mask])),
    ]:
        h.create_dataset(
            key,
            data=value,
            maxshape=(None, *value.shape[1:]),
            chunks=(1, *value.shape[1:]),
            compression="gzip",
        )
    h["particle_id"] = ids
    h["radii_m"] = r
    h.attrs["complete"] = False
    h.flush()
    dtvalues = []
    error = None
    status = "INCOMPLETE"
    last_report = time.perf_counter()
    io = 0.0
    tracker.translation_log = (out / "TRANSLATION_SOLVES.jsonl").open("x")
    monitor = ResourceMonitor().__enter__()

    def alarm(*args):
        raise TimeoutError("REFERENCE_NATIVE_WATCHDOG_4500S")

    prior = signal.signal(signal.SIGALRM, alarm)
    signal.setitimer(signal.ITIMER_REAL, 4500)
    try:
        f = tracker.evaluate(x, rr, t)
        with (out / "ACCEPTED.jsonl").open("x") as log:
            while elapsed < duration - 1e-14:
                if time.perf_counter() - start > 4500 or len(dtvalues) >= 100000:
                    raise TimeoutError("REFERENCE_BUDGET")
                end = min((np.floor((elapsed + 1e-14) / 1e-5) + 1) * 1e-5, duration)
                dt = tracker.select_timestep(x, rr, f, end - elapsed)
                if dt <= 0:
                    raise RuntimeError("TIMESTEP_COLLAPSE")
                x, wall, contacts, correction, iterations = tracker.step(x, rr, t, dt, f)
                elapsed = min(elapsed + dt, end)
                t = t0 + elapsed
                f = tracker.evaluate(x, rr, t)
                full[mask] = x
                dtvalues.append(dt)
                row = dict(
                    time_s=t,
                    elapsed_s=elapsed,
                    dt_s=dt,
                    wall_s=time.perf_counter() - start,
                    projection_max_m=tracker.projection_diagnostics["max_correction_m"],
                    final_projection_m=correction,
                    **tracker.last_translation.diagnostics,
                )
                log.write(json.dumps(row, allow_nan=False, default=lambda v: v.item()) + "\n")
                begin = time.perf_counter()
                n = len(h["time_s"])
                for key, value in [("time_s", t), ("positions_m", full), ("active", mask)]:
                    h[key].resize(n + 1, axis=0)
                    h[key][n] = value
                h.flush()
                log.flush()
                io += time.perf_counter() - begin
                if time.perf_counter() - last_report > 10:
                    progress = dict(
                        run=out.name,
                        elapsed_s=elapsed,
                        target_s=duration,
                        accepted=len(dtvalues),
                        wall_s=time.perf_counter() - start,
                    )
                    write(out / "PROGRESS.json", progress)
                    print(json.dumps(progress), flush=True)
                    last_report = time.perf_counter()
            status = "COMPLETE"
    except (Exception, KeyboardInterrupt) as exc:
        error = repr(exc)
        (out / "ERROR.txt").write_text(traceback.format_exc())
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, prior)
        h.attrs["complete"] = status == "COMPLETE"
        h.close()
        tracker.translation_log.close()
        monitor.__exit__()
        write(out / "RESOURCE_SAMPLES.json", monitor.rows)
    summary = dict(
        status=status,
        error=error,
        case=case,
        target_s=duration,
        completed_duration_s=elapsed,
        accepted_steps=len(dtvalues),
        wall_s=time.perf_counter() - start,
        initialization_s=initialization,
        process_cpu_s=time.process_time() - cpu,
        io_s=io,
        pair_backend=pair,
        repeat=repeat,
        run=str(out),
        actual_dt_quantiles_s=np.quantile(dtvalues, [0, 0.5, 1]).tolist() if dtvalues else None,
        translation=tracker.translation_summary,
        projection=tracker.projection_diagnostics,
        tracker_timings_inclusive_s=tracker.timings,
        monitoring_inclusive_s=monitor.overhead_s,
        input_unchanged=all(sha(p) == s for p, s in inputs.items()),
    )
    write(out / "SUMMARY.json", summary)
    print(json.dumps(summary), flush=True)
    return summary


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--case", required=True)
    p.add_argument("--duration", type=float, required=True)
    p.add_argument("--pair", default="cupy")
    p.add_argument("--repeat", type=int, default=0)
    a = p.parse_args()
    run(a.case, a.duration, a.pair, a.repeat)
