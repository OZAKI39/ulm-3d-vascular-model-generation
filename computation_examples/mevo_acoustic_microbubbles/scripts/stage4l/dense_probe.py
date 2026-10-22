"""Finite, uniquely identified Stage 4L checkpoint experiments. Gold is read-only."""

from __future__ import annotations

import argparse
import copy
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import time
import traceback
import zipfile

import h5py
import numpy as np

from mevo_mb.bubbles.translation_backends import (
    LEGACY_LABEL,
    LUBRICATION_PHYSICS,
    INTERFACE_SCOPE,
    NOT_INCLUDED,
)

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "results/stage4/lubrication_upgrade"


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False, default=lambda x: x.item()) + "\n")


def load_case(legacy, case):
    from mevo_mb.acoustics.basis import ComplexAcousticField
    from mevo_mb.flow.fields import VTKFlowField

    gold = legacy / "results/stage4/performance_refactor/gold_baseline"
    identity = gold / "cases" / case / "RUN_IDENTITY.json"
    config = copy.deepcopy(json.loads(identity.read_text())["config"])
    inputs = {}

    def frozen(path):
        rel = str(path).removeprefix("/workspace/mevo_acoustic_microbubbles/")
        p = gold / "inputs" / rel
        if not p.is_file():
            raise FileNotFoundError(p)
        inputs[str(p)] = sha(p)
        return str(p)

    config["geometry"]["mesh_surface"] = frozen(config["geometry"]["mesh_surface"])
    config["bubble_response"]["path"] = frozen(config["bubble_response"]["path"])
    config["effective_acoustic_field"]["path"] = frozen(config["effective_acoustic_field"]["path"])
    config["flow"]["snapshot"] = frozen(config["flow"]["snapshot"])
    benchmark = legacy / "results/stage4/performance_refactor/protocol/PAIR_BACKEND_STATE.json"
    config["particles"]["pair_benchmark_path"] = str(benchmark)
    inputs[str(benchmark)] = sha(benchmark)
    field = ComplexAcousticField.load(config["effective_acoustic_field"]["path"])
    f = config["flow"]
    flow = VTKFlowField(
        f["snapshot"],
        velocity_array=f["velocity_array"],
        coordinate_scale_to_m=f["coordinate_scale_to_m"],
        velocity_scale_to_m_s=f["velocity_scale_to_m_s"],
        spatial_method=f["spatial_method"],
    )
    return config, field, flow, gold, inputs


def run(legacy, case, backend, duration, max_wall_s=1800, max_steps=100000):
    from mevo_mb.bubbles.ordered_contact_sweep import OrderedContactSweep
    from mevo_mb.coupling.bifurcation import BranchClassifier
    from mevo_mb.coupling.integrators import GlobalSynchronousIntegrator

    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    out = BASE / "dense_probe" / f"{stamp}_{case}_{backend}_{duration * 1000:g}ms"
    out.mkdir(parents=True, exist_ok=False)
    start = time.perf_counter()
    config, field, flow, gold, inputs = load_case(legacy, case)
    config["particles"]["translation_backend"] = backend
    config["particles"]["lubrication"] = dict(activation=0.1, floor=0.001, residual_tolerance=1e-11)
    trajectory = gold / "cases" / case / "particles/trajectories.h5"
    with h5py.File(trajectory, "r") as h:
        frame = int(np.argmin(np.abs(h["time_s"][:] - 3.0)))
        t0 = float(h["time_s"][frame])
        xfull = h["positions_m"][frame]
        radii = h["radii_m"][:]
        ids = h["particle_id"][:]
        mask = h["active"][frame]
    if abs(t0 - 3) > 1e-9:
        raise ValueError("Missing 3 s checkpoint")
    np.savez_compressed(
        out / "INPUT_SNAPSHOT.npz", time_s=t0, positions_m=xfull, radii_m=radii, particle_id=ids, active=mask
    )
    (out / "INPUT_SNAPSHOT.npz").chmod(0o444)
    inputs[str(trajectory)] = sha(trajectory)
    inputs[str(gold / "GOLD_BASELINE_MANIFEST.json")] = sha(gold / "GOLD_BASELINE_MANIFEST.json")
    sources = {
        str(p.relative_to(ROOT)): sha(p)
        for folder in ["src", "scripts/stage4l"]
        for p in (ROOT / folder).rglob("*.py")
    }
    with zipfile.ZipFile(out / "source_snapshot.zip", "x", zipfile.ZIP_DEFLATED) as z:
        for p in sources:
            z.write(ROOT / p, p)
    write(
        out / "IDENTITY.json",
        dict(
            case=case,
            backend=backend,
            duration_s=duration,
            input_sha256=inputs,
            snapshot_sha256=sha(out / "INPUT_SNAPSHOT.npz"),
            source_sha256=sources,
            legacy_label=LEGACY_LABEL,
            git_commit=subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
            base_snapshot="aaf2549fc506933df74229be8e0ec9717f97f7ed",
            source_base="a6c3153f484113d37ffb5e9b4b57f0c6594ef396",
            python=platform.python_version(),
            config=config,
            initial_time_s=t0,
            finite_budget=dict(max_wall_s=max_wall_s, max_steps=max_steps),
            scope="Only translation resistance/contact metric/gap guard upgraded; identical acoustic and radial inputs",
        ),
    )
    write(
        out / "PHYSICS_BANNER.json",
        dict(
            translation_backend=backend,
            production_physics=[LUBRICATION_PHYSICS, INTERFACE_SCOPE]
            if backend == "lubrication_coupled"
            else [LEGACY_LABEL],
            not_included=NOT_INCLUDED,
            full_near_contact_bubble_validation=False,
        ),
    )
    tracker = GlobalSynchronousIntegrator(config, field, flow)
    # Exact original ordered projection algorithm, own cache and same backend for both arms.
    import mevo_mb.bubbles.tracker as tracking

    original = tracking.project_contacts
    sweep = OrderedContactSweep(BASE / "unit_tests/ordered_contact_cache")
    tracking.project_contacts = sweep.project
    write(out / "PROJECTION_RUNTIME.json", sweep.identity)
    active = ids[mask].copy()
    x = xfull[mask].copy()
    classifier = BranchClassifier(tracker.cylinder, ids)
    old_branches = json.loads((gold / "cases" / case / "BRANCH_METRICS.json").read_text())
    for event in old_branches["crossings"]:
        if event.get("crossing_time_s", np.inf) <= t0:
            classifier.assignment[int(event["particle_id"])] = event["branch"]
    t = t0
    final = t0 + duration
    dts = []
    speeds = []
    times = [t0]
    positions = [xfull.copy()]
    masks = [mask.copy()]
    crossings = []
    status = "RUNNING"
    error = None
    tracker.translation_log = (out / "TRANSLATION_SOLVES.jsonl").open("x")
    try:
        state = tracker.evaluate(x, radii[active], t)
        np.savez_compressed(out / "INITIAL_FORCES_AND_VELOCITIES.npz", **state)
        if backend == "lubrication_coupled":
            a = tracker.last_translation.matrix.toarray() / tracker.last_translation.scale
            eigen = np.linalg.eigvalsh(a)
            write(
                out / "INITIAL_MATRIX.json",
                dict(
                    size=len(a),
                    min_eigenvalue=float(eigen[0]),
                    max_eigenvalue=float(eigen[-1]),
                    condition_number_2=float(eigen[-1] / eigen[0]),
                    symmetry_max_abs=float(np.max(np.abs(a - a.T))),
                    **tracker.last_translation.diagnostics,
                ),
            )
        last_report = time.perf_counter()
        nominal = config["particles"]["dt_s"]
        with (out / "steps.jsonl").open("x") as log:
            while t < final - 1e-14 and len(active):
                if time.perf_counter() - start > max_wall_s or len(dts) >= max_steps:
                    status = "DEVELOPMENT_BUDGET_REACHED"
                    break
                target = min(t0 + (int(np.floor((t - t0 + 1e-14) / nominal)) + 1) * nominal, final)
                dt = tracker.select_timestep(x, radii[active], state, target - t)
                if dt <= 0 or t + dt == t:
                    raise RuntimeError("TIMESTEP_COLLAPSE")
                previous = x.copy()
                oldstate = state
                old_projection = tracker.projection_diagnostics.copy()
                x, wall, contacts, correction, iterations = tracker.step(x, radii[active], t, dt, state)
                tnew = min(t + dt, target)
                state = tracker.evaluate(x, radii[active], tnew)
                events = classifier.observe(
                    active, previous, x, oldstate["total_velocity_m_s"], state["total_velocity_m_s"], t, tnew
                )
                crossings.extend(events)
                xfull[active] = x
                leaving = {e["particle_id"] for e in events}
                for event in events:
                    xfull[event["particle_id"]] = event["position_m"]
                keep = np.array([int(pid) not in leaving for pid in active], bool)
                speeds.extend(np.linalg.norm(state["total_velocity_m_s"], axis=1).tolist())
                diag = tracker.last_translation.diagnostics.copy()
                if backend == "lubrication_coupled":
                    diag.update(tracker.translation_diagnostics[-1])
                diag.pop("time_s", None)
                log.write(
                    json.dumps(
                        dict(
                            time_s=tnew,
                            dt_s=dt,
                            active=len(active),
                            projection_calls=tracker.projection_diagnostics["calls"]
                            - old_projection["calls"],
                            projection_max_m=tracker.projection_diagnostics["max_correction_m"],
                            final_projection_m=correction,
                            projection_iterations=iterations,
                            secondary_force_max_N=float(
                                np.max(np.linalg.norm(state["secondary_force_N"], axis=1))
                            ),
                            velocity_max_m_s=float(
                                np.max(np.linalg.norm(state["total_velocity_m_s"], axis=1))
                            ),
                            **diag,
                        ),
                        default=lambda a: a.item(),
                    )
                    + "\n"
                )
                active = active[keep]
                x = x[keep]
                state = {k: v[keep] for k, v in state.items()}
                t = tnew
                if leaving and len(active):
                    state = tracker.evaluate(x, radii[active], t)
                tracker.monitor.record_accepted_interval(dt)
                mask = np.isin(ids, active)
                times.append(t)
                dts.append(dt)
                positions.append(xfull.copy())
                masks.append(mask)
                if time.perf_counter() - last_report > 15:
                    progress = dict(
                        case=case,
                        backend=backend,
                        time_s=t,
                        elapsed_probe_s=t - t0,
                        target_probe_s=duration,
                        steps=len(dts),
                        wall_s=time.perf_counter() - start,
                        latest_dt_s=dt,
                    )
                    write(out / "PROGRESS.json", progress)
                    print(json.dumps(progress), flush=True)
                    log.flush()
                    last_report = time.perf_counter()
            else:
                status = "COMPLETE"
    except Exception as exc:
        status = "FAILED"
        error = str(exc)
        (out / "ERROR.txt").write_text(traceback.format_exc())
        if hasattr(tracker, "_last_query"):
            np.savez_compressed(out / "LAST_QUERY.npz", **tracker._last_query)
    finally:
        tracker.translation_log.close()
        tracker.translation_log = None
        write(
            out / "TRANSLATION_DIAGNOSTICS.json",
            dict(**tracker.translation_summary, failure=getattr(tracker.translation, "last_failure", None)),
        )
        tracking.project_contacts = original
        with h5py.File(out / "trajectories.h5", "x") as h:
            h.attrs.update(
                complete=status == "COMPLETE",
                translation_backend=backend,
                checkpoint_time_s=t0,
                all_accepted_substeps_recorded=True,
            )
            h.create_dataset("time_s", data=times)
            h.create_dataset("positions_m", data=positions, compression="gzip")
            h.create_dataset("active", data=masks, compression="gzip")
            h.create_dataset("particle_id", data=ids)
            h.create_dataset("radii_m", data=radii)
        p = tracker.projection_diagnostics.copy()
        p["rms_correction_m"] = float(np.sqrt(p["squared_correction_m2"] / max(1, p["particle_queries"])))
        p["state"] = "PROJECTION_FALLBACK"
        summary = dict(
            status=status,
            case=case,
            backend=backend,
            requested_duration_s=duration,
            completed_duration_s=t - t0,
            accepted_steps=len(dts),
            wall_s=time.perf_counter() - start,
            error=error,
            min_dt_s=min(dts) if dts else None,
            median_dt_s=float(np.median(dts)) if dts else None,
            max_velocity_m_s=max(speeds) if speeds else None,
            median_velocity_m_s=float(np.median(speeds)) if speeds else None,
            saved_Reynolds_max=(
                2
                * float(np.max(radii))
                * config["medium"]["density_kg_m3"]
                * max(speeds)
                / config["medium"]["viscosity_Pa_s"]
            )
            if speeds
            else None,
            projection=p,
            physics_validity=tracker.monitor.summary(),
            crossings=crossings,
            source_inputs_unchanged=all(sha(path) == value for path, value in inputs.items()),
            run=str(out.relative_to(ROOT)),
        )
        write(out / "SUMMARY.json", summary)
        print(json.dumps(summary), flush=True)
    return out, summary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--legacy-root", type=Path, required=True)
    parser.add_argument("--case", choices=["B2", "B3"], required=True)
    parser.add_argument("--backend", choices=["legacy_stokes", "lubrication_coupled"], required=True)
    parser.add_argument("--duration", type=float, default=0.001)
    parser.add_argument("--max-wall-s", type=float, default=1800)
    args = parser.parse_args()
    if args.duration not in (0.001, 0.005):
        raise ValueError("Dense probe supports only 1 or 5 ms")
    _, summary = run(args.legacy_root, args.case, args.backend, args.duration, args.max_wall_s)
    if summary["status"] != "COMPLETE":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
