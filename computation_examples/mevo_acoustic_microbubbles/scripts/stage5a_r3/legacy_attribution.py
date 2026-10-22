"""Controlled short legacy replays: source, macro clock, pair device and dt refinement."""

from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from time import perf_counter
import signal
import h5py
import numpy as np
from scripts.stage5a_r3.common import BASE, load, read, write, identity
from scripts.stage5a_r3.cross_replay import legacy_module


def job(task):
    case, variant, clock, pair, cap, label = task
    begin = perf_counter()
    tracker, _, cfg, state, _ = load(case, pair)
    lineage = read(BASE / "diagnosis" / f"{case}_LINEAGE_V1.json")
    module, path = legacy_module(variant, Path(lineage["paths"]["archived"]).parent)
    tracker.translation = module.LubricationCoupledTranslation(
        cfg["medium"]["viscosity_Pa_s"], **cfg["particles"]["lubrication"]
    )
    tracker.motion_constrained = False
    tracker.lubrication_enabled = True
    import mevo_mb.bubbles.tracker as tracking
    from mevo_mb.bubbles.ordered_contact_sweep import OrderedContactSweep

    out = BASE / "diagnosis/legacy_attribution_v1" / (case + "_" + label)
    out.mkdir(parents=True)
    ordered = OrderedContactSweep(out / "ordered_cache")
    original = tracking.project_contacts
    tracking.project_contacts = ordered.project
    mask, ids, rr = state["active"], state["particle_id"], state["radii_m"]
    full = state["positions_m"].copy()
    x, r = full[mask].copy(), rr[mask]
    t0 = float(state["time_s"])
    elapsed, t = 0.0, t0
    duration = 40e-6 if cap is None else 10e-6
    h = h5py.File(out / "trajectories.h5", "x")
    for key, val in [
        ("time_s", np.array([t])),
        ("positions_m", np.array([full])),
        ("active", np.array([mask])),
    ]:
        h.create_dataset(
            key, data=val, maxshape=(None, *val.shape[1:]), chunks=(1, *val.shape[1:]), compression="gzip"
        )
    h["particle_id"], h["radii_m"] = ids, rr
    h.attrs["complete"] = False
    rows = []
    error = None

    def alarm(*args):
        raise TimeoutError("DEVELOPMENT_300S_BUDGET")

    prior = signal.signal(signal.SIGALRM, alarm)
    signal.setitimer(signal.ITIMER_REAL, 300)
    try:
        f = tracker.evaluate(x, r, t)
        while elapsed < duration - 1e-14:
            if perf_counter() - begin > 300 or len(rows) > 100000:
                raise TimeoutError("DEVELOPMENT_300S_BUDGET")
            end = min((np.floor((elapsed + 1e-14) / clock) + 1) * clock, duration)
            dt = tracker.select_timestep(x, r, f, end - elapsed)
            if cap is not None:
                dt = min(dt, cap)
            if dt <= 0 or t + dt == t:
                raise RuntimeError("TIMESTEP_COLLAPSE")
            x, wall, contacts, correction, iterations = tracker.step(x, r, t, dt, f)
            elapsed = min(elapsed + dt, end)
            t = t0 + elapsed
            f = tracker.evaluate(x, r, t)
            full[mask] = x
            index = len(h["time_s"])
            for key, val in [("time_s", t), ("positions_m", full), ("active", mask)]:
                h[key].resize(index + 1, axis=0)
                h[key][index] = val
            rows.append(
                dict(
                    time_s=t,
                    dt_s=dt,
                    projection_m=correction,
                    kkt=tracker.translation_summary["max_contact_kkt_residual"],
                )
            )
        status = "COMPLETE"
    except Exception as exc:
        status, error = "INCOMPLETE", repr(exc)
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, prior)
        tracking.project_contacts = original
        h.attrs["complete"] = error is None
        h.close()
    report = dict(
        case=case,
        variant=variant,
        clock_s=clock,
        pair_setting=pair,
        max_dt_s=cap,
        label=label,
        target_s=duration,
        completed_duration_s=elapsed,
        status=status,
        error=error,
        source=str(path),
        wall_s=perf_counter() - begin,
        steps=len(rows),
        projection=tracker.projection_diagnostics,
        translation=tracker.translation_summary,
        scope="Controlled legacy-branch adapter; full frozen forcing, all active IDs, relative elapsed clock in every arm",
    )
    write(out / "STEPS.json", rows)
    write(out / "RESULT.json", report)
    print(report, flush=True)
    return report


def main():
    tasks = []
    for case in ["B2", "B3"]:
        tasks += [
            (case, "archived", 5e-4, "auto", None, "archived_500us_auto"),
            (case, "fresh", 5e-4, "auto", None, "fresh_500us_auto"),
            (case, "fresh", 1e-5, "auto", None, "fresh_10us_auto"),
            (case, "fresh", 1e-5, "cupy", None, "fresh_10us_cupy"),
        ]
        tasks += [
            (case, "fresh", 1e-5, "numpy", 1e-7 / 2**level, f"precision_L{level}") for level in range(3)
        ]
    with ProcessPoolExecutor(max_workers=2) as pool:
        rows = list(pool.map(job, tasks))
    write(
        BASE / "diagnosis/LEGACY_ATTRIBUTION_V1.json",
        dict(rows=rows, identity=identity(), charged_experiment_wall_s=sum(r["wall_s"] for r in rows)),
    )


if __name__ == "__main__":
    main()
