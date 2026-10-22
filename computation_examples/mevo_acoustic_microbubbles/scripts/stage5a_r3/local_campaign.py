"""Bounded, full-population development refinement; not formal holdout evidence."""

import cProfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import h5py
import numpy as np
from scripts.stage5a_r3.common import BASE, R2, B2_RUNS, read, write, sha
from scripts.stage5a_r3.run import run, checkpoint
from mevo_mb.bubbles.constrained_motion import AcceptedCache


def actual_fixture(case, at):
    if case == "B2":
        source = R2 / "runs" / B2_RUNS["nominal"]
    else:
        ledger = read(R2 / "final/CAMPAIGN_V4A/RESULT.json")["ledger"]
        source = Path(
            next(
                x["run"]
                for x in ledger
                if x["case"] == case
                and x["duration"] == 0.005
                and x["kind"] == "nominal"
                and x.get("repeat") == 0
            )
        )
    path = BASE / "fixtures" / f"{case}_actual_{at:.8f}.npz"
    path.parent.mkdir(exist_ok=True)
    with h5py.File(source / "trajectories.h5") as h:
        index = int(np.searchsorted(h["time_s"][:], at, side="right") - 1)
        x, r, ids, mask, t = (
            h["positions_m"][index],
            h["radii_m"][:],
            h["particle_id"][:],
            h["active"][index],
            float(h["time_s"][index]),
        )
    metadata = dict(
        cache_provenance="EXPLICIT_COLD_START; original intermediate accepted cache not stored",
        source=str(source / "trajectories.h5"),
        source_sha256=sha(source / "trajectories.h5"),
        actual_index=index,
        requested_time_s=at,
        original_time_s=t,
    )
    checkpoint(path, x, np.zeros_like(x), r, ids, mask, t, metadata, AcceptedCache(), metadata)
    return str(path)


def job(name):
    case = "B3" if name.startswith("B3") else "B2"
    terminal = BASE / "diagnosis/reproduction_v1/nominal/terminal.npz"
    if name == "B2_event":
        fixture = BASE / "diagnosis/reproduction_v1/nominal/before_persistent_sweep.npz"
    elif name == "B2_S0":
        fixture = terminal
    elif name == "B2_stable":
        fixture = actual_fixture(case, 3.003)
    elif name == "B3_late":
        fixture = actual_fixture(case, 3.0049)
    else:
        fixture = None
    rows = []
    for level in range(3):
        method = "S0" if name == "B2_S0" else "S1"
        duration = 400e-9 if method == "S0" else 20e-6
        cap = (50e-9 if method == "S0" else 5e-6) / 2**level
        options = dict(
            case=case,
            duration=duration,
            method=method,
            max_dt=cap,
            tol_factor=0.25 / 2**level,
            pair="numpy",
            budget=300,
            initial_dt=min(50e-9, cap),
            nonlinear_strategy="quasi_newton",
            nonsmooth_recovery=method == "S1",
            fixture=fixture,
            s0_dt=50e-9 / 2**level,
            label=f"local_v1_{name}_L{level}",
        )
        profiler = cProfile.Profile() if level == 0 else None
        if profiler:
            profiler.enable()
        out, summary = run(**options)
        if profiler:
            profiler.disable()
            profiler.dump_stats(str(out / "DYNAMIC_PROFILE.pstats"))
        rows.append(dict(name=name, level=level, run=str(out), summary=summary))
        write(BASE / "diagnosis" / f"LOCAL_REFINEMENT_{name}.json", rows)
        if summary["status"] != "COMPLETE":
            break
    return rows


def main():
    names = ["B2_event", "B3_late", "B2_S0", "B2_stable", "B2_initial", "B3_initial"]
    with ProcessPoolExecutor(max_workers=2) as pool:
        rows = [row for batch in pool.map(job, names) for row in batch]
    write(
        BASE / "diagnosis/LOCAL_REFINEMENT_V1.json",
        dict(
            rows=rows,
            resource_scope="Two independent CPU workers, BLAS=1; development timing only",
            charged_experiment_wall_s=sum(row["summary"]["wall_s"] for row in rows),
        ),
    )


if __name__ == "__main__":
    main()
