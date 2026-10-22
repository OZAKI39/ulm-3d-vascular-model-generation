"""Bounded three-level all-active replay around first formal B2 ID divergence."""

from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import h5py
import numpy as np
from scripts.stage5a_r3.common import BASE, read, write, sha
from scripts.stage5a_r3.run import run, checkpoint
from mevo_mb.bubbles.constrained_motion import AcceptedCache


def job(kind):
    ledger = read(BASE / "final/CAMPAIGN_R3_V1/CLOSED_RESULT.json")["ledger"]
    entry = next(r for r in ledger if r["case"] == "B2" and r["kind"] == kind and r["duration_s"] == 0.005)
    source = Path(entry["run"]) / "trajectories.h5"
    with h5py.File(source) as h:
        times = h["time_s"][:]
        index = int(np.argmin(abs(times - 3.00003)))
        assert abs(times[index] - 3.00003) < 1e-14
        x, r, ids, active, t = (
            h["positions_m"][index],
            h["radii_m"][:],
            h["particle_id"][:],
            h["active"][index],
            float(times[index]),
        )
    fixture = BASE / "fixtures" / f"B2_early_{kind}_3p00003.npz"
    metadata = dict(
        source=str(source),
        source_sha256=sha(source),
        actual_index=index,
        time_s=t,
        cache="Explicit cold solver initialization at common actual state; no intermediate canonical cache was stored",
        scope="Full 186 active particles; diagnostic independent starts, not new original-3-s formal completion",
    )
    checkpoint(fixture, x, np.zeros_like(x), r, ids, active, t, metadata, AcceptedCache(), metadata)
    rows = []
    for level in range(3):
        out, summary = run(
            "B2",
            20e-6,
            max_dt=5e-6 / 2**level,
            tol_factor=0.25 / 4**level,
            pair="cupy",
            budget=300,
            label=f"early_precision_{kind}_L{level}",
            initial_dt=1e-9,
            components=False,
            nonlinear_strategy="quasi_newton",
            nonsmooth_recovery=True,
            fixture=fixture,
            s0_dt=50e-9 / 2**level,
            optimized=True,
        )
        rows.append(dict(kind=kind, level=level, run=str(out), summary=summary))
        write(BASE / "diagnosis" / f"EARLY_PRECISION_{kind}.json", rows)
    return rows


if __name__ == "__main__":
    with ProcessPoolExecutor(max_workers=2) as pool:
        rows = [row for batch in pool.map(job, ["nominal", "refined"]) for row in batch]
    write(
        BASE / "diagnosis/EARLY_PRECISION.json",
        dict(
            rows=rows,
            charged_experiment_wall_s=sum(r["summary"]["wall_s"] for r in rows),
            resource_scope="Two independent development workers, BLAS=1, no formal timing concurrently",
            hypothesis="Local truncation and contact/constraint transitions around first +30..40 us same-ID separation; no new global-accuracy claim",
        ),
    )
