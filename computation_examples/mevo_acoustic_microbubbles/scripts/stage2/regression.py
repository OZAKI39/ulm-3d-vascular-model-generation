"""Retest cached Stage1 acoustics and deterministic migration without recomputing bases."""

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np

from mevo_mb.acoustics.basis import AcousticBasisDatabase, ComplexAcousticField
from mevo_mb.acoustics.convergence import compare_grid
from mevo_mb.bubbles.population import BubblePopulation
from mevo_mb.bubbles.tracker import ParticleTracker
from mevo_mb.coupling.experiment import flow_from_config
from mevo_mb.coupling.geometry import Cylinder
from mevo_mb.io.artifacts import new_run
from mevo_mb.utils.config import load_config


def main():
    root = Path.cwd()
    out = (
        root
        / "results/stage2/audit"
        / ("regression_" + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ"))
    )
    out.mkdir(parents=True)
    manifest = json.loads((root / "docs/DELIVERY_MANIFEST.json").read_text())
    missing = []
    changed = []
    count = 0
    for name, meta in manifest["files"].items():
        # Scientific baseline outputs are immutable; source changes are tracked
        # in git and baseline tag, rather than expected to have identical hashes.
        if not name.startswith("results/"):
            continue
        path = root / name
        count += 1
        if not path.exists():
            missing.append(name)
        elif hashlib.file_digest(path.open("rb"), "sha256").hexdigest() != meta["sha256"]:
            changed.append(name)
    basis_checks = []
    grids = {}
    for f, scenario in ((1, "aggregation"), (2, "migration")):
        state = json.loads((root / f"results/acoustic{f}_STATE.json").read_text())
        db = AcousticBasisDatabase.load_basis(state["basis"])
        v = state["validation"]
        meta = v["direct_metadata"]
        u = np.array(meta["source_control_real_Pa"]) + 1j * np.array(meta["source_control_imag_Pa"])
        predicted = db.combine(u).pressure_Pa
        with np.load(Path(state["validation_run"]) / "acoustic/direct/field.npz") as z:
            direct = z["P_complex"]
        c = load_config(f"configs/demo_{scenario}.yaml")
        mask = Cylinder(**c["vessel"]).mask(db.axes_m)
        error = float(np.linalg.norm((predicted - direct)[mask]) / np.linalg.norm(direct[mask]))
        basis_checks.append(
            dict(
                frequency_Hz=f * 1e6,
                relative_L2=error,
                baseline_relative_L2=v["relative_L2"],
                unchanged=bool(np.isclose(error, v["relative_L2"], rtol=1e-10, atol=1e-15)),
            )
        )
        grid = json.loads((root / f"results/grid{f}_STATE.json").read_text())
        run = new_run(out, f"grid{f}", c)
        grids[scenario] = compare_grid(c, Path(state["basis"]), Path(grid["refined_field"]), run)
    old = root / "results/20261001T140059_607305Z_migration_baseline"
    c = load_config(old / "config.yaml")
    field = ComplexAcousticField.load(old / "acoustic/combined_field.h5")
    pop = BubblePopulation.initialize(c)
    run = new_run(out, "migration_baseline_repeat", c)
    ParticleTracker(c, field, flow_from_config(c)).run(pop, run)
    with h5py.File(old / "particles/trajectories.h5") as a, h5py.File(run / "particles/trajectories.h5") as b:
        max_error = float(np.max(np.abs(a["positions_m"][:] - b["positions_m"][:])))
        same_times = bool(np.array_equal(a["time_s"][:], b["time_s"][:]))
    report = dict(
        status="PASS"
        if not missing
        and not changed
        and all(x["unchanged"] for x in basis_checks)
        and max_error < 1e-14
        and same_times
        else "FAIL",
        immutable_result_files=count,
        missing=missing,
        changed=changed,
        basis_direct=basis_checks,
        grid_diagnostics=grids,
        migration_max_coordinate_error_m=max_error,
        migration_same_times=same_times,
        migration_run=str(run.relative_to(root)),
        legacy_demo_only=True,
        note="All cached fields reused. Diagnostic grid differences remain, not evidence of grid independence.",
    )
    (out / "SUMMARY.json").write_text(json.dumps(report, indent=2) + "\n")
    (out.parent / "REGRESSION_STATE.json").write_text(
        json.dumps(dict(report, run=str(out.relative_to(root))), indent=2) + "\n"
    )
    print(report["status"], max_error, flush=True)


if __name__ == "__main__":
    main()
