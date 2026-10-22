"""Resolve nominal-dt ambiguity with a smaller effective step and onset-window probe."""
from pathlib import Path
import json

import h5py
import numpy as np
import yaml

from mevo_mb.bubbles.metrics import aggregation_metrics
from mevo_mb.coupling.experiment import run_particles
from mevo_mb.io.artifacts import new_run, write_json
from mevo_mb.utils.config import load_config


def main():
    root = Path("results")
    campaign = json.loads((root/"DEMO_CAMPAIGN.json").read_text())
    checks = yaml.safe_load(Path("configs/numerical_checks.yaml").read_text())
    results = {}
    for scenario in ["migration", "aggregation"]:
        original = next(c for c in campaign["scenarios"] if c["scenario"] == scenario and c["pattern"] == "baseline")
        coarse = Path(original["run"])
        c = load_config(coarse/"config.yaml")
        c["particles"]["dt_s"] *= checks[scenario]["nominal_dt_factor"]
        c["particles"]["max_displacement_per_radius"] *= checks[scenario]["displacement_limit_factor"]
        if scenario == "aggregation":
            c["particles"]["duration_s"] = checks[scenario]["onset_horizon_s"]
        run = new_run(root, scenario+"_effective_step_check", c)
        summary = run_particles(c, Path(original["summary"]["basis"]), run)
        t = c["particles"]["duration_s"]
        with h5py.File(coarse/"particles/trajectories.h5") as a, h5py.File(run/"particles/trajectories.h5") as b:
            times = a["time_s"][:]
            i = np.argmin(abs(times-t))
            assert abs(times[i]-t) < 1e-14
            assert np.array_equal(a["positions_m"][0], b["positions_m"][0])
            x0, x1, x2, radii = a["positions_m"][0], a["positions_m"][i], b["positions_m"][-1], a["radii_m"][:]
        m1 = aggregation_metrics(x1, radii, c["metrics"])[0]
        m2 = summary["final"]
        results[scenario] = {"status": "MEASURED", "run": str(run), "coarse_run": str(coarse), "horizon_s": t,
            "coarse_steps_to_horizon": int(i), "refined_steps_to_horizon": summary["accepted_substeps"],
            "coarse_max_displacement_per_radius": c["particles"]["max_displacement_per_radius"]/checks[scenario]["displacement_limit_factor"],
            "refined_max_displacement_per_radius": c["particles"]["max_displacement_per_radius"],
            "endpoint_RMS_difference_m": float(np.sqrt(np.mean(np.sum((x1-x2)**2, axis=1)))),
            "coarse_mean_displacement_m": float(np.linalg.norm(x1-x0, axis=1).mean()),
            "refined_mean_displacement_m": m2["mean_displacement_m"],
            "coarse_mean_NN_m": m1["mean_nearest_neighbor_m"], "refined_mean_NN_m": m2["mean_nearest_neighbor_m"],
            "coarse_largest_cluster": m1["largest_cluster_size"], "refined_largest_cluster": m2["largest_cluster_size"],
            "limitation": "Two-level diagnostic; aggregation probe covers onset through 2 ms only, not the full packed-state evolution"}
        write_json(root/"EFFECTIVE_STEP_PARTIAL.json", results)
    write_json(root/"EFFECTIVE_STEP_STATE.json", {"status": "COMPLETE", "results": results})


if __name__ == "__main__":
    main()
