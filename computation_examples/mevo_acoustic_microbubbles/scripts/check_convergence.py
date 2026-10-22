"""Execute dt/2 particle runs and compare already computed dx/1.25 acoustic fields."""
from __future__ import annotations

import copy
import json
from pathlib import Path

import h5py
import numpy as np

from mevo_mb.acoustics.convergence import compare_grid, load_refined
from mevo_mb.coupling.experiment import run_particles
from mevo_mb.io.artifacts import new_run, write_json
from mevo_mb.utils.config import load_config


def compare_particles(coarse: Path, fine: Path) -> dict:
    with h5py.File(coarse/"particles/trajectories.h5") as a, h5py.File(fine/"particles/trajectories.h5") as b:
        if not np.array_equal(a["positions_m"][0], b["positions_m"][0]) or not np.array_equal(a["radii_m"][:], b["radii_m"][:]):
            raise ValueError("Convergence requires identical initial particles")
        delta = a["positions_m"][-1]-b["positions_m"][-1]
        rms = float(np.sqrt(np.mean(np.sum(delta**2, axis=1))))
        maximum = float(np.linalg.norm(delta, axis=1).max())
    a = json.loads((coarse/"particles/SUMMARY.json").read_text())
    b = json.loads((fine/"particles/SUMMARY.json").read_text())
    measures = {}
    for key in ["mean_displacement_m", "mean_lateral_displacement_m", "mean_nearest_neighbor_m", "median_nearest_neighbor_m", "largest_cluster_fraction", "near_contact_fraction"]:
        v0, v1 = a["final"][key], b["final"][key]
        measures[key] = {"coarse": v0, "fine": v1, "absolute_difference": abs(v0-v1),
                         "relative_difference": abs(v0-v1)/abs(v1) if v1 != 0 else None}
    return {"status": "MEASURED_TWO_LEVEL_DIAGNOSTIC", "coarse_run": str(coarse), "fine_run": str(fine),
            "endpoint_RMS_difference_m": rms, "endpoint_max_difference_m": maximum,
            "coarse_nominal_dt_s": a["nominal_dt_s"], "fine_nominal_dt_s": b["nominal_dt_s"],
            "coarse_substeps": a["accepted_substeps"], "fine_substeps": b["accepted_substeps"],
            "coarse_minimum_dt_s": a["minimum_accepted_dt_s"], "fine_minimum_dt_s": b["minimum_accepted_dt_s"],
            "metrics": measures, "limitation": "Geometry guard may also reduce steps; endpoint agreement is not validation of the radial or contact model"}


def main():
    root = Path("results")
    campaign = json.loads((root/"DEMO_CAMPAIGN.json").read_text())
    report = new_run(root, "convergence", load_config("configs/demo.yaml"))
    results = {}
    for scenario, frequency in [("migration", 2), ("aggregation", 1)]:
        base = next(x for x in campaign["scenarios"] if x["scenario"] == scenario and x["pattern"] == "baseline")
        coarse = Path(base["run"])
        c = load_config(coarse/"config.yaml")
        basis = Path(base["summary"]["basis"])
        refined = copy.deepcopy(c)
        refined["particles"]["dt_s"] /= 2
        cached_dt = root/f"{scenario}_dt_half_STATE.json"
        if cached_dt.exists():
            state = json.loads(cached_dt.read_text())
            assert state["status"] == "COMPLETE"
            target = Path(state["run"])
            assert load_config(target/"config.yaml") == refined
        else:
            target = new_run(root, scenario+"_dt_half", refined)
            print(f"DT_CONVERGENCE_RUN={target}", flush=True)
            run_particles(refined, basis, target)
        dt_result = compare_particles(coarse, target)
        results[scenario] = {"timestep": dt_result}
        write_json(report/"PARTIAL.json", results)
        grid_state = json.loads((root/f"grid{frequency}_STATE.json").read_text())
        assert grid_state["status"] == "COMPLETE"
        cached_field = root/f"{scenario}_grid_diagnostic_STATE.json"
        if cached_field.exists():
            state = json.loads(cached_field.read_text())
            assert state["status"] == "COMPLETE"
            grid_report = Path(state["run"])
            assert load_config(grid_report/"config.yaml") == c
            field_result = state["result"]
        else:
            grid_report = new_run(root, scenario+"_grid_comparison", c)
            field_result = compare_grid(c, basis, Path(grid_state["refined_field"]), grid_report)
        # Use identical physical source controls for particle grid comparison.
        control = json.loads((coarse/"acoustic/CONTROL.json").read_text())
        scale = control["aggregation_source_scale"]
        fine_field = load_refined(Path(grid_state["refined_field"]), c["acoustics"]["frequency_Hz"], scale)
        cached_grid = root/f"{scenario}_refined_acoustic_particles_STATE.json"
        if cached_grid.exists():
            state = json.loads(cached_grid.read_text())
            assert state["status"] == "COMPLETE"
            fine_particles = Path(state["run"])
            assert load_config(fine_particles/"config.yaml") == c
        else:
            fine_particles = new_run(root, scenario+"_refined_acoustic_particles", c)
            run_particles(c, basis, fine_particles, field_override=fine_field)
        transport_result = compare_particles(coarse, fine_particles)
        results[scenario].update(grid_field=field_result, grid_particle=transport_result, grid_report=str(grid_report))
        write_json(report/"PARTIAL.json", results)
    write_json(report/"CONVERGENCE.json", results)
    from mevo_mb.visualization.plots import save
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    for ax, metric, title in [(axs[0], "mean_displacement_m", "Mean displacement"), (axs[1], "mean_nearest_neighbor_m", "Mean nearest neighbor")]:
        labels, values = [], []
        for scenario in results:
            row = results[scenario]["timestep"]["metrics"][metric]
            labels += [scenario+" dt", scenario+" dt/2"]
            values += [row["coarse"]*1e6, row["fine"]*1e6]
        ax.bar(labels, values)
        ax.set(ylabel=title+" (µm)")
        ax.tick_params(axis="x", rotation=20)
    save(fig, report/"figures/timestep_convergence")
    write_json(report/"metadata.json", {"status": "ANALYSIS_COMPLETE", "results": results})
    write_json(root/"CONVERGENCE_STATE.json", {"status": "COMPLETE", "run": str(report), "results": results})


if __name__ == "__main__":
    main()
