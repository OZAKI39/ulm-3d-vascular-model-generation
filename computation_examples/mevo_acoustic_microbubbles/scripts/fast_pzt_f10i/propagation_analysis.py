"""Contrast the one same-dt tolerance control, without changing old admission."""

from scripts.fast_pzt_f10h.analyze import compare_paths
from .common import ROOT, BASE, TRIALS, js, write, sha, identity, gz


def main():
    r = js(BASE / "tolerance_v1/RESULTS.json")["rows"][0]
    assert r["status"] == "COMPLETE" and r["accepted_steps"] == 200
    path = ROOT / r["path"]
    a = compare_paths(TRIALS["DT10"], path)
    b = compare_paths(TRIALS["DT5"], path)
    steps = gz(path / "SOFT_STEPS.json.gz")
    assert all(
        s["accepted"] and s["result"]["converged"] and s["result"]["mobility_equation_residual_m"] <= 1.01e-14
        for s in steps
    )
    cfg = js(path / "CONFIG.json")
    old = js(TRIALS["DT10"] / "CONFIG.json")
    for key in ["initial_positions_m", "events"]:
        assert cfg[key] == old[key], key
    evidence = dict(
        source=identity(),
        same_dt_tolerance_comparison=a,
        tighter_DT10_vs_original_DT5=b,
        trajectory=r,
        max_equation_residual_m=max(s["result"]["mobility_equation_residual_m"] for s in steps),
        scopes="One original-initial-state 2ms control at same dt and stiffness; compare only actual common accepted points. 2.5us NOT_RUN.",
    )
    write(BASE / "analysis/TOLERANCE_PROPAGATION.json", evidence)
    admitted = (
        all(
            a["gates"][k]
            for k in [
                "centroid_m",
                "fixed_ID_m",
                "Kabsch_m",
                "Rg_relative",
                "aspect_relative",
                "tensor_relative",
                "largest_members",
                "all_members",
            ]
        )
        and a["matched_edge_max_time_difference_s"] <= 1e-5 + 1e-14
    )
    write(
        BASE / "analysis/GPU_START_GATE.json",
        dict(
            admitted=admitted,
            source_sha256=sha(BASE / "analysis/TOLERANCE_PROPAGATION.json"),
            reason="Archived 10/5us trajectories and single same-dt tolerance control separate tiny iterative effects from trajectory timestep differences; the slow near-threshold excursion magnifies small spatial shifts. Original strict edge FAIL retained. Only fixed-state GPU feasibility allowed.",
            allowed="Two/four followed conditionally by one saved32 state; no new physics trajectory, no B64",
            training_ready=False,
        ),
    )
    print("Same dt tolerance:", a["maxima"], a["matched_edge_max_time_difference_s"], a["event_counts"])
    print(
        "Tight10 vs original5:",
        b["matched_edge_max_time_difference_s"],
        b["event_counts"],
        "GPU fixed-state admitted",
        admitted,
    )


if __name__ == "__main__":
    main()
