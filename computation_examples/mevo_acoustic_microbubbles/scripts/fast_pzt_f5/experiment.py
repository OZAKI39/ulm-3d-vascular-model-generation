"""Bounded true-state F5 trajectories and morphology, never visualization motion."""

from time import perf_counter
import json
import h5py
import numpy as np
from mevo_mb.fast_sim.f5_engine import F5FlowSimulator
from mevo_mb.fast_sim.nearfield import F2Physics
from mevo_mb.fast_sim.storage import save_run
from mevo_mb.fast_sim.morphodynamics import morphology_series, topology
from scripts.fast_pzt_f2.common import quality_pass
from .common import write, identity


def run_case(out, cfg, flow, geometry, basis, response):
    out.mkdir(parents=True, exist_ok=False)
    begin = perf_counter()
    p = F2Physics()
    event = cfg["events"][0]
    acoustic = bool(
        np.any(
            np.array(event.get("amplitude_Pa", [event.get("pressure_Pa", 0)]))
            * np.array(event.get("enabled", [True] * len(event.get("amplitude_Pa", [0]))))
        )
    )
    sim = F5FlowSimulator(
        np.asarray(cfg["positions_m"]),
        p,
        response,
        basis,
        flow=flow,
        acoustic_enabled=acoustic,
        flow_enabled=event["flow_enabled"],
        vessel=geometry,
        vessel_enabled=True,
        execution_mode="stable_reduction",
    )
    run_id = dict(
        identity(),
        flow=flow.metadata,
        geometry=geometry.metadata,
        basis=None if basis is None else basis.metadata,
        response=response.metadata,
    )
    write(out / "CONFIG.json", cfg)
    write(out / "IDENTITY.json", run_id)
    result = sim.run(
        cfg["duration_s"], cfg["dt_s"], cfg["events"], cfg["sample_interval_s"], cfg.get("wall_budget_s", 240)
    )
    metrics, _ = save_run(out / "trajectory.h5", result, p, run_id, cfg)
    s = result["summary"]
    q = s["quality"]
    active = bool(result["active"].all())
    passed = (
        s["status"] == "COMPLETE"
        and quality_pass(q, p)
        and q["wall"]["minimum_queried_gap_m"] >= -0.01 * p.radius_m
        and q["wall"]["max_pass_correction_m"] <= 0.1 * p.radius_m
    )
    s.update(
        numerical_quality="PASS" if passed else "FAIL",
        fixed_ID_metrics_valid=active,
        configuration=cfg,
        minimum_saved_wall_gap_m=float(result["wall_gap_m"].min()),
        maximum_saved_pair_overlap_m=float(max(0, -metrics["minimum_gap_m"].min())),
    )
    with h5py.File(out / "trajectory.h5", "r+") as h:
        h.attrs.update(
            schema="fast_pzt_f5_v1",
            flow_physics="PRESCRIBED_STEADY_FLOW",
            wall_physics="GEOMETRIC_WALL_ONLY",
            numerical_quality=s["numerical_quality"],
            fixed_ID_metrics_valid=active,
        )
        h.attrs["model_json"] = json.dumps(s["model"])
        for key in ["active", "wall_gap_m", "wall_stats"]:
            h.create_dataset(key, data=result[key], compression="gzip")
        h.create_dataset(
            "quality_samples/raw_F1_diagonal", data=result["quality_samples"], compression="gzip"
        )
        for key, array in result["flow"].items():
            h.create_dataset("flow/" + key, data=array, compression="gzip")
        h["boundary_events_json"] = json.dumps(result["boundary_events"])
        if "flow_domain_failure" in result:
            write(out / "FLOW_DOMAIN_FAILURE.json", result["flow_domain_failure"])
        if not active:
            del h["metrics"]
            h.attrs["metrics_status"] = "WITHHELD_AFTER_EXIT; positions and actual active IDs retained"
    if active:
        x = result["positions_m"][:, 0]
        morph = morphology_series(x, np.arange(sim.n), metrics["labels"][:, 0], 0.1)
        np.savez_compressed(out / "MORPHOLOGY.npz", time_s=result["time_s"], **morph)
        top, _, _ = topology(
            x[0],
            result["time_s"],
            metrics["labels"][:, 0],
            result["edge_event_time_s"],
            result["edge_event_ids"],
            result["time_s"][-1],
            0.002,
        )
        write(out / "TOPOLOGY.json", top)
        s["topology_counts"] = top["counts"]
    s["process_body_wall_s"] = perf_counter() - begin
    write(out / "SUMMARY.json", s)
    print(
        json.dumps(
            dict(
                run=out.name,
                status=s["status"],
                quality=s["numerical_quality"],
                wall_s=s["wall_s"],
                largest=s["final_largest_size"],
            )
        ),
        flush=True,
    )
    return sim, result, metrics
