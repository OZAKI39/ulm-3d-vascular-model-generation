"""One bounded F2 trajectory, immutable source/data identity and quality sidecars."""

import argparse
import json
from time import perf_counter
from dataclasses import asdict
import h5py
import numpy as np
from mevo_mb.fast_sim.nearfield import F2Physics
from mevo_mb.fast_sim.models import ResponseSection, Basis, random_cloud
from mevo_mb.fast_sim.f2_engine import F2Simulator
from mevo_mb.fast_sim.storage import save_run
from .common import BASE, ROOT, PROTOCOL, source_identity, write, quality_pass


def run_case(out, cfg):
    out.mkdir(parents=True, exist_ok=False)
    begin = perf_counter()
    p = F2Physics(**cfg.get("physics", {}))
    table = ResponseSection(BASE / "inputs/response.h5", p)
    basis = Basis(BASE / "inputs/basis.h5", p.frequency_Hz) if cfg.get("basis") else None
    x = (
        np.asarray(cfg["positions_m"])
        if "positions_m" in cfg
        else np.array(
            [
                random_cloud(
                    cfg["count"], p.radius_m, cfg["cloud_radius_m"], seed, cfg.get("center_m", [0, 0, 0])
                )
                for seed in cfg.get("seeds", [cfg["seed"]])
            ]
        )
    )
    if x.ndim == 2:
        x = x[None]
    identity = dict(
        source_sha256=source_identity(),
        execution_commit=(ROOT / "SOURCE_F2_COMMIT.txt").read_text().strip()
        if (ROOT / "SOURCE_F2_COMMIT.txt").exists()
        else None,
        f1_archive=PROTOCOL["f1_archive"],
        response=table.metadata,
        basis=None if basis is None else basis.metadata,
    )
    write(out / "CONFIG.json", cfg)
    write(out / "IDENTITY.json", identity)
    sim = F2Simulator(x, p, table, basis, particle_block_dim=cfg.get("particle_block_dim", 32))
    result = sim.run(
        cfg["duration_s"],
        cfg.get("dt_s", PROTOCOL["dt_s"]),
        cfg["events"],
        cfg.get("sample_interval_s", 0.0002),
        cfg.get("wall_budget_s", PROTOCOL["wall_budget_s"]),
        substeps=1,
    )
    metrics, _ = save_run(out / "trajectory.h5", result, p, identity, cfg)
    s = result["summary"]
    s["sample_numerical_quality"] = [
        "PASS" if s["status"] == "COMPLETE" and quality_pass(q, p) else "FAIL"
        for q in s["quality_per_sample"]
    ]
    s["numerical_quality"] = "PASS" if all(x == "PASS" for x in s["sample_numerical_quality"]) else "FAIL"
    with h5py.File(out / "trajectory.h5", "r+") as h:
        h.attrs["schema"] = "fast_pzt_f2_v1"
        h.attrs["numerical_quality"] = s["numerical_quality"]
        h.attrs["sample_numerical_quality_json"] = json.dumps(s["sample_numerical_quality"])
        h.create_dataset(
            "quality_samples/raw_F1_diagonal", data=result["quality_samples"], compression="gzip"
        )
        h["quality_samples"].attrs["columns"] = (
            "pressure_min,pressure_max,max_predicted_overlap,max_accepted_overlap,cumulative_correction,cumulative_mobility_path,max_pass_correction,max_speed,reserved,max_substep_displacement,max_balance_velocity"
        )
    s.update(
        configuration=cfg,
        physics_parameters=asdict(p),
        formation_time_s=[],
        edge_events_recorded=len(result["edge_event_time_s"]),
    )
    for b in range(len(x)):
        idx = np.flatnonzero(metrics["largest_size"][:, b] >= 16)
        s["formation_time_s"].append(float(result["time_s"][idx[0]]) if len(idx) else None)
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
                overlap_nm=s["quality"]["max_accepted_penetration_m"] * 1e9,
            )
        ),
        flush=True,
    )
    return s


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    ap.add_argument("--config", required=True)
    a = ap.parse_args()
    run_case(BASE / a.label, json.loads((ROOT / a.config).read_text()))
