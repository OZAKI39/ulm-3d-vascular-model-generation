"""One immutable, bounded F4 acoustic trajectory and geometry observations."""

from dataclasses import asdict
from time import perf_counter
import json
import h5py
import numpy as np
from mevo_mb.fast_sim.f4_engine import F4VesselSimulator, GPUVessel
from mevo_mb.fast_sim.models import ResponseSection, Basis
from mevo_mb.fast_sim.nearfield import F2Physics
from mevo_mb.fast_sim.storage import save_run
from mevo_mb.fast_sim.morphodynamics import morphology_series, topology
from scripts.fast_pzt_f2.common import quality_pass
from .common import BASE, write, identity


def run_case(out, cfg, geometry=None):
    out.mkdir(parents=True, exist_ok=False)
    started = perf_counter()
    p = F2Physics()
    response = ResponseSection(BASE / "inputs/response.h5", p)
    basis = Basis(BASE / "inputs/basis.h5", p.frequency_Hz) if cfg.get("basis", True) else None
    if geometry is None:
        geometry = GPUVessel(BASE / "inputs/prepared_geometry.npz", BASE / "inputs/labelled_lumen_m.vtp")
    run_id = dict(
        identity(),
        response=response.metadata,
        basis=None if basis is None else basis.metadata,
        geometry=geometry.metadata,
    )
    write(out / "CONFIG.json", cfg)
    write(out / "IDENTITY.json", run_id)
    wall = cfg.get("wall_enabled", True)
    sim = F4VesselSimulator(
        np.asarray(cfg["positions_m"]),
        p,
        response,
        basis,
        vessel=geometry,
        vessel_enabled=wall,
        execution_mode=cfg.get("execution_mode", "stable_reduction"),
    )
    result = sim.run(
        cfg["duration_s"],
        cfg.get("dt_s", 1e-5),
        cfg["events"],
        cfg.get("sample_interval_s", 0.0002),
        cfg.get("wall_budget_s", 180),
    )
    if not wall:
        begin = perf_counter()
        x = result["positions_m"]
        dist = geometry.query(x.reshape(-1, 3))["distance_m"]
        result["wall_gap_m"] = dist.reshape(x.shape[:-1]) - p.radius_m
        result["active"] = np.ones(x.shape[:-1], np.int32)
        result["wall_stats"] = np.zeros((*x.shape[:-1], 12))
        result["boundary_events"] = []
        result["summary"]["wall_off_observation_query_s"] = perf_counter() - begin
        result["summary"]["vessel_enabled"] = False
    # Scientific morphology comparison is defined on retained actual IDs. Refuse
    # silent inclusion of exited particles in the inherited fixed-N metric schema.
    all_active = bool(result["active"].all())
    metrics, _ = save_run(out / "trajectory.h5", result, p, run_id, cfg)
    s = result["summary"]
    q = s["quality"]
    passed = s["status"] == "COMPLETE" and quality_pass(q, p)
    if wall:
        passed = (
            passed
            and q["wall"]["minimum_queried_gap_m"] >= -0.01 * p.radius_m
            and q["wall"]["max_pass_correction_m"] <= 0.1 * p.radius_m
        )
    s.update(
        numerical_quality="PASS" if passed else "FAIL",
        physics_parameters=asdict(p),
        wall_physics="GEOMETRIC_WALL_ONLY",
        wall_contact_observed=bool(
            np.any(result["wall_gap_m"][result["active"].astype(bool)] <= 0.01 * p.radius_m)
        ),
        fixed_ID_metrics_valid=all_active,
        configuration=cfg,
    )
    with h5py.File(out / "trajectory.h5", "r+") as h:
        h.attrs.update(
            schema="fast_pzt_f4_v1",
            wall_physics="GEOMETRIC_WALL_ONLY",
            numerical_quality=s["numerical_quality"],
            fixed_ID_metrics_valid=all_active,
        )
        h.create_dataset("active", data=result["active"], compression="gzip")
        h.create_dataset("wall_gap_m", data=result["wall_gap_m"], compression="gzip")
        h.create_dataset("wall_stats", data=result["wall_stats"], compression="gzip")
        h.create_dataset(
            "quality_samples/raw_F1_diagonal", data=result["quality_samples"], compression="gzip"
        )
        h.create_dataset(
            "boundary_events_json", data=json.dumps(result["boundary_events"]), dtype=h5py.string_dtype()
        )
        if not all_active:
            # Positions remain untouched; fixed-N metrics cannot be interpreted as
            # active-population metrics after an exit and are explicitly withheld.
            del h["metrics"]
            h.attrs["metrics_status"] = "WITHHELD_AFTER_BOUNDARY_EXIT; use active mask"
    if all_active:
        x = result["positions_m"][:, 0]
        m = morphology_series(x, np.arange(len(x[0])), metrics["labels"][:, 0], 0.1)
        np.savez_compressed(out / "MORPHOLOGY.npz", time_s=result["time_s"], **m)
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
    s["process_body_wall_s"] = perf_counter() - started
    s["minimum_saved_wall_gap_m"] = float(result["wall_gap_m"].min())
    write(out / "SUMMARY.json", s)
    print(
        json.dumps(
            dict(
                run=out.name,
                status=s["status"],
                quality=s["numerical_quality"],
                wall_s=s["wall_s"],
                largest=s["final_largest_size"],
                gap_um=s["minimum_saved_wall_gap_m"] * 1e6,
            )
        ),
        flush=True,
    )
    return s
