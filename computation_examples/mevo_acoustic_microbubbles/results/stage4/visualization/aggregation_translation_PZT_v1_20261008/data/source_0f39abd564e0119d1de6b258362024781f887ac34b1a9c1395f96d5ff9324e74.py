"""Bounded real-PZT Secondary-OFF counterfactual for the archived F4 parent cloud."""

import argparse
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--engine-root", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    root = args.engine_root.resolve()
    sys.path[:0] = [str(root / "src"), str(root)]
    import numpy as np
    import h5py
    from mevo_mb.fast_sim.f4_engine import F4VesselSimulator, GPUVessel
    from mevo_mb.fast_sim.nearfield import F2Physics
    from mevo_mb.fast_sim.models import ResponseSection, Basis
    from mevo_mb.fast_sim.storage import save_run
    from scripts.fast_pzt_f2.common import quality_pass, write

    base = root / "results/fast_pzt/f4"
    source = base / "aggregation_v1/parent_wall_on"
    cfg = json.loads((source / "CONFIG.json").read_text())
    assert (root / "SOURCE_F4_COMMIT.txt").read_text().strip() == "e247fa0c61f668a9b0b2e26c909dd687eedde389"
    p = F2Physics()
    response = ResponseSection(base / "inputs/response.h5", p)
    basis = Basis(base / "inputs/basis.h5", p.frequency_Hz)
    vessel = GPUVessel(base / "inputs/prepared_geometry.npz", base / "inputs/labelled_lumen_m.vtp")

    class SecondaryToggle(F4VesselSimulator):
        def pairs(self):
            if self.p.secondary:
                return super().pairs()

    # F3 intentionally validates the frozen baseline during construction. This
    # standalone, explicitly labelled ablation changes only the Secondary toggle;
    # immutable F1/F2/F3/F4 engine files are not edited.
    sim = SecondaryToggle(
        np.array(cfg["positions_m"]), p, response, basis, vessel=vessel, execution_mode="stable_reduction"
    )
    sim.p = replace(p, secondary=False)
    sim.control(cfg["events"][0])
    sim.fields()
    primary = sim.force.numpy().copy()
    sim.pairs()
    np.testing.assert_array_equal(primary, sim.force.numpy())
    assert np.linalg.norm(primary) > 0
    out = args.output
    out.mkdir(parents=True, exist_ok=False)
    source_sha = hashlib.sha256((source / "trajectory.h5").read_bytes()).hexdigest()
    identity = dict(
        f4_archive="e247fa0c61f668a9b0b2e26c909dd687eedde389",
        matched_secondary_on=str(source / "trajectory.h5"),
        matched_secondary_on_sha256=source_sha,
        runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        engine_sha256={
            f.name: hashlib.sha256(f.read_bytes()).hexdigest()
            for f in sorted((root / "src/mevo_mb/fast_sim").glob("*.py"))
        },
        basis=basis.metadata,
        response=response.metadata,
        geometry=vessel.metadata,
    )
    cfg = dict(
        cfg,
        secondary_force_enabled=False,
        physics_override={"secondary": False},
        scenario="MATCHED_SECONDARY_FORCE_ABLATION; actual PZT and original F4 WALL remain enabled",
    )
    write(out / "CONFIG.json", cfg)
    write(out / "IDENTITY.json", identity)
    result = sim.run(cfg["duration_s"], cfg["dt_s"], cfg["events"], cfg["sample_interval_s"], 180)
    result["summary"]["numerical_semantics"] = (
        "F4 stable baseline with ONLY Secondary pair-force addition disabled; Primary, local mobility, pair contact and WALL retained"
    )
    metrics, _ = save_run(out / "trajectory.h5", result, sim.p, identity, cfg)
    with h5py.File(out / "trajectory.h5", "r+") as h:
        for key in ["active", "wall_gap_m", "wall_stats"]:
            h.create_dataset(key, data=result[key], compression="gzip")
        h["boundary_events_json"] = json.dumps(result["boundary_events"])
        h.create_dataset(
            "quality_samples/raw_F1_diagonal", data=result["quality_samples"], compression="gzip"
        )
        h.attrs.update(schema="f4_visual_secondary_ablation_v1", secondary=False)
    with h5py.File(source / "trajectory.h5") as h:
        np.testing.assert_array_equal(h["positions_m"][0], result["positions_m"][0])
        np.testing.assert_array_equal(h["time_s"][:], result["time_s"])
    s = result["summary"]
    q = s["quality"]
    passed = (
        s["status"] == "COMPLETE"
        and quality_pass(q, sim.p)
        and q["wall"]["minimum_queried_gap_m"] >= -0.01 * p.radius_m
        and q["wall"]["max_pass_correction_m"] <= 0.1 * p.radius_m
        and result["active"].all()
    )
    s.update(
        numerical_quality="PASS" if passed else "FAIL",
        primary_retained_nonzero=True,
        secondary_addition_exactly_zero=True,
        matched_initial_positions_bitwise=True,
        matched_saved_times_bitwise=True,
        configuration=cfg,
        conclusion="No aggregation in this matched finite observation"
        if metrics["largest_size"][-1, 0] == 1
        else "Observed membership change; report actual result",
    )
    write(out / "SUMMARY.json", s)
    print(
        json.dumps(
            dict(
                status=s["status"],
                quality=s["numerical_quality"],
                largest=int(metrics["largest_size"][-1, 0]),
                rg_um=float(metrics["global_rg_m"][-1, 0] * 1e6),
                wall_s=s["wall_s"],
            )
        ),
        flush=True,
    )
    assert passed


if __name__ == "__main__":
    main()
