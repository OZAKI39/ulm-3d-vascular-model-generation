"""One continuous initial -> Secondary ON -> Secondary OFF F4 trajectory."""

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

    class SequentialToggle(F4VesselSimulator):
        def pairs(self):
            if self.p.secondary:
                return super().pairs()

        def control(self, event):
            enabled = bool(event["secondary_force_enabled"])
            if enabled != self.p.secondary:
                before = self.x.numpy().copy()
                cached = len(self.graphs)
                self.p = replace(self.p, secondary=enabled)
                # Captured graphs contain the actual Python-selected kernel set.
                # ON graphs must never be reused for OFF evolution.
                self.graphs.clear()
                np.testing.assert_array_equal(before, self.x.numpy())
                self.switch_records.append(
                    dict(
                        time_s=event["time_s"],
                        secondary=enabled,
                        discarded_graphs=cached,
                        positions_unchanged=True,
                    )
                )
            record = super().control(event)
            if not enabled:
                self.fields()
                primary = self.force.numpy().copy()
                self.pairs()
                np.testing.assert_array_equal(primary, self.force.numpy())
                assert np.linalg.norm(primary) > 0
                self.switch_records[-1].update(primary_nonzero=True, secondary_addition_exactly_zero=True)
            return record

    sim = SequentialToggle(
        np.array(cfg["positions_m"]), p, response, basis, vessel=vessel, execution_mode="stable_reduction"
    )
    sim.switch_records = []
    event = cfg["events"][0]
    cfg.update(
        duration_s=0.6,
        sample_interval_s=0.0005,
        events=[
            dict(event, time_s=0.0, secondary_force_enabled=True),
            dict(event, time_s=0.3, secondary_force_enabled=False),
        ],
        scenario="CONTINUOUS_SECONDARY_ON_THEN_OFF; same PZT, same particles, no state reset at switch",
    )
    out = args.output
    out.mkdir(parents=True, exist_ok=False)

    def digest(f):
        return hashlib.sha256(Path(f).read_bytes()).hexdigest()

    identity = dict(
        f4_archive="e247fa0c61f668a9b0b2e26c909dd687eedde389",
        reference_on_sha256=digest(source / "trajectory.h5"),
        runner_sha256=digest(__file__),
        engine_sha256={f.name: digest(f) for f in sorted((root / "src/mevo_mb/fast_sim").glob("*.py"))},
        basis=basis.metadata,
        response=response.metadata,
        geometry=vessel.metadata,
    )
    write(out / "CONFIG.json", cfg)
    write(out / "IDENTITY.json", identity)
    result = sim.run(cfg["duration_s"], cfg["dt_s"], cfg["events"], cfg["sample_interval_s"], 180)
    s = result["summary"]
    s["numerical_semantics"] = (
        "Uninterrupted F4 evolution; Secondary force switched OFF at 300 ms, cached CUDA graphs discarded; all other physics retained"
    )
    metrics, _ = save_run(out / "trajectory.h5", result, p, identity, cfg)
    with h5py.File(out / "trajectory.h5", "r+") as h:
        for key in ["active", "wall_gap_m", "wall_stats"]:
            h.create_dataset(key, data=result[key], compression="gzip")
        h["boundary_events_json"] = json.dumps(result["boundary_events"])
        h["secondary_force_enabled"] = h["control_index"][:] == 0
        h.attrs.update(
            schema="f4_sequential_secondary_toggle_v1",
            secondary="time-dependent; see secondary_force_enabled",
        )
        h.attrs["model_json"] = json.dumps(
            dict(p.metadata(), secondary="time-dependent", secondary_schedule=cfg["events"])
        )
    common = np.arange(301)
    with h5py.File(source / "trajectory.h5") as h:
        # Common 1 ms samples through the switch must reproduce the archived ON run.
        np.testing.assert_array_equal(result["positions_m"][2 * common], h["positions_m"][:1501:5])
    assert len(sim.switch_records) == 1
    assert sim.switch_records[0]["discarded_graphs"] > 0
    q = s["quality"]
    passed = (
        s["status"] == "COMPLETE"
        and quality_pass(q, p)
        and q["wall"]["minimum_queried_gap_m"] >= -0.01 * p.radius_m
        and q["wall"]["max_pass_correction_m"] <= 0.1 * p.radius_m
        and result["active"].all()
    )
    s.update(
        numerical_quality="PASS" if passed else "FAIL",
        switch_records=sim.switch_records,
        on_reference_common_frames_bitwise_equal=301,
        configuration=cfg,
        switch_largest=int(metrics["largest_size"][600, 0]),
        final_largest=int(metrics["largest_size"][-1, 0]),
        switch_rg_um=float(metrics["global_rg_m"][600, 0] * 1e6),
        final_rg_um=float(metrics["global_rg_m"][-1, 0] * 1e6),
    )
    s["model"]["secondary"] = "time-dependent; ON before 300 ms, OFF from 300 ms"
    write(out / "SUMMARY.json", s)
    print(
        json.dumps(
            {
                k: s[k]
                for k in [
                    "status",
                    "numerical_quality",
                    "wall_s",
                    "switch_largest",
                    "final_largest",
                    "switch_rg_um",
                    "final_rg_um",
                ]
            }
        ),
        flush=True,
    )
    assert passed


if __name__ == "__main__":
    main()
