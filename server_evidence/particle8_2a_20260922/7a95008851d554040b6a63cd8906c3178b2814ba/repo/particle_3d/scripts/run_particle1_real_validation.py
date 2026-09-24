#!/usr/bin/env python3
"""Fixed-before-run three-dt real-field validation, no timestep tuning."""
from pathlib import Path
import argparse
import csv
import json
import os
import sys
import time

PACKAGE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PACKAGE / "src"))
from particle_3d.audit import read_frozen
from particle_3d.field import FrozenFEMField
from particle_3d.particle1_audit import VALIDATION_SEED, scope_record
from particle_3d.sonovue_adapter import sample_single_validation_size
from particle_3d.particle1_cases import select_inlet_centroid, real_trajectory
from particle_3d.validation_boundary import ValidationBoundaryClassifier


def write_json(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n")


def write_csv(path, rows):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=PACKAGE / "reports/particle1/data")
    parser.add_argument("--fem-root", type=Path, default=PACKAGE.parent / "formal_3D_flow_solver/FEM_SimVascular")
    parser.add_argument("--sonovue-root", type=Path, default=Path(os.environ.get("SONOVUE_ROOT", "/home/lzy/projects/sonovue_size_distribution_v0")))
    args = parser.parse_args(); args.output.mkdir(parents=True, exist_ok=True)
    audit, mesh, flow, boundaries = read_frozen(args.fem_root)
    write_json(args.output / "00_particle1_scope_and_provenance.json", scope_record(PACKAGE.parent, args.sonovue_root, audit))
    field = FrozenFEMField.from_grids(mesh, flow)
    size = sample_single_validation_size(args.sonovue_root, seed=VALIDATION_SEED)
    write_json(args.output / "01_single_mb_size_provenance.json", size.to_dict())
    # Mandatory adjacent CSV-bound metadata for the sole DEMO sample.
    write_csv(args.output / "01_single_mb_sample.csv", [size.to_dict()])
    from particle_3d.audit import sha256
    write_json(args.output / "01_single_mb_sample.metadata.json", dict(**size.to_dict(), population_file="01_single_mb_sample.csv",
        population_sha256=sha256(args.output / "01_single_mb_sample.csv"), schema="SONOVUE_SINGLE_MB_ADAPTER_PROVENANCE_V0"))
    initial, candidates = select_inlet_centroid(field, boundaries["INLET"])
    write_json(args.output / "05_real_initialization.json", initial)
    write_csv(args.output / "05_inlet_initialization_candidates.csv", candidates)
    print(json.dumps(dict(size=size.to_dict(), initialization=initial)), flush=True)
    classifier = ValidationBoundaryClassifier(boundaries)
    summaries = []
    for index, dt in enumerate(initial["validation_timesteps_s"]):
        start = time.perf_counter()
        rows, summary = real_trajectory(field, classifier, initial, size.radius_m, dt,
                                       progress=lambda step, t: print(f"dt case {index}: step={step}, t={t:.6g}s", flush=True))
        summary["elapsed_wall_time_s"] = time.perf_counter() - start
        summary["case_id"] = index
        write_csv(args.output / f"05_real_trajectory_dt{index}.csv", rows)
        write_json(args.output / f"05_real_trajectory_dt{index}.json", summary)
        summaries.append(summary)
        print(json.dumps(summary), flush=True)
        # Failure evidence is retained. No rerouting/reflection/retuned retry.
    write_json(args.output / "05_real_trajectory_sweep.json", summaries)
    write_csv(args.output / "07_validation_timestep_comparison.csv", [
        {k: s[k] for k in ["case_id", "validation_dt_s", "timestep_role", "exit_boundary", "exit_time_s", "trajectory_length_m",
                          "max_trial_step_length_m", "row_count", "wall_crossing", "passed"]} for s in summaries])
    if not all(s["passed"] for s in summaries): raise SystemExit(1)


if __name__ == "__main__": main()
