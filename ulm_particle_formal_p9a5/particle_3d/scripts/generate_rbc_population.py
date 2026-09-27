#!/usr/bin/env python3
"""Generate only RBC validation geometry; explicit seed; no dynamics or injection."""
import argparse
import csv
import json
from pathlib import Path
import sys
import numpy as np

PACKAGE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PACKAGE / "src"))
from particle_3d.rbc_distribution import sample_rbc_geometries, write_population, statistics, quantile_indices, stratified_indices


def write_json(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n", type=int, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    pop = sample_rbc_geometries(args.n, args.seed)
    path = args.output / f"C57BL6_RBC_GEOMETRY_VALIDATION_{args.n}.csv"
    meta = write_population(pop, path)
    with (args.output / "candidate_ledger.csv").open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(pop.candidates.dtype.names)
        for row in pop.candidates:
            writer.writerow([int(row["candidate_id"]), format(row["D_raw_um"], ".17g"),
                             format(row["V_raw_fL"], ".17g"), row["status"]])
    summary = {name: statistics(pop.samples[name]) for name in ["D_um","V_fL","a_um","c_um","r","jeffery_lambda"]}
    summary["full_thickness_um"] = statistics(2*pop.samples["c_um"])
    summary["latent_D_V_assumption"] = "independent"
    summary["accepted_D_V_sample_correlation"] = float(np.corrcoef(pop.samples["D_um"],pop.samples["V_fL"])[0,1])
    write_json(args.output / "accepted_statistics.json", summary)
    selected = []
    for fraction, index in zip([.05,.25,.5,.75,.95], quantile_indices(pop.samples)):
        selected.append(dict(quantile=fraction, **{n: pop.samples[index][n].item() for n in pop.samples.dtype.names}))
    write_json(args.output / "selected_geometries.json", selected)
    if args.n >= 64:
        indices = stratified_indices(pop.samples, seed=2026092064)
        write_json(args.output / "orientation_sweep_selection.json",
                   dict(seed=2026092064, role="VALIDATION_GEOMETRY_SELECTION_ONLY", method="one random original sample from each of 64 equal-count r strata",
                        samples=[{n: pop.samples[i][n].item() for n in pop.samples.dtype.names} for i in indices]))
    print(json.dumps(meta, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
