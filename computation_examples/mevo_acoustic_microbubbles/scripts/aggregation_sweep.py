"""Declared pressure/spacing/exposure sensitivity; no force-law modification."""
import copy
from pathlib import Path

import pandas as pd

from mevo_mb.coupling.campaign import find_verified_basis
from mevo_mb.coupling.experiment import run_particles
from mevo_mb.io.artifacts import new_run, write_json
from mevo_mb.utils.config import load_config


def main():
    c = load_config("configs/demo_aggregation.yaml")
    root = Path("results")
    basis = find_verified_basis(root, c)
    rows = []
    variants = [("pressure", p) for p in c["aggregation"]["pressure_sweep_Pa"]]
    variants += [("spacing", radius) for radius in c["aggregation"]["cloud_radius_sweep_m"]]
    for kind, value in variants:
        config = copy.deepcopy(c)
        if kind == "pressure":
            config["aggregation"]["target_initial_local_PNP_Pa"] = value
        else:
            config["particles"]["initialization"]["cloud_radius_m"] = value
        run = new_run(root, f"aggregation_{kind}_sweep", config)
        run_particles(config, basis, run)
        metrics = pd.read_csv(run/"metrics.csv")
        for exposure in c["aggregation"]["duration_sweep_s"]:
            row = metrics.iloc[(metrics.time_s-exposure).abs().argmin()].to_dict()
            rows.append({"kind": kind, "value": value, "run": str(run), "requested_exposure_s": exposure, **row})
    pd.DataFrame(rows).to_csv(root/"AGGREGATION_SWEEP.csv", index=False)
    write_json(root/"AGGREGATION_SWEEP_STATE.json", {"status": "COMPLETE", "rows": rows})


if __name__ == "__main__":
    main()
