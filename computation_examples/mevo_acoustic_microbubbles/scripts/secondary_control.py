"""Matched secondary-force-off control, keeping initial cloud and source controls."""
from pathlib import Path

from mevo_mb.coupling.campaign import find_verified_basis
from mevo_mb.coupling.experiment import run_particles
from mevo_mb.io.artifacts import new_run, write_json
from mevo_mb.utils.config import load_config


def main():
    config = load_config("configs/demo_aggregation.yaml")
    config["particles"]["secondary_enabled"] = False
    root = Path("results")
    basis = find_verified_basis(root, config)
    run = new_run(root, "aggregation_secondary_off", config)
    summary = run_particles(config, basis, run)
    write_json(root/"SECONDARY_OFF_STATE.json", {"status": "COMPLETE", "run": str(run), "summary": summary})


if __name__ == "__main__":
    main()
