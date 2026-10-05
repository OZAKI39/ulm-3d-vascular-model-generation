"""Run independent two-level particle comparisons on two bounded CPU workers."""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import copy
import json

import numpy as np

from mevo_mb.acoustics.basis import AcousticBasisDatabase
from mevo_mb.acoustics.convergence import load_refined
from mevo_mb.bubbles.population import BubblePopulation
from mevo_mb.control.pzt import PZTArray
from mevo_mb.coupling.campaign import find_verified_basis
from mevo_mb.coupling.experiment import run_particles
from mevo_mb.io.artifacts import new_run, write_json
from mevo_mb.utils.config import load_config


def execute(job):
    scenario, kind = job
    root = Path("results")
    c = load_config(f"configs/demo_{scenario}.yaml")
    basis = find_verified_basis(root, c)
    target = copy.deepcopy(c)
    field = None
    if kind == "dt_half":
        target["particles"]["dt_s"] /= 2
    else:
        frequency = 2 if scenario == "migration" else 1
        state = json.loads((root/f"grid{frequency}_STATE.json").read_text())
        assert state["status"] == "COMPLETE"
        scale = 1.
        if "aggregation" in c:
            coarse = AcousticBasisDatabase.load_basis(basis).combine(PZTArray.from_config(c).get_control_vector())
            cloud = BubblePopulation.initialize(c)
            scale = c["aggregation"]["target_initial_local_PNP_Pa"] / np.median(np.abs(coarse.interpolate_complex_pressure(cloud.positions_m)))
        field = load_refined(Path(state["refined_field"]), c["acoustics"]["frequency_Hz"], scale)
    run = new_run(root, scenario+"_"+kind, target)
    print("CONVERGENCE_PARTICLES", run, flush=True)
    summary = run_particles(target, basis, run, field_override=field)
    result = {"scenario": scenario, "kind": kind, "run": str(run), "summary": summary}
    write_json(root/f"{scenario}_{kind}_STATE.json", {"status": "COMPLETE", **result})
    return result


def main():
    jobs = [(s, k) for s in ["migration", "aggregation"] for k in ["dt_half", "refined_acoustic_particles"]]
    with ProcessPoolExecutor(max_workers=2) as executor:
        results = list(executor.map(execute, jobs))
    write_json(Path("results/CONVERGENCE_PARTICLES_STATE.json"), {"status": "COMPLETE", "runs": results})


if __name__ == "__main__":
    main()
