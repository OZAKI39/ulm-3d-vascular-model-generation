"""Supervisor-managed finite stage runner; no automatic repetition of failed tasks."""
import json
import time
import traceback
from pathlib import Path

from mevo_mb.acoustics.kwave_backend import precompute
from mevo_mb.acoustics.validation import validate_superposition
from mevo_mb.io.artifacts import new_run, write_json
from mevo_mb.utils.config import load_config


def main():
    job = json.loads(Path("scripts/JOB.json").read_text())
    stage = job["stage"]
    start = time.perf_counter()
    state = {"stage": stage, "status": "RUNNING", "job": job}
    write_json(Path("results/WORKER_STATE.json"), state)
    try:
        c = load_config(job["config"])
        run = new_run(Path("results"), stage, c)
        state["run"] = str(run)
        write_json(Path("results/WORKER_STATE.json"), state)
        if stage.startswith("acoustic"):
            basis = precompute(c, run)
            check = new_run(Path("results"), stage+"_superposition", c)
            validation = validate_superposition(c, basis, check)
            state.update(basis=str(basis), validation_run=str(check), validation=validation)
        elif stage.startswith("grid"):
            from mevo_mb.acoustics.convergence import run_refined
            state.update(run_refined(c, run))
        else:
            raise ValueError(f"Unknown stage {stage}")
        state.update(status="COMPLETE", elapsed_s=time.perf_counter()-start)
    except BaseException as error:
        state.update(status="FAILED", error=repr(error), traceback=traceback.format_exc(), elapsed_s=time.perf_counter()-start)
        raise
    finally:
        write_json(Path("results/WORKER_STATE.json"), state)
        write_json(Path("results")/(stage+"_STATE.json"), state)


if __name__ == "__main__":
    main()
