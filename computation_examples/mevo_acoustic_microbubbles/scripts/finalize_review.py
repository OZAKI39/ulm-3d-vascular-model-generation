"""Supervisor-owned finite completion gate and closed-data report pipeline."""
from pathlib import Path
import json
import subprocess
import sys
import time
import traceback

from mevo_mb.io.artifacts import write_json


def main():
    start = time.monotonic()
    state_path = Path("results/REVIEW_STATE.json")
    write_json(state_path, {"status": "WAITING_FOR_PARTICLE_CONVERGENCE"})
    try:
        ready = Path("results/CONVERGENCE_PARTICLES_STATE.json")
        while not ready.exists():
            if time.monotonic()-start > 600:
                raise TimeoutError("Finite report gate exceeded 10 minutes")
            time.sleep(5)
        assert json.loads(ready.read_text())["status"] == "COMPLETE"
        for script in ["check_convergence.py", "review_artifacts.py", "build_report.py"]:
            write_json(state_path, {"status": "RUNNING", "script": script})
            print("REVIEW_STAGE", script, flush=True)
            subprocess.run([sys.executable, "scripts/"+script], check=True)
        with Path("docs/RUFF_RESULTS.txt").open("w") as output:
            subprocess.run([sys.executable, "-m", "ruff", "check", "src", "tests", "scripts"], check=True, stdout=output, stderr=subprocess.STDOUT)
        write_json(state_path, {"status": "COMPLETE", "elapsed_s": time.monotonic()-start})
    except BaseException:
        write_json(state_path, {"status": "FAILED", "traceback": traceback.format_exc()})
        raise


if __name__ == "__main__":
    main()
