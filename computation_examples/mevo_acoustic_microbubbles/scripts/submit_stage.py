"""Safely submit a new finite acoustic stage to this project's supervisor only."""
import argparse
import json
import subprocess
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["acoustic1", "acoustic2", "grid1", "grid2"])
    ap.add_argument("--config", required=True)
    args = ap.parse_args()
    result = subprocess.run(["supervisorctl", "status", "mevo_mb_worker"], capture_output=True, text=True, check=False)
    # supervisorctl returns 3 for normal EXITED/STOPPED. Text state is authoritative.
    if any(state in result.stdout for state in ["RUNNING", "STARTING", "STOPPING", "BACKOFF"]):
        raise RuntimeError("Previous stage still active")
    if (Path("results")/(args.stage+"_STATE.json")).exists():
        raise FileExistsError("Stage already has evidence; use a new isolated run, never overwrite its ledger")
    Path("scripts/JOB.json").write_text(json.dumps({"stage": args.stage, "config": args.config}))
    subprocess.run(["supervisorctl", "start", "mevo_mb_worker"], check=True)


if __name__ == "__main__":
    main()
