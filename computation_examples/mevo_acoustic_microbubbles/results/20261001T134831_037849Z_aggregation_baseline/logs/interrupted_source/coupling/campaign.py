"""One command builds and validates two acoustic frequencies and runs both demos."""
from pathlib import Path
import json

from mevo_mb.acoustics.basis import AcousticBasisDatabase
from mevo_mb.acoustics.kwave_backend import acoustic_identity, precompute
from mevo_mb.acoustics.validation import validate_superposition
from mevo_mb.coupling.experiment import run_particles
from mevo_mb.io.artifacts import new_run, write_json
from mevo_mb.utils.config import load_config


def find_verified_basis(root: Path, config: dict) -> Path:
    expected = acoustic_identity(config)
    for state_path in sorted(root.glob("acoustic*_STATE.json")):
        state = json.loads(state_path.read_text())
        if state.get("status") != "COMPLETE" or state.get("validation", {}).get("status") != "PASS":
            continue
        path = Path(state["basis"])
        db = AcousticBasisDatabase.load_basis(path)
        if db.metadata["acoustic_identity"] == expected:
            return path
    raise FileNotFoundError("No matching closed and superposition-validated basis")


def run_demo(root: Path, reuse_acoustics: bool = False):
    root.mkdir(exist_ok=True, parents=True)
    result = {"status": "RUNNING", "scenarios": []}
    campaign = new_run(root, "demo_campaign", load_config("configs/demo.yaml"))
    write_json(campaign/"metadata.json", result)
    for scenario, stage, patterns in [("migration", "acoustic2", ["baseline", "shifted"]), ("aggregation", "acoustic1", ["baseline"])]:
        c = load_config(f"configs/demo_{scenario}.yaml")
        if reuse_acoustics:
            basis = find_verified_basis(root, c)
        else:
            run = new_run(root, stage, c)
            basis = precompute(c, run)
            check = new_run(root, stage+"_superposition", c)
            validation = validate_superposition(c, basis, check)
            write_json(root/(stage+"_STATE.json"), {"status": "COMPLETE", "basis": str(basis), "validation": validation, "validation_run": str(check), "run": str(run)})
        for pattern in patterns:
            run = new_run(root, scenario+"_"+pattern, c)
            print(f"PARTICLE_RUN={run}", flush=True)
            summary = run_particles(c, basis, run, pattern)
            result["scenarios"].append({"scenario": scenario, "pattern": pattern, "run": str(run), "summary": summary})
            write_json(campaign/"metadata.json", result)
    result["status"] = "COMPLETE"
    write_json(campaign/"metadata.json", result)
    write_json(root/"DEMO_CAMPAIGN.json", result)
    return result
