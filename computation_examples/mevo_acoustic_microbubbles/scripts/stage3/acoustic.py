"""One new medium/frequency identity; all subsequent controls reuse this basis."""

from pathlib import Path
from mevo_mb.acoustics.kwave_backend import precompute
from mevo_mb.acoustics.validation import validate_superposition
from mevo_mb.io.artifacts import new_run
from mevo_mb.utils.config import load_config
from scripts.stage3.common import atomic_json

c = load_config("configs/stage3/cylinder.yaml")
root = Path("results/stage3/acoustic")
run = new_run(root, "water_2p5MHz", c)
state = dict(status="RUNNING", run=str(run))
atomic_json(root / "STATE.json", state)
try:
    basis = precompute(c, run)
    check = new_run(root, "superposition", c)
    validation = validate_superposition(c, basis, check)
    state.update(status="COMPLETE", basis=str(basis), validation_run=str(check), validation=validation)
except Exception as e:
    state.update(status="FAILED", error=repr(e))
    raise
finally:
    atomic_json(root / "STATE.json", state)
