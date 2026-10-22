"""Two-grid acoustic diagnostic using the actual nonlinear table at demo points."""

import json
from pathlib import Path
import numpy as np
from mevo_mb.acoustics.basis import AcousticBasisDatabase
from mevo_mb.acoustics.convergence import run_refined, load_refined
from mevo_mb.bubbles.population import BubblePopulation
from mevo_mb.io.artifacts import new_run
from mevo_mb.control.pzt import PZTArray
from mevo_mb.utils.config import load_config
from scripts.stage3.common import atomic_json

c = load_config("configs/stage3/cylinder.yaml")
run = new_run(Path("results/stage3/acoustic_grid"), "refinement", c)
fine = run_refined(c, run)
a = AcousticBasisDatabase.load_basis(
    json.loads(Path("results/stage3/acoustic/STATE.json").read_text())["basis"]
).combine(PZTArray.from_config(c).get_control_vector())
b = load_refined(Path(fine["refined_field"]), a.frequency_Hz)
x = BubblePopulation.initialize(c).positions_m
result = dict(
    status="TWO_GRID_DIAGNOSTIC_NOT_GRID_INDEPENDENCE", run=str(run), refined_field=fine["refined_field"]
)
for key, method in [
    ("pressure", "interpolate_complex_pressure"),
    ("gradient", "interpolate_complex_gradient"),
]:
    av = getattr(a, method)(x)
    bv = getattr(b, method)(x)
    result[key + "_relative_L2"] = float(np.linalg.norm(av - bv) / np.linalg.norm(bv))
result["initial_pressure_Pa"] = abs(a.interpolate_complex_pressure(x)).tolist()
atomic_json(run / "GRID_DIAGNOSTIC.json", result)
atomic_json(Path("results/stage3/acoustic_grid/STATE.json"), result)
print(result)
