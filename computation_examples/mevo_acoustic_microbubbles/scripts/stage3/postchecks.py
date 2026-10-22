"""Accepted cylinder: full brute contact regression, timestep and grid diagnostics."""

import copy
import json
from pathlib import Path
import h5py
import numpy as np
from mevo_mb.acoustics.basis import ComplexAcousticField
from mevo_mb.acoustics.convergence import load_refined
from mevo_mb.bubbles.adaptive_table import AdaptiveResponseTableBubbleModel
from mevo_mb.bubbles.forces import primary_bjerknes
from mevo_mb.bubbles.population import BubblePopulation
from mevo_mb.bubbles.tracker import ParticleTracker
from mevo_mb.coupling.experiment import flow_from_config
from mevo_mb.io.artifacts import new_run as particle_run
from mevo_mb.utils.config import load_config
from scripts.stage3.common import new_run, finish

s = json.loads(Path("results/stage3/cylinder/STATE.json").read_text())
assert s["status"] == "PASS"
old = Path(next(r for r in s["rows"] if r["scenario"] == "C3")["run"])
c = load_config(old / "config.yaml")
field = ComplexAcousticField.load(old / "acoustic/combined_field.h5")
out = new_run(
    "postchecks",
    dict(contact_trajectory_absolute_tolerance_m=1e-14, early_timestep_relative_limit=0.02),
    ["Exact candidate regression and temporal refinement"],
)
with h5py.File(old / "particles/trajectories.h5") as h:
    xo = h["positions_m"][:]
    to = h["time_s"][:]
results = {}
for name, backend, dt in [
    ("brute", "brute", c["particles"]["dt_s"]),
    ("fine", "ckdtree", c["particles"]["dt_s"] / 2),
]:
    cc = copy.deepcopy(c)
    cc["particles"].update(contact_search_backend=backend, dt_s=dt)
    run = particle_run(out, name, cc)
    ParticleTracker(cc, field, flow_from_config(cc)).run(BubblePopulation.initialize(cc), run)
    with h5py.File(run / "particles/trajectories.h5") as h:
        x = h["positions_m"][:]
        t = h["time_s"][:]
    if name == "brute":
        results["contact_regression"] = dict(
            run=str(run),
            same_times=bool(np.array_equal(t, to)),
            bitwise_equal_positions=bool(np.array_equal(x, xo)),
            maximum_coordinate_error_m=float(np.max(abs(x - xo))),
        )
    else:
        i = np.argmin(abs(t - 1e-4))
        j = np.argmin(abs(to - 1e-4))
        results["timestep"] = dict(
            run=str(run),
            early_time_s=float(t[i]),
            early_displacement_relative_difference=float(
                np.linalg.norm(x[i] - xo[j]) / np.linalg.norm(x[i] - x[0])
            ),
            final_RMS_position_difference_m=float(np.sqrt(np.mean(np.sum((x[-1] - xo[-1]) ** 2, axis=1)))),
        )
# Evaluate nonlinear table on both acoustic grids, at identical initial particles.
grid = json.loads(Path("results/stage3/acoustic_grid/STATE.json").read_text())
fine = load_refined(Path(grid["refined_field"]), field.frequency_Hz)
model = AdaptiveResponseTableBubbleModel(c["bubble_response"]["path"])
pop = BubblePopulation.initialize(c)
values = []
for f in [field, fine]:
    p = f.interpolate_complex_pressure(pop.positions_m)
    v = model.response(pop.radii_m, f.frequency_Hz, p, c["medium"])
    values.append(primary_bjerknes(v, f.interpolate_complex_gradient(pop.positions_m)))
results["acoustic_primary_grid_relative_L2"] = float(
    np.linalg.norm(values[0] - values[1]) / np.linalg.norm(values[1])
)
results["grid_scope"] = "Two grids at initial cloud positions only; does not establish grid independence"
passed = (
    results["contact_regression"]["same_times"]
    and results["contact_regression"]["maximum_coordinate_error_m"] < 1e-14
    and results["timestep"]["early_displacement_relative_difference"] < 0.02
)
finish(out, dict(status="PASS" if passed else "FAIL", **results))
