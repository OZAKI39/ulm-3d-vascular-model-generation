"""Publication exports from saved reference curves, without rerunning physics."""

import json
from pathlib import Path
import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scripts.stage3.common import watermark

s = json.loads(Path("results/stage3/history/STATE.json").read_text())
root = Path(s["run"])
out = root / "figures"
out.mkdir(exist_ok=True)
curves = {}
for name in ["H1", "H2", "H3"]:
    with h5py.File(root / (name + ".h5")) as h:
        curves[name] = {key: h[key][:] for key in ["time_s", "radius_m", "instantaneous_force_N"]}
        for model in ["Aprod", "Aresolved", "B", "C"]:
            curves[name][model] = {key: h[model][key][:] for key in h[model]}


def save(fig, name):
    watermark(fig, "REFERENCE SENSITIVITY; NO EXPERIMENTAL CALIBRATION")
    for suffix in ["png", "pdf"]:
        fig.savefig(out / (name + "." + suffix), dpi=240)
    plt.close(fig)


fig, axs = plt.subplots(2, 3, figsize=(13, 7.5), layout="constrained")
for col, (name, d) in enumerate(curves.items()):
    for row, end in enumerate([8e-6, 1e-4]):
        ax = axs[row, col]
        mask = d["time_s"] <= end + 1e-15
        for model, label in [("Aprod", "Production"), ("B", "Added mass"), ("C", "Added mass + history")]:
            ax.plot(d["time_s"][mask] * 1e6, d[model]["displacement_m"][mask] * 1e6, label=label, lw=1.1)
        ax.set(
            xlabel="Time (us)",
            ylabel="Signed displacement (um)",
            title=f"{name}: " + ("20 cycles" if row == 0 else "0.1 ms"),
        )
        ax.legend(fontsize=8)
save(fig, "HF1")
fig, ax = plt.subplots(figsize=(7, 4), layout="constrained")
for offset, window in [(-0.2, "short"), (0.2, "long")]:
    ax.bar(
        np.arange(3) + offset,
        [100 * r["total_production_vs_C_relative_error"] for r in s["rows"] if r["window"] == window],
        0.4,
        label=window,
    )
ax.set(
    xticks=np.arange(3),
    xticklabels=["H1", "H2", "H3"],
    ylabel="Production vs history-reference displacement error (%)",
)
ax.legend()
save(fig, "HF2")
d = curves["H2"]
mask = d["time_s"] <= 4e-6
fig, ax = plt.subplots(figsize=(8, 4), layout="constrained")
for name, v in [
    ("Acoustic", d["instantaneous_force_N"]),
    ("Stokes drag", d["C"]["quasisteady_drag_N"]),
    ("Added mass", d["C"]["added_mass_N"]),
    ("History", d["C"]["history_force_N"]),
]:
    ax.plot(d["time_s"][mask] * 1e6, v[mask] * 1e9, label=name, lw=1)
ax.set(xlabel="Time (us)", ylabel="Force (nN)", title="H2: matched instantaneous forces")
ax.legend()
save(fig, "HF3")
fig, ax = plt.subplots(figsize=(7, 4), layout="constrained")
for window in ["short", "long"]:
    rows = [r for r in s["rows"] if r["window"] == window]
    x = [max(r["max_R_over_R0"] - 1, 1 - r["min_R_over_R0"]) for r in rows]
    y = [100 * r["total_production_vs_C_relative_error"] for r in rows]
    ax.scatter(x, y, label=window)
    for xx, yy, r in zip(x, y, rows):
        ax.annotate(r["case"], (xx, yy), xytext=(5, 4), textcoords="offset points")
ax.set(
    xlabel="Maximum absolute radial excursion / R0",
    ylabel="Production vs history-reference displacement error (%)",
)
ax.legend()
save(fig, "HF4")
