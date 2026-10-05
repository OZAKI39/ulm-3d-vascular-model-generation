"""Export unsmoothed statistical comparisons; no edits to scientific trajectories."""

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from scripts.stage4l.dense_probe import ROOT, write, sha


def local_path(value):
    p = Path(value)
    return (
        p if p.exists() else ROOT / str(value).removeprefix("/workspace/mevo_acoustic_microbubbles_stage4l/")
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("comparison", type=Path)
    args = parser.parse_args()
    out = args.comparison / "figures"
    out.mkdir(exist_ok=True)
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8,
            "axes.labelsize": 8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "svg.fonttype": "none",
        }
    )
    fig, axes = plt.subplots(2, 3, figsize=(7.3, 4.2), layout="constrained")
    for row, case in enumerate(["B2", "B3"]):
        report = json.loads((args.comparison / (case + ".json")).read_text())
        for mode, color, label in [
            ("legacy", "#6e6e6e", "Legacy Stokes"),
            ("lubrication", "#483d8b", "Normal lubrication"),
        ]:
            folder = local_path(report[mode])
            artifact = local_path(
                report["metrics"]["10us"][mode].get("metric_artifact", folder / "cluster_metrics_10us.npz")
            )
            with np.load(artifact) as d:
                elapsed = (d["time_s"] - d["time_s"][0]) * 1000
                axes[row, 0].plot(elapsed, d["rigid_fit_residual_m"] * 1e6, color=color, lw=1.1, label=label)
                axes[row, 1].plot(elapsed[1:], d["turnover_fraction"] * 100, color=color, lw=0.9)
                axes[row, 2].plot(
                    elapsed, d["radius_of_gyration_m"] / d["radius_of_gyration_m"][0], color=color, lw=1.1
                )
        for col in range(3):
            axes[row, col].set_xlabel("Time from dense checkpoint (ms)")
            axes[row, col].grid(color="#dfdfdf", lw=0.45)
            axes[row, col].set_title(
                f"{case}  " + ["Internal rearrangement", "Neighbor turnover", "Cluster envelope"][col],
                loc="left",
                fontsize=9,
            )
        axes[row, 0].set_ylabel("Rigid-fit residual (µm)")
        axes[row, 1].set_ylabel("Changed edges / union (%)")
        axes[row, 2].set_ylabel(r"$R_g(t)/R_g(0)$")
    axes[0, 0].legend(frameon=False, fontsize=7, loc="center right")
    for suffix in ["png", "pdf", "svg"]:
        fig.savefig(out / ("dense_cluster_comparison." + suffix), dpi=300, facecolor="white")
    plt.close(fig)
    write(
        out / "FIGURE_IDENTITY.json",
        dict(
            analysis="Initial largest cluster; same stable IDs; 10 µs common physical sampling",
            interpolation=report.get("sampling_mode", "legacy_linear_statistics"),
            script_sha256=sha(Path(__file__)),
            report_sha256={c: sha(args.comparison / (c + ".json")) for c in ["B2", "B3"]},
        ),
    )


if __name__ == "__main__":
    main()
