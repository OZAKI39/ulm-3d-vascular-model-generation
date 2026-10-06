"""Paper-style diagnostics from new measured R2 evidence (PNG and vector PDF)."""

import json
from pathlib import Path
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scripts.stage5a_r2.seal import BASE
from scripts.stage4l.dense_probe import write, sha

plt.rcParams.update(
    {
        "font.size": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.facecolor": "white",
    }
)


def qp_figure():
    raw = BASE / "diagnosis/native_v1/SUMMARY.json"
    rows = json.loads(raw.read_text())
    replay_path = BASE / "diagnosis/reconstruction_v1/REPLAY.json"
    replay = json.loads(replay_path.read_text())
    fixtures = sorted({r["fixture"] for r in rows})
    fig, axes = plt.subplots(3, 3, figsize=(13.5, 9.2), layout="constrained")
    labels = ["OSQP", "ProxQP dense", "ProxQP sparse", "Clarabel"]
    names = ["osqp", "prox_dense", "prox_sparse", "clarabel"]
    for k, fixture in enumerate(fixtures):
        selected = [
            next(r for r in rows if r["fixture"] == fixture and r["solver"] == name and r["eps"] == 1e-12)
            for name in names
        ]
        residual = [r["kkt"] for r in selected]
        seconds = [r["seconds"] for r in selected]
        colors = ["#365f91" if v <= 1e-9 else "#ad4d49" for v in residual]
        ax = axes[k, 0]
        ax.bar(np.arange(4), residual, color=colors)
        ax.set_yscale("log")
        ax.axhline(1e-9, color="black", ls="--", lw=1, label="KKT gate")
        ax.set_xticks(range(4), labels, rotation=22, ha="right")
        ax.set_ylabel("Original KKT residual")
        ax.legend(fontsize=8)
        ax = axes[k, 1]
        ax.bar(np.arange(4), np.array(seconds) * 1e3, color=colors)
        ax.set_yscale("log")
        ax.set_xticks(range(4), labels, rotation=22, ha="right")
        ax.set_ylabel("Setup + solve wall time (ms)")
        record = next(r for r in replay if r["fixture"] == fixture and "spectrum" in r)
        near = record["spectrum"]["near_rows"]
        ax = axes[k, 2]
        ax.scatter(
            [max(1e-17, r["row_difference_l2"]) for r in near],
            [max(1e-17, abs(r["lower_difference"])) for r in near],
            s=14,
            alpha=0.7,
            color="#704f95",
        )
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Row difference, L2")
        ax.set_ylabel("Absolute lower-bound difference")
        title = [
            "Earlier active set (114 variables)",
            "B2 OSQP failure (51 variables)",
            "B3 OSQP failure (18 variables)",
        ][k]
        axes[k, 0].set_title(title, loc="left", fontweight="bold")
        axes[k, 2].set_title(f"{len(near)} near-parallel row pairs", loc="left")
    fig.suptitle(
        "R3 | Original failed QPs: solver recovery and near-dependent inequalities", fontweight="bold"
    )
    fig.text(
        0.02,
        0.005,
        "Original matrices retained; internal tolerance 1e−12. Red = independent certificate rejected. Parallel diagnostic timing; not trajectory speedup.",
        fontsize=8,
    )
    out = BASE / "figures"
    out.mkdir(exist_ok=True)
    for suffix in ["png", "pdf"]:
        fig.savefig(out / f"R3_QP_DIAGNOSTIC.{suffix}", dpi=240, bbox_inches="tight")
    plt.close(fig)
    # Separate KKT terms retain native and failed polishing candidates.
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.3), layout="constrained")
    terms = ["stationarity", "feasibility", "dual_feasibility", "active_slack", "complementarity"]
    for ax, fixture in zip(axes, fixtures[1:]):
        rec = next(r for r in replay if r["fixture"] == fixture and "reconstruction" in r)
        attempt = rec["diagnostics"]["attempts"][-1]
        original = [max(1e-17, attempt["kkt_terms"][k]) for k in terms]
        recovered = next(
            r for r in rows if r["fixture"] == fixture and r["solver"] == "prox_dense" and r["eps"] == 1e-12
        )
        ax.plot(range(5), original, "o-", label="Old ADMM candidate", color="#a65c48")
        ax.plot(
            range(5),
            [max(1e-17, recovered["terms"][k]) for k in terms],
            "o-",
            label="R2 ProxQP",
            color="#315f89",
        )
        ax.axhline(
            attempt["external_sparse_polish"]["kkt"], color="#af3737", ls=":", label="Old polish: full KKT"
        )
        ax.axhline(1e-9, color="black", ls="--", label="Unchanged gate")
        ax.set_yscale("log")
        ax.set_xticks(range(5), terms, rotation=25, ha="right")
        ax.set_ylabel("Independent residual")
        ax.legend(fontsize=8)
        ax.set_title(fixture.split("_")[-2])
    for suffix in ["png", "pdf"]:
        fig.savefig(out / f"R3b_KKT_TERMS.{suffix}", dpi=240, bbox_inches="tight")
    plt.close(fig)
    write(
        out / "R3_DATA_SOURCES.json",
        dict(
            raw=str(raw),
            raw_sha256=sha(raw),
            reconstruction=str(replay_path),
            reconstruction_sha256=sha(replay_path),
            script_sha256=sha(Path(__file__)),
            definitions="Fixed Stage5A normalization, every original inequality. Timing is parallel native API diagnostic, never end-to-end speedup.",
        ),
    )


if __name__ == "__main__":
    qp_figure()
