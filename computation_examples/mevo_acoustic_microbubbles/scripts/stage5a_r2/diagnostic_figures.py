"""Paper-style diagnostics from new measured R2 evidence (PNG and vector PDF)."""

import json
from pathlib import Path
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scripts.stage5a_r2.seal import BASE, ROOT
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
    polish_path = BASE / "diagnosis/polish_parts/SUMMARY.json"
    polish = json.loads(polish_path.read_text())
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
    fig.get_layout_engine().set(rect=(0, 0.045, 1, 0.955))
    out = BASE / "figures"
    out.mkdir(exist_ok=True)
    for suffix in ["png", "pdf"]:
        fig.savefig(out / f"R3_QP_DIAGNOSTIC.{suffix}", dpi=240, bbox_inches="tight")
    plt.close(fig)
    # Separate KKT terms retain native and failed polishing candidates.
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.3), layout="constrained")
    terms = ["stationarity", "feasibility", "dual_feasibility", "active_slack", "complementarity"]
    for ax, fixture in zip(axes, fixtures[1:]):
        parts = next(r for r in polish["rows"] if r["fixture"] == fixture)
        original = [max(1e-17, parts["before"]["kkt_terms"][k]) for k in terms]
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
        ax.plot(
            range(5),
            [max(1e-17, parts["after"]["kkt_terms"][k]) for k in terms],
            "s:",
            color="#af3737",
            label="Old one-shot polish: replayed",
        )
        ax.axhline(1e-9, color="black", ls="--", label="Unchanged gate")
        ax.set_yscale("log")
        ax.set_xticks(range(5), terms, rotation=25, ha="right")
        ax.set_ylabel("Independent residual")
        ax.legend(fontsize=8)
        ax.set_title(fixture.split("_")[-2])
    fig.suptitle("R3b | Immutable candidates before and after old saddle polish")
    fig.text(
        0.01,
        -0.06,
        "Original Stage5A base 7c39828d9350; reconstruction_v1 / polish_parts. Zeros drawn at 1e-17 only on this log plot.\n"
        "Supplementary local SciPy polish replay; all original residual values and server results remain in source data.",
        fontsize=7,
    )
    for suffix in ["png", "pdf"]:
        fig.savefig(out / f"R3b_KKT_TERMS.{suffix}", dpi=240, bbox_inches="tight")
    plt.close(fig)
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.5), layout="constrained")
    for ax, fixture in zip(axes, fixtures):
        selected = [
            next(r for r in rows if r["fixture"] == fixture and r["solver"] == name and r["eps"] == 1e-12)
            for name in names
        ]
        ax.bar(range(4), [r["iterations"] for r in selected], color="#43678f")
        ax.set_xticks(range(4), labels, rotation=25, ha="right")
        ax.set_yscale("log")
        ax.set_ylabel("Native reported iterations")
        ax.set_title(fixture.split("_")[-2])
        for i, r in enumerate(selected):
            ax.text(i, r["iterations"] * 1.08, str(r["iterations"]), ha="center", fontsize=8)
    fig.suptitle("R3c | Native iteration counts on the three unchanged failed QPs")
    fig.text(
        0.01,
        -0.07,
        "native_v1 / original matrices from Stage5A base 7c39828d9350. Counts describe different algorithms and are not equal work units.\n"
        "Native statuses are retained in SUMMARY.json; full internal iteration histories are unavailable and are not reconstructed.",
        fontsize=7,
    )
    for suffix in ["png", "pdf"]:
        fig.savefig(out / f"R3c_NATIVE_ITERATIONS.{suffix}", dpi=240, bbox_inches="tight")
    plt.close(fig)
    write(
        out / "R3_DATA_SOURCES.json",
        dict(
            raw=str(raw),
            raw_sha256=sha(raw),
            reconstruction=str(replay_path),
            reconstruction_sha256=sha(replay_path),
            polish_parts=str(polish_path),
            polish_parts_sha256=sha(polish_path),
            original_failed_code_commit="7c39828d93509fff5551b7a16537ca3aa4ee9c23",
            native_replay_script_sha256=sha(ROOT / "scripts/stage5a_r2/qp_explore.py"),
            script_sha256=sha(Path(__file__)),
            definitions="Fixed Stage5A normalization, every original inequality. Timing is parallel native API diagnostic, never end-to-end speedup.",
        ),
    )


def convergence_figure():
    folder = BASE / "diagnosis/manufactured_v2"
    rows = json.loads((folder / "SUMMARY.json").read_text())
    fig, axes = plt.subplots(2, 3, figsize=(13, 7.5), layout="constrained")
    colors = ["#72869b", "#43678f", "#2d3c69"]
    for i, case in enumerate(["pair", "wall"]):
        items = sorted(
            [r for r in rows if r["case"] == case and r["tolerance_factor"] == 1],
            key=lambda r: -r["max_dt_s"],
        )
        for color, row in zip(colors, items):
            grid = row["max_dt_s"] / 1e-4
            d = np.load(folder / f"{case}_tol1_grid{grid:g}.npz")
            label = f"Max macro dt = {row['max_dt_s'] * 1e6:g} µs ({row['accepted_steps']} steps)"
            axes[i, 0].plot(d["time_s"] * 1e3, d["gap_m"] * 1e9, color=color, label=label)
            axes[i, 1].plot(
                d["time_s"] * 1e3, abs(d["gap_m"] - d["reference_gap_m"]) * 1e9, color=color, label=label
            )
        axes[i, 0].plot(
            d["time_s"] * 1e3,
            d["reference_gap_m"] * 1e9,
            "--",
            color="black",
            lw=1.2,
            label="Independent scalar ODE",
        )
        axes[i, 0].scatter(
            [1.0],
            [items[0]["raw_large_step_gap_m"] * 1e9],
            marker="^",
            color="#b05353",
            s=55,
            label="Single 1 ms step: safe, inaccurate",
        )
        axes[i, 0].set_ylabel("Surface gap (nm)")
        axes[i, 1].set_ylabel("Absolute gap error (nm)")
        for ax in axes[i, :2]:
            ax.set_xlabel("Physical time (ms)")
        axes[i, 0].set_title(
            "Stiff pair squeeze" if case == "pair" else "Stiff wall approach", loc="left", fontweight="bold"
        )
        axes[i, 0].legend(fontsize=7)
        h = np.array([r["max_dt_s"] / 2 for r in items]) * 1e6
        errors = np.array([r["endpoint_error_m"] for r in items]) * 1e9
        axes[i, 2].loglog(h, errors, "o-", color="#315f89", label="Actual terminal error")
        axes[i, 2].loglog(h, errors[0] * h / h[0], "--", color="gray", label="First-order slope")
        axes[i, 2].set_xlabel("Maximum accepted half-step dt (µs)")
        axes[i, 2].set_ylabel("Terminal gap error (nm)")
        axes[i, 2].legend(fontsize=8)
    fig.suptitle("R4 | S1 V4: continuous stiff motion and actual grid refinement", fontweight="bold")
    out = BASE / "figures"
    for suffix in ["png", "pdf"]:
        fig.savefig(out / f"R4_CONVERGENCE.{suffix}", dpi=240, bbox_inches="tight")
    plt.close(fig)
    write(
        out / "R4_DATA_SOURCES.json",
        dict(
            summary=str(folder / "SUMMARY.json"),
            summary_sha256=sha(folder / "SUMMARY.json"),
            script_sha256=sha(Path(__file__)),
            model="Frozen Stage4L self drag + normal lubrication; constant imposed flow; float64",
            reference="Independent scalar solve_ivp, rtol=2e-12, atol=1e-20; no coordinate interpolation in simulated trajectory",
            source_note="Development manufactured_v2 diagnostic; formal B2/B3 gates are separate. Explicit max-dt refinement uses 20/40/80 actual accepted steps.",
        ),
    )


if __name__ == "__main__":
    qp_figure()
    convergence_figure()
