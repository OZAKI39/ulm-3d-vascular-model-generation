"""Four English PNG/PDF groups from archived measured arrays; no simulated artwork."""

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import h5py
from .common import ROOT, BASE, js, sha, write

COLORS = ["#B64E48", "#737D90", "#287C90", "#3A815D"]
NAMES = ["F10I v2\nFAIL", "F10J corrected\nPASS", "F10J optimized\nPASS", "CPU reference\nPASS"]


def main():
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "legend.fontsize": 8.5,
            "pdf.fonttype": 42,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "savefig.facecolor": "white",
        }
    )
    out = BASE / "figures"
    out.mkdir(exist_ok=True)
    old = js(BASE / "analysis/ORIGINAL_AUDIT.json")
    v2 = js(BASE / "optimized_v2/RESULTS.json")
    oldraw = js(BASE / "audit_v1/PLAIN.json")["result"]
    new = v2["rows"][-1]["gpu"]
    metrics = v2["rows"][-1]["metrics"]
    perf = js(BASE / "performance_v1/RESULTS.json")
    cold = js(BASE / "process_v1/RESULTS.json")
    screen = js(BASE / "analysis/TOPOLOGY_SCREEN.json")
    provenance = []

    def save(fig, name, sources, caption):
        fig.subplots_adjust(
            top=0.89,
            bottom=0.155 if name.endswith("sync_cost") else 0.13,
            left=0.09,
            right=0.96,
            hspace=0.48,
            wspace=0.32,
        )
        fig.text(0.09, 0.03, caption, fontsize=8.5, va="bottom", color="#30343B")
        for ext in ["png", "pdf"]:
            fig.savefig(out / (name + "." + ext), dpi=220)
        plt.close(fig)
        provenance.append(
            dict(
                figure=name,
                source_sha256={p: sha(BASE / p) for p in sources},
                caption=caption,
                files={ext: sha(out / (name + "." + ext)) for ext in ["png", "pdf"]},
            )
        )

    fig, aa = plt.subplots(2, 2, figsize=(12, 8.5))
    fig.suptitle("PCG residuals and independent fixed-state accuracy", fontsize=16, y=0.97)
    ax = aa[0, 0]
    trace = np.array(old["linear"][0]["trace"])
    for col, label, color, style in [
        (0, "Recursive residual", COLORS[0], "-"),
        (1, "Recomputed true residual", "#252C36", "--"),
        (4, "Preconditioned true residual", COLORS[1], ":"),
    ]:
        ax.semilogy(np.arange(len(trace)), trace[:, col], style, label=label, color=color)
    ax.axhline(old["linear"][0]["tolerance"], color="black", lw=0.9, ls="-.", label="Original PCG target")
    ax.scatter(
        [132],
        [new["linear_audit"][0][1]],
        s=65,
        marker="*",
        color=COLORS[2],
        zorder=5,
        label="F10J true residual at 132",
    )
    ax.set(
        xlabel="Base PCG iteration",
        ylabel="Dimensionless residual norm",
        title="a  Original 128-iteration base solve fails",
    )
    ax.legend(loc="lower left")
    ax = aa[0, 1]
    for kind, rows, color, offset in [
        ("F10I", old["linear"][1:], COLORS[0], -0.06),
        ("F10J", new["linear_audit"][1:], COLORS[2], 0.06),
    ]:
        if kind == "F10I":
            rec = [r["recursive_final"] / r["tolerance"] for r in rows]
            true = [r["true_final"] / r["tolerance"] for r in rows]
        else:
            rec = [r[0] / r[2] for r in rows]
            true = [r[1] / r[2] for r in rows]
        ax.semilogy(
            np.arange(1, len(rows) + 1) + offset, true, "o-", color=color, label=kind + " true / target"
        )
        ax.scatter(
            np.arange(1, len(rows) + 1) + offset,
            rec,
            marker="x",
            color=color,
            s=55,
            label=kind + " recursive / target",
        )
    ax.axhline(1, color="black", ls="--", lw=0.9)
    ax.set(
        xlabel="Contact direction solve",
        ylabel="Residual / original PCG target",
        title="b  First F10I direction also misses its target",
    )
    ax.legend(loc="upper right", ncol=2, fontsize=8)
    ax = aa[1, 0]
    olderr = js(ROOT / "results/fast_pzt/f10i/gpu_v2/RESULTS.json")["rows"][-1]["metrics"]
    labels = ["Gradient bound", "Position equation\nresidual", "Position vs.\nindependent CPU"]
    for idx, (vals, color, label) in enumerate(
        [
            (
                [
                    oldraw["conservative_equation_bound_m"],
                    olderr["max_equation_residual_m"],
                    olderr["max_position_vs_independent_m"],
                ],
                COLORS[0],
                "F10I: rejected",
            ),
            (
                [
                    new["conservative_equation_bound_m"],
                    metrics["max_equation_residual_m"],
                    metrics["max_position_vs_independent_m"],
                ],
                COLORS[2],
                "F10J: fixed-state pass",
            ),
        ]
    ):
        ax.bar(np.arange(3) + (idx - 0.5) * 0.32, vals, width=0.3, color=color, label=label)
    ax.set_yscale("log")
    ax.set_xticks(range(3), labels)
    ax.axhline(1e-14, color="black", lw=0.9, ls="--", label="Unchanged 1e-14 m target")
    ax.set(ylabel="Magnitude (m)", title="c  Residual bounds are not position-error identities")
    ax.legend(loc="upper right", fontsize=8)
    ax = aa[1, 1]
    curve = old["linear"][0]["position_equivalent_error_curve_m"]
    ax.semilogy(
        np.arange(len(curve)),
        np.maximum(curve, 1e-30),
        color=COLORS[1],
        label="R times max. per-bubble linear solution error",
    )
    ax.set(
        xlabel="Original base PCG iteration",
        ylabel="Position-equivalent linear error (m)",
        title="d  Independent dense reference: same operator",
    )
    ax.text(
        0.04,
        0.10,
        "A is SPD: minimum eigenvalue ≈ 1\nNo nonpositive PCG curvature observed\nBase true residual / target = 2.416 at iteration 128\nLinear solution error is distinct from nonlinear error",
        transform=ax.transAxes,
        fontsize=9,
    )
    save(
        fig,
        "Fig_F10J_1_pcg_convergence",
        ["analysis/ORIGINAL_AUDIT.json", "corrected_v1/RESULTS.json", "optimized_v2/RESULTS.json"],
        "Same 32-bubble saved state, double precision, unchanged B2 model and stiffness. F10I failure is retained.\nNew residual marker is a measured endpoint; no unrecorded F10J convergence curve is interpolated.",
    )

    fig, aa = plt.subplots(2, 2, figsize=(12, 8.5))
    fig.suptitle("Nonlinear convergence and cancellation in Armijo acceptance", fontsize=16, y=0.97)
    ax = aa[0, 0]
    for result, color, label in [
        (oldraw, COLORS[0], "F10I: line-search limit"),
        (new, COLORS[2], "F10J: converged"),
    ]:
        vals = [r["gradient_bound_before_m"] for r in result["iterations"]] + [
            result["conservative_equation_bound_m"]
        ]
        ax.semilogy(
            range(len(vals)),
            vals,
            "x--" if result is oldraw else "o-",
            color=color,
            label=label,
            zorder=4 if result is oldraw else 3,
            markersize=8 if result is oldraw else 5,
        )
    ax.axhline(1e-14, color="black", ls="--", lw=0.9, label="Original stopping target")
    ax.set(
        xlabel="Accepted outer updates",
        ylabel="Conservative gradient bound (m)",
        title="a  One additional accepted update resolves the failure",
    )
    ax.legend()
    last = [r for r in old["evaluations"] if r["trial"] and r["outer"] == 4]
    ax = aa[0, 1]
    steps = np.arange(20)
    ax.plot(
        steps,
        [r["original_energy_change"] for r in last],
        "o-",
        ms=3,
        color=COLORS[0],
        label="Subtracted total energies: reject all",
    )
    ax.plot(
        steps,
        [r["stable_energy_change"] for r in last],
        "o-",
        ms=3,
        color=COLORS[2],
        label="Stable independent increment: accept all",
    )
    ax.plot(
        steps, [r["armijo_change_rhs"] for r in last], "--", color="black", label="Armijo increment threshold"
    )
    ax.set_yscale("symlog", linthresh=1e-23)
    ax.set(
        xlabel="Backtrack index k (step = 2⁻ᵏ)",
        ylabel="Dimensionless energy change",
        title="b  Original failed outer update: all 20 candidates",
    )
    ax.legend(fontsize=8, loc="upper right", bbox_to_anchor=(1, 0.80))
    ax = aa[1, 0]
    ax.semilogy(
        steps,
        [r["equation_residual_m"] for r in last],
        "o-",
        color=COLORS[1],
        label="True position equation residual",
    )
    ax.semilogy(
        steps, [r["gradient_bound_m"] for r in last], "s--", color=COLORS[0], label="Original gradient bound"
    )
    ax.axhline(1e-14, color="black", ls="--", lw=0.9)
    ax.set(
        xlabel="Original backtrack index k",
        ylabel="Magnitude (m)",
        title="c  The rejected full step already meets the target",
    )
    ax.legend(fontsize=8)
    ax = aa[1, 1]
    correction = np.array([r["max_correction_m"] for r in new["iterations"]]) * 1e9
    ax.semilogy(
        np.arange(1, len(correction) + 1),
        correction,
        "o-",
        color=COLORS[2],
        label="F10J accepted maximum correction",
    )
    ax.axhline(178, color="black", ls="--", lw=0.9, label="Unchanged 178 nm budget")
    sets = [set(map(tuple, r["active_pairs"])) for r in last]
    ax.text(
        0.25,
        0.97,
        f"All {len(new['candidate_trace'])} F10J trials accepted at full step\nOriginal failed candidates: {len(sets[0])} active contacts\nActive set unchanged across all 20: {all(v == sets[0] for v in sets)}\nFinal overlap: {metrics['max_overlap_m'] * 1e9:.6f} nm",
        transform=ax.transAxes,
        fontsize=8.5,
        va="top",
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.9),
    )
    ax.set(
        xlabel="Accepted F10J outer update",
        ylabel="Maximum single correction (nm)",
        title="d  Geometry and original correction budget retained",
    )
    ax.legend(loc="lower left", fontsize=8)
    save(
        fig,
        "Fig_F10J_2_nonlinear_linesearch",
        ["analysis/ORIGINAL_AUDIT.json", "audit_v1/PLAIN.json", "optimized_v2/RESULTS.json"],
        "Stable ΔE evaluates the same unilateral-contact energy difference; Armijo c = 1e-4, maximum 20 backtracks, maximum 24 updates.\nNo force, stiffness, geometric graph threshold, overlap budget, or stopping tolerance was relaxed.",
    )

    names = ["original", "corrected", "optimized", "cpu"]
    vals = {
        n: np.array(
            [
                next(v for v in s["measurements"] if v["name"] == n)["matched_full_wall_s"]
                for s in perf["samples"]
            ]
        )
        * 1e3
        for n in names
    }
    fig, aa = plt.subplots(2, 2, figsize=(12, 9))
    fig.suptitle("Numerically admitted GPU solve has no B = 1 cost advantage", fontsize=16, y=0.97)
    ax = aa[0, 0]
    for i, n in enumerate(names):
        ax.bar(i, np.median(vals[n]), color=COLORS[i], alpha=0.8, hatch="//" if i == 0 else None)
        ax.scatter(i + np.array([-0.07, 0, 0.07]), vals[n], color="black", s=16, zorder=5)
        ax.text(i, max(vals[n]) + 1, f"{np.median(vals[n]):.2f}", ha="center", fontsize=9)
    ax.set_xticks(range(4), NAMES, fontsize=8)
    ax.set(
        ylabel="Matched full wall time (ms)",
        ylim=(0, 50),
        title="a  Three actual pairs; dots show each measurement",
    )
    ax = aa[0, 1]
    decisions = [perf["diagnostics"][n]["result"]["timing"]["scalar_reads"] for n in names[:3]]
    ax.bar(range(3), decisions, color=COLORS[:3], label="Intermediate decision readbacks")
    ax.bar(
        range(3), [12, 14, 14], bottom=decisions, color="#CFD4DD", label="Final vector / residual readbacks"
    )
    for i, v in enumerate(decisions):
        ax.text(i, v / 2, str(v), color="white", ha="center", weight="bold")
    ax.set_xticks(range(3), NAMES[:3], fontsize=8)
    ax.set(ylabel="Host array readback calls", title="b  Readbacks reduced; final evidence is still read")
    ax.legend(fontsize=8)
    ax = aa[1, 0]
    for k, (label, extract, color) in enumerate(
        [
            (
                "PCG graph device span",
                lambda d: sum(e["graph_cuda_span_s"] for e in d["pcg_graphs"]),
                COLORS[1],
            ),
            (
                "Scalar wait + read wall",
                lambda d: d["result"]["timing"]["scalar_sync_read_wall_s"],
                COLORS[2],
            ),
            (
                "Armijo evaluation wall",
                lambda d: sum(e["wall_s"] for e in d["evaluations"] if e["trial"]),
                COLORS[0],
            ),
        ]
    ):
        ax.bar(
            np.arange(3) + (k - 1) * 0.25,
            [extract(perf["diagnostics"][n]) * 1e3 for n in names[:3]],
            width=0.23,
            label=label,
            color=color,
        )
    ax.set_xticks(range(3), NAMES[:3], fontsize=8)
    ax.set(ylabel="One diagnostic profile (ms)", title="c  Overlapping categories: do not add these bars")
    ax.legend(fontsize=8)
    ax = aa[1, 1]
    builds = perf["construction"]
    ax.bar(
        range(3),
        [r["original_graph_build_s"] * 1e3 for r in builds],
        color=COLORS[1],
        label="Inherited graph construction",
    )
    ax.bar(
        range(3),
        [r["corrected_graph_build_s"] * 1e3 for r in builds],
        bottom=[r["original_graph_build_s"] * 1e3 for r in builds],
        color=COLORS[2],
        label="Opt-in replacement graph construction",
    )
    ax.set_xticks(range(3), NAMES[:3], fontsize=8)
    ax.set(ylabel="Graph construction wall time (ms)", title="d  Graph costs are outside the resident timing")
    ax.legend(fontsize=8)
    startup = ", ".join(f"{r['method']}: {r['process_wall_s']:.3f} s" for r in cold["rows"])
    prep = " / ".join(
        f"{perf['diagnostics'][n]['result']['timing']['preparation_wall_s'] * 1e3:.3f}" for n in names[:3]
    )
    save(
        fig,
        "Fig_F10J_3_gpu_sync_cost",
        ["performance_v1/RESULTS.json", "process_v1/RESULTS.json"],
        "Full time includes preparation, solve, identical independent quality checks and matched serialized output. CPU uses the same accuracy targets.\nIndependent new processes (one each; warm compiler cache): "
        + startup
        + "\nGPU buffer/input preparation in diagnostic profiles (original / corrected / optimized): "
        + prep
        + " ms."
        + "\nF10I is a failed diagnostic, not a production speed baseline. No long trajectory or batch throughput was measured.",
    )

    fig = plt.figure(figsize=(12, 9))
    fig.suptitle("Fixed-state numerical admission and limits of scientific coverage", fontsize=16, y=0.97)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.1, 1.0])
    ax = fig.add_subplot(gs[0, :])
    ax.axis("off")
    table_rows = [
        ["Two-bubble references", "12 / 12 matched", "PASS", "Original equation and forces"],
        [
            "Four-bubble references",
            "4 accepted; 2 budget negatives",
            "PASS",
            "Both over-budget states rejected",
        ],
        [
            "Real 32-bubble state",
            f"{metrics['max_equation_residual_m']:.3e} m residual",
            "PASS",
            "Old 1e-14 m gradient criterion met",
        ],
        [
            "Contact budgets",
            f"{metrics['max_overlap_m'] * 1e9:.3f} nm; {metrics['max_single_iteration_m'] * 1e9:.3f} nm",
            "PASS",
            "Overlap ≤17.8 nm; one correction ≤178 nm",
        ],
        ["GPU implementation", "Corrected / optimized positions identical", "PASS", "Fixed-state tests only"],
        [
            "Fair B = 1 performance",
            f"GPU {np.median(vals['optimized']):.2f} / CPU {np.median(vals['cpu']):.2f} ms",
            "NO ADVANTAGE",
            "Stop production integration",
        ],
        [
            "Contact / near-field physics",
            "Exploratory stiffness and approximations",
            "UNCALIBRATED",
            "Not experimental truth",
        ],
        ["AI data admission", "training_ready = false", "NOT ADMITTED", "No new dataset or training"],
    ]
    table = ax.table(
        cellText=table_rows,
        colLabels=["Scope", "Measured evidence", "Decision", "Limit"],
        cellLoc="left",
        colLoc="left",
        colWidths=[0.20, 0.32, 0.15, 0.33],
        bbox=[0, 0, 1, 1],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    for (i, j), cell in table.get_celld().items():
        cell.set_edgecolor("#DEE1E7")
        if i == 0:
            cell.set_facecolor("#E6EBF1")
            cell.set_text_props(weight="bold")
        elif j == 2:
            cell.set_text_props(
                color=COLORS[3] if table_rows[i - 1][2] == "PASS" else COLORS[0], weight="bold"
            )
    ax = fig.add_subplot(gs[1, 0])
    with h5py.File(ROOT / "results/fast_pzt/f5/science_v1/main/trajectory.h5") as f:
        ax.step(
            f["time_s"][:] * 1e3,
            f["metrics/cluster_count"][:, 0],
            where="post",
            color=COLORS[1],
            label="Historical F5 M0 (saved frames)",
        )
    ax.axvline(300, color="black", ls="--", lw=0.9, label="Ultrasound switched off")
    ax.axvline(417.5, color=COLORS[2], ls=":", label="Candidate saved state: 417.5 ms")
    ax.set(
        xlabel="Historical F5 time (ms)",
        ylabel="Connected components",
        title="Historical risk candidate; not B2 ground truth",
        xlim=(300, 600),
        ylim=(0.4, 12),
    )
    ax.legend(fontsize=8, loc="upper left")
    ax = fig.add_subplot(gs[1, 1])
    ax.axis("off")
    ax.text(
        0, 1, "PERSISTENT TOPOLOGY VALIDATION NOT RUN", color=COLORS[0], weight="bold", fontsize=11, va="top"
    )
    text = (
        f"{screen['reconstruction_status_counts']['ACCEPTED_LEDGER_RECONSTRUCTION_PASS']} historical ledgers reconstructed; original 2 ms dwell.\n"
        "F5 release: 11 persistent splits, 2 merges.\n"
        "Candidate needs B2 restart and near-wall validation.\n\n"
        "F10H/I geometric B2: one connected group throughout.\n"
        "No observed persistent split or merge to validate.\n"
        "Original strict edge-time failure remains: 12 / 29.\n\n"
        "Group position / shape evidence is retained.\n"
        "No new topology trajectory, B64 run, or AI training.\n"
        "Decision: preserve CPU reference; stop GPU production route."
    )
    ax.text(0, 0.83, text, fontsize=10, va="top", linespacing=1.5)
    save(
        fig,
        "Fig_F10J_4_solver_admission_scope",
        [
            "optimized_v2/RESULTS.json",
            "performance_v1/RESULTS.json",
            "analysis/TOPOLOGY_SCREEN.json",
            "analysis/TOPOLOGY_CANDIDATE.json",
        ],
        "EXPLORATORY / PHYSICS_UNCALIBRATED. Numerical fixed-state agreement does not calibrate contact physics or validate persistent topology.\nHistorical runs are correlated evidence, not independent scientific samples. Original F10I failures and all quality budgets are retained.",
    )
    write(
        out / "FIGURE_PROVENANCE.json",
        dict(
            script_sha256=sha(__file__),
            figures=provenance,
            old_F10I_reference_sha256=sha(ROOT / "results/fast_pzt/f10i/gpu_v2/RESULTS.json"),
            historical_F5_trajectory_sha256=sha(ROOT / "results/fast_pzt/f5/science_v1/main/trajectory.h5"),
        ),
    )


if __name__ == "__main__":
    main()
