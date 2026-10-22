"""Four English publication figures from saved trajectories and quality records only."""

from pathlib import Path
import numpy as np
import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist
from scripts.straight_tube_exp1.f10o.analyze import load, metrics
from .common import BASE, OLD, ROOT, CONFIG, read, write, sha, identity

PURPLE = "#654a86"
TEAL = "#237b80"
RED = "#b65050"
GRAY = "#65717e"
plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.labelcolor": "#20232a",
        "text.color": "#20232a",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,
        "savefig.dpi": 220,
    }
)


def raw(folder):
    t, x = load(folder)
    with h5py.File(folder / "trajectory.h5") as h:
        m = h["step_metrics"][:] if "step_metrics" in h else None
    return t, x, m


def main():
    out = BASE / "figures"
    out.mkdir(exist_ok=True)
    a = BASE / "analysis"
    runs = read(a / "RUN_ANALYSIS.json")
    comps = read(a / "COMPARISONS.json")
    gate = read(a / "FINAL_GATE.json")
    cost = read(a / "COSTS.json")
    val = BASE / "validation_v1/runs"
    b2 = ROOT / "results/straight_tube_exp1/f10o/validation_v1/runs"
    names = []

    def save(fig, name):
        fig.text(
            0.5,
            0.012,
            "EXPLORATORY  |  HOMOGENEOUS_WATER_CONTROL_PROXY  |  PHYSICS_UNCALIBRATED  |  training_ready=false",
            ha="center",
            fontsize=8,
            color="#41464d",
        )
        for ext in ["png", "pdf"]:
            fig.savefig(out / (name + "." + ext), bbox_inches="tight", facecolor="white")
        names.append(name)
        plt.close(fig)

    def panel(ax, title, xlabel=None, ylabel=None):
        ax.set_title(title, loc="left", fontweight="bold", fontsize=11, pad=10)
        if xlabel:
            ax.set_xlabel(xlabel)
        if ylabel:
            ax.set_ylabel(ylabel)
        ax.grid(alpha=0.15)

    def saved(name):
        return np.load(a / (name + ".npz"))

    fig, axs = plt.subplots(2, 2, figsize=(12, 8))
    fig.subplots_adjust(hspace=0.5, wspace=0.35, bottom=0.11, top=0.92)
    fig.suptitle("Hard-sphere contact: geometry, complementarity and cost", fontweight="bold", fontsize=15)
    ax = axs[0, 0]
    for folder, label, color in [
        (OLD / "campaign_v1/runs/OPTIMIZED_REPLAY", "M0 / PBD (saved frames)", GRAY),
        (b2 / "STRONG_DT1", "B2 soft contact", RED),
        (val / "STRONG_DT1", "Hard spheres (legal prefix)", PURPLE),
    ]:
        t, x, m = raw(folder)
        y = np.array([max(0.0, -np.min(pdist(row) - 3.56e-6)) for row in x]) * 1e9
        mask = t <= 0.0005
        ax.plot(t[mask] * 1e3, y[mask], label=label, color=color)
    ax.axhline(17.8, color=RED, ls="--", label="Original overlap budget: 17.8 nm")
    ax.set_xlim(0, 0.5)
    panel(ax, "a  Actual overlap", "Time (ms)", "Maximum overlap (nm)")
    ax.legend(fontsize=8)
    ax = axs[0, 1]
    t, x, m = raw(val / "STRONG_DT1")
    steps = np.arange(1, len(m) + 1) * 1e-3
    ax.plot(steps, m[:, 15] * 1e9, color=TEAL, label="Local-plane complementarity")
    ax.plot(steps, m[:, 17] * 1e9, color=PURPLE, label="True-geometry complementarity")
    ax.axhline(0.01, color=RED, ls="--", label="True-geometry limit: 0.01 nm")
    ax.axhline(0.001, color=TEAL, ls=":", label="Local-plane limit: 0.001 nm")
    ax.set_yscale("symlog", linthresh=1e-10)
    ax.set_ylim(0, .03)
    panel(ax, "b  Independent reaction-gap checks", "Attempted time (ms)", "Force-normalized gap (nm)")
    ax.legend(fontsize=8, loc="upper left")
    ax = axs[1, 0]
    ax.axis("off")
    f = read(BASE / "local_foundations_v2.json")
    table = ax.table(
        cellText=[
            [c["name"].replace("_", " "), c["status"], f"{c['independent_max_error_m'] * 1e9:.2g}"]
            for c in f["cases"]
        ],
        colLabels=["Reference case", "Gate", "Error (nm)"],
        cellLoc="center",
        colWidths=[0.55, 0.18, 0.27],
        bbox=[0, 0.03, 1, 0.84],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    panel(ax, "c  Two- and four-sphere independent references")
    ax = axs[1, 1]
    bv = read(b2 / "STRONG_DT1/RESULT.json")["mean_step_s"] * 1e3
    hs = runs["STRONG_DT1"]["mean_step_s"] * 1e3
    ax.bar(["B2\n0–0.175 ms", "Hard spheres\n0–0.464 ms"], [bv, hs], color=[RED, PURPLE], width=0.5)
    for i, v in enumerate([bv, hs]):
        ax.text(i, v + 0.1, f"{v:.3f}", ha="center")
    ax.set_ylim(0, max(bv, hs) * 1.28)
    panel(
        ax,
        "d  Measured whole-step cost",
        "Different valid windows; no matched-time speedup claim",
        "Wall time per step (ms)",
    )
    save(fig, "Fig_F10P_1_hard_sphere_contact")

    fig = plt.figure(figsize=(14, 8))
    gs = fig.add_gridspec(2, 3, hspace=0.5, wspace=0.4, bottom=0.1, top=0.90)
    fig.suptitle(
        "Strong control: actual legal prefixes and numerical sensitivity", fontweight="bold", fontsize=15
    )
    ax = fig.add_subplot(gs[0, 0])
    d = saved("STRONG_DT1")
    t = d["time_s"] * 1e3
    for j, lab, col in [(0, "X", PURPLE), (1, "Y", TEAL), (2, "Z", GRAY)]:
        ax.plot(t, (d["centroid_m"][:, j] - d["centroid_m"][0, j]) * 1e6, label=lab, color=col)
    ax.axvline(0.175, color=RED, ls="--", label="Old B2 rejection")
    panel(ax, "a  Center of the original 32 IDs", "Time (ms)", "Displacement (µm)")
    ax.legend(fontsize=8)
    ax = fig.add_subplot(gs[0, 1])
    ax.plot(t, (d["rg_m"] - d["rg_m"][0]) * 1e9, color=PURPLE)
    panel(ax, "b  Centered radius of gyration", "Time (ms)", "ΔRg (nm)")
    ax = fig.add_subplot(gs[0, 2])
    for suffix, col in [("DT05", TEAL), ("FINE", PURPLE)]:
        file = a / ("compare_STRONG_PREFIX_" + suffix + ".npz")
        if file.exists():
            q = np.load(file)
            ax.plot(q["time_s"] * 1e3, q["fixed_id_m"] * 1e9, color=col, label=suffix)
    ax.axhline(50, color=RED, ls="--", label="ID limit: 50 nm")
    panel(ax, "c  Same-ID difference, common legal time", "Time (ms)", "Maximum difference (nm)")
    ax.legend(fontsize=8)
    ax = fig.add_subplot(gs[1, 0], projection="3d")
    _, x, _ = raw(val / "STRONG_DT1")
    centered = (x[-1] - x[-1].mean(0)) * 1e6
    ax.scatter(*centered.T, c=np.arange(32), cmap="viridis", s=28, depthshade=True)
    for k in range(32):
        ax.text(*centered[k], str(k), fontsize=6)
    ax.set(xlabel="X (µm)", ylabel="Y (µm)", zlabel="Z (µm)")
    ax.set_title("d  Actual particle centers and IDs", loc="left", fontsize=11, fontweight="bold")
    ax.set_box_aspect(np.ptp(centered, axis=0))
    ax.view_init(24, -50)
    ax = fig.add_subplot(gs[1, 1])
    ax.axis("off")
    r = runs["STRONG_DT1"]
    last = read(val / "STRONG_DT1/LAST_STEP.json")["result"]
    text = f"Requested: 2.000 ms\nLegal prefix: {r['completed_s'] * 1e3:.3f} ms\nRejected proposal: 0.464 ms\n\nLocal QP equations: PASS\nSphere nonpenetration: PASS\nTrue-gap complementarity: FAIL\nMeasured: {last['actual_complementarity_equivalent_gap_m'] * 1e9:.5f} nm\nRegistered limit: 0.01000 nm\n\nNo full 2 ms strong admission"
    ax.text(0, 0.92, text, va="top", linespacing=1.65, fontsize=10)
    panel(ax, "e  Separate quality decisions")
    ax = fig.add_subplot(gs[1, 2])
    tags = ["STRONG_DT1", "STRONG_PREFIX_DT05", "STRONG_PREFIX_FINE"]
    present = [k for k in tags if k in runs]
    ax.barh(
        [k.replace("STRONG_", "").replace("PREFIX_", "prefix ") for k in present],
        [runs[k]["completed_s"] * 1e3 for k in present],
        color=[PURPLE, TEAL, GRAY][: len(present)],
    )
    ax.axvline(0.175, color=RED, ls="--")
    panel(ax, "f  Actually covered time", "Legal physical time (ms)")
    save(fig, "Fig_F10P_2_strong_control_validation")

    fig = plt.figure(figsize=(13, 8))
    gs = fig.add_gridspec(2, 2, hspace=0.52, wspace=0.32, bottom=0.1, top=0.89, height_ratios=[1, 1.15])
    fig.suptitle(
        "Group and particle admission are independent scientific statements", fontweight="bold", fontsize=15
    )
    ax = fig.add_subplot(gs[0, :])
    ax.axis("off")
    cols = [
        "BASIC_NUMERICS_PASS",
        "GROUP_POSITION_NUMERICALLY_ADMITTED",
        "GROUP_SHAPE_NUMERICALLY_ADMITTED",
        "GROUP_TOPOLOGY_NUMERICALLY_ADMITTED",
        "PARTICLE_ID_NUMERICALLY_ADMITTED",
        "EXPERIMENTALLY_CALIBRATED",
    ]
    rownames = ["STRONG", "OFF", "WEAK"]
    values = [[gate["groups"][tag][k] for k in cols] for tag in rownames]
    table = ax.table(
        cellText=values,
        rowLabels=rownames,
        colLabels=["Basic", "Group position", "Group shape", "2 ms network", "Particle ID", "Calibration"],
        cellLoc="center",
        bbox=[0.02, 0.08, 0.98, 0.77],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    for i, row in enumerate(values, 1):
        for j, value in enumerate(row):
            table[i, j].set_facecolor(
                {"PASS": "#dceae7", "FAIL": "#f1dedb", "NOT_RUN": "#ececf0", "NOT_VERIFIED": "#ececf0"}[value]
            )
    panel(ax, "a  No admission can be inferred from another column")
    ax = fig.add_subplot(gs[1, 0])
    vals = []
    for folder in [OLD / "campaign_v1/runs/SHAPE_REPLAY", b2 / "WEAK_DT1", val / "WEAK_DT1"]:
        tt, xx = load(folder)
        mm = metrics(tt, xx)
        vals.append((mm["rg_m"][-1] - mm["rg_m"][0]) * 1e9)
    ax.bar(["M0 / PBD", "B2 soft", "Hard sphere"], vals, color=[GRAY, RED, PURPLE])
    ax.axhline(0, color="#444", lw=0.8)
    for j, v in enumerate(vals):
        ax.annotate(
            f"{v:+.4f}", (j, v), xytext=(0, 5 if v >= 0 else -15), textcoords="offset points", ha="center"
        )
    ax.set_ylim(min(vals) - 1.5, max(vals) + 2)
    panel(
        ax,
        "b  Weak-control shape is contact-model dependent",
        "Same input, full 2 ms; different models",
        "ΔRg (nm)",
    )
    ax = fig.add_subplot(gs[1, 1])
    ax.axis("off")
    dt = comps["WEAK_DT05"]
    ca = comps["WEAK_FINE"]
    text = f"Weak control: full 2 ms numerical comparison\n\nCOM difference: {dt['COM_m'] * 1e9:.3g} / {ca['COM_m'] * 1e9:.3g} nm\nID difference: {dt['fixed_id_m'] * 1e9:.3g} / {ca['fixed_id_m'] * 1e9:.3g} nm\nRg difference: {dt['Rg_m'] * 1e9:.3g} / {ca['Rg_m'] * 1e9:.3g} nm\nNormalized Q: {dt['Q_relative']:.3g} / {ca['Q_relative']:.3g}\n(values: step / cache)\n\nOne initial network survives 2 ms.\nNo new persistent merger or physical fusion.\nNo observed group-pass / ID-fail hard-sphere case.\ntraining_ready = false"
    ax.text(0, 0.96, text, va="top", linespacing=1.7, fontsize=10)
    panel(ax, "c  Scope and measured differences")
    save(fig, "Fig_F10P_3_group_vs_particle_admission")

    fig, axs = plt.subplots(2, 2, figsize=(12, 8))
    fig.subplots_adjust(hspace=0.5, wspace=0.38, bottom=0.1, top=0.9)
    fig.suptitle(
        "XY navigation: registered targets, actual control and stopping evidence",
        fontweight="bold",
        fontsize=14,
    )
    ax = axs[0, 0]
    for name, label, col in [
        ("STRONG_DT1", "Strong: legal prefix only", PURPLE),
        ("WEAK_DT1", "Weak: 2 ms", TEAL),
        ("OFF_DT1", "OFF: stationary", GRAY),
    ]:
        d = saved(name)
        v = (d["centroid_m"] - d["centroid_m"][0]) * 1e6
        ax.plot(v[:, 0], v[:, 1], label=label, color=col, marker="o", markevery=[-1], markersize=4)
    ax.scatter([0.2, 0], [0, 0.2], marker="x", c=RED, s=70, label="Registered X / Y targets")
    ax.set_xlim(-0.025, 0.23)
    ax.set_ylim(-0.025, 0.23)
    ax.set_aspect("equal")
    panel(ax, "a  Actual COM; XY search NOT_RUN", "X (µm)", "Y (µm)")
    ax.legend(fontsize=8, loc="upper right")
    ax = axs[0, 1]
    for name, col in [("STRONG_DT1", PURPLE), ("WEAK_DT1", TEAL)]:
        d = saved(name)
        v = d["centroid_m"] - d["centroid_m"][0]
        for target, style in [(np.array([0.2e-6, 0, 0]), "-"), (np.array([0, 0.2e-6, 0]), "--")]:
            ax.plot(
                d["time_s"] * 1e3,
                np.linalg.norm(v - target, axis=1) * 1e6,
                style,
                color=col,
                label=name.replace("_DT1", "") + (" → X" if style == "-" else " → Y"),
            )
    ax.axhline(0.075, color=RED, ls=":", label="Target tolerance")
    panel(ax, "b  Error is not target admission", "Time (ms)", "Target error (µm)")
    ax.legend(fontsize=8, ncol=2)
    ax = axs[1, 0]
    v = [
        cost["CPU_solver_s"],
        cost["CPU_resistance_s"],
        cost["integration_s"] - cost["CPU_solver_s"] - cost["CPU_resistance_s"],
    ]
    ax.bar(["CPU contact", "CPU resistance", "Other integration*"], v, color=[PURPLE, TEAL, GRAY])
    panel(
        ax,
        "c  Measured integration cost",
        "*Includes GPU calls, transfers and evidence; not pure GPU time",
        "Wall time (s)",
    )
    ax = axs[1, 1]
    ax.axis("off")
    ax.text(
        0,
        0.95,
        f"Scientific numeric processes: {cost['total_numerical_s']:.2f} s\nRegistered budget: approximately 600 s\nPeak process RSS: {cost['peak_RSS_MiB']:.1f} MiB\n\nStop: true-gap complementarity at 0.464 ms\nStrong full 2 ms: FAIL\nNew XY optimization: NOT_RUN\n5 / 10 / 30 ms extensions: NOT_RUN\n\nNext: geometric consistency of contact reactions;\nthen repeat full-window group admission.\nNo spring tuning or quality-budget relaxation.",
        va="top",
        linespacing=1.75,
        fontsize=10,
    )
    panel(ax, "d  Stop on numerical evidence, not exhausted budget")
    save(fig, "Fig_F10P_4_navigation_and_cost")
    inputs = [
        CONFIG,
        *a.glob("*.json"),
        *a.glob("*.npz"),
        BASE / "local_foundations_v2.json",
        val / "STRONG_DT1/LAST_STEP.json",
    ]
    inputs += list(BASE.glob("*/runs/*/trajectory.h5"))
    inputs += [OLD / "campaign_v1/runs" / n / "trajectory.h5" for n in ["OPTIMIZED_REPLAY", "SHAPE_REPLAY"]]
    inputs += [b2 / n / f for n in ["STRONG_DT1", "WEAK_DT1"] for f in ["trajectory.h5", "RESULT.json"]]
    write(
        out / "FIGURE_PROVENANCE.json",
        dict(
            identity=identity(),
            source_script_sha256=sha(Path(__file__)),
            input_sha256={str(f.relative_to(ROOT)): sha(f) for f in inputs},
            scientific_execution_commits=sorted(
                {read(f)["identity"]["execution_commit"] for f in BASE.glob("*/CAMPAIGN.json")}
            ),
            figures=[
                dict(figure=n, file_sha256={ext: sha(out / (n + "." + ext)) for ext in ["png", "pdf"]})
                for n in names
            ],
        ),
    )


if __name__ == "__main__":
    main()
