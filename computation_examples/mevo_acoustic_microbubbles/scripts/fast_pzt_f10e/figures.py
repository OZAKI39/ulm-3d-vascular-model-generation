"""Publication figures from retained F10D/F10E arrays; no synthetic trajectories."""

import numpy as np
import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from .common import BASE, ROOT, CAMPAIGN, js, sha, write, identity
from .analyze import trace

COLORS = ["#3d536e", "#21827a", "#ad702f", "#8b619f"]
LABELS = ["C0: 10 µs / 12 rounds", "C1: 5 µs / 12 rounds", "C2: 10 µs / 24 rounds", "C3: 5 µs / 24 rounds"]
NAMES = [
    "Fig_F10E_1_batch_failure_isolation",
    "Fig_F10E_2_contact_convergence",
    "Fig_F10E_3_trajectory_morphology",
    "Fig_F10E_4_valid_yield_cost",
]
OUT = BASE / "figures"


def style():
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "legend.fontsize": 8.5,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "pdf.fonttype": 42,
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
            "axes.edgecolor": "#444444",
            "text.color": "#222222",
            "axes.labelcolor": "#222222",
        }
    )


def save(fig, index, note):
    fig.text(0.5, 0.016, note, ha="center", fontsize=8, color="#444444")
    fig.savefig(OUT / (NAMES[index] + ".png"), dpi=240)
    fig.savefig(OUT / (NAMES[index] + ".pdf"))
    plt.close(fig)


def figure1():
    fig, axs = plt.subplots(1, 2, figsize=(13.8, 6.2), gridspec_kw={"width_ratios": [1.25, 1]})
    fig.subplots_adjust(left=0.085, right=0.975, bottom=0.17, top=0.83, wspace=0.34)
    fig.suptitle("F10E | Batch isolation reveals four independent S/O failures", fontsize=15, y=0.96)
    names = ["S / seed1101", "S / seed1102", "O / seed1101", "O / seed1102"]
    cases = ["S_1000_1101", "S_1000_1102", "O_1000_1101", "O_1000_1102"]
    intervals = []
    for case in cases:
        with h5py.File(BASE / "isolation_v1" / case / "trajectory.h5") as h:
            intervals.append([h["rejected_interval"].attrs[k] * 1e3 for k in ["start_s", "attempted_end_s"]])
    ax = axs[0]
    for i, (case, (lo, hi)) in enumerate(zip(cases, intervals)):
        parent_end = 5.0 if case.startswith("S") else 4.85
        ax.plot([0, lo], [i, i], color="#bc7845", lw=8, solid_capstyle="butt")
        ax.plot([0, parent_end], [i + 0.16, i + 0.16], color="#8b98a5", lw=3)
        ax.plot(hi, i, marker="x", color="#923f3f", ms=9, mew=2)
        emitter = i in (0, 3)
        ax.text(hi + 0.15, i, f"({lo:.2f}, {hi:.2f}] ms", va="center", fontsize=9)
        ax.text(
            0.12,
            i - 0.18,
            "Original batch: " + ("overlap observed" if emitter else "stopped by peer"),
            fontsize=8,
            color="#4f5e70",
        )
    ax.set(
        yticks=range(4), yticklabels=names, xlim=(0, 9), ylim=(3.6, -0.8), xlabel="Branch physical time (ms)"
    )
    ax.set_title("(a) Accepted prefixes and independent failure intervals", loc="left")
    ax.grid(axis="x", alpha=0.2)
    ax.legend(
        handles=[
            plt.Line2D([0], [0], color="#8b98a5", lw=3, label="Original batch accepted prefix"),
            plt.Line2D([0], [0], color="#bc7845", lw=6, label="Independent retry accepted prefix"),
        ],
        loc="upper left",
        bbox_to_anchor=(0, -0.08),
        frameon=False,
    )
    ax = axs[1]
    states = np.array([[1, 1, 0, 0], [1, 1, 0, 0], [0, -1, -1, -1], [0, -1, -1, -1]])
    colors = {1: "#d9ebe4", 0: "#eed9d3", -1: "#eeeeee"}
    texts = [
        ["PASS", "PASS", "FAIL\n(5.00, 5.05]", "FAIL\n(6.60, 6.65]"],
        ["PASS", "PASS", "FAIL\n(6.15, 6.20]", "FAIL\n4.900 ms"],
        ["FAIL\n(9.45, 9.50]", "Not run", "Not run", "Not run"],
        ["FAIL\n2.240 ms", "Not run", "Not run", "Not run"],
    ]
    for i in range(4):
        for j in range(4):
            ax.add_patch(
                plt.Rectangle(
                    (j - 0.48, i - 0.48), 0.96, 0.96, facecolor=colors[states[i, j]], edgecolor="white"
                )
            )
            ax.text(j, i, texts[i][j], ha="center", va="center", fontsize=8.5)
    ax.set(
        xlim=(-0.6, 3.6),
        ylim=(3.6, -0.6),
        xticks=range(4),
        xticklabels=["Preparation", "OFF", "S", "O"],
        yticks=range(4),
        yticklabels=["seed1101", "seed1102", "seed1103", "seed1104"],
    )
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title("(b) Original 10 µs / 12-round outcomes\nFailure times in ms", loc="left")
    ax.text(
        0.5,
        -0.17,
        "Recovered complete trajectories: 0\nComplete OFF/S/O seed families: 0 / 4",
        transform=ax.transAxes,
        ha="center",
        fontsize=10,
    )
    save(
        fig,
        0,
        "Intervals are verified 50-µs brackets; 2.240 and 4.900 ms are step-resolved diagnostics. All four retry prefixes match the original batch.",
    )


def sphere_cloud(ax, positions, center, title, highlight=(28, 31), edges=False):
    x = (np.asarray(positions) - center) * 1e6
    radius = 1.78
    u = np.linspace(0, 2 * np.pi, 21)
    v = np.linspace(0, np.pi, 13)
    sx = radius * np.outer(np.cos(u), np.sin(v))
    sy = radius * np.outer(np.sin(u), np.sin(v))
    sz = radius * np.outer(np.ones_like(u), np.cos(v))
    for i, p in enumerate(x):
        color = "#b9723f" if i in highlight else "#7194b5"
        ax.plot_surface(
            p[0] + sx,
            p[1] + sy,
            p[2] + sz,
            color=color,
            alpha=0.88 if i in highlight else 0.66,
            linewidth=0,
            antialiased=True,
            shade=True,
            rasterized=True,
        )
        if i in highlight:
            label_position = np.array([-9.0, -3.0, 9.0]) if i == highlight[0] else np.array([7.0, 4.0, 10.0])
            ax.plot(*np.vstack([p, label_position]).T, color="#444444", lw=0.7, zorder=999)
            ax.text(
                *label_position,
                "ID " + str(i),
                fontsize=9,
                color="#161616",
                zorder=1000,
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=0.3),
            )
    if edges:
        d = np.linalg.norm(x[:, None] - x[None, :], axis=-1)
        for i, j in zip(*np.where(np.triu(d <= 2.1 * radius, 1))):
            ax.plot(*x[[i, j]].T, color="#555555", lw=0.7, alpha=0.75)
    ax.set_box_aspect([1, 1, 1])
    ax.view_init(elev=24, azim=-53)
    for axis in "xyz":
        getattr(ax, "set_" + axis + "lim")(-12, 12)
        getattr(ax, "set_" + axis + "label")(axis.upper() + " (µm)", labelpad=0)
        getattr(ax, "set_" + axis + "ticks")([-10, 0, 10])
    ax.tick_params(pad=0, labelsize=8)
    ax.set_title(title, pad=7)


def mode_traces():
    folders = [
        BASE / "contact_v1/O_1000_1102_C0",
        *[BASE / "alternatives_v1" / ("O_1000_1102_C" + str(i)) for i in range(1, 4)],
    ]
    return [trace(p) for p in folders]


def figure2():
    data = mode_traces()
    # Each scheme reached t=4.900 ms from the SAME original initial state;
    # the current internal arrangements have diverged and are not held fixed.
    selected = [min(d["traces"], key=lambda t: abs(t["attempted_end_s"] - 0.0049)) for d in data]
    assert all(abs(t["attempted_end_s"] - 0.0049) < 1e-12 for t in selected)
    fig = plt.figure(figsize=(14.5, 9.5))
    gs = fig.add_gridspec(2, 3, left=0.075, right=0.97, bottom=0.12, top=0.85, wspace=0.38, hspace=0.42)
    fig.suptitle(
        "F10E | Residual overlap and single-pass correction are different budgets", fontsize=15, y=0.965
    )
    fig.text(
        0.5,
        0.918,
        "O / seed1102 at 4.900 ms: same original state and controls, scheme-specific replay histories",
        ha="center",
        fontsize=10,
    )
    ax = fig.add_subplot(gs[0, 0])
    for i, t in enumerate(selected):
        values = [t["prediction"]["max_overlap_m"] * 1e9, *[r["max_overlap_m"] * 1e9 for r in t["rounds"]]]
        ax.plot(range(len(values)), values, "o-", ms=3, lw=1.6, color=COLORS[i], label=LABELS[i])
    ax.axhline(17.8, color="#943f43", ls="--", lw=1.2, label="Residual budget: 17.8 nm")
    ax.set(
        title="(a) Prediction (round 0) and PBD residual",
        xlabel="PBD round",
        ylabel="Maximum overlap (nm)",
        ylim=(0, 135),
        xticks=[0, 6, 12, 18, 24],
    )
    ax.grid(alpha=0.18)
    ax = fig.add_subplot(gs[0, 1])
    for i, t in enumerate(selected):
        ax.plot(
            [r["round"] for r in t["rounds"]],
            [r["max_correction_m"] * 1e9 for r in t["rounds"]],
            "-o",
            ms=3,
            color=COLORS[i],
        )
    ax.axhline(178, color="#943f43", ls="--", lw=1.2)
    ax.text(23, 182, "178 nm budget", ha="right", color="#943f43", fontsize=9)
    ax.set(
        title="(b) Single-pass correction",
        xlabel="PBD round",
        ylabel="Maximum correction (nm)",
        ylim=(0, 205),
        xticks=[0, 6, 12, 18, 24],
    )
    ax.grid(alpha=0.18)
    ax = fig.add_subplot(gs[0, 2])
    prep = trace(BASE / "contact_v1/prepare_1104_C0")["traces"][-1]
    ax.plot(
        range(13),
        [prep["prediction"]["max_overlap_m"] * 1e9, *[r["max_overlap_m"] * 1e9 for r in prep["rounds"]]],
        "o-",
        color=COLORS[0],
        ms=4,
    )
    ax.axhline(17.8, color="#943f43", ls="--")
    ax.set(
        title="(c) Preparation seed1104: step 224",
        xlabel="PBD round",
        ylabel="Maximum overlap (nm)",
        ylim=(0, 75),
        xticks=[0, 3, 6, 9, 12],
    )
    ax.text(
        0.46,
        0.92,
        "2.230 → 2.240 ms\nRejected: pair 24–28\nFinal residual: 17.831 nm",
        transform=ax.transAxes,
        va="top",
        fontsize=9,
    )
    ax.grid(alpha=0.18)
    ax = fig.add_subplot(gs[1, 0])
    t = selected[0]
    x = np.array(t["proposal_positions_m"])
    xx = (x - x.mean(0)) * 1e6
    for i, p in enumerate(xx):
        ax.add_patch(Circle(p[:2], 1.78, facecolor="#7194b5", edgecolor="#7194b5", alpha=0.10))
    dist = np.linalg.norm(xx[:, None] - xx[None, :], axis=-1)
    for i, j in zip(*np.where(np.triu(dist <= 2.1 * 1.78, 1))):
        ax.plot(*xx[[i, j], :2].T, lw=0.6, color="#707b85", alpha=0.45)
    for i, offset in [(28, (-7, -7)), (31, (5, 7))]:
        ax.add_patch(Circle(xx[i, :2], 1.78, facecolor="#bc7845", edgecolor="#733c17", alpha=0.8, zorder=5))
        ax.annotate(
            f"ID {i}",
            xy=xx[i, :2],
            xytext=offset,
            arrowprops=dict(arrowstyle="->", color="#333333"),
            fontsize=9,
            zorder=10,
        )
    ax.set(
        title="(d) Rejected contact graph: XY projection",
        xlabel="X (µm)",
        ylabel="Y (µm)",
        xlim=(-10, 10),
        ylim=(-10, 10),
    )
    ax.set_aspect("equal")
    ax.grid(alpha=0.16)
    ax = fig.add_subplot(gs[1, 1])
    analysis = js(BASE / "analysis/CONTACT_ANALYSIS.json")[0]
    values = np.array(analysis["correction_path_per_id_nm"])
    ax.bar(range(32), values, color=[COLORS[2] if i in (28, 31) else "#849caf" for i in range(32)])
    ax.set(
        title="(e) PBD correction path in the rejected step",
        xlabel="Fixed particle ID",
        ylabel="Sum of round displacements (nm)",
        xticks=[0, 8, 16, 24, 31],
    )
    ax.grid(axis="y", alpha=0.18)
    ax = fig.add_subplot(gs[1, 2])
    for i, d in enumerate(data):
        ts = np.array([t["attempted_end_s"] for t in d["traces"]]) * 1e3
        ax.plot(ts, [t["final"]["max_overlap_m"] * 1e9 for t in d["traces"]], color=COLORS[i], lw=1.8)
    ax.axhline(17.8, color="#943f43", ls="--")
    ax.plot(4.9, selected[0]["final"]["max_overlap_m"] * 1e9, "x", color="#943f43", ms=9, mew=2)
    ax.set(
        title="(f) End-of-step overlap near rejection",
        xlabel="Physical time (ms)",
        ylabel="Final residual (nm)",
        xlim=(4.84, 5.005),
        ylim=(0, 21),
    )
    ax.grid(alpha=0.18)
    fig.legend(
        handles=[plt.Line2D([0], [0], color=c, lw=2, label=label) for c, label in zip(COLORS, LABELS)],
        loc="lower center",
        bbox_to_anchor=(0.5, 0.051),
        ncol=4,
        frameon=False,
    )
    save(
        fig,
        1,
        "True sphere radius: 1.78 µm; links denote geometric proximity, not fusion. The rejected proposal is explicitly marked and is not an accepted trajectory.",
    )


def load_series(path):
    with h5py.File(path / "trajectory.h5") as h:
        return dict(
            t=h["time_s"][:] * 1e3,
            x=h["positions_m"][:, 0],
            rg=h["metrics/global_rg_m"][:, 0] * 1e6,
            aspect=h["morphology/fixed_aspect_ratio"][:, 0],
            rmsd=h["morphology/aligned_internal_rmsd_m"][:, 0] * 1e6,
            labels=h["metrics/labels"][:, 0],
            center=h["metrics/global_centroid_m"][:, 0],
        )


def figure3():
    paths = [
        BASE / "contact_v1/O_1000_1102_C0",
        BASE / "extension_v1/O_1000_1102_C1",
        BASE / "alternatives_v1/O_1000_1102_C2",
        BASE / "alternatives_v1/O_1000_1102_C3",
    ]
    series = [load_series(p) for p in paths]
    fig = plt.figure(figsize=(14, 11))
    gs = fig.add_gridspec(3, 3, left=0.07, right=0.965, bottom=0.13, top=0.87, hspace=0.42, wspace=0.36)
    fig.suptitle(
        "F10E | Contact-budget success does not establish particle-level shape fidelity", fontsize=15, y=0.963
    )
    fig.text(
        0.5,
        0.922,
        "O / seed1102 · fixed IDs · common accepted prefix ends at 4.850 ms · C0 has no reference beyond its failure",
        ha="center",
        fontsize=10,
    )
    common = len(series[0]["t"]) - 1
    center = series[0]["x"][common].mean(0)
    for j, i in enumerate([0, 1, 3]):
        ax = fig.add_subplot(gs[0, j], projection="3d")
        sphere_cloud(ax, series[i]["x"][common], center, f"({chr(97 + j)}) {LABELS[i]}\nat 4.850 ms")
    axes = [fig.add_subplot(gs[r, c]) for r in [1, 2] for c in range(3)]
    for i, s in enumerate(series):
        c = COLORS[i]
        axes[0].plot(s["t"], (s["center"][:, 2] - s["center"][0, 2]) * 1e6, color=c)
        axes[1].plot(s["t"], s["rg"], color=c)
        axes[2].plot(s["t"], s["aspect"], color=c)
        axes[3].plot(s["t"], s["rmsd"], color=c)
        axes[4].plot(s["t"], [len(np.unique(v)) for v in s["labels"]], color=c, ls=["-", "--", ":", "-."][i])
        if i:
            delta = np.linalg.norm(s["x"][: common + 1] - series[0]["x"], axis=-1).max(1) * 1e6
            axes[5].plot(s["t"][: common + 1], delta, color=c)
    desc = [
        ("(d) Centroid Z displacement", "Displacement (µm)"),
        ("(e) Radius of gyration", "Rg (µm)"),
        ("(f) Fixed-ID aspect ratio", "Major / minor axis"),
        ("(g) Internal Kabsch RMSD", "RMSD to initial state (µm)"),
        ("(h) Geometric components", "Component count"),
        ("(i) Fixed-ID sensitivity to C0", "Largest position difference (µm)"),
    ]
    for ax, (title, ylabel) in zip(axes, desc):
        ax.set(title=title, xlabel="Physical time (ms)", ylabel=ylabel)
        ax.axvline(4.85, color="#777777", ls=":", lw=1)
        ax.grid(alpha=0.18)
    axes[4].set(ylim=(0.8, 1.2), yticks=[1])
    axes[5].set_xlim(0, 4.85)
    fig.legend(
        handles=[plt.Line2D([0], [0], color=c, lw=2, label=label) for c, label in zip(COLORS, LABELS)],
        loc="lower center",
        bbox_to_anchor=(0.5, 0.057),
        ncol=4,
        frameon=False,
    )
    save(
        fig,
        2,
        "All spheres have their physical radius (1.78 µm); equal spatial axes and one shared origin. C2/C3 end at the registered 5 ms window; only C1 was extended to 10 ms.",
    )


def figure4():
    d = js(BASE / "analysis/EFFECTIVE_YIELD.json")
    rates = d["body_rates"]
    fig, axs = plt.subplots(2, 2, figsize=(13, 9))
    fig.subplots_adjust(left=0.085, right=0.965, bottom=0.14, top=0.855, hspace=0.47, wspace=0.3)
    fig.suptitle("F10E | Effective yield includes every failed attempt and replay", fontsize=15, y=0.96)
    fig.text(
        0.5,
        0.91,
        "Four independent initial seed families; controls and numerical variants never create new independent samples",
        ha="center",
        fontsize=10,
    )
    names = ["F10D\noriginal", "+ isolation\n(original mode)", "+ numerical\ndiagnostics / C1"]
    xx = np.arange(3)
    ax = axs[0, 0]
    ax.bar(xx, [2, 2, 2], color=COLORS[0], width=0.55, label="Original-mode complete controls")
    ax.bar(
        xx,
        [0, 0, d["alternative_only_complete_trajectories"]],
        bottom=2,
        color=COLORS[1],
        width=0.55,
        label="C1-only complete control",
    )
    for x, row in zip(xx, rates):
        ax.text(
            x,
            row["valid_completed_trajectories"] + 0.08,
            str(row["valid_completed_trajectories"]),
            ha="center",
        )
    ax.set(
        title="(a) Complete 10 ms control trajectories",
        xticks=xx,
        xticklabels=names,
        ylabel="Numerically qualified trajectories",
        ylim=(0, 4.3),
        yticks=range(5),
    )
    ax.legend(frameon=False, loc="upper left", fontsize=8)
    ax = axs[0, 1]
    ax.bar(xx, [0, 0, 0], color=COLORS[1])
    for x in xx:
        ax.text(x, 0.08, "0 / 4", ha="center", fontsize=14)
    ax.set(
        title="(b) Complete OFF/S/O initial-state families",
        xticks=xx,
        xticklabels=names,
        ylabel="Same-mode complete triplets",
        ylim=(0, 1),
        yticks=[0, 1],
    )
    ax.text(
        0.5,
        0.67,
        "OFF remains at 10 µs / 12 rounds.\nOne C1 O branch is not a matched triplet.",
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=10,
    )
    ax = axs[1, 0]
    original = np.full(3, d["original_generation_s"])
    retries = np.array([0, d["isolation_production_replay_s"], d["isolation_production_replay_s"]])
    extra = np.array([0, 0, d["extension_production_s"]])
    other = np.array([0, 0, d["all_attempts_body_s"] - original[2] - retries[2] - extra[2]])
    bottom = np.zeros(3)
    for vals, col, label in [
        (original, COLORS[0], "F10D: all preparation / branch attempts"),
        (retries, "#bc7845", "Four failed singleton replays + static input"),
        (extra, COLORS[1], "One full C1 replay + static input"),
        (other, "#b6a7c2", "Diagnostics, checks and observer output"),
    ]:
        ax.bar(xx, vals, bottom=bottom, color=col, width=0.55, label=label)
        bottom += vals
    for x, v in zip(xx, bottom):
        ax.text(x, v + 0.4, f"{v:.2f} s", ha="center")
    ax.set(
        title="(c) Retained workflow-body cost",
        xticks=xx,
        xticklabels=names,
        ylabel="Elapsed time (s)",
        ylim=(0, 24),
    )
    ax.legend(fontsize=7.5, frameon=False, loc="upper left")
    ax = axs[1, 1]
    values = [r["valid_completed_trajectories_per_gpu_hour"] for r in rates]
    ax.bar(xx, values, color=[COLORS[0], "#bc7845", COLORS[1]], width=0.55)
    for x, v in zip(xx, values):
        ax.text(x, v + 45, f"{v:.0f}", ha="center")
    ax.set(
        title="(d) Valid controls per reserved GPU hour",
        xticks=xx,
        xticklabels=names,
        ylabel="10 ms trajectories / hour",
        ylim=(0, 2250),
    )
    lo, hi = d["all_process_rate_interval_per_gpu_hour"]
    ax.text(
        0.03,
        0.94,
        f"All costs + process startup: {lo:.0f}–{hi:.0f} / h\nTriplets: 0 / h in every workflow",
        transform=ax.transAxes,
        va="top",
        fontsize=9,
    )
    save(
        fig,
        3,
        "Historical cost accounting, not a paired speed benchmark or sustained production estimate. Device reservation includes CPU work; CUDA spans are not pure kernel time.",
    )


def main():
    OUT.mkdir(exist_ok=True)
    style()
    figure1()
    figure2()
    figure3()
    figure4()
    inputs = {
        str(p.relative_to(ROOT)): sha(p)
        for folder in [
            CAMPAIGN,
            BASE / "isolation_v1",
            BASE / "contact_v1",
            BASE / "alternatives_v1",
            BASE / "extension_v1",
            BASE / "analysis",
        ]
        for p in sorted(folder.rglob("*"))
        if p.is_file() and p.suffix in [".h5", ".json", ".gz"]
    }
    write(
        OUT / "FIGURE_PROVENANCE.json",
        dict(
            source=identity(),
            inputs=inputs,
            figures={p.name: sha(p) for p in OUT.iterdir() if p.suffix in [".png", ".pdf"]},
            physical_sphere_radius_m=1.78e-6,
            synthetic_positions=False,
            training_ready=False,
        ),
    )


if __name__ == "__main__":
    main()
