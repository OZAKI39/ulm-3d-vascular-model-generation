"""Four English PNG/PDF groups from archived physical steps; no synthetic frames."""

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from .common import ROOT, BASE, js, sha, write, identity
from .analyze import gz, read_steps
from mevo_mb.fast_sim.f10g_contact import shape

OUT = BASE / "figures"
C = ["#405a77", "#91629c", "#2d8b83", "#b47932", "#a84246"]
plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.labelcolor": "#222733",
        "text.color": "#222733",
        "axes.edgecolor": "#434953",
        "xtick.color": "#30343b",
        "ytick.color": "#30343b",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,
        "savefig.facecolor": "white",
    }
)


def title(ax, s):
    ax.set_title(s, loc="left", fontweight="bold", fontsize=11, pad=10)


def finish(fig, name, caption):
    fig.get_layout_engine().set(rect=(0, 0.065, 1, 0.96))
    fig.text(0.025, 0.012, caption, fontsize=9, color="#454953")
    fig.savefig(OUT / (name + ".png"), dpi=300, bbox_inches="tight")
    fig.savefig(OUT / (name + ".pdf"), bbox_inches="tight")
    plt.close(fig)


def lubrication(L, T):
    fig, aa = plt.subplots(2, 3, figsize=(15, 9), layout="constrained")
    for j, c in enumerate(L["cases"]):
        rows = c["rows"]
        name = ["O / 1102", "Preparation / 1104"][j]
        aa[0, 0].plot(
            range(4),
            [100 * r["residuals"]["balanced"]["normalized_preconditioned"] for r in rows],
            "o-",
            color=C[j],
            label=name,
        )
        aa[0, 1].plot(
            range(4),
            [-min(r["relative_normal_velocity_m_s"]) * 1e3 for r in rows],
            "o-",
            color=C[j],
            label=name,
        )
        aa[0, 2].plot(
            range(4),
            [r["predicted_overlap_m"] * 1e9 for r in rows],
            "o-",
            color=C[j],
            label=name + " prediction",
        )
        aa[0, 2].plot(
            range(4), [r["final_overlap_m"] * 1e9 for r in rows], "s--", color=C[j], label=name + " PBD12"
        )
    for ax in aa[0]:
        ax.set(
            xticks=range(4),
            xticklabels=["2", "4", "8", "Reference"],
            xlabel="Jacobi iterations / direct reference",
        )
        ax.legend(fontsize=8)
        ax.grid(alpha=0.18)
    title(aa[0, 0], "a  Frozen-state algebraic residual (F10G reused)")
    aa[0, 0].set_ylabel("Block-preconditioned residual / free speed (%)")
    title(aa[0, 1], "b  Frozen-state maximum approach speed")
    aa[0, 1].set_ylabel("Normal approach speed (mm/s)")
    title(aa[0, 2], "c  Frozen-state overlap and original budget")
    aa[0, 2].set_ylabel("Maximum overlap (nm)")
    aa[0, 2].axhline(17.8, color=C[4], ls=":", label="17.8 nm")
    aa[0, 2].legend(fontsize=8)
    for row in T["rows"]:
        if "/short_v1/" not in row["path"]:
            continue
        j = 0 if "O_1000" in row["path"] else 1
        rr = row["records"]
        time = np.array([s["time_s"] for s in rr]) * 1e3
        aa[1, 0].plot(
            time,
            [s["normalized_preconditioned"] * 100 for s in rr],
            color=C[j],
            label=["O / 1102", "Preparation / 1104"][j],
        )
        aa[1, 1].plot(time, [s["normalized_velocity_error"] * 100 for s in rr], color=C[j])
        aa[1, 2].plot(time, [s["max_approach_speed_m_s"] * 1e3 for s in rr], color=C[j])
    for ax in aa[1]:
        ax.set_xlabel("Actual accepted physical time (ms)")
        ax.grid(alpha=0.18)
    aa[1, 0].legend(fontsize=8)
    title(aa[1, 0], "d  Actual LUB4 trajectory residual")
    aa[1, 0].set_ylabel("Block-preconditioned residual (%)")
    title(aa[1, 1], "e  Actual LUB4 discrete-reference difference")
    aa[1, 1].set_ylabel("Max. ID velocity error / free speed (%)")
    title(aa[1, 2], "f  Actual LUB4 normal approach speed")
    aa[1, 2].set_ylabel("Maximum normal approach (mm/s)")
    finish(
        fig,
        "Fig_F10H_1_lubrication_contact",
        "Top: identical saved positions, forces and 10 µs step; reused 2/4/8/direct resistance references. Bottom: genuinely evolved LUB4 trajectories.\nAlgebraic residuals refer to the existing approximate resistance model; they are not measured physical fluid-velocity errors.",
    )


def windows(T):
    fig, aa = plt.subplots(2, 3, figsize=(15, 9), layout="constrained")
    for j, case in enumerate(["O_1000_1102", "prepare_1104"]):
        old = (
            gz(ROOT / f"results/fast_pzt/f10e/contact_v1/{case}_C0/SUBSTEP_LEDGER.json.gz")["steps"]
            if (ROOT / f"results/fast_pzt/f10e/contact_v1/{case}_C0/SUBSTEP_LEDGER.json.gz").exists()
            else None
        )
        configs = []
        # Old last 32 physical steps and original rejection are exact saved evidence.
        old = gz(ROOT / f"results/fast_pzt/f10e/contact_v1/{case}_C0/CONTACT_TRACE.json.gz")["traces"]
        configs.append(("Original LUB2", old, C[0], "--", True))
        for row in T["trials"]:
            if row["case"] != case:
                continue
            if case == "O_1000_1102" and row["mode"].endswith("DT10") and row["target_s"] == 0.005:
                continue
            d = gz(ROOT / row["path"] / "SUBSTEP_LEDGER.json.gz")["steps"]
            configs.append(
                (
                    "LUB4 / " + ("5 µs" if row["mode"].endswith("DT5") else "10 µs"),
                    d,
                    C[2] if row["mode"].endswith("DT5") else C[1],
                    "-",
                    False,
                )
            )
        for name, steps, color, ls, isold in configs:
            times = np.array([s["attempted_end_s"] for s in steps]) * 1e3
            ov = [(s["final"]["max_overlap_m"] if isold else s["final_overlap_m"]) * 1e9 for s in steps]
            corr = [
                max(r["max_correction_m"] for r in s["rounds"]) * 1e9
                if isold
                else max(s["round_max_correction_m"]) * 1e9
                for s in steps
            ]
            aa[j, 0].plot(times, ov, ls, color=color, label=name)
            aa[j, 1].plot(times, corr, ls, color=color, label=name)
            if not steps[-1]["accepted"]:
                aa[j, 0].plot(times[-1], ov[-1], "x", color=color, ms=9, mew=2)
        aa[j, 0].axhline(17.8, color=C[4], ls=":", label="17.8 nm budget")
        aa[j, 1].axhline(178, color=C[4], ls=":", label="178 nm budget")
        aa[j, 1].set_ylim(0, 190)
        for row in T["rows"]:
            if case not in row["path"] or "/short_v1/" not in row["path"]:
                continue
            sh = row["per_step_shapes"]
            tm = np.arange(len(sh)) * 1e-5 * 1e3
            center = np.array([s["centroid_m"] for s in sh])
            disp = np.linalg.norm(center - center[0], axis=1) * 1e6
            aa[j, 2].plot(tm, disp, color=C[1], label="Centroid displacement")
            ax = aa[j, 2].twinx()
            ax.plot(tm, [s["Rg_m"] * 1e6 for s in sh], color=C[2], label="Rg")
            ax.set_ylabel("Rg (µm)", color=C[2])
            ax.spines["top"].set_visible(False)
        name = ["O / 1102", "Preparation / 1104"][j]
        for i, sub in enumerate(
            ["Contact overlap", "Single PBD-round correction", "Accepted LUB4 / 10 µs motion"]
        ):
            title(aa[j, i], f"{chr(97 + 3 * j + i)}  {name}: {sub}")
            aa[j, i].set_xlabel("Physical time (ms)")
            aa[j, i].grid(alpha=0.18)
        aa[j, 0].set_ylabel("Maximum overlap (nm)")
        aa[j, 1].set_ylabel("Maximum correction (nm)")
        aa[j, 2].set_ylabel("Centroid displacement (µm)")
        aa[j, 0].legend(fontsize=8)
        aa[j, 1].legend(fontsize=8)
        aa[j, 2].legend(fontsize=8)
    finish(
        fig,
        "Fig_F10H_2_failure_window_replay",
        "Crosses are rejected proposals, not accepted states. Original failures: O 4.900 ms; preparation 2.240 ms.\nLUB4 / 10 µs passes both short windows, then rejects O at 6.370 ms; no accepted continuation to 10 ms is fabricated.",
    )


def spheres(ax, x, origin, lim):
    y = (x - origin) * 1e6
    r = 1.78
    u = np.linspace(0, 2 * np.pi, 18)
    v = np.linspace(0, np.pi, 13)
    sx = np.outer(np.cos(u), np.sin(v))
    sy = np.outer(np.sin(u), np.sin(v))
    sz = np.outer(np.ones_like(u), np.cos(v))
    for i, j in shape(x, 1.78e-6)["edges"]:
        ax.plot(*y[[i, j]].T, color="#485566", alpha=0.4, lw=0.6)
    for i, p in enumerate(y):
        ax.plot_surface(
            p[0] + r * sx,
            p[1] + r * sy,
            p[2] + r * sz,
            color=plt.cm.viridis(i / 31),
            linewidth=0,
            antialiased=True,
            shade=True,
            alpha=0.95,
        )
    ax.set(
        xlim=(-lim, lim),
        ylim=(-lim, lim),
        zlim=(-lim, lim),
        xlabel="x (µm)",
        ylabel="y (µm)",
        zlabel="z (µm)",
    )
    ax.zaxis.labelpad = 0
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=23, azim=-53)
    ax.tick_params(labelsize=8, pad=0)


def geometry(B):
    E = js(BASE / "geometry_v1/GEOMETRY_EVENTS.json")
    S = js(BASE / "geometry_v1/SMALL_SYSTEMS.json")
    D = js(BASE / "geometry_v1/DENSE_STATIC.json")
    fig = plt.figure(figsize=(16, 10.5), layout="constrained")
    fig.get_layout_engine().set(wspace=0.14, hspace=0.1)
    gs = fig.add_gridspec(2, 3)
    ax = fig.add_subplot(gs[0, 0])
    events = E["events"]
    ax.scatter(
        [e["actual_gap_m"] * 1e9 for e in events],
        [e["old_force_N"] * 1e12 for e in events],
        color=C[1],
        label="Frozen normals: separated pair",
    )
    ax.scatter(
        [e["new_gap_m"] * 1e9 for e in events],
        [e["new_force_N"] * 1e12 for e in events],
        color=C[2],
        marker="x",
        label="Geometric re-solve, same input",
    )
    title(ax, "a  Archived separated-pair reactions")
    ax.set(xlabel="True surface gap of separated pair (nm)", ylabel="Pair reaction (pN)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.18)
    ax = fig.add_subplot(gs[0, 1])
    rows = S["rows"]
    ax.bar(
        np.arange(len(rows)),
        [s["reference_error_m"] * 1e12 for s in rows],
        color=[C[0] if s["bubble_count"] == 2 else C[2] for s in rows],
    )
    title(ax, "b  Independent nonlinear reference")
    ax.set(xlabel="12 two-bubble + 6 four-bubble checks", ylabel="Max. position difference (pm)")
    ax.grid(axis="y", alpha=0.18)
    ax = fig.add_subplot(gs[0, 2])
    xx = np.arange(2)
    old = js(ROOT / "results/fast_pzt/f10g/prototype_v1/DENSE_STATIC.json")["rows"]
    old = [
        next(v for v in old if v["case"] == s["name"] and v["mode"] == "B2" and v["multiplier"] == 2)
        for s in D["rows"]
    ]
    ax.bar(
        xx - 0.18,
        [max(s["result"]["geometric_overlap_m"]) * 1e9 for s in old],
        width=0.36,
        color=C[1],
        label="Frozen normals",
    )
    ax.bar(
        xx + 0.18,
        [max(s["result"]["geometric_overlap_m"]) * 1e9 for s in D["rows"]],
        width=0.36,
        color=C[2],
        label="Updated geometry",
    )
    ax.axhline(17.8, color=C[4], ls=":", label="17.8 nm budget")
    ax.set(
        xticks=xx,
        xticklabels=["O / 1102", "Preparation / 1104"],
        ylabel="Same-state final overlap (nm)",
        ylim=(16.8, 18.1),
    )
    title(ax, "c  Saved32 static states, fixed 2kc")
    ax.legend(fontsize=8)
    states = [read_steps(ROOT / t["path"])[1][-1] for t in B["trials"]]
    origin = states[0].mean(0)
    lim = max(abs(s - origin).max() for s in states) * 1e6 + 2.4
    for i, x in enumerate(states):
        ax = fig.add_subplot(gs[1, i], projection="3d")
        spheres(ax, x, origin, lim)
        from matplotlib.colors import Normalize

        barax = ax.inset_axes([0.06, 0.89, 0.50, 0.022])
        cb = fig.colorbar(
            plt.cm.ScalarMappable(norm=Normalize(0, 31), cmap="viridis"),
            cax=barax,
            orientation="horizontal",
            ticks=[0, 15, 31],
        )
        cb.ax.text(1.07, 0.5, "Bubble ID", fontsize=8, va="center", transform=cb.ax.transAxes)
        cb.ax.tick_params(labelsize=7, pad=1)
        title(ax, f"{chr(100 + i)}  Geometric B2 / {10 if i == 0 else 5} µs, t = 5 ms")
    ax = fig.add_subplot(gs[1, 2])
    graphs = [np.linalg.norm(x[:, None] - x[None, :], axis=-1) <= 2.1 * 1.78e-6 for x in states]
    code = graphs[0].astype(int) + 2 * graphs[1].astype(int)
    np.fill_diagonal(code, 0)
    cmap = ListedColormap(["#ffffff", C[1], C[2], "#c0c8d0"])
    ax.imshow(code, cmap=cmap, vmin=0, vmax=3, interpolation="nearest", origin="lower")
    title(ax, "f  Fixed-ID connectivity at 5 ms")
    ax.set(xlabel="Bubble ID", ylabel="Bubble ID", xticks=[0, 8, 16, 24, 31], yticks=[0, 8, 16, 24, 31])
    from matplotlib.patches import Patch

    ax.legend(
        handles=[
            Patch(color=C[1], label="10 µs only"),
            Patch(color=C[2], label="5 µs only"),
            Patch(color="#c0c8d0", label="Both"),
        ],
        fontsize=8,
        loc="upper left",
        bbox_to_anchor=(1.01, 1),
    )
    finish(
        fig,
        "Fig_F10H_3_geometry_morphology",
        "Fixed kc = 0.0464132 N/m; an uncalibrated numerical candidate. Geometric spring reactions use current distance and normal.\nTrue radius = 1.78 µm; equal XYZ scales and fixed ID colors. Connections use d <= 2.1R; matrix preserves every fixed-ID edge.\nSeventeen same-input pair-event re-solves and evolved 5 ms states are distinct evidence. No sphere shrinkage, smoothing or fabricated trajectory.",
    )


def admission(M, B, K):
    fig, aa = plt.subplots(2, 3, figsize=(15, 9.5), layout="constrained")
    comparisons = M["comparisons"] + B["comparisons"]
    names = ["LUB4 O", "LUB4 prep.", "Geometric B2 O"]
    for j, c in enumerate(comparisons):
        ss = c["series"]
        tm = np.array([s["time_s"] for s in ss]) * 1e3
        aa[0, 0].plot(tm, [s["max_fixed_ID_m"] * 1e6 for s in ss], color=C[j], label=names[j] + " max ID")
        aa[0, 0].plot(tm, [s["Kabsch_RMSD_m"] * 1e6 for s in ss], "--", color=C[j])
        aa[0, 1].bar(
            np.arange(5) + (j - 1) * 0.25,
            [
                c["maxima"][v] / lim
                for v, lim in [
                    ("max_fixed_ID_m", 178e-9),
                    ("Kabsch_RMSD_m", 89e-9),
                    ("Rg_relative", 0.01),
                    ("aspect_relative", 0.02),
                    ("tensor_relative", 0.02),
                ]
            ],
            width=0.25,
            color=C[j],
            label=names[j],
        )
    title(aa[0, 0], "a  Same-model 10 / 5 µs shape sensitivity")
    aa[0, 0].set(xlabel="Common accepted physical time (ms)", ylabel="Difference (µm)")
    aa[0, 0].legend(fontsize=7)
    aa[0, 0].grid(alpha=0.18)
    title(aa[0, 1], "b  Preregistered shape thresholds")
    aa[0, 1].axhline(1, color=C[4], ls=":")
    aa[0, 1].set(
        xticks=range(5),
        xticklabels=["ID", "Kabsch", "Rg", "Aspect", "Tensor"],
        ylabel="Maximum difference / threshold",
        yscale="log",
    )
    aa[0, 1].legend(fontsize=7)
    E = js(BASE / "analysis/EDGE_AUDIT.json")
    for j, w in enumerate(E["worst_pair_window"]):
        aa[0, 2].plot(
            np.array(w["time_s"]) * 1e3,
            np.array(w["surface_gap_m"]) * 1e9,
            "o-" if j == 0 else ".-",
            ms=3,
            color=C[j],
            label="B2 / " + ("10 µs" if j == 0 else "5 µs"),
        )
    aa[0, 2].axhline(178, color=C[4], ls=":", label="Connectivity threshold")
    for key, col in [("time_DT10_s", C[0]), ("time_DT5_s", C[1])]:
        aa[0, 2].axvline(E["worst"][key] * 1e3, color=col, ls="--", lw=0.8)
    title(aa[0, 2], "c  B2 edge timing: fixed IDs 7 and 24")
    aa[0, 2].set(xlabel="Physical time (ms)", ylabel="True surface gap (nm)")
    aa[0, 2].text(
        0.97, 0.06, "Event lag: 130 µs\nGate: 10 µs", transform=aa[0, 2].transAxes, ha="right", fontsize=9
    )
    aa[0, 2].legend(fontsize=7)
    aa[0, 2].grid(alpha=0.18)
    med = K["matched_diagnostic_medians"]
    labels = ["LUB2", "LUB4"]
    modes = ["M0_LUB2_PBD12_DT10", "M0_LUB4_PBD12_DT10"]
    bottom = np.zeros(2)
    for key, color, name in [
        ("integration_quality_wall_s", C[0], "Integration + quality"),
        ("observer_output_s", C[1], "Diagnostic output"),
        ("hdf5_output_s", C[2], "HDF5 output"),
    ]:
        vals = np.array([med[m][key] for m in modes])
        aa[1, 0].bar(range(2), vals, bottom=bottom, color=color, label=name)
        bottom += vals
    title(aa[1, 0], "d  Matched 0.2 ms diagnostic cost")
    aa[1, 0].set(xticks=range(2), xticklabels=labels, ylabel="Median wall time (s), three repeats")
    aa[1, 0].legend(fontsize=8)
    rr = [
        next(
            r
            for r in K["actual_science"]
            if r["case"] == "O_1000_1102" and r["mode"] == m and r["target_s"] == 0.005
        )
        for m in ["M0_LUB4_PBD12_DT10", "B2_GEOM_K2_DT10", "B2_GEOM_K2_DT5"]
    ]
    bottom = np.zeros(3)
    for key, col, label in [
        ("integration_quality_wall_s", C[0], "Integration + quality"),
        ("observer_output_s", C[1], "Detailed output"),
    ]:
        vals = [r["timing"][key] for r in rr]
        aa[1, 1].bar(range(3), vals, bottom=bottom, color=col, label=label)
        bottom += vals
    title(aa[1, 1], "e  Actual 5 ms diagnostic cost")
    aa[1, 1].set(
        xticks=range(3),
        xticklabels=["LUB4\n10 µs", "B2 geom.\n10 µs", "B2 geom.\n5 µs"],
        ylabel="Wall time (s)",
    )
    aa[1, 1].legend(fontsize=8)
    ax = aa[1, 2]
    ax.axis("off")
    title(ax, "f  Separate admission decisions")
    lines = [
        "LUB4: short windows PASS; 10 ms FAIL",
        "LUB4 lubrication: NOT FULLY CONVERGED",
        "LUB4 morphology: NOT CONVERGED",
        "B2 geometric unilateral spring: PASS",
        "B2: 5 ms contact PASS (10 / 5 µs)",
        "B2 positions / shape: PASS; edge timing: FAIL",
        "Production batch: NOT ADMITTED",
        "Physics: UNCALIBRATED",
        "AI training_ready = false",
    ]
    for i, s in enumerate(lines):
        ax.text(0, 0.96 - i * 0.1, s, transform=ax.transAxes, va="top", fontsize=9.5, color="#303743")
    finish(
        fig,
        "Fig_F10H_4_accuracy_cost_admission",
        "Shape thresholds are engineering comparisons, not experimental accuracy. Solid: fixed-ID maximum; dashed: Kabsch RMSD.\nMatched short timings share outputs and checks; B2 detailed outputs differ from PBD. CPU/diagnostic costs and CUDA spans are not pure-kernel production benchmarks.",
    )


def main():
    OUT.mkdir(exist_ok=True)
    L = js(ROOT / "results/fast_pzt/f10g/lubrication_v1/RESULTS.json")
    T = js(BASE / "analysis/LUBRICATION_TRAJECTORIES.json")
    M = js(BASE / "analysis/MORPHOLOGY.json")
    B = js(BASE / "analysis/B2_GEOMETRY_TRAJECTORIES.json")
    K = js(BASE / "analysis/COST_AND_YIELD.json")
    lubrication(L, T)
    windows(T)
    geometry(B)
    admission(M, B, K)
    sources = [
        ROOT / "results/fast_pzt/f10g/lubrication_v1/RESULTS.json",
        ROOT / "results/fast_pzt/f10g/prototype_v1/DENSE_STATIC.json",
        *[
            BASE / "analysis" / name
            for name in [
                "LUBRICATION_TRAJECTORIES.json",
                "MORPHOLOGY.json",
                "B2_GEOMETRY_TRAJECTORIES.json",
                "BETWEEN_MODELS.json",
                "COST_AND_YIELD.json",
                "FINAL_GATE.json",
                "EDGE_AUDIT.json",
            ]
        ],
        *sorted((BASE / "geometry_v1").glob("*.json")),
        *sorted(BASE.glob("*/*/trajectory.h5")),
        *sorted(BASE.glob("*/*/ALL_STEPS.json.gz")),
        *sorted(BASE.glob("*/*/SUBSTEP_LEDGER.json.gz")),
        *sorted(ROOT.glob("results/fast_pzt/f10e/contact_v1/*_C0/CONTACT_TRACE.json.gz")),
    ]
    write(
        OUT / "FIGURE_PROVENANCE.json",
        dict(
            source=identity(),
            files={str(p.relative_to(ROOT)): sha(p) for p in sources},
            radius_m=1.78e-6,
            interpolation="none",
            PNG_dpi=300,
            PDF_fonttype=42,
            training_ready=False,
        ),
    )
    print("Four PNG/PDF groups generated")


if __name__ == "__main__":
    main()
