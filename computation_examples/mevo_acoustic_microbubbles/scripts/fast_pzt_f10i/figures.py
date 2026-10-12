"""Four English publication diagnostics, entirely from archived accepted states."""

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from .common import BASE, js, sha, write, identity

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


def title(ax, value):
    ax.set_title(value, loc="left", fontsize=11, fontweight="bold", pad=10)


def finish(fig, name, caption):
    # The constrained-layout rect is (left, bottom, width, height).
    # Reserve independent space for the overall title and the caption.
    fig._suptitle.set_in_layout(False)
    fig._suptitle.set_position((0.5, 0.985))
    fig.get_layout_engine().set(rect=(0, 0.07, 1, 0.87))
    caption_artist = fig.text(0.025, 0.012, caption, fontsize=9, color="#454953", linespacing=1.5)
    extras = (fig._suptitle, caption_artist)
    fig.savefig(OUT / (name + ".png"), dpi=300, bbox_inches="tight", bbox_extra_artists=extras)
    fig.savefig(OUT / (name + ".pdf"), bbox_inches="tight", bbox_extra_artists=extras)
    plt.close(fig)


def edge(pair, events, threshold):
    fig, aa = plt.subplots(2, 2, figsize=(15, 10), layout="constrained")
    fig.suptitle("B2 edge timing: accepted-state threshold crossings", fontsize=15)
    for k, s in enumerate(pair["series"]):
        time = np.array(s["times_s"]) * 1e3
        gap = np.array(s["surface_gap_m"]) * 1e9
        mask = (time >= 1.35) & (time <= 2.00)
        label = ["10 µs", "5 µs"][k]
        aa[0, 0].plot(time[mask], gap[mask], ".-", color=C[k], ms=4, lw=1.3, label=label)
        last = s["events"][-1]
        aa[0, 0].axvline(last["time_s"] * 1e3, color=C[k], ls="--", lw=1)
        aa[0, 0].annotate(
            f"{last['time_s'] * 1e3:.3f} ms",
            (last["time_s"] * 1e3, 179.5 + 0.8 * k),
            color=C[k],
            xytext=(4, 0),
            textcoords="offset points",
        )
        speed = np.array(s["accepted_interval_normal_speed_m_s"]) * 1e6
        t = time[1:]
        aa[0, 1].plot(
            t[(t >= 1.35) & (t <= 2)], speed[(t >= 1.35) & (t <= 2)], ".-", color=C[k], ms=3, label=label
        )
        aa[0, 1].plot(last["time_s"] * 1e3, last["secant_normal_speed_m_s"] * 1e6, "o", color=C[k], ms=7)
    ax = aa[0, 0]
    ax.axhline(178, color=C[4], ls=":", label="178 nm graph threshold")
    ax.annotate("130 µs", (1.725, 178.75), ha="center", color=C[4])
    ax.annotate("", (1.790, 178.55), xytext=(1.660, 178.55), arrowprops=dict(arrowstyle="<->", color=C[4]))
    ax.set(
        xlim=(1.35, 2),
        ylim=(176.3, 182),
        xlabel="Accepted physical time (ms)",
        ylabel="ID 7–24 surface gap (nm)",
    )
    title(ax, "a  Small gap differences near a turning point")
    ax.legend(fontsize=9, loc="upper left")
    ax = aa[0, 1]
    ax.axhline(0, color="#777777", lw=0.7)
    ax.set(xlim=(1.35, 2), xlabel="Accepted interval endpoint (ms)", ylabel="Normal gap speed (µm/s)")
    title(ax, "b  Final disconnect: +4.73 / +15.58 µm/s")
    ax.legend(fontsize=9)
    ax = aa[1, 0]
    rows = events["events"]
    lag = np.array([e["time_difference_s"] for e in rows]) * 1e6
    ax.bar(np.arange(29), lag, color=[C[4] if e["timing_gate"] == "FAIL" else C[0] for e in rows])
    ax.axhline(10, color=C[4], ls=":", label="Original 10 µs limit")
    ax.set(
        xticks=np.arange(29),
        xticklabels=[f"{e['pair_id'][0]}–{e['pair_id'][1]}\n{'+' if e['on'] else '−'}" for e in rows],
        ylabel="Absolute accepted-event lag (µs)",
        xlabel="Matched pair ID and switch (+ on, − off)",
    )
    ax.tick_params(axis="x", rotation=70, labelsize=7)
    ax.legend(fontsize=9)
    title(ax, "c  All 29 events matched; 12 timing failures")
    ax = aa[1, 1]
    values = threshold["thresholds"]
    ax.bar(
        range(3),
        [q["comparison"]["max_event_difference_s"] * 1e6 for q in values],
        color=[C[0], C[1], C[2]],
        width=0.55,
    )
    ax.axhline(10, color=C[4], ls=":", label="Unchanged strict limit")
    ax.set(
        xticks=range(3),
        xticklabels=["2.100 R\n178.0 nm (original)", "2.095 R\n169.1 nm", "2.105 R\n186.9 nm"],
        ylabel="Maximum event lag (µs)",
        xlabel="Offline geometric definition only",
    )
    for i, q in enumerate(values):
        ax.text(
            i,
            q["comparison"]["max_event_difference_s"] * 1e6 + 3,
            f"{len(q['comparison']['events'])} matched events",
            ha="center",
            fontsize=9,
        )
    ax.set_ylim(0, 151)
    ax.legend(fontsize=9)
    title(ax, "d  Threshold sensitivity; no changed trajectory")
    for ax in aa.flat:
        ax.grid(axis="y", alpha=0.16)
    finish(
        fig,
        "Fig_F10I_1_edge_crossing",
        "O / seed 1102, R = 1.78 µm. Markers are actual accepted steps; connecting lines guide the eye, not new physical frames.\nGraph connectivity at a 178 nm surface gap is not sphere contact. No interpolated crossing replaces the original strict event gate.",
    )


def topology(events, group):
    fig, aa = plt.subplots(3, 2, figsize=(15, 12), layout="constrained")
    fig.suptitle("Instantaneous edges and persistent group topology are separate labels", fontsize=15)
    pairs = sorted({tuple(e["pair_id"]) for e in events["events"]})
    for k, mode in enumerate(["DT10", "DT5"]):
        z = np.load(BASE / f"analysis/ACCEPTED_GEOMETRY_{mode}.npz")
        lookup = {tuple(p): i for i, p in enumerate(z["pairs"])}
        ii = [lookup[p] for p in pairs]
        ax = aa[0, k]
        ax.pcolormesh(
            z["time_s"] * 1e3,
            np.arange(len(ii) + 1),
            z["linked"][:-1, ii].T,
            cmap=ListedColormap(["#f0f1f3", C[k]]),
            vmin=0,
            vmax=1,
            rasterized=True,
        )
        ax.set(
            yticks=np.arange(len(ii)) + 0.5,
            yticklabels=[f"{i}–{j}" for i, j in pairs],
            xlabel="Accepted physical time (ms)",
            ylabel="27 pairs with at least one changed edge",
        )
        ax.tick_params(axis="y", labelsize=7)
        title(ax, f"{chr(97 + k)}  {10 if k == 0 else 5} µs: colored = connected")
        m = group["modes"][mode]
        t = np.array(m["time_s"]) * 1e3
        for ax, key in [(aa[1, 0], "largest_fraction"), (aa[2, 0], "Rg_m"), (aa[2, 1], "aspect")]:
            v = np.array(m[key]) * (1e6 if key == "Rg_m" else 1)
            ax.plot(t, v, "-" if k == 0 else "--", color=C[k], label=f"{10 if k == 0 else 5} µs", lw=1.5)
            ax.set_xlabel("Actual accepted physical time (ms)")
            ax.grid(alpha=0.16)
            ax.legend(fontsize=9)
    aa[1, 0].set(ylim=(0.98, 1.02), ylabel="Largest connected group fraction")
    title(aa[1, 0], "c  Identical members: all IDs 0–31, all frames")
    aa[1, 0].text(0.15, 0.985, "1 component / 32 members throughout 0–5 ms", fontsize=10)
    ax = aa[1, 1]
    ax.axis("off")
    title(ax, "d  Layered evaluation (original gates retained)")
    for y, text, color in [
        (0.89, "Edge identities and per-pair switch sequences: PASS", C[0]),
        (0.73, "Edge timing: FAIL (12/29); cross-pair order: FAIL", C[4]),
        (0.57, "Group position and shape: PASS in this window", C[2]),
        (0.41, "2 ms persistence: no split or merge observed", C[0]),
        (0.25, "Right-censored topology events: 0", C[0]),
        (0.09, "Unseen persistent-event accuracy: NOT EVALUATED", C[4]),
    ]:
        ax.text(0.01, y, text, color=color, transform=ax.transAxes, fontsize=10)
    title(aa[2, 0], "e  Radius of gyration; max. relative difference 0.00636%")
    aa[2, 0].set_ylabel("Rg (µm)")
    title(aa[2, 1], "f  Aspect ratio; max. relative difference 0.03765%")
    aa[2, 1].set_ylabel("Principal-axis aspect ratio")
    finish(
        fig,
        "Fig_F10I_2_topology_significance",
        "Raster shows the last accepted graph state on each interval; all final accepted states are included in group statistics and saved member lists.\nNo multi-component event was observed. Stable membership here does not validate unobserved fragmentation, physical fusion or AI training readiness.",
    )


def solver(s, group, prop):
    fig, aa = plt.subplots(2, 3, figsize=(17, 10), layout="constrained")
    fig.suptitle("B2 solver tolerance and time discretization: distinct numerical tests", fontsize=15)
    for k, mode in enumerate(["DT10", "DT5"]):
        z = np.load(BASE / f"analysis/ACCEPTED_GEOMETRY_{mode}.npz")
        aa[0, 0].semilogy(
            z["time_s"][1:] * 1e3,
            z["solver_residual_m"],
            color=C[k],
            lw=0.8,
            label=f"{10 if k == 0 else 5} µs",
        )
    aa[0, 0].axhline(1e-13, color=C[4], ls=":", label="Original tolerance")
    aa[0, 0].set(xlabel="Actual accepted time (ms)", ylabel="Position equation residual (m)")
    aa[0, 0].legend(fontsize=8)
    title(aa[0, 0], "a  Archived nonlinear residuals")
    for k, tol in enumerate([1e-13, 1e-14, 1e-15]):
        vs = [next(v for v in row["variants"] if v["tolerance_m"] == tol) for row in s["rows"]]
        x = np.arange(14) + 0.14 * (k - 1)
        aa[0, 1].semilogy(
            x,
            [v["result"]["mobility_equation_residual_m"] for v in vs],
            "o",
            color=C[k],
            ms=4,
            label=f"tol {tol:.0e} m",
        )
        for i, v in enumerate(vs):
            if not v["result"]["converged"]:
                aa[0, 1].plot(
                    x[i], v["result"]["mobility_equation_residual_m"], "x", color=C[4], ms=11, mew=2
                )
        if k:
            aa[0, 2].plot(
                range(14),
                [v["vs_archived"]["max_position_difference_m"] * 1e9 for v in vs],
                "o-",
                color=C[k],
                ms=4,
                label=f"tol {tol:.0e} m",
            )
            aa[1, 0].plot(
                range(14),
                [v["vs_archived"]["max_force_difference_N"] * 1e15 for v in vs],
                "o-",
                color=C[k],
                ms=4,
                label=f"tol {tol:.0e} m",
            )
    aa[0, 1].semilogy(
        range(14),
        [row["reference"]["residual_m"] for row in s["rows"]],
        "+",
        color="#777777",
        ms=5,
        label="Independent root",
    )
    for ax in [aa[0, 1], aa[0, 2], aa[1, 0]]:
        ax.set(
            xticks=range(14),
            xticklabels=[q["name"].replace("DT", "") for q in s["rows"]],
            xlabel="Saved state: step size (µs) / accepted step",
        )
        ax.tick_params(axis="x", labelsize=7, rotation=65)
        ax.legend(fontsize=8)
    aa[0, 1].set_ylabel("Position equation residual (m)")
    title(aa[0, 1], "b  14 same-input CPU reference states")
    aa[0, 2].set_ylabel("Max. fixed-ID change from archive (nm)")
    aa[0, 2].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    title(aa[0, 2], "c  Position changes after tightening")
    aa[1, 0].set_ylabel("Max. bubble contact-force change (fN)")
    title(aa[1, 0], "d  Active contacts and edges unchanged")
    ax = aa[1, 1]
    for rows, label, color in [
        (group["comparison"], "10 / 5 µs; original tolerance", C[0]),
        (prop["same_dt_tolerance_comparison"]["series"], "10 µs; tolerance 1e−13 / 1e−14 m", C[2]),
    ]:
        rr = [q for q in rows if 0 < q["time_s"] <= 0.002]
        ax.semilogy(
            [q["time_s"] * 1e3 for q in rr], [q["max_fixed_ID_m"] * 1e9 for q in rr], color=color, label=label
        )
    ax.set(xlabel="Common actual accepted time (ms)", ylabel="Max. fixed-ID trajectory difference (nm)")
    ax.legend(fontsize=8, loc="lower right")
    title(ax, "e  One original-initial-state 2 ms control")
    ax = aa[1, 2]
    ax.axis("off")
    title(ax, "f  Limits and observed event sensitivity")
    ax.text(
        0.01,
        0.92,
        "1e−14 m tolerance: 14/14 states converge\n1e−15 m: 13/14 (red × marks failure)\nFailed residual: 1.0026e−15 m\n\n2 ms tighter-tolerance replay:\n200 accepted steps; max. ID change 0.00559 nm\nAll 23 prefix event times unchanged\nWorst 10 / 5 µs edge lag still 130 µs\n\nNear-threshold slope amplifies small gap changes.\n2.5 µs reference: NOT_RUN\nNo absolute temporal-order claim.",
        va="top",
        transform=ax.transAxes,
        fontsize=10,
        linespacing=1.7,
    )
    for ax in aa.flat:
        if ax.axison:
            ax.grid(axis="y", alpha=0.16)
    finish(
        fig,
        "Fig_F10I_3_solver_vs_timestep",
        "Same stiffness, acoustic force, CFD and discrete resistance law; no additional contact model. Nonconverged CPU iterate is retained as a failure.\nTighter local solves alone do not bound trajectory error; the single 2 ms replay tests propagation without an unverified checkpoint restart.",
    )


def gpu(raw):
    fig, aa = plt.subplots(2, 3, figsize=(17, 10), layout="constrained")
    fig.suptitle("GPU B2 fixed-state prototype: GPU_B2_NOT_ADMITTED", fontsize=15)
    rows = raw["rows"]
    indices = np.arange(len(rows)) + 1
    metrics = [
        ("max_position_vs_cpu_m", 1e-12, "Position error vs. CPU (m)"),
        ("max_force_vs_cpu_N", 1e-15, "Contact-force error vs. CPU (N)"),
        ("max_equation_residual_m", 1e-13, "Actual position equation residual (m)"),
    ]
    for k, (key, limit, label) in enumerate(metrics):
        ax = aa[0, k]
        vals = np.array([r["metrics"][key] for r in rows])
        positive = vals > 0
        ax.semilogy(indices[positive], vals[positive], "o", color=C[k], ms=5, label="Stored nonzero values")
        ax.axhline(limit, color=C[4], ls=":", label="Original comparison limit")
        ax.axvspan(18.6, 19.4, color=C[4], alpha=0.08)
        ax.set(
            xlabel="Case index: 1–12 two / 13–18 four / 19 saved32",
            ylabel=label,
            xticks=[1, 6, 12, 16, 18, 19],
        )
        ax.text(
            0.98,
            0.04,
            f"{sum(~positive)} exact zeros omitted",
            ha="right",
            fontsize=8,
            transform=ax.transAxes,
        )
        ax.legend(fontsize=8, loc="upper left")
        title(
            ax,
            f"{chr(97 + k)}  "
            + ["Position agreement", "Force agreement", "Independently checked residual"][k],
        )
        ax.grid(axis="y", alpha=0.16)
    r = rows[-1]
    t = r["gpu"]["timing"]
    ax = aa[1, 0]
    parts = [
        ("Prepare", t["preparation_wall_s"], C[0]),
        ("Scalar sync/read", t["scalar_sync_read_wall_s"], C[1]),
        ("Other solve wall", t["solve_wall_s"] - t["scalar_sync_read_wall_s"], C[2]),
        ("Final audit read", t["final_sync_read_wall_s"], C[3]),
    ]
    left = 0
    for label, value, color in parts:
        ax.barh(1, value * 1e3, left=left, color=color, label=label)
        left += value * 1e3
    ax.barh(0, r["cpu"]["total_s"] * 1e3, color="#8993a0", label="CPU contact (assembly excluded)")
    ax.set(
        yticks=[0, 1],
        yticklabels=["CPU 12.80 ms", "GPU 35.46 ms"],
        xlabel="Single saved32 attempt wall time (ms)",
        xlim=(0, 40),
        ylim=(-0.6, 2.1),
    )
    ax.legend(fontsize=8, loc="upper right")
    title(ax, "d  Measured attempt costs; not a speed benchmark")
    ax = aa[1, 1]
    ax.axis("off")
    title(ax, "e  Why the saved32 solver was rejected")
    ax.text(
        0.01,
        0.93,
        "Geometry / active contacts / edges: agree\nOverlap 17.498 nm < 17.8 nm budget\nCorrection 1.932 nm < 178 nm budget\n\nConservative bound: 1.850e−14 m\nRegistered stopping target: 1e−14 m\nArmijo exhausted 20 line-search trials\nBase PCG missed its recursive-residual target\n\n18 small cases match the CPU reference;\n2 are expected overlap-budget rejections.\nNo changed tolerance; no production admission.",
        va="top",
        transform=ax.transAxes,
        fontsize=10,
        linespacing=1.7,
    )
    ax = aa[1, 2]
    ax.axis("off")
    title(ax, "f  Full cost and scope retained")
    ax.text(
        0.01,
        0.93,
        "saved32 allocation: 0.534 ms\nReusable PCG graph build: 40.560 ms\n35 scalar reads: 22.242 ms within solve\nExplicit buffers: 78,376 bytes (lower bound)\nCG workspace / context / graphs: extra\n\nFive matched timing repeats: NOT_RUN\nStopped after the solver admission failure\nGPU v1 engineering failure also archived\n\nB = 1 fixed-state experiment only\nNo new 10 ms / B64 / training set\nPHYSICS_UNCALIBRATED; training_ready=false",
        va="top",
        transform=ax.transAxes,
        fontsize=10,
        linespacing=1.7,
    )
    finish(
        fig,
        "Fig_F10I_4_gpu_feasibility_admission",
        "Original budgets retained. Small cases 17–18 intentionally verify the same CPU/GPU overlap rejection; they are not positive contact-quality passes.\nGPU comparison agreement does not override nonconvergence. CPU assembly is outside the shown CPU timing; graph construction is outside the GPU step.",
    )


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    a = BASE / "analysis"
    pair = js(a / "PAIR_7_24.json")
    events = js(a / "EVENT_MATCHES.json")
    group = js(a / "GROUP_TOPOLOGY.json")
    edge(pair, events, js(a / "THRESHOLD_SENSITIVITY.json"))
    topology(events, group)
    solver(js(a / "SOLVER_TOLERANCE.json"), group, js(a / "TOLERANCE_PROPAGATION.json"))
    gpu(js(BASE / "gpu_v2/RESULTS.json"))
    inputs = [
        a / n
        for n in [
            "PAIR_7_24.json",
            "EVENT_MATCHES.json",
            "GROUP_TOPOLOGY.json",
            "THRESHOLD_SENSITIVITY.json",
            "SOLVER_TOLERANCE.json",
            "TOLERANCE_PROPAGATION.json",
            "ACCEPTED_GEOMETRY_DT10.npz",
            "ACCEPTED_GEOMETRY_DT5.npz",
        ]
    ] + [BASE / "gpu_v2/RESULTS.json"]
    write(
        OUT / "FIGURE_IDENTITY.json",
        dict(
            source=identity(),
            inputs={str(p.relative_to(BASE)): sha(p) for p in inputs},
            outputs={p.name: sha(p) for p in sorted(OUT.glob("Fig*"))},
            dpi=300,
            pdf_fonts="Embedded TrueType / English",
            new_physical_steps=0,
        ),
    )
    print("Four PNG/PDF groups written from archived data.")


if __name__ == "__main__":
    main()
