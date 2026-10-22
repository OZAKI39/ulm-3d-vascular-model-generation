"""Four reproducible publication figures from closed F9/F8 numerical evidence."""

import subprocess
import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
from .common import ROOT, BASE, F5, js, sha, write

COLORS = ["#416d9c", "#8855a2", "#dca44b", "#40877d", "#8d8b8a", "#bd4e4a", "#b7ac8a", "#182936"]
NAMES = ["CFD", "Primary", "Secondary", "Mobility balance", "Pair PBD", "Wall", "Boundary", "Accepted"]
CC = {"Y0": "#68747b", "Y1": "#bc514a", "Y2": "#356caa"}


def style(ax, title, xlabel=None, ylabel=None):
    ax.set_title(title, loc="left", fontsize=12, fontweight="bold", pad=10)
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if not hasattr(ax, "zaxis"):
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(alpha=0.15)


def clip_polygon(poly, axis, bound, keep_greater):
    out = []
    for a, b in zip(poly, np.roll(poly, -1, axis=0)):
        ia = a[axis] >= bound if keep_greater else a[axis] <= bound
        ib = b[axis] >= bound if keep_greater else b[axis] <= bound
        if ia:
            out.append(a)
        if ia != ib:
            out.append(a + (bound - a[axis]) / (b[axis] - a[axis]) * (b - a))
    return np.asarray(out)


def main():
    out = BASE / "figures"
    out.mkdir(exist_ok=True)
    r = js(BASE / "analysis/RESULTS.json")
    d = np.load(BASE / "analysis/plot_data.npz")
    cases = r["cases"]
    # Original target comes from the actual registered run, never inferred from its endpoint.
    cfg = js(ROOT / cases["Y1_10"]["path"] / "CONFIG.json")
    target = np.asarray(cfg["environment"]["target_m"])
    origin = d["Y1_10_c"][0]
    paths = {k: ROOT / v["path"] / "trajectory.h5" for k, v in cases.items()}
    raw = {}
    for key, path in paths.items():
        with h5py.File(path) as h:
            if "motion_ledger" in h:
                raw[key] = {
                    n: h["motion_ledger/" + n][:]
                    for n in ["displacement_m", "step_meta", "wall_contact_events"]
                }
    source_paths = [
        BASE / "analysis/RESULTS.json",
        BASE / "analysis/plot_data.npz",
        ROOT / "configs/fast_pzt/f9/protocol.json",
        ROOT / "scripts/fast_pzt_f9/figures.py",
        F5 / "inputs/f4/prepared_geometry.npz",
        *paths.values(),
    ]
    sources = {str(p.relative_to(ROOT)): sha(p) for p in source_paths}
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.labelsize": 10,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "pdf.fonttype": 42,
            "savefig.facecolor": "white",
            "axes.labelcolor": "#20262b",
            "text.color": "#20262b",
        }
    )

    def save(fig, name, notes):
        fig.savefig(out / (name + ".png"), dpi=220)
        fig.savefig(out / (name + ".pdf"))
        write(
            out / (name + ".json"),
            dict(
                sources=sources,
                input_hashes=js(F5 / "INPUT_MANIFEST.json"),
                plotting_commit=commit,
                execution_commits={k: v["execution_commit"] for k, v in r["provenance"].items()},
                interpretation=notes,
                no_smoothing=True,
                radius_exaggeration=1,
            ),
        )
        plt.close(fig)

    fig, axs = plt.subplots(2, 2, figsize=(12, 8), layout="constrained")
    for ax, key, letter in zip(axs[0], ["Y1_10", "Y1_5"], ["a", "b"]):
        cum = d[key + "_ledger_cumulative"]
        for i, (color, name) in enumerate(zip(COLORS, NAMES)):
            ax.plot(
                d[key + "_t"] * 1e3,
                cum[:, i, 2] * 1e6,
                color=color,
                lw=2 if i in [5, 7] else 1.3,
                ls="--" if i == 7 else "-",
                label=name,
            )
        style(ax, f"({letter}) Y1, {key.split('_')[1]} µs", "Time (ms)", "Mean cumulative signed Z (µm)")
    axs[0, 0].legend(ncol=2, fontsize=8.5, loc="upper left")
    for ax in axs[0]:
        ax.set_ylim(-0.19, 0.70)
    ax = axs[1, 0]
    for j, key in enumerate(["Y1_10", "Y1_5"]):
        ax.bar(
            np.arange(8) + (j - 0.5) * 0.36,
            d[key + "_ledger_cumulative"][-1, :8, 2] * 1e6,
            width=0.36,
            color=["#9dafbb", "#38566d"][j],
            label=f"{key.split('_')[1]} µs",
        )
    ax.set_xticks(
        range(8),
        ["CFD", "Primary", "Secondary", "Mobility", "PBD", "Wall", "Boundary", "Accepted"],
        rotation=35,
        ha="right",
    )
    ax.axhline(0, color="black", lw=0.6)
    ax.legend()
    style(ax, "(c) Signed endpoint accounting", ylabel="Mean Z displacement at 300 ms (µm)")
    ax = axs[1, 1]
    for key, color in [("Y1_10", CC["Y1"]), ("Y1_5", "#294e77")]:
        a = raw[key]
        ax.semilogy(
            (a["step_meta"][:, 0] + a["step_meta"][:, 1]) * 1e3,
            np.maximum(a["step_meta"][:, 3], 1e-25),
            lw=0.7,
            color=color,
            label=f"{key.split('_')[1]} µs",
        )
    ax.axhline(5e-17, color="black", ls="--", lw=1, label="Registered limit")
    ax.set_ylim(1e-22, 1e-16)
    style(ax, "(d) Closure at every accepted step", "Time (ms)", "Maximum particle residual (m)")
    ax.legend(fontsize=9, loc="lower right")
    fig.suptitle("Exact displacement accounting on the unchanged Y1 trajectory", fontsize=15)
    save(
        fig,
        "Fig_F9_1_displacement_ledger",
        "Signed fixed-32 means. Lines use exact per-step sums at save times; closure uses every physical step. Not a no-wall counterfactual.",
    )

    fig = plt.figure(figsize=(12, 9.2), layout="constrained")
    gs = fig.add_gridspec(2, 2, height_ratios=[1.1, 1])
    local = gs[0, 0].subgridspec(3, 1, height_ratios=[0.20, 0.73, 0.07])
    header = fig.add_subplot(local[0, 0])
    header.axis("off")
    header.set_title(
        "(a) Original facets, physical bubble radius", loc="left", fontsize=12, fontweight="bold"
    )
    caption = fig.add_subplot(local[2, 0])
    caption.axis("off")
    ax = fig.add_subplot(local[1, 0], projection="3d", computed_zorder=False)
    e = raw["Y1_10"]["wall_contact_events"]
    event = e[np.argmax(e[:, 14])]
    center = event[5:8] * 1e6
    mesh = np.load(F5 / "inputs/f4/prepared_geometry.npz")
    faces = mesh["faces"][mesh["labels"] == 1]
    polys = []
    face_ids = np.unique(e[:, 19]).astype(int)
    for face in face_ids:
        poly = mesh["vertices_m"][faces[face]] * 1e6 - center
        for axis, lo, hi in [(0, -8, 8), (1, -10, 10)]:
            poly = clip_polygon(poly, axis, lo, True)
            if len(poly):
                poly = clip_polygon(poly, axis, hi, False)
        polys.append(poly)
    ax.add_collection3d(
        Poly3DCollection(
            polys,
            facecolors=["#bba8a0", "#d4c4b9"],
            edgecolors="#6f6560",
            linewidths=1.2,
            alpha=0.65,
            zorder=1,
        )
    )
    u, v = np.meshgrid(np.linspace(0, 2 * np.pi, 50), np.linspace(0, np.pi, 25))
    ax.plot_surface(
        1.78 * np.cos(u) * np.sin(v),
        1.78 * np.sin(u) * np.sin(v),
        1.78 * np.cos(v),
        color="#3b78a5",
        linewidth=0,
        alpha=0.92,
        shade=True,
        zorder=2,
    )
    n = event[16:19]
    surface = (1.78 + event[15] * 1e6) * n
    ax.quiver(*surface, *(n * 3), color="#49434b", arrow_length_ratio=0.22, linewidth=2)
    ax.quiver(0, 0, 0, *(event[11:14] * 1e6 * 2000), color="#c13e38", arrow_length_ratio=0.22, linewidth=2.5)
    ax.set(
        xlim=(-8, 8),
        ylim=(-10, 10),
        zlim=(-5, 5),
        xlabel="Local X (µm)",
        ylabel="Local Y (µm)",
        zlabel="Local Z (µm)",
    )
    ax.set_box_aspect((16, 20, 10))
    ax.view_init(elev=24, azim=-63)
    header.text(
        0.02,
        0.85,
        f"ID {int(event[3])}; t = {event[1] * 1e3:.2f} ms\nR = 1.78 µm; correction arrow ×2000\nNormal arrow: 3 µm display length",
        transform=header.transAxes,
        fontsize=8.5,
        va="top",
    )
    caption.text(
        0.02,
        0.6,
        "Facets 25305 / 271 clipped to viewing box;\norigin = actual pre-projection bubble center",
        transform=caption.transAxes,
        fontsize=8,
    )
    ax = fig.add_subplot(gs[0, 1])
    for j, key in enumerate(["Y1_10", "Y1_5"]):
        vals = np.asarray(cases[key]["ledger"]["per_ID_signed_wall_m"])
        ax.bar(
            np.arange(32) + (j - 0.5) * 0.38,
            vals[:, 2] * 1e6,
            width=0.38,
            color=["#b6574d", "#335f83"][j],
            label=key.replace("_", ", ") + " µs",
        )
    style(
        ax,
        "(b) Individual signed wall contributions",
        "Fixed bubble ID",
        "Cumulative signed Z per bubble (µm)",
    )
    ax.set_xticks([0, 2, 8, 16, 18, 24, 31])
    ax.legend()
    edges = np.linspace(0, 0.3, 151)
    heats = []
    for key in ["Y1_10", "Y1_5"]:
        ev = raw[key]["wall_contact_events"]
        a = np.zeros((32, len(edges) - 1))
        bins = np.minimum(np.searchsorted(edges, ev[:, 1], side="right") - 1, len(edges) - 2)
        np.add.at(a, (ev[:, 3].astype(int), bins), ev[:, 13] * 1e9)
        heats.append(a)
    vmax = max(abs(a).max() for a in heats)
    for j, (key, a) in enumerate(zip(["Y1_10", "Y1_5"], heats)):
        ax = fig.add_subplot(gs[1, j])
        im = ax.imshow(
            a,
            extent=(0, 300, -0.5, 31.5),
            origin="lower",
            aspect="auto",
            cmap="RdBu_r",
            vmin=-vmax,
            vmax=vmax,
            interpolation="nearest",
        )
        style(ax, f"({'cd'[j]}) Y1, {key.split('_')[1]} µs: wall events", "Time (ms)", "Fixed bubble ID")
        ax.set_yticks([0, 2, 8, 16, 18, 24, 31])
        ax.grid(False)
        fig.colorbar(im, ax=ax, label="Signed Z sum per 2 ms bin (nm)", shrink=0.85)
    fig.suptitle("Where the geometric wall correction enters the trajectory", fontsize=15)
    save(
        fig,
        "Fig_F9_2_wall_correction_mechanism",
        "Original mesh, physical 1.78µm sphere, actual maximum wall event before/after vector magnified x2000. Display clipping preserves facet planes; normal arrow shows orientation only. Heat maps are sums of raw events, not interpolated trajectories.",
    )

    fig, axs = plt.subplots(2, 3, figsize=(13.3, 8.1), layout="constrained")
    ax = axs[0, 0]
    target_rel = (target - origin) * 1e6
    ax.add_patch(Circle(target_rel[1:], 0.15, color="#c5d5bd", alpha=0.45))
    ax.scatter(*target_rel[1:], marker="*", s=130, color="black", label="Original target")
    for c in ["Y0", "Y1", "Y2"]:
        pts = np.array([cases[c + "_" + dt]["gate"]["derived_centroid_m"] for dt in ["10", "5"]])
        pts = (pts - origin) * 1e6
        ax.plot(pts[:, 1], pts[:, 2], color=CC[c], lw=0.8, label=c)
        for i, marker in enumerate(["o", "s"]):
            ax.scatter(
                pts[i, 1],
                pts[i, 2],
                edgecolors=CC[c],
                marker=marker,
                s=42,
                facecolors=CC[c] if i == 0 else "none",
            )
    ax.set_aspect("equal", adjustable="box")
    ax.legend(fontsize=8, ncol=2, loc="lower left")
    ax.set_ylim(-0.11, 0.46)
    ax.text(0.02, 0.98, "Circle: 10 µs; square: 5 µs", transform=ax.transAxes, va="top", fontsize=8)
    style(ax, "(a) Same downstream X gate", "Y − initial centroid (µm)", "Z − initial centroid (µm)")
    ax = axs[0, 1]
    for c in ["Y0", "Y1", "Y2"]:
        ax.plot(
            [10, 5],
            [cases[c + "_" + dt]["gate"]["transverse_error_m"] * 1e6 for dt in ["10", "5"]],
            "o-",
            color=CC[c],
            label=c,
        )
    ax.axhline(0.15, color="black", ls="--", lw=1, label="Original tolerance")
    ax.set_xticks([5, 10])
    ax.legend(fontsize=8)
    style(ax, "(b) Target error at the gate", "Physical timestep (µs)", "YZ target error (µm)")
    ax = axs[0, 2]
    for c in ["Y0", "Y1", "Y2"]:
        ax.plot(
            [10, 5],
            [cases[c + "_" + dt]["gate"]["largest_fraction"] for dt in ["10", "5"]],
            "o-",
            color=CC[c],
            label=c,
        )
    ax.axhline(0.9, color="black", ls="--", lw=1, label="Original minimum")
    ax.set_xticks([5, 10])
    ax.set_ylim(0, 1.08)
    ax.legend(fontsize=8)
    style(ax, "(c) Geometric connectivity", "Physical timestep (µs)", "Largest connected fraction at gate")
    ax = axs[1, 0]
    for c in ["Y0", "Y1", "Y2"]:
        ax.plot(
            d[c + "_10_t"] * 1e3,
            np.linalg.norm(d[c + "_step_centroid_delta"][:, 1:], axis=1) * 1e6,
            color=CC[c],
            label=c,
        )
    ax.axhline(0.05, color="black", ls="--", lw=1, label="50 nm task scale")
    ax.axvline(cases["Y1_10"]["ledger"]["first_wall_step_bracket_s"][0] * 1e3, color=CC["Y1"], ls=":", lw=1)
    ax.legend(fontsize=8)
    style(ax, "(d) Same-time 10 / 5 µs sensitivity", "Time (ms)", "YZ centroid difference (µm)")
    ax = axs[1, 1]
    delta = d["Y1_step_ledger_delta"]
    for i in [0, 1, 4, 5, 7]:
        ax.plot(
            d["Y1_10_t"] * 1e3,
            delta[:, i, 2] * 1e6,
            color=COLORS[i],
            ls="--" if i == 7 else "-",
            label=NAMES[i],
        )
    ax.legend(fontsize=8)
    style(ax, "(e) Y1 difference by signed source", "Time (ms)", "Fine − coarse cumulative Z (µm)")
    ax = axs[1, 2]
    for dt, ls in [("10", "-"), ("5", "--")]:
        for c in ["Y0", "Y1"]:
            z = (
                np.asarray(cases["Y2_" + dt]["gate"]["derived_centroid_m"])
                - np.asarray(cases[c + "_" + dt]["gate"]["derived_centroid_m"])
            )[2] * 1e6
            ax.plot(
                [0 if c == "Y0" else 1],
                [z],
                marker="o" if dt == "10" else "s",
                color=CC[c],
                ms=7,
                label=f"vs {c}, {dt} µs",
            )
    ax.axhline(0, color="black", lw=0.8)
    ax.axhline(0.05, color="black", ls=":", lw=0.8)
    ax.set_xticks([0, 1], ["Y2 − Y0", "Y2 − Y1"])
    ax.set_xlim(-0.5, 1.5)
    ax.legend(fontsize=8)
    style(ax, "(f) Extra target-direction progress", ylabel="Gate Z increment (µm)")
    fig.suptitle("Timestep sensitivity at fixed time and the original gate", fontsize=15)
    save(
        fig,
        "Fig_F9_3_timestep_and_gate",
        "Accepted-step linear gate intersections, not saved frames. Target tolerance unchanged. Differences are deterministic two-step comparisons, not confidence intervals. First Y1 wall event dotted in (d); Y1 controls never switch.",
    )

    fig, axs = plt.subplots(
        1, 3, figsize=(13.4, 4.8), layout="constrained", gridspec_kw={"width_ratios": [1.2, 1.1, 1]}
    )
    ax = axs[0]
    keys = ["Y0_5", "Y1_10", "Y1_5", "Y2_10"]
    for j, (idx, label, color) in enumerate(
        [(1, "Acoustic free", COLORS[1]), (5, "Wall", COLORS[5]), (7, "Accepted", COLORS[7])]
    ):
        vals = []
        for key in keys:
            cum = d[key + "_ledger_cumulative"][-1]
            vals.append((cum[1, 2] + cum[2, 2] if idx == 1 else cum[idx, 2]) * 1e6)
        ax.bar(np.arange(4) + (j - 1) * 0.24, vals, width=0.24, label=label, color=color)
    ax.set_xticks(range(4), ["Y0\n5 µs", "Y1\n10 µs", "Y1\n5 µs", "Y2\n10 µs"])
    ax.axhline(0, color="black", lw=0.6)
    ax.legend(fontsize=8)
    style(ax, "(a) Audited signed movement", ylabel="Mean Z displacement at 300 ms (µm)")
    ax = axs[1]
    for j, (field, label, color) in enumerate(
        [
            ("gate_target_direction_increment_m", "Measured Y2 increment", "#356caa"),
            ("numerical_pair_bound_m", "Paired step sensitivity", "#d5ad63"),
            ("required_increment_m", "Required increment", "#9a737e"),
        ]
    ):
        ax.bar(
            np.arange(2) + (j - 1) * 0.23,
            [r["effects"][c][field] * 1e6 for c in ["Y0", "Y1"]],
            width=0.23,
            label=label,
            color=color,
        )
    ax.set_xticks([0, 1], ["Against Y0", "Against Y1"])
    ax.axhline(0, color="black", lw=0.6)
    ax.legend(fontsize=8)
    style(ax, "(b) Resolving active control", ylabel="Gate Z increment / threshold (µm)")
    ax = axs[2]
    ax.axis("off")
    ax.set_title("(c) Separate evidence gates", loc="left", fontweight="bold", fontsize=12, pad=10)
    texts = [
        ("NUMERICAL PASS", "Exact ledger; original budgets pass.", "#42694f"),
        ("WALL SENSITIVE", "Y1 wall term changes by 0.209 µm\nwhen the step is halved.", "#a24c42"),
        (
            "NAVIGATION NOT RESOLVED",
            "Y2 does not exceed both controls\nby the registered numerical margin.",
            "#a24c42",
        ),
        ("COHESION ONLY", "Acoustic groups remain connected;\nconnectivity is not fusion.", "#536476"),
        (
            "PHYSICS UNCALIBRATED",
            "Geometric wall projection is not\nvalidated near-wall propulsion.",
            "#536476",
        ),
    ]
    for y, (title, body, color) in zip([0.93, 0.75, 0.55, 0.34, 0.14], texts):
        ax.text(0, y, title, transform=ax.transAxes, fontsize=10, weight="bold", color=color)
        ax.text(0, y - 0.055, body, transform=ax.transAxes, fontsize=9, va="top", linespacing=1.4)
    fig.suptitle("Near-wall control remains unresolved at the original task scale", fontsize=15)
    save(
        fig,
        "Fig_F9_4_wall_vs_acoustic_authority",
        "Free acoustic = primary plus secondary before mobility. Panel(a) omits CFD and other terms only for comparison; full identity in Figure1 and raw data. Y0 10µs and Y2 5µs exact ledgers were not computed. Bounds are sensitivity thresholds, not statistical uncertainty.",
    )


if __name__ == "__main__":
    main()
