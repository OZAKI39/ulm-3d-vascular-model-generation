"""Four data-traceable F8 review figures. No altered coordinates or smoothed paths."""

import subprocess
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Circle
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from .common import ROOT, BASE, js, write, sha

OUT = BASE / "figures"
COL = {"Y0": "#777777", "Y1": "#aa7537", "Y2": "#326a9c", "A0": "#aa7537", "A1": "#326a9c"}
LABEL = {
    "Y0": "CFD only",
    "Y1": "Fixed PZT",
    "Y2": "Target feedback",
    "A0": "Stokes prediction",
    "A1": "RPY prediction",
}
plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "axes.edgecolor": "#333333",
        "xtick.color": "#333333",
        "ytick.color": "#333333",
        "pdf.fonttype": 42,
        "savefig.facecolor": "white",
        "figure.facecolor": "white",
    }
)


def save(fig, name, caption, files):
    fig.savefig(OUT / (name + ".png"), dpi=220)
    fig.savefig(OUT / (name + ".pdf"))
    plt.close(fig)
    sources = {str(p.relative_to(ROOT)): sha(p) for p in files}
    write(
        OUT / (name + ".json"),
        dict(
            caption=caption,
            plot_commit=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
            sources=sources,
            input_hashes=js(ROOT / "results/fast_pzt/f5/INPUT_MANIFEST.json"),
            execution_identities=js(BASE / "analysis/RESULTS.json")["provenance"],
            notes="Actual saved SI coordinates. Derived gate intersections explicitly labelled; bubble centers are markers, not enlarged physical spheres. No temporal smoothing or synthetic trajectory frames.",
        ),
    )


def section(axis, coordinate):
    data = np.load(ROOT / "results/fast_pzt/f5/inputs/f4/prepared_geometry.npz")
    tri = data["vertices_m"][data["faces"]]
    other = [i for i in range(3) if i != axis]
    segments = []
    for t in tri[(tri[:, :, axis].min(1) <= coordinate) & (tri[:, :, axis].max(1) >= coordinate)]:
        points = []
        for i, j in [(0, 1), (1, 2), (2, 0)]:
            a, b = t[i], t[j]
            if (a[axis] - coordinate) * (b[axis] - coordinate) < 0:
                points.append(a + (b - a) * (coordinate - a[axis]) / (b[axis] - a[axis]))
        if len(points) == 2:
            segments.append(np.array(points)[:, other] * 1e6)
    return segments


def plain(ax):
    ax.ticklabel_format(useOffset=False, style="plain")
    ax.grid(alpha=0.16)


def figure1():
    root = BASE / "y_static_v1"
    reg = js(root / "TARGETS_BEFORE_TRAJECTORIES.json")
    checks = js(root / "CHECKS.json")
    data = np.load(ROOT / "results/fast_pzt/f5/inputs/f4/prepared_geometry.npz")
    x = np.array(reg["initial_positions_m"]) * 1e6
    c = x.mean(0)
    target = np.array(reg["target_m"]) * 1e6
    support_file = ROOT / "results/fast_pzt/f7/reachability_v1/REACHABILITY.json"
    support = js(support_file)["neighborhood"]
    lower = np.array(support["field_lower_m"]) * 1e3
    upper = np.array(support["field_upper_m"]) * 1e3
    fig = plt.figure(figsize=(15, 10), layout="constrained")
    gs = fig.add_gridspec(2, 2)
    ax = fig.add_subplot(gs[0, 0], projection="3d")
    v = data["vertices_m"] * 1e3
    tri = v[data["faces"][data["labels"] == 1]]
    surf = Poly3DCollection(tri, facecolors="#b66b70", edgecolors="none", alpha=0.16, rasterized=True)
    ax.add_collection3d(surf)
    import itertools

    corners = np.array(list(itertools.product(*zip(lower, upper))))
    for i in range(8):
        for j in range(i):
            if np.count_nonzero(corners[i] != corners[j]) == 1:
                segment = corners[[i, j]]
                ax.plot(*segment.T, c="#826393", ls=":", lw=0.7, alpha=0.7)
    ax.plot([], [], [], c="#826393", ls=":", label="Stored acoustic support")
    bounds = np.column_stack([v.min(0), v.max(0)])
    for axis, setter in enumerate([ax.set_xlim, ax.set_ylim, ax.set_zlim]):
        setter(bounds[axis])
    ax.set_box_aspect(np.ptp(v, axis=0))
    ax.view_init(elev=22, azim=-62)
    ax.scatter(*(c / 1000), s=38, c=COL["Y2"], depthshade=False, label="Initial group / local study")
    ax.scatter(*(target / 1000), s=50, c="black", marker="*", depthshade=False)
    ax.set(xlabel="X (mm)", ylabel="Y (mm)", zlabel="Z (mm)", title="(a) Original Y lumen: location overview")
    ax.legend(loc="upper left", fontsize=8)
    ax.text2D(
        0.02,
        0.02,
        "Local parent-vessel transport; no bifurcation crossing",
        transform=ax.transAxes,
        fontsize=9,
    )
    for axis in [ax.xaxis, ax.yaxis, ax.zaxis]:
        axis.set_major_locator(MaxNLocator(4))
    ax = fig.add_subplot(gs[0, 1])
    plain(ax)
    ax.set_facecolor("#f5f1f8")
    ax.add_collection(LineCollection(section(1, c[1] / 1e6), colors="#b26870", lw=2))
    ax.scatter(x[:, 0], x[:, 2], s=9, c=COL["Y2"], label="Actual bubble centers")
    flow = np.load(root / "static_data.npz")["flow_velocity_m_s"] * 1e6
    q = ax.quiver(
        x[:, 0],
        x[:, 2],
        flow[:, 0],
        flow[:, 2],
        angles="xy",
        scale_units="xy",
        scale=20,
        color="#777777",
        width=0.004,
    )
    ax.quiverkey(q, 0.68, 0.90, 80, "80 µm/s", coordinates="axes", labelpos="E")
    ax.axvline(target[0], c="black", ls="--", lw=1, label="Registered X gate")
    ax.scatter(target[0], target[2], marker="*", s=100, c="black")
    ax.set(
        xlim=(c[0] - 12, target[0] + 7),
        ylim=(c[2] - 12, c[2] + 14),
        xlabel="Original X (µm)",
        ylabel="Original Z (µm)",
        title="(b) Actual-scale local XZ projection",
    )
    ax.set_aspect("equal")
    ax.legend(fontsize=8, loc="lower right")
    ax.text(
        0.02,
        0.98,
        "Entire window is within stored acoustic support",
        transform=ax.transAxes,
        va="top",
        fontsize=8,
    )
    ax = fig.add_subplot(gs[1, 0])
    plain(ax)
    seg = section(0, target[0] / 1e6)
    ax.add_collection(LineCollection(seg, colors="#b26870", lw=2, label="True YZ wall intersection"))
    ax.scatter(x[:, 1], x[:, 2], s=10, c=COL["Y2"], label="Initial centers projected to YZ")
    ax.scatter(target[1], target[2], s=90, c="black", marker="*", label="Registered centroid target")
    ax.add_patch(
        Circle(target[1:], reg["design"]["target_tolerance_m"] * 1e6, fill=False, color="black", lw=1)
    )
    ax.set(
        xlim=(c[1] - 12, c[1] + 12),
        ylim=(c[2] - 12, c[2] + 12),
        xlabel="Original Y (µm)",
        ylabel="Original Z (µm)",
        title="(c) Downstream plane and finite-radius clearance",
    )
    ax.set_aspect("equal")
    ax.legend(fontsize=8, loc="upper left")
    ax.text(
        0.02,
        0.02,
        f"Target surface gap for R = 1.78 µm: {reg['target_wall_gap_m'] * 1e6:.3f} µm\nCentroid tolerance: 0.15 µm; target offset: +Z 0.25 µm",
        transform=ax.transAxes,
        fontsize=8,
    )
    ax = fig.add_subplot(gs[1, 1])
    plain(ax)
    for row in checks["rows"]:
        if row["production_mode"] != "M0":
            continue
        vec = np.array(row["acoustic_mean_m_s"]) * 1e6
        name = row["condition"]
        ax.scatter(vec[1], vec[2], s=55, c=COL[name], label=LABEL[name])
        ax.annotate("", vec[1:], (0, 0), arrowprops=dict(arrowstyle="->", color=COL[name]))
    ax.axhline(0, c="#aaaaaa", lw=0.7)
    ax.axvline(0, c="#aaaaaa", lw=0.7)
    ax.set(
        xlabel="Acoustic mean drift Y (µm/s)",
        ylabel="Acoustic mean drift Z (µm/s)",
        title="(d) Initial-state M0 control authority",
    )
    ax.legend(fontsize=9)
    ax.set_aspect("equal")
    ax.set_xlim(-0.55, 0.45)
    ax.set_ylim(-0.80, 0.50)
    ax.text(
        0.03,
        0.04,
        "Static velocities precede PBD and wall correction.\nEstimated CFD flight: 247.6 ms; registered horizon: 300 ms.\nThe finite-candidate estimate is not a future trajectory.",
        transform=ax.transAxes,
        fontsize=8,
    )
    fig.suptitle("Original Y geometry, CFD transport and preregistered transverse target", fontsize=16)
    save(
        fig,
        "Fig_F8_1_Y_flow_and_target",
        "Original full lumen, local metric XZ and true gate YZ wall sections, actual CFD arrows and preregistered target. Static directions are not integrated trajectories.",
        [
            root / "TARGETS_BEFORE_TRAJECTORIES.json",
            root / "CHECKS.json",
            root / "static_data.npz",
            support_file,
            ROOT / "results/fast_pzt/f5/inputs/f4/prepared_geometry.npz",
        ],
    )


def figure2(d, r):
    reg = js(BASE / "y_static_v1/TARGETS_BEFORE_TRAJECTORIES.json")
    target = np.array(reg["target_m"]) * 1e6
    fig, axes = plt.subplots(2, 2, figsize=(14, 10), layout="constrained")
    for col, axis in enumerate([1, 2]):
        ax = axes[0, col]
        plain(ax)
        for name in ["Y0", "Y1", "Y2"]:
            c = d[name + "_c"] * 1e6
            g = r["cases"][name]["gate"]
            ax.plot(c[:, 0], c[:, axis], color=COL[name], label=LABEL[name])
            ax.scatter(
                d[name + "_cross_c"][0] * 1e6, d[name + "_cross_c"][axis] * 1e6, c=COL[name], s=36, marker="D"
            )
            for ms in [0, 100, 200, 300]:
                i = np.argmin(abs(d[name + "_t"] - ms / 1000))
                ax.plot(c[i, 0], c[i, axis], ".", c=COL[name], ms=4)
        ax.axvline(target[0], c="black", ls="--", lw=1)
        ax.axhline(target[axis], c="black", ls=":", lw=1)
        ax.set(
            xlabel="Original X (µm)",
            ylabel=f"Original {'XYZ'[axis]} (µm)",
            title=f"({'ab'[col]}) Downstream path: X versus {'XYZ'[axis]}",
        )
        ax.legend(fontsize=8)
    ax = axes[1, 0]
    plain(ax)
    ax.add_patch(
        Circle(
            target[1:],
            0.15,
            facecolor="#eeeeee",
            edgecolor="black",
            lw=1,
            label="Centroid tolerance: 0.15 µm",
        )
    )
    ax.scatter(target[1], target[2], marker="*", s=100, c="black", label="Registered target")
    for name in ["Y0", "Y1", "Y2"]:
        c = d[name + "_cross_c"] * 1e6
        g = r["cases"][name]["gate"]
        ax.scatter(
            c[1],
            c[2],
            s=60,
            c=COL[name],
            marker="D",
            label=f"{LABEL[name]}: {g['transverse_error_m'] * 1e6:.4f} µm",
        )
        ax.annotate(
            f"{g['time_s'] * 1e3:.2f} ms",
            c[1:],
            xytext=(-105, 12) if name == "Y1" else (10, -9) if name == "Y2" else (7, 7),
            arrowprops=dict(arrowstyle="-", color=COL[name], lw=0.6),
            textcoords="offset points",
            fontsize=8,
            c=COL[name],
        )
        if name + "_half_cross_c" in d:
            h = d[name + "_half_cross_c"] * 1e6
            ax.plot([c[1], h[1]], [c[2], h[2]], c=COL[name], lw=1)
            ax.scatter(h[1], h[2], s=45, facecolors="none", edgecolors=COL[name], marker="o")
    ax.set_aspect("equal")
    ax.set(
        xlim=(target[1] - 0.20, target[1] + 0.30),
        ylim=(target[2] - 0.32, target[2] + 0.20),
        xlabel="Original Y (µm)",
        ylabel="Original Z (µm)",
        title="(c) Same X gate: fixed-32-ID centroid",
    )
    ax.legend(fontsize=8, loc="upper left")
    ax.text(
        0.02,
        0.02,
        "Diamonds: 10 µs; open circles: 5 µs\nDerived step intersections, not saved sample frames",
        transform=ax.transAxes,
        fontsize=8,
    )
    ax = axes[1, 1]
    plain(ax)
    for name in ["Y0", "Y1", "Y2"]:
        x = d[name + "_cross_x"] * 1e6
        g = r["cases"][name]["gate"]
        ax.scatter(
            x[:, 1],
            x[:, 2],
            s=20,
            alpha=0.65,
            facecolors="none" if name == "Y1" else COL[name],
            edgecolors=COL[name],
            label=f"{LABEL[name]}: connected {g['largest_fraction'] * 32:.0f}/32",
        )
    ax.scatter(target[1], target[2], marker="*", s=90, c="black")
    ax.set_aspect("equal")
    ax.set(
        xlabel="Original Y (µm)",
        ylabel="Original Z (µm)",
        title="(d) All 32 centers at derived gate intersections",
    )
    ax.legend(fontsize=8)
    ax.text(
        0.02,
        0.02,
        "Marker size is not bubble radius. Connectivity is geometric.",
        transform=ax.transAxes,
        fontsize=8,
    )
    fig.suptitle("Matched downstream transport: sound off, fixed PZT and target feedback", fontsize=16)
    save(
        fig,
        "Fig_F8_2_transverse_navigation",
        "True saved centroid paths, same-gate linear intersections, preregistered centroid tolerance and all32 particle centers. Initial coordinates and physical models are matched.",
        [BASE / "analysis/plot_data.npz", BASE / "analysis/RESULTS.json"],
    )


def figure3(d, r):
    fig, axes = plt.subplots(3, 3, figsize=(16, 12), layout="constrained")
    ax = axes.ravel()
    for name in ["Y0", "Y1", "Y2"]:
        t = d[name + "_t"] * 1e3
        col = COL[name]
        ax[0].plot(t, d[name + "_error"] * 1e6, c=col, label=LABEL[name])
        ax[3].plot(t, d[name + "_fraction"], c=col, label=LABEL[name])
        ax[4].plot(t, d[name + "_rg"] * 1e6, c=col)
        ax[7].plot(t, d[name + "_quality"][:, :, 3].max(1) * 1e9, c=col)
        ax[8].plot(t, d[name + "_wall"][:, :, 3].sum(1) / 32 * 1e9, c=col)
        ax[8].plot(t, d[name + "_sampled_wall_residual_m"][:, 2] * 1e9, c=col, ls=":", lw=1)
        for a, ls in enumerate(["-", "--", ":"]):
            ax[2].plot(t, (d[name + "_c"][:, a] - d[name + "_c"][0, a]) * 1e6, c=col, ls=ls, lw=1)
    ax[0].axhline(0.15, c="black", ls=":", label="0.15 µm gate tolerance")
    ax[0].legend(fontsize=8)
    for control, ls in [("Y0", "-"), ("Y1", "--")]:
        ax[1].plot(
            d["Y2_t"] * 1e3,
            d["Y2_minus_" + control][:, 2] * 1e6,
            c=COL[control],
            ls=ls,
            label="Feedback minus " + LABEL[control],
        )
    ax[1].axhline(0.05, c="black", ls=":", label="Minimum positive effect: 0.05 µm")
    ax[1].axhline(0, c="#777777", lw=0.6)
    ax[1].legend(fontsize=7)
    ax[2].text(
        0.03,
        0.97,
        "Solid: X; dashed: Y; dotted: Z\nSame condition colors as panel (a)",
        va="top",
        transform=ax[2].transAxes,
        fontsize=8,
    )
    ax[3].axhline(0.9, c="black", ls=":", lw=1)
    ax[3].set_ylim(0, 1.06)
    for i, col in enumerate(["#356a93", "#a46d3b", "#7b5b8c", "#518476"]):
        times = np.r_[d["Y2_control_t"], d["Y2_t"][-1]] * 1e3
        ax[5].step(
            times,
            np.r_[d["Y2_amplitudes"][:, i], d["Y2_amplitudes"][-1, i]],
            where="post",
            c=col,
            label=f"PZT {i + 1}",
        )
        ax[6].step(times, np.r_[d["Y2_phases"][:, i], d["Y2_phases"][-1, i]], where="post", c=col)
    ax[5].legend(fontsize=7, ncol=2)
    ax[5].axhline(291.237687721856, c="black", ls=":", lw=0.7)
    if "Y1_half_wall" in d:
        ax[8].plot(
            d["Y1_half_t"] * 1e3,
            d["Y1_half_wall"][:, :, 3].sum(1) / 32 * 1e9,
            c=COL["Y1"],
            ls="--",
            label="Fixed PZT, 5 µs",
        )
        ax[8].legend(fontsize=7, loc="lower right")
    ax[7].axhline(17.8, c="black", ls=":", label="0.01R")
    ax[7].legend(fontsize=8)
    ax[8].text(
        0.03,
        0.97,
        "Solid: summed wall correction / 32 (norm bound)\nDotted: +Z drift minus sampled CFD + primary",
        va="top",
        transform=ax[8].transAxes,
        fontsize=8,
    )
    titles = [
        "Transverse target error",
        "Same-time +Z control increments",
        "Fixed-ID displacement",
        "Largest connected fraction",
        "Fixed-ID radius of gyration",
        "Feedback amplitudes",
        "Feedback phases",
        "Cumulative max bubble overlap",
        "Wall correction and group-drift diagnostic",
    ]
    units = [
        "Error (µm)",
        "Increment (µm)",
        "Displacement (µm)",
        "Fraction",
        "Rg (µm)",
        "Source amplitude (Pa)",
        "Phase (rad)",
        "Overlap (nm)",
        "Correction / residual (nm)",
    ]
    for i, a in enumerate(ax):
        a.set(title=f"({chr(97 + i)}) {titles[i]}", xlabel="Time (ms)", ylabel=units[i], xlim=(0, 300))
        a.grid(alpha=0.15)
        for name in ["Y0", "Y1", "Y2"]:
            a.axvline(r["cases"][name]["gate"]["time_s"] * 1e3, color=COL[name], lw=0.65, ls=":", alpha=0.6)
    fig.suptitle(
        "Control events, transverse motion and numerical quality\nVertical dotted lines: each condition crosses the same downstream gate",
        fontsize=15,
    )
    save(
        fig,
        "Fig_F8_3_control_error_morphology",
        "Actual0.5ms saved samples, recorded25ms control events and actual-step derived gate times. Geometric connectivity is not fusion.",
        [BASE / "analysis/plot_data.npz", BASE / "analysis/RESULTS.json"],
    )


def figure4(d, r):
    a = js(BASE / "a_static_v1/CHECKS.json")
    contact = js(BASE / "contact_v1/CHECKS.json")
    cd = np.load(BASE / "contact_v1/contact_states.npz")
    fig, axes = plt.subplots(3, 2, figsize=(14, 13), layout="constrained")
    ax = axes[0, 0]
    values = []
    labels = []
    for row in a["rows"]:
        for name in ["A0", "A1"]:
            values.append(row["modes"][name]["event"]["amplitude_Pa"])
            labels.append(f"t={row['time_s'] * 1e3:.0f} ms / {LABEL[name]}")
    im = ax.imshow(values, aspect="auto", cmap="Purples", vmin=0, vmax=291.237687721856)
    ax.set(
        xticks=range(4),
        xticklabels=["PZT 1", "PZT 2", "PZT 3", "PZT 4"],
        yticks=range(6),
        yticklabels=labels,
        title="(a) Matched-state selected source amplitudes",
    )
    fig.colorbar(im, ax=ax, label="Source amplitude (Pa)", shrink=0.8)
    ax = axes[1, 0]
    xx = np.arange(3)
    for j, name in enumerate(["A0", "A1"]):
        vals = [row["modes"][name]["mean_prediction_residual_m_s"] * 1e6 for row in a["rows"]]
        ax.bar(xx + (j - 0.5) * 0.3, vals, width=0.3, color=COL[name], label=LABEL[name])
    ax.set(
        xticks=xx,
        xticklabels=["Initial", "500 ms", "900 ms"],
        ylabel="Mean velocity residual (µm/s)",
        title="(c) Prediction versus instantaneous production RPY",
    )
    ax.legend(fontsize=8)
    ax.text(0.97, 0.93, "Before contact projection", ha="right", transform=ax.transAxes, fontsize=8)
    ax = axes[2, 0]
    for name in ["A0", "A1"]:
        ax.plot(d[name + "_t"] * 1e3, d[name + "_error"] * 1e6, color=COL[name], label=LABEL[name])
    ax.axhline(2, c="black", ls=":", label="2 µm")
    ax.axvspan(893.5, 1000, color=COL["A1"], alpha=0.12)
    ax.set(
        xlabel="Time (ms)",
        ylabel="Target error (µm)",
        title="(e) Same RPY integration; only predictor differs",
    )
    ax.legend(fontsize=8)
    ax.text(0.03, 0.03, "RPY: 106.5 ms sampled dwell; requirement ≥50 ms", transform=ax.transAxes, fontsize=8)
    ax = axes[0, 1]
    names = ["two", "four", "actual32"]
    ix = np.arange(3)
    ax.bar(
        ix - 0.17,
        [contact["rows"][n]["overlap_PBD_m"] * 1e9 for n in names],
        width=0.34,
        color="#aa7537",
        label="Original PBD",
    )
    ax.bar(
        ix + 0.17,
        [contact["rows"][n]["overlap_candidate_m"] * 1e9 for n in names],
        width=0.34,
        color="#326a9c",
        label="RPY candidate: 16 sweeps",
    )
    ax.set(
        xticks=ix,
        xticklabels=["Two bubbles", "Four asymmetric", "Actual 32"],
        ylabel="Residual overlap (nm)",
        title="(b) Contact geometry checks",
    )
    ax.legend(fontsize=8)
    ax = axes[1, 1]
    delta = np.array(contact["rows"]["actual32"]["mean_candidate_correction_m"]) * 1e9
    pbd_delta = np.array(contact["rows"]["actual32"]["mean_PBD_correction_m"]) * 1e9
    ax.bar(ix - 0.17, pbd_delta, width=0.34, color="#aa7537", label="Original PBD")
    ax.bar(ix + 0.17, delta, width=0.34, color="#326a9c", label="RPY candidate")
    ax.axhline(0, c="black", lw=0.6)
    ax.set_ylim(-2.6, 1.8)
    ax.legend(fontsize=8, loc="upper right")
    ax.set(
        xticks=ix,
        xticklabels=["X", "Y", "Z"],
        ylabel="Mean correction (nm)",
        title="(d) Actual 32: one 10 µs free-step proposal",
    )
    ax.text(
        0.03,
        0.03,
        "RPY mean correction follows M Jᵀλ;\nno imposed group shift; equal/opposite normal impulses.",
        transform=ax.transAxes,
        fontsize=8,
    )
    ax = axes[2, 1]
    ax.semilogy(
        np.arange(1, 17),
        np.maximum(cd["actual32_trace"][:, 0] * 1e9, 1e-8),
        c="#326a9c",
        label="16-sweep candidate",
    )
    ax.semilogy(
        np.arange(1, 65),
        np.maximum(cd["actual32_refined_trace"][:, 0] * 1e9, 1e-8),
        c="#326a9c",
        ls=":",
        label="64-sweep diagnostic",
    )
    ax.axhline(1.78, c="black", ls=":", label="Registered residual budget")
    ax.set(
        xlabel="Contact sweep (not physical time)",
        ylabel="Projected residual (nm)",
        title="(f) NOT ADMITTED: current cost exceeds budget",
    )
    cr = contact["rows"]["actual32"]
    ax.text(
        0.98,
        0.95,
        f"Setup + solve / PBD: {cr['setup_inclusive_ratio']:.1f}×\nAllowed ratio: 2×\nNo contact-candidate navigation run",
        ha="right",
        va="top",
        transform=ax.transAxes,
        fontsize=8,
    )
    ax.legend(fontsize=8, loc="lower left")
    for ax in axes.ravel():
        ax.grid(alpha=0.12)
    fig.suptitle("Single-factor controller ablation and an unpromoted RPY contact prototype", fontsize=16)
    save(
        fig,
        "Fig_F8_4_ablation_contact",
        "Matched optimizer/constraints/production with only candidate mobility changed; actual small-state contact overlap, mean displacement, residual and measured cost. Candidate not admitted to navigation.",
        [
            BASE / "analysis/plot_data.npz",
            BASE / "analysis/RESULTS.json",
            BASE / "a_static_v1/CHECKS.json",
            BASE / "contact_v1/CHECKS.json",
            BASE / "contact_v1/contact_states.npz",
        ],
    )


def main():
    OUT.mkdir(exist_ok=True)
    d = np.load(BASE / "analysis/plot_data.npz")
    r = js(BASE / "analysis/RESULTS.json")
    figure1()
    figure2(d, r)
    figure3(d, r)
    figure4(d, r)
    print("Four F8 PNG/PDF figures generated; visual inspection still required.")


if __name__ == "__main__":
    main()
