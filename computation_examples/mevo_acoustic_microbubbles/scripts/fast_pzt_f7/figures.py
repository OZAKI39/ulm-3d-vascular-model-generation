"""Four English review figures from saved F7 data; no smoothing or fake frames."""

import subprocess
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.collections import LineCollection
from .common import ROOT, BASE, js, write, sha

COLORS = dict(
    M0="#35679a",
    M1="#bf8737",
    M2="#a54259",
    old="#785891",
    new="#167d7e",
    failed="#858585",
    half="#167d7e",
    slow="#406382",
)
LABELS = dict(
    M0="M0: Stokes + lubrication",
    M1="M1: Stokes",
    M2="M2: RPY",
    old="RPY / F6 feedback",
    new="RPY / F7 constrained feedback",
    failed="RPY / F7 first attempt",
    half="F7, half time step",
    slow="F7, slow affine flow",
)
plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.labelcolor": "#222222",
        "text.color": "#222222",
        "xtick.color": "#222222",
        "ytick.color": "#222222",
        "axes.edgecolor": "#555555",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
    }
)


def style3(ax, lim):
    for axis, label in zip("xyz", ["X", "Y", "Z"]):
        getattr(ax, "set_" + axis + "lim")(*lim)
        getattr(ax, "set_" + axis + "label")(label + " offset (µm)", labelpad=6)
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=24, azim=-52)
    for a in [ax.xaxis, ax.yaxis, ax.zaxis]:
        a.pane.set_alpha(0.05)


def save(fig, name, caption, files):
    out = BASE / "figures"
    out.mkdir(exist_ok=True)
    fig.savefig(out / (name + ".png"), dpi=240, bbox_inches="tight")
    fig.savefig(out / (name + ".pdf"), bbox_inches="tight")
    plt.close(fig)
    write(
        out / (name + ".json"),
        dict(
            caption=caption,
            plot_commit=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
            sources={str(p.relative_to(ROOT)): sha(p) for p in files},
            input_hashes=js(ROOT / "results/fast_pzt/f5/INPUT_MANIFEST.json"),
            notes="Actual saved positions and times. Center offsets are coordinate origins only. Center markers are not radius-scaled spheres. No smoothing, interpolated trajectory frames, bubble-size inflation or gas fusion.",
        ),
    )


def figure1():
    fixed = BASE / "fixed_v1"
    d = np.load(fixed / "fixed_states.npz")
    r = js(fixed / "CHECKS.json")
    fig = plt.figure(figsize=(14, 8))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.2, 0.8], hspace=0.30, wspace=0.28)
    x = d["mature_positions_m"]
    center = x.mean(axis=0)
    q = (x - center) * 1e6
    for i, m in enumerate(["M0", "M1", "M2"]):
        ax = fig.add_subplot(gs[0, i], projection="3d")
        v = d["mature_" + m + "_velocity_m_s"] * 1e6
        ax.scatter(*q.T, s=10, c="#454545", depthshade=False)
        ax.quiver(*q.T, *(v * 0.0005).T, color=COLORS[m], arrow_length_ratio=0.2, linewidth=0.8)
        style3(ax, (-13, 13))
        ax.set_title(LABELS[m], fontsize=12, pad=5)
        ax.text2D(0.02, 0.95, "Mature cluster; same forces", transform=ax.transAxes, fontsize=9)
    for i, state in enumerate(["mature", "dispersed", "near_wall"]):
        ax = fig.add_subplot(gs[1, i])
        models = r["states"][state]["models"]
        vp = [models[m]["target_component_m_s"] * 1e6 for m in ["M0", "M1", "M2"]]
        vt = [models[m]["transverse_speed_m_s"] * 1e6 for m in ["M0", "M1", "M2"]]
        xx = np.arange(3)
        ax.bar(
            xx - 0.18, vp, 0.35, color=[COLORS[m] for m in ["M0", "M1", "M2"]], label="Toward target (signed)"
        )
        ax.bar(
            xx + 0.18,
            vt,
            0.35,
            facecolor="none",
            edgecolor=[COLORS[m] for m in ["M0", "M1", "M2"]],
            hatch="///",
            label="Transverse magnitude",
        )
        ax.axhline(0, c="#777777", lw=0.7)
        ax.set_xticks(xx, ["M0", "M1", "M2"])
        ax.set_ylabel("Group velocity (µm/s)")
        ax.set_title(state.replace("_", " ").capitalize())
        ax.set_ylim(-90, 145)
        ax.grid(axis="y", alpha=0.16)
        if i == 0:
            ax.legend(fontsize=8, loc="upper left")
    fig.suptitle("Same-state comparison isolates lubrication and collective mobility", fontsize=16, y=0.99)
    fig.text(
        0.5,
        0.025,
        "Top: actual bubble centers, common velocity scale (arrow displacement = velocity × 0.5 ms).\nBottom: group means at three saved states. Velocities precede PBD; near-wall RPY is uncalibrated.",
        ha="center",
        fontsize=10,
    )
    fig.subplots_adjust(bottom=0.14, top=0.89)
    save(
        fig,
        "Fig_F7_1_mobility_velocity_comparison",
        "M0/M1 have the same group mean by component balance but different relative motion; M2 changes drift and transverse response. Same-state, same force, pre-contact velocities.",
        [fixed / "fixed_states.npz", fixed / "CHECKS.json"],
    )


def target_sphere(ax, c, r=2):
    a = np.linspace(0, 2 * np.pi, 30)
    b = np.linspace(0, np.pi, 18)
    ax.plot_wireframe(
        c[0] + r * np.outer(np.cos(a), np.sin(b)),
        c[1] + r * np.outer(np.sin(a), np.sin(b)),
        c[2] + r * np.outer(np.ones_like(a), np.cos(b)),
        color="#777777",
        alpha=0.35,
        lw=0.5,
        rstride=3,
        cstride=3,
    )
    ax.scatter(*c, s=70, marker="*", c="black", label="Fixed target; 2 µm radius")


def figure2(d):
    fig = plt.figure(figsize=(15, 7.2))
    for panel, keys, title in [
        (1, ["M0", "M1", "M2"], "(a) Identical 40-event waveform replay"),
        (2, ["old", "failed", "new"], "(b) Feedback on the same RPY model"),
    ]:
        ax = fig.add_subplot(1, 2, panel, projection="3d")
        origin = d[keys[0] + "_c"][0] * 1e6
        allc = []
        for key in keys:
            c = d[key + "_c"] * 1e6 - origin
            allc.append(c)
            ax.plot(
                *c.T,
                color=COLORS[key],
                label=LABELS[key],
                lw=2.7 if key == "M0" else 1.8,
                ls="--" if key in ["M1", "failed"] else "-",
            )
            for t, marker in [(0.25, "^"), (0.5, "s"), (0.75, "D"), (1.0, "o")]:
                ii = np.argmin(abs(d[key + "_t"] - t))
                ax.scatter(*c[ii], marker=marker, s=23, c=COLORS[key], depthshade=False)
        target = d[keys[0] + "_target"] * 1e6 - origin
        target_sphere(ax, target)
        ax.scatter(0, 0, 0, c="black", marker="x", s=40, label="Common start")
        allc = np.vstack(allc + [target[None], np.zeros((1, 3))])
        lo = allc.min(axis=0) - 4
        hi = allc.max(axis=0) + 4
        for i, axis in enumerate("xyz"):
            getattr(ax, "set_" + axis + "lim")(lo[i], hi[i])
            getattr(ax, "set_" + axis + "label")("XYZ"[i] + " offset (µm)", labelpad=8)
        ax.set_box_aspect(hi - lo)
        ax.xaxis.set_major_locator(MaxNLocator(3))
        ax.yaxis.set_major_locator(MaxNLocator(5))
        ax.zaxis.set_major_locator(MaxNLocator(4))
        ax.view_init(elev=23, azim=-46)
        ax.set_title(title, pad=15, fontsize=12)
        ax.legend(loc="upper left", fontsize=8, framealpha=0.95)
    fig.suptitle(
        "Prescribed target navigation: open-loop physics and closed-loop control", fontsize=15, y=0.99
    )
    fig.text(
        0.5,
        0.025,
        "Coordinates relative to the common initial centroid. Time markers: ▲ 0.25 s, ■ 0.50 s, ◆ 0.75 s, ● 1.00 s.\nPanels use separate extents with metric aspect ratios. Geometric connectivity is not gas fusion.",
        ha="center",
        fontsize=10,
    )
    fig.subplots_adjust(left=0.03, right=0.97, bottom=0.13, top=0.90, wspace=0.07)
    save(
        fig,
        "Fig_F7_2_target_navigation",
        "Saved 3-D centroid paths with fixed 2µm target sphere. Separate causal replay and feedback comparisons; all paths cover1s.",
        [BASE / "analysis/plot_data.npz", BASE / "analysis/RESULTS.json"],
    )


def figure3(d, r):
    fig, axes = plt.subplots(3, 2, figsize=(14, 11), sharex=True, layout="constrained")
    ax = axes[0, 0]
    for key in ["old", "failed", "new", "half"]:
        ax.plot(
            d[key + "_t"] * 1e3,
            d[key + "_error"] * 1e6,
            c=COLORS[key],
            ls="--" if key in ["failed", "half"] else "-",
            lw=1.5,
            label=LABELS[key],
        )
    if "slow_t" in d:
        ax.plot(
            d["slow_t"] * 1e3,
            d["slow_error"] * 1e6,
            c=COLORS["slow"],
            ls=":",
            label="Slow affine: full validation",
        )
        ax.plot(
            d["slow_failed_t"][-1] * 1e3,
            d["slow_failed_error"][-1] * 1e6,
            marker="x",
            ms=8,
            c="black",
            label="Earlier slow attempt stopped",
        )
    ax.axhline(2, c="black", ls=":", label="Arrival threshold: 2 µm")
    nav = r["cases"]["new"]["navigation"]
    text = (
        "F7: no sustained arrival"
        if not nav["arrived_and_held"]
        else f"F7 dwell: {nav['longest_sampled_dwell_s'] * 1e3:.1f} ms"
    )
    ax.set_title("(a) Target error; required sampled dwell ≥ 50 ms")
    ax.set_ylabel("Error (µm)")
    ax.legend(fontsize=7, ncols=2)
    ax.text(0.98, 0.22, text, ha="right", transform=ax.transAxes, fontsize=9)
    ax = axes[0, 1]
    c = (d["new_c"] - d["new_c"][0]) * 1e6
    for i, col in enumerate(["#a54259", "#35679a", "#167d7e"]):
        ax.plot(d["new_t"] * 1e3, c[:, i], c=col, label="XYZ"[i])
    ax.axhline(-10, color="#777777", ls=":", lw=1)
    ax.set_title("(b) F7 displacement; target ΔY = −10 µm")
    ax.set_ylabel("Centroid displacement (µm)")
    ax.legend(fontsize=9)
    ax = axes[1, 0]
    for key in ["old", "new", "half"]:
        ax.plot(
            d[key + "_t"] * 1e3,
            d[key + "_fraction"],
            c=COLORS[key],
            ls="--" if key == "half" else "-",
            label=LABELS[key],
        )
    ax.axhline(0.9, c="black", ls=":")
    ax.set_ylim(0.85, 1.03)
    ax.set_title("(c) Largest geometric connected component")
    ax.set_ylabel("Fraction of 32 bubbles")
    ax.legend(fontsize=8)
    ax = axes[1, 1]
    for i, col in enumerate(["#35679a", "#a54259", "#8b7841", "#167d7e"]):
        ax.step(
            np.r_[d["new_control_t"], 1.0] * 1e3,
            np.r_[d["new_amplitudes"][:, i], d["new_amplitudes"][-1, i]],
            where="post",
            c=col,
            label=f"PZT {i + 1}",
        )
    ax.axhline(291.237687721856, c="black", ls=":", lw=1)
    ax.set_ylim(-5, 315)
    ax.set_title("(d) Actual 25 ms control events")
    ax.set_ylabel("Source amplitude (Pa)")
    ax.legend(ncols=2, fontsize=8)
    ax = axes[2, 0]
    for key in ["old", "new", "half"]:
        q = d[key + "_quality"]
        overlap = q[:, :, 3].max(axis=1) * 1e9
        ax.plot(
            d[key + "_t"] * 1e3, overlap, c=COLORS[key], ls="--" if key == "half" else "-", label=LABELS[key]
        )
    ax.axhline(17.8, c="black", ls=":", label="Unchanged 0.01R budget")
    ax.set_ylabel("Cumulative max overlap (nm)")
    ax.set_title("(e) Original contact quality gate")
    ax.legend(fontsize=8)
    ax = axes[2, 1]
    pred = d["new_predicted_velocity"] * 1e6
    actual = d["new_actual_segment_velocity"] * 1e6
    t = d["new_control_t"][: len(pred)] * 1e3
    for i, col in enumerate(["#a54259", "#35679a", "#167d7e"]):
        ax.plot(t, pred[:, i], c=col, lw=1, ls="--")
        ax.plot(t, actual[:, i], c=col, marker=".", ms=3, label="XYZ"[i])
    ax.set_ylabel("Group velocity (µm/s)")
    ax.set_title("(f) Dashed: instantaneous prediction; solid: next interval")
    ax.legend(fontsize=8)
    for ax in axes.ravel():
        ax.grid(alpha=0.15)
        ax.set_xlim(0, 1000)
    for ax in axes[-1]:
        ax.set_xlabel("Time (ms)")
    fig.suptitle("Feedback response, cohesion and numerical diagnostics", fontsize=16)
    save(
        fig,
        "Fig_F7_3_feedback_and_quality",
        "Real saved error, displacement, connectivity, amplitude events, contact budget and next-segment prediction residual. Dwell requires contiguous sampled arrival >=50ms, not just entering once.",
        [BASE / "analysis/plot_data.npz", BASE / "analysis/RESULTS.json"],
    )


def wall_intersection(y):
    data = np.load(ROOT / "results/fast_pzt/f5/inputs/f4/prepared_geometry.npz")
    vertices = data["vertices_m"]
    tri = vertices[data["faces"]]
    segments = []
    # True triangle intersections with this one named XZ plane; no inferred tube.
    for t in tri[(tri[:, :, 1].min(axis=1) <= y) & (tri[:, :, 1].max(axis=1) >= y)]:
        points = []
        for i, j in [(0, 1), (1, 2), (2, 0)]:
            a, b = t[i], t[j]
            if (a[1] - y) * (b[1] - y) < 0:
                points.append(a + (b - a) * ((y - a[1]) / (b[1] - a[1])))
        if len(points) == 2:
            segments.append(np.array(points)[:, [0, 2]] * 1e6)
    return segments


def figure4():
    p = BASE / "reachability_v1"
    r = js(p / "REACHABILITY.json")
    d = np.load(p / "Y_state.npz")
    x = d["positions_m"] * 1e6
    flow = d["flow_velocity_m_s"] * 1e6
    c = np.array(r["initial_centroid_m"]) * 1e6
    target = np.array(r["target_m"]) * 1e6
    fig, axes = plt.subplots(2, 2, figsize=(14, 10), layout="constrained")
    ax = axes[0, 0]
    ax.set_facecolor("#f4f0f8")
    ax.add_collection(
        LineCollection(
            wall_intersection(c[1] / 1e6),
            colors="#a55459",
            linewidths=2,
            label="Original wall / named XZ plane",
        )
    )
    ax.scatter(x[:, 0], x[:, 2], s=10, c="#3c587c", label="Actual bubble centers (projection)")
    ax.quiver(
        x[:, 0],
        x[:, 2],
        flow[:, 0],
        flow[:, 2],
        angles="xy",
        scale_units="xy",
        scale=20,
        width=0.0035,
        color="#406382",
    )
    ax.quiverkey(ax.collections[-1], 0.72, 0.90, 100, "100 µm/s", labelpos="E", coordinates="axes")
    ax.scatter(target[0], target[2], marker="*", s=120, c="black", label="Prescribed +Z target")
    ax.annotate("", target[[0, 2]], c[[0, 2]], arrowprops=dict(arrowstyle="->", color="black", lw=1.2))
    ax.set_xlim(c[0] - 18, c[0] + 23)
    ax.set_ylim(c[2] - 12, c[2] + 21)
    ax.set_aspect("equal")
    ax.ticklabel_format(useOffset=False, style="plain")
    ax.set_xlabel("Original X (µm)")
    ax.set_ylabel("Original Z (µm)")
    ax.set_title(f"(a) Original Y local region; wall plane Y = {c[1]:.3f} µm")
    ax.legend(loc="lower right", fontsize=7)
    ax.text(
        0.02,
        0.98,
        "Entire displayed window lies inside stored acoustic support",
        va="top",
        transform=ax.transAxes,
        fontsize=8,
    )
    ax = axes[0, 1]
    names = ["off", "unchanged", "target_phase", "X-", "Z+", "hold_fit"]
    rows = [q for q in r["rows"] if q["name"] in names]
    for row, color in zip(rows, ["#555555", "#bf8737", "#35679a", "#a54259", "#785891", "#167d7e"]):
        drift = np.array(row["models"]["M2"]["acoustic_mean_m_s"]) * 1e6
        ax.scatter(drift[0], drift[2], s=40, c=color, label=row["name"])
        ax.annotate(
            "",
            (drift[0], drift[2]),
            (0, 0),
            arrowprops=dict(arrowstyle="->", color="#785891", alpha=0.6, lw=1),
        )
    allv = (
        np.array([v["mean_velocity_m_s"] for v in r["candidates"]]) * 1e6 - np.array(r["flow_mean_m_s"]) * 1e6
    )
    ax.scatter(allv[:, 0], allv[:, 2], s=4, alpha=0.18, c="#785891", label="Tested optimizer candidates")
    uf = np.array(r["flow_mean_m_s"]) * 1e6
    ax.text(
        0.03,
        0.03,
        f"To cancel CFD: X = {-uf[0]:.2f} µm/s (outside displayed range)",
        transform=ax.transAxes,
        fontsize=9,
    )
    ax.set_xlim(-1, 4.5)
    ax.set_ylim(-9, 3)
    ax.axhline(0, c="#999999", lw=0.7)
    ax.axvline(0, c="#999999", lw=0.7)
    ax.set_xlabel("Acoustic group drift X (µm/s)")
    ax.set_ylabel("Acoustic group drift Z (µm/s)")
    ax.set_title("(b) RPY acoustic authority: XZ projection")
    ax.legend(fontsize=8, loc="best")
    ax.grid(alpha=0.15)
    ax = axes[1, 0]
    lo = np.array(r["tested_total_velocity_min_m_s"]) * 1e6
    hi = np.array(r["tested_total_velocity_max_m_s"]) * 1e6
    for i in range(3):
        ax.plot([i, i], [lo[i], hi[i]], lw=9, c="#785891", solid_capstyle="round")
        ax.scatter(i, uf[i], marker="x", s=70, c="black")
    ax.axhline(0, c="black", ls=":")
    ax.set_xticks(range(3), ["X", "Y", "Z"])
    ax.set_ylabel("Total group velocity (µm/s)")
    ax.set_title("(c) Observed candidate ranges; × = CFD alone")
    ax.grid(axis="y", alpha=0.15)
    ax = axes[1, 1]
    ax.axis("off")
    nb = r["neighborhood"]
    geom = nb["geometry"]
    status = "No full Y run admitted" if not r["full_Y_admitted"] else "Static admission passed"
    lines = [
        status,
        "",
        f"Mean CFD: ({uf[0]:.2f}, {uf[1]:.2f}, {uf[2]:.2f}) µm/s",
        f"Best tested total X: {lo[0]:.2f} µm/s",
        f"Best tested +Z progress: {hi[2]:.2f} µm/s",
        f"Target wall distance: {geom['distance_m'][1] * 1e6:.2f} µm",
        f"Target and six neighbors: field/CFD valid = {all(nb['field_valid']) and all(nb['flow']['valid'])}",
        "",
        "Finite candidates at the original F5 300 ms state.",
        "Ranges are not global controllability bounds.",
        "One XZ view is not a 3-D reachability proof.",
        "Free-space RPY has no calibrated wall mobility.",
    ]
    ax.text(0.04, 0.96, "\n".join(lines), va="top", fontsize=11, linespacing=1.55)
    fig.suptitle("Original Y-vessel: fixed-target reachability diagnosis", fontsize=16)
    save(
        fig,
        "Fig_F7_4_Y_target_reachability",
        "Original actual Y bubble coordinates and CFD, finite bounded acoustic candidates, XYZ drift ranges and target-support checks. Local XZ projection explicitly limited, not global3D proof.",
        [
            p / "REACHABILITY.json",
            p / "Y_state.npz",
            ROOT / "results/fast_pzt/f5/inputs/f4/prepared_geometry.npz",
        ],
    )


def main():
    d = np.load(BASE / "analysis/plot_data.npz")
    r = js(BASE / "analysis/RESULTS.json")
    figure1()
    figure2(d)
    figure3(d, r)
    figure4()
    print("Saved four PNG/PDF figures; visual inspection still required.")


if __name__ == "__main__":
    main()
