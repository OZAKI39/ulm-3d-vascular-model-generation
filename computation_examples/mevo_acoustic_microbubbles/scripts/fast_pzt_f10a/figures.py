"""Four scientific figures from unmodified actual results; no trajectory smoothing."""

import time
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import h5py
from .common import BASE, ROOT, js, write, sha, identity

COLORS = ["#38598b", "#398979", "#b47835", "#954b6d"]
plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "legend.fontsize": 8.5,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.18,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
    }
)


def label(ax, title, x=None, y=None):
    ax.set_title(title, loc="left", fontweight="bold", pad=10)
    if x:
        ax.set_xlabel(x)
    if y:
        ax.set_ylabel(y)


def trajectory(campaign, condition):
    with h5py.File(BASE / campaign / condition / "trajectory.h5") as h:
        t = h["time_s"][:]
        x = h["positions_m"][:, 0]
        d = h["motion_ledger/displacement_m"][:]
        meta = h["motion_ledger/step_meta"][:]
        return dict(
            t=t,
            c=x.mean(1),
            rg=h["metrics/global_rg_m"][:, 0],
            fraction=h["metrics/largest_size"][:, 0] / 32,
            d=d,
            meta=meta,
        )


def main():
    start = time.perf_counter()
    out = BASE / "figures"
    out.mkdir(exist_ok=True)
    rad = js(BASE / "radial_v1/RESULTS.json")
    stat = js(BASE / "static_v1/RESULTS.json")
    res = js(BASE / "analysis/RESULTS.json")
    sources = [
        BASE / "radial_v1/RESULTS.json",
        BASE / "static_v1/RESULTS.json",
        BASE / "analysis/RESULTS.json",
        BASE / "support_v1/CFD_cross_section.npz",
        ROOT / "results/fast_pzt/f5/inputs/f4/prepared_geometry.npz",
    ]
    sources += list((BASE / "radial_v1").glob("*/cycles.npz"))
    sources += (
        list((BASE / "pilot_v3").glob("*/trajectory.h5"))
        + list((BASE / "extended_10_v1").glob("*/trajectory.h5"))
        + list((BASE / "extended_5_v1").glob("*/trajectory.h5"))
    )
    source_record = dict(
        plotting_identity=identity(),
        inputs={str(p.relative_to(ROOT)): sha(p) for p in sources},
        execution_commits={
            k: (
                js(BASE / k / "CHECKS.json")["identity"]
                if k == "support_v1"
                else js(BASE / k / "IDENTITY.json")
            )["execution_commit"]
            for k in ["radial_v1", "static_v1", "pilot_v3", "extended_10_v1", "extended_5_v1", "support_v1"]
        },
        trajectory_rendering="Unsmoothed stored samples, straight connecting segments; no invented intermediate particle coordinates",
    )
    write(out / "FIGURE_SOURCES.json", source_record)
    produced = []

    def save(fig, name, scope):
        fig.text(
            0.01,
            0.014,
            scope + "  |  EXPLORATORY; no hardware / survival calibration",
            fontsize=8,
            color="#404040",
        )
        fig.tight_layout(rect=(0, 0.045, 1, 0.97), pad=1.5)
        for suffix in ["png", "pdf"]:
            path = out / (name + "." + suffix)
            fig.savefig(
                path, dpi=210, facecolor="white", metadata={"Title": name} if suffix == "pdf" else None
            )
            produced.append(path.name)
        plt.close(fig)

    prod = [r for r in rad["rows"] if r["name"].startswith("production_") and r["samples_per_cycle"] == 512]
    p = np.array([r["pressure_Pa"] for r in prod]) / 1000
    fig, ax = plt.subplots(2, 2, figsize=(11.5, 8))
    fig.suptitle("F10A | Pressure response and shell assumptions", fontsize=15, fontweight="bold")
    for color, r in zip(COLORS, prod):
        z = np.load(BASE / "radial_v1" / r["name"] / "cycles.npz")
        t = (z["time_s"] - z["time_s"][0]) * 1e6
        take = t <= 0.80001
        ax[0, 0].plot(
            t[take], z["radius_m"][take] * 1e6, color=color, label=f"{r['pressure_Pa'] / 1000:g} kPa"
        )
    for value, style, text in [
        (prod[0]["Rb_m"], "--", r"$R_b$"),
        (1.78e-6, ":", r"$R_0$"),
        (prod[0]["Rr_m"], "--", r"$R_r$"),
    ]:
        ax[0, 0].axhline(value * 1e6, color="#555555", ls=style, lw=1)
        ax[0, 0].text(0.82, value * 1e6, text, va="center")
    ax[0, 0].set_xlim(0, 0.91)
    ax[0, 0].legend(loc="upper center", ncol=4)
    label(ax[0, 0], "a  Direct production-shell waveforms", "Time within steady window (µs)", "Radius (µm)")
    for family, style, col in [("production", "-", COLORS[0]), ("TuTable2", "--", COLORS[2])]:
        rr = [
            r
            for r in rad["rows"]
            if r["name"].startswith(family + "_")
            and r["sigma_initial_N_m"] == 0.036
            and r["samples_per_cycle"] == 512
        ]
        amp = np.array([(r["Rmax_m"] - r["Rmin_m"]) / 2 for r in rr]) * 1e9
        ax[0, 1].plot(
            p,
            amp,
            marker="o",
            ls=style,
            color=col,
            label="Production Fig. 4" if family == "production" else "Separate Tu Table 2",
        )
    ax[0, 1].set_xscale("log")
    ax[0, 1].set_xticks(p, [str(v) for v in p])
    ax[0, 1].legend()
    label(
        ax[0, 1],
        "b  Radius half-range; two distinct parameter sets",
        "Local peak pressure (kPa)",
        "(Rmax − Rmin) / 2 (nm)",
    )
    cases = [
        next(
            r
            for r in rad["rows"]
            if r["shell_family"] == "TuTable2" and r["pressure_Pa"] == pp and r["sigma_initial_N_m"] == sigma
        )
        for pp in [3000, 10000]
        for sigma in [0, 0.036, 0.072]
    ]
    bottom = np.zeros(6)
    for branch, col in zip(["buckled", "elastic", "ruptured"], [COLORS[0], COLORS[1], COLORS[3]]):
        v = np.array([r["branch_fractions"][branch] for r in cases])
        ax[1, 0].bar(np.arange(6), v, bottom=bottom, color=col, label=branch)
        bottom += v
    ax[1, 0].set_xticks(
        np.arange(6),
        [f"{r['pressure_Pa'] / 1000:g} / {r['sigma_initial_N_m']:g}" for r in cases],
        rotation=25,
        ha="right",
    )
    ax[1, 0].set_ylim(0, 1.23)
    ax[1, 0].legend(ncol=3, loc="upper center")
    label(
        ax[1, 0],
        "c  Separate initial-tension assumptions (Tu Table 2)",
        "Pressure (kPa) / initial tension (N/m)",
        "Fraction of steady cycles",
    )
    ax[1, 1].plot(
        p,
        [100 * r["table_fundamental_complex_relative_error"] for r in prod],
        "o-",
        color=COLORS[0],
        label="Complex fundamental",
    )
    ax[1, 1].plot(
        p,
        [100 * r["table_dotV_energy_relative_error"] for r in prod],
        "s-",
        color=COLORS[2],
        label="Volume-rate energy",
    )
    ax[1, 1].axhline(5, ls="--", color="#555", label="Registered 5% admission")
    ax[1, 1].set_xscale("log")
    ax[1, 1].set_ylim(0, 5.6)
    ax[1, 1].set_xticks(p, [str(v) for v in p])
    ax[1, 1].legend(loc="center left")
    label(
        ax[1, 1],
        "d  Production table versus direct APECSS",
        "Local peak pressure (kPa)",
        "Relative difference (%)",
    )
    save(
        fig,
        "Fig1_radial_response",
        "R0=1.78 µm; f=2.5 MHz; main σinit=0.036 N/m; main cases remain elastic. APECSS "
        + source_record["execution_commits"]["radial_v1"][:9],
    )
    fig, axes = plt.subplots(2, 3, figsize=(13, 8))
    fig.suptitle("F10A | More pressure does not select the right direction", fontsize=15, fontweight="bold")
    axes = axes.ravel()
    selected = {
        mode: [
            next(
                r
                for r in stat["rows"]
                if r["state"] == "interior" and r["mode"] == mode and r["pressure_Pa"] == int(pp * 1000)
            )
            for pp in p
        ]
        for mode in ["S", "O"]
    }
    for mode, col, mark in [("S", COLORS[2], "s"), ("O", COLORS[0], "o")]:
        rows = selected[mode]
        name = "S: original phases" if mode == "S" else "O: bounded direction search"
        axes[0].plot(p, [r["primary_mean_N"][2] * 1e12 for r in rows], marker=mark, color=col, label=name)
        axes[1].loglog(p, [r["secondary_rms_N"] * 1e12 for r in rows], marker=mark, color=col, label=mode)
        axes[2].plot(
            p, [r["mean_acoustic_M0_m_s"][2] * 1e6 for r in rows], marker=mark, color=col, label=mode
        )
        axes[4].loglog(
            p, [max(r["event"]["amplitude_Pa"]) / 1000 for r in rows], marker=mark, color=col, label=mode
        )
        y = np.array([r["pressure_statistics_Pa"]["median"] for r in rows]) / 1000
        lo = np.array([r["pressure_statistics_Pa"]["min"] for r in rows]) / 1000
        hi = np.array([r["pressure_statistics_Pa"]["max"] for r in rows]) / 1000
        axes[5].errorbar(
            p * (0.97 if mode == "S" else 1.03),
            y,
            yerr=[y - lo, hi - y],
            fmt=mark,
            color=col,
            capsize=4,
            label=mode + " bubble min–max",
        )
    rr = selected["S"]
    for key, col, lbl in [
        ("primary_rms_N", COLORS[2], "S primary RMS"),
        ("secondary_rms_N", COLORS[0], "S secondary RMS"),
    ]:
        ratios = np.array([r[key] / rr[0][key] / (pp / p[0]) ** 2 for r, pp in zip(rr, p)])
        axes[3].plot(p, ratios, marker="o", color=col, label=lbl)
    axes[3].plot(
        p,
        [v["equal_in_phase_far_pair_over_quadratic"] for v in res["direct_radial_scaling"]],
        "^--",
        color=COLORS[1],
        label="Direct radial pair response",
    )
    axes[3].axhline(1, color="#444", ls=":")
    axes[3].set_ylim(0.995, 1.005)
    axes[4].axhline(0.2912377, color="#555", ls=":", label="Historical per-source limit")
    axes[4].axhline(20, color="#555", ls="--", label="Exploratory per-source cap")
    for k in [0, 2]:
        axes[k].axhline(0, color="#666", lw=0.8)
        axes[k].set_xscale("log")
        axes[k].set_yscale("symlog", linthresh=1 if k == 2 else 0.01)
        axes[k].margins(y=0.15)
    for k in [3, 5]:
        axes[k].set_xscale("log")
    axes[5].set_yscale("log")
    titles = [
        "a  Signed primary force",
        "b  Internal secondary force",
        "c  M0 acoustic group drift",
        "d  Quadratic scaling check",
        "e  Source commands, not voltage",
        "f  Actual local pressure spread",
    ]
    yl = [
        "Mean primary Fz (pN)",
        "Secondary force RMS (pN)",
        "Mean acoustic Uz (µm/s)",
        "Response / 300 Pa quadratic reference",
        "Maximum source amplitude (kPa)",
        "Bubble peak pressure (kPa)",
    ]
    for k, ax in enumerate(axes):
        label(ax, titles[k], "Initial median pressure (kPa)", yl[k])
        ax.legend(fontsize=7.8)
        ax.set_xticks(p, ["0.3", "1", "3", "10"])
    budget = [r["source_budget_Pa2"] / r["S_budget_Pa2"] for r in selected["O"]]
    axes[5].text(
        0.04,
        0.95,
        "O/S command² budget: " + ", ".join(f"{v:.3f}" for v in budget),
        transform=axes[5].transAxes,
        va="top",
        fontsize=7.5,
    )
    save(
        fig,
        "Fig2_pressure_force_control",
        "Same interior 32-ID state; O median matched, source²≤S, max pressure≤1.01×S. Static "
        + source_record["execution_commits"]["static_v1"][:9],
    )
    # Exact triangle-plane intersection of the original wall: no ideal circular tube substitute.
    geom = np.load(ROOT / "results/fast_pzt/f5/inputs/f4/prepared_geometry.npz")
    verts = geom["vertices_m"]
    faces = geom["faces"][geom["labels"] == 1]
    plane = np.array(stat["states"]["near"]["positions_m"]).mean(0)[0]
    segments = []
    for triangle in verts[faces]:
        delta = triangle[:, 0] - plane
        if delta.min() > 0 or delta.max() < 0:
            continue
        points = []
        for i, j in [(0, 1), (1, 2), (2, 0)]:
            if delta[i] * delta[j] < 0:
                points.append(triangle[i] + (triangle[j] - triangle[i]) * (-delta[i]) / (delta[j] - delta[i]))
        if len(points) == 2:
            segments.append(np.array(points)[:, 1:] * 1000)
    fig = plt.figure(figsize=(13, 8))
    fig.suptitle("F10A | Wall clearance versus flow and acoustic authority", fontsize=15, fontweight="bold")
    grid = fig.add_gridspec(2, 3, width_ratios=[1.25, 1, 1])
    ax = fig.add_subplot(grid[:, 0])
    q = np.load(BASE / "support_v1/CFD_cross_section.npz")
    shape = q["grid_shape"]
    v = q["velocity_m_s"][:, 0] * 1000
    mask = q["valid"].astype(bool) & q["inside"].astype(bool)
    v = np.ma.array(v, mask=~mask).reshape(shape)
    xyz = q["positions_m"].reshape(*shape, 3)
    img = ax.pcolormesh(
        xyz[:, :, 1] * 1000, xyz[:, :, 2] * 1000, v, cmap="Purples", shading="auto", rasterized=True
    )
    ax.add_collection(LineCollection(segments, colors="#852f42", linewidths=1.5))
    fig.colorbar(img, ax=ax, orientation="horizontal", pad=0.1, label="Actual CFD Ux (mm/s)", shrink=0.9)
    names = ["near", "middle", "interior"]
    centers = []
    for i, name in enumerate(names):
        xx = np.array(stat["states"][name]["positions_m"])
        centers.append(xx.mean(0))
        ax.scatter(xx[:, 1] * 1000, xx[:, 2] * 1000, s=3, color=COLORS[i], zorder=4)
        ax.annotate(
            name.capitalize(),
            xy=xx.mean(0)[1:] * 1000,
            xytext=(0.30, xx.mean(0)[2] * 1000 + 0.12),
            arrowprops=dict(arrowstyle="-", color=COLORS[i]),
            color=COLORS[i],
            fontsize=10,
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=2),
        )
    ax.set_aspect("equal")
    ax.set_xlim(-1.45, 1.45)
    ax.set_ylim(-1.45, 1.45)
    label(ax, "a  Original wall / equal spatial scale", "Y (mm)", "Z (mm)")
    ax.text(
        0.02,
        0.98,
        f"Cross-section X={plane * 1000:.3f} mm\nDots: actual 32 bubble centers",
        transform=ax.transAxes,
        va="top",
        fontsize=8,
    )
    axes = [
        fig.add_subplot(grid[0, 1]),
        fig.add_subplot(grid[0, 2]),
        fig.add_subplot(grid[1, 1]),
        fig.add_subplot(grid[1, 2]),
    ]
    vals = [
        [stat["states"][n]["minimum_wall_gap_m"] * 1e6 for n in names],
        [stat["states"][n]["mean_flow_m_s"][0] * 1000 for n in names],
        [np.mean(stat["states"][n]["shear_rate_s_inv"]) for n in names],
        [
            next(
                r["mean_acoustic_M0_m_s"][2] * 1e6
                for r in stat["rows"]
                if r["state"] == n and r["mode"] == "O" and r["pressure_Pa"] == 3000
            )
            for n in names
        ],
    ]
    for i, (aa, vv, title, unit) in enumerate(
        zip(
            axes,
            vals,
            [
                "b  Minimum finite-radius gap",
                "c  Downstream advection",
                "d  Local strain-rate magnitude",
                "e  Optimized acoustic drift (3 kPa)",
            ],
            ["Wall gap (µm)", "Mean CFD Ux (mm/s)", "sqrt(2 D:D) (1/s)", "Acoustic Uz (µm/s)"],
        )
    ):
        aa.bar(names, vv, color=COLORS[:3])
        label(aa, title, None, unit)
        for j, value in enumerate(vv):
            aa.text(j, value, f"{value:.2f}", ha="center", va="bottom" if value >= 0 else "top", fontsize=9)
        aa.margins(y=0.22)
    times = [20e-6 / stat["states"][n]["mean_flow_m_s"][0] * 1000 for n in names]
    axes[1].text(
        0.02,
        0.94,
        "20 µm transit: " + ", ".join(f"{t:.2f}" for t in times) + " ms",
        transform=axes[1].transAxes,
        va="top",
        fontsize=8,
    )
    save(
        fig,
        "Fig3_wall_flow_tradeoff",
        "Rigidly translated F9 cloud; all current 32-bubble supports legal. Interior1/3 kPa static admission; 3 kPa trajectory later rejected.",
    )
    fig, axes = plt.subplots(2, 3, figsize=(13, 8))
    axes = axes.ravel()
    fig.suptitle(
        "F10A | Resolved group drift with substantial internal contact", fontsize=15, fontweight="bold"
    )
    data = {key: trajectory("extended_10_v1", key) for key in ["OFF", "S_1000", "O_1000"]}
    fine = {key: trajectory("extended_5_v1", key) for key in data}
    for key, col in zip(data, ["#666666", COLORS[2], COLORS[0]]):
        d = data[key]
        f = fine[key]
        name = key.replace("_1000", " (1 kPa)")
        axes[0].plot(d["t"] * 1000, (d["c"][:, 2] - d["c"][0, 2]) * 1e6, color=col, label=name)
        axes[0].plot(f["t"] * 1000, (f["c"][:, 2] - f["c"][0, 2]) * 1e6, color=col, ls="--", lw=1)
        axes[2].plot(d["t"] * 1000, d["rg"] * 1e6, color=col, label=name)
        axes[2].plot(f["t"] * 1000, f["rg"] * 1e6, color=col, ls="--", lw=1)
        axes[3].plot(d["t"] * 1000, abs(d["c"][:, 2] - f["c"][:, 2]) * 1e9, color=col, label=name)
        if len(d["meta"]):
            axes[4].plot(
                (d["meta"][:, 0] + d["meta"][:, 1]) * 1000,
                np.cumsum(d["meta"][:, 4]) * 1e6,
                color=col,
                label=name,
            )
    for other, col in [("OFF", COLORS[0]), ("S_1000", COLORS[1])]:
        axes[1].plot(
            data["O_1000"]["t"] * 1000,
            (data["O_1000"]["c"][:, 2] - data[other]["c"][:, 2]) * 1e6,
            color=col,
            label="O − " + other.replace("_1000", ""),
        )
    weak = trajectory("pilot_v3", "O_300")
    weak_off = trajectory("pilot_v3", "OFF")
    axes[1].plot(
        weak["t"] * 1000,
        (weak["c"][:, 2] - weak_off["c"][:, 2]) * 1e6,
        ls=":",
        color=COLORS[2],
        label="0.3 kPa: O − OFF (2 ms)",
    )
    axes[1].axhline(0.05, color="#555", ls="--", label="Registered resolution floor")
    axes[2].text(
        0.04,
        0.95,
        "Largest connected fraction = 1.0\nfor every saved sample shown",
        transform=axes[2].transAxes,
        va="top",
        fontsize=8,
    )
    axes[4].text(
        0.04,
        0.95,
        "Wall correction = 0\nNet centroid contact ≈ 0\nLarge individual contact paths persist",
        transform=axes[4].transAxes,
        va="top",
        fontsize=8,
    )
    for i, title, y in [
        (0, "a  Actual fixed-ID centroid", "Z displacement (µm)"),
        (1, "b  Paired acoustic increments", "Additional Z displacement (µm)"),
        (2, "c  Shape and connectivity", "Radius of gyration (µm)"),
        (3, "d  New 10 / 5 µs sensitivity", "Absolute centroid Z difference (nm)"),
        (4, "e  Unsigned individual PBD path", "Mean cumulative contact path (µm)"),
    ]:
        label(axes[i], title, "Time (ms)", y)
        axes[i].legend(loc="best", fontsize=7.8)
    axes[0].text(0.02, 0.94, "Solid: 10 µs; dashed: 5 µs", transform=axes[0].transAxes, va="top", fontsize=8)
    axes[5].axis("off")
    label(axes[5], "f  Quality and interpretation")
    effect = res["effects"]["extended_5_v1"]["OFF"] * 1e6
    sensitivity = res["timestep"]["O_1000"]["max_abs_Z_difference_m"] * 1e9
    lines = [
        f"1 kPa paired increment: +{effect:.3f} µm",
        f"O step sensitivity: {sensitivity:.2f} nm",
        "1 kPa: all six 10 ms runs PASS",
        "3 kPa: S and O FAIL in first 50 µs",
        "Reason: residual overlap >0.01 R",
        "No accepted 3 kPa trajectory fabricated",
        f"Individual step difference: {res['timestep']['O_1000']['max_particle_norm_difference_m'] * 1e6:.2f} µm",
        f"ACOUSTIC_EFFECT_RESOLVED = {str(res['flags']['ACOUSTIC_EFFECT_RESOLVED']).lower()}",
        f"CONTACT_DOMINATED = {str(res['flags']['CONTACT_DOMINATED']).lower()}",
        "WALL_DOMINATED = false",
        "NAVIGATION_NOT_ESTABLISHED",
        "TRAINING_READY = false",
    ]
    axes[5].text(
        0, 0.93, "\n".join(lines), transform=axes[5].transAxes, va="top", fontsize=9, linespacing=1.7
    )
    save(
        fig,
        "Fig4_short_trajectories_quality",
        "Same initial 32 IDs, interior CFD; fixed controls, initial median1 kPa. Trajectory "
        + source_record["execution_commits"]["extended_10_v1"][:9],
    )
    write(
        out / "PLOT_RUN.json",
        dict(
            identity=identity(),
            files=produced,
            wall_s=time.perf_counter() - start,
            source_manifest_sha256=sha(out / "FIGURE_SOURCES.json"),
        ),
    )


if __name__ == "__main__":
    main()
