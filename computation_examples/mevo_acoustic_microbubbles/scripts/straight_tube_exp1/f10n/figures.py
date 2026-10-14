"""Four data-driven English figures; no simulated success illustrations."""

import shutil
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from mevo_mb.acoustics.water_disc_proxy import CachedWaterBasis
from mevo_mb.fast_sim.models import ResponseSection
from mevo_mb.fast_sim.nearfield import F2Physics
from .common import BASE, ROOT, configs, read, write, sha, identity
from .analyze import load_run

COLORS = {"OFF": "#60666b", "BASELINE": "#be7b30", "OPTIMIZED_REPLAY": "#247ab4", "SHAPE_REPLAY": "#774c9e"}
LABELS = {
    "OFF": "OFF",
    "BASELINE": "BASELINE",
    "OPTIMIZED_REPLAY": "OPTIMIZED",
    "SHAPE_REPLAY": "SHAPE CONTROL",
}


def style():
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.titlesize": 12,
            "axes.labelsize": 11,
            "axes.edgecolor": "#33383d",
            "axes.labelcolor": "#23272b",
            "xtick.color": "#33383d",
            "ytick.color": "#33383d",
            "text.color": "#23272b",
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


def spheres(ax, x, r, title):
    x = (x - x.mean(0)) * 1e6
    r *= 1e6
    u = np.linspace(0, 2 * np.pi, 25)
    v = np.linspace(0, np.pi, 17)
    sphere = (
        np.stack(
            np.broadcast_arrays(
                np.cos(u)[:, None] * np.sin(v),
                np.sin(u)[:, None] * np.sin(v),
                np.ones((len(u), 1)) * np.cos(v),
            ),
            axis=-1,
        )
        * r
    )
    for point in x:
        q = sphere + point
        ax.plot_surface(
            q[:, :, 0],
            q[:, :, 1],
            q[:, :, 2],
            color="#368fc6",
            edgecolor="none",
            shade=True,
            antialiased=True,
            rstride=1,
            cstride=1,
        )
    ax.set(
        xlim=(-8, 8),
        ylim=(-8, 8),
        zlim=(-8, 8),
        xlabel="X (µm)",
        ylabel="Y (µm)",
        zlabel="Z (µm)",
        title=title,
    )
    ax.set_box_aspect((1, 1, 1), zoom=0.78)
    ax.view_init(25, -55)
    ax.tick_params(labelsize=8, pad=0)
    for axis in [ax.xaxis, ax.yaxis, ax.zaxis]:
        axis.pane.fill = False
        axis.labelpad = 1


def main():
    style()
    p = configs()["protocol"]
    accepted = read(BASE / "ACCEPTED_CAMPAIGN.json")
    out = BASE / accepted["path"]
    folder = BASE / "figures"
    folder.mkdir(exist_ok=True)
    gate = read(BASE / "analysis/FINAL_GATE.json")
    cmp = read(BASE / "analysis/NUMERICAL_SENSITIVITY.json")
    campaign = read(out / "CAMPAIGN.json")
    rows = read(BASE / "analysis/CONTROL_COMPARISON.json")
    records = {r["name"]: r for r in rows}
    runs = {n: load_run(out, n) for n in ["OFF", "BASELINE", "OPTIMIZED_REPLAY", "SHAPE_REPLAY"]}
    provenance = []

    def save(fig, name, inputs):
        fig.savefig(folder / (name + ".png"), dpi=240)
        fig.savefig(folder / (name + ".pdf"))
        plt.close(fig)
        provenance.append(
            dict(
                figure=name,
                file_sha256={ext: sha(folder / (name + "." + ext)) for ext in ["png", "pdf"]},
                source_sha256={str(path.relative_to(ROOT)): sha(path) for path in inputs},
                identity=identity(),
            )
        )

    shared = [
        BASE / "CONFIG_SNAPSHOT.json",
        BASE / "PREREGISTRATION.json",
        BASE / "analysis/FINAL_GATE.json",
        BASE / "analysis/CONTROL_COMPARISON.json",
        BASE / "analysis/NUMERICAL_SENSITIVITY.json",
    ]
    shared += [out / "runs" / n / "metrics.npz" for n in [*runs, "OPTIMIZED_HALF_DT", "OPTIMIZED_FINE_CACHE"]]
    shared += [BASE / "inputs/response.h5"]

    fig = plt.figure(figsize=(14, 9), layout="constrained")
    gs = fig.add_gridspec(2, 2, height_ratios=[0.95, 1.05])
    ax = fig.add_subplot(gs[0, :])
    ax.axis("off")
    apparatus = ROOT / "results/straight_tube_exp1/f10m_r3/figures/apparatus_render.png"
    shutil.copy2(apparatus, folder / "apparatus_geometry_reference.png")
    ax.imshow(plt.imread(apparatus))
    ax.set_title("A  Four 7 mm sources on the sealed-X, dual-port geometry", loc="left")
    ax.text(
        0.99,
        0.04,
        "Homogeneous-water control proxy\n32.5 × 9 × 9 mm; bore ID 1.6 mm",
        ha="right",
        va="bottom",
        transform=ax.transAxes,
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.85),
    )
    ax = fig.add_subplot(gs[1, 0])
    basis = CachedWaterBasis(BASE / "cache/coarse.h5")
    ev = read(out / "runs/OPTIMIZED_REPLAY/INPUT.json")["events"][0]
    u = basis.vector(ev["amplitude_Pa"], ev["phase_rad"])
    q = np.linspace(-40e-6, 40e-6, 81)
    xx, yy = np.meshgrid(q, q)
    points = np.c_[xx.ravel(), yy.ravel(), np.zeros(xx.size)]
    pb, gb = basis.query_basis(points)
    P = u @ pb
    G = np.einsum("i,inj->nj", u, gb)
    response = ResponseSection(
        BASE / "inputs/response.h5", F2Physics(radius_m=p["radius_m"], frequency_Hz=p["frequency_Hz"])
    )
    v = response.query(P)
    F = -0.5 * np.real(v[:, 0, None] * G.conj())
    im = ax.pcolormesh(
        xx * 1e6, yy * 1e6, abs(P).reshape(xx.shape) / 1000, cmap="magma", shading="auto", rasterized=True
    )
    Q = F.reshape((*xx.shape, 3))
    step = 12
    ax.quiver(
        xx[::step, ::step] * 1e6,
        yy[::step, ::step] * 1e6,
        Q[::step, ::step, 0],
        Q[::step, ::step, 1],
        color="white",
        width=0.004,
    )
    initial = runs["OPTIMIZED_REPLAY"][1]["positions_m"][0, 0] * 1e6
    ax.scatter(
        initial[:, 0], initial[:, 1], s=7, c="#56b7ec", edgecolors="none", label="Initial bubble centres"
    )
    ax.set(
        xlabel="X (µm)",
        ylabel="Y (µm)",
        aspect="equal",
        title="B  Local XY field; arrows: primary-force direction",
    )
    ax.legend(loc="upper right", fontsize=8, facecolor="white", framealpha=0.9)
    cb = fig.colorbar(im, ax=ax, shrink=0.8)
    cb.set_label(r"$|\,P\,|$ (kPa)")
    ax = fig.add_subplot(gs[1, 1])
    idx = np.arange(4)
    ax.bar(idx, ev["amplitude_Pa"], color="#819aaa", width=0.58)
    ax.set(
        xticks=idx,
        xticklabels=["Y−", "Y+", "X−", "X+"],
        ylabel="Equivalent source amplitude (Pa)",
        title="C  Selected fixed control: 2.5 MHz, 2 ms",
    )
    for i, (a, ph) in enumerate(zip(ev["amplitude_Pa"], ev["phase_rad"])):
        ax.text(i, a + 45, f"{a:.1f} Pa\n{np.degrees(ph):.1f}°", ha="center", va="bottom", fontsize=9)
    ax.set_ylim(0, max(ev["amplitude_Pa"]) * 1.35 + 100)
    ax.text(
        0.02,
        0.96,
        "Amplitude = ρc vₙ; phase convention exp(+iωt)\nIndependent amplitudes, phases and enable flags\nNo PDMS reflection or voltage calibration",
        transform=ax.transAxes,
        va="top",
        fontsize=9,
    )
    fig.suptitle("Four-PZT control field | HOMOGENEOUS_WATER_CONTROL_PROXY", fontsize=15)
    save(
        fig,
        "Fig_F10N_1_four_pzt_control_field",
        shared + [BASE / "cache/coarse.h5", apparatus, out / "runs/OPTIMIZED_REPLAY/INPUT.json"],
    )

    fig, axs = plt.subplots(2, 2, figsize=(13.5, 9), layout="constrained")
    target = np.array(p["target_position_m"]) * 1e6
    for name in ["OFF", "BASELINE", "OPTIMIZED_REPLAY"]:
        rec, a, m = runs[name]
        c = m["shape_centroid_m"] * 1e6
        t = a["time_s"] * 1000
        err = np.linalg.norm(c - target, axis=-1)
        for ax in axs[0]:
            ax.plot(c[:, 0], c[:, 1], color=COLORS[name], label=LABELS[name], lw=2)
        axs[0, 1].scatter(c[-1, 0], c[-1, 1], c=COLORS[name], s=25)
        axs[1, 0].plot(t, err, color=COLORS[name], label=LABELS[name], lw=2)
    axs[0, 0].scatter(*target[:2], marker="*", s=170, color="#a6273c", label="Registered target")
    axs[0, 0].plot([0, target[0]], [0, target[1]], "--", color="#aaaaaa", lw=1)
    axs[0, 0].set(
        title="A  Goal fixed before all trajectories",
        xlabel="COM X (µm)",
        ylabel="COM Y (µm)",
        aspect="equal",
        xlim=(-0.25, 3.4),
        ylim=(-0.3, 2.5),
    )
    axs[0, 0].text(
        0.04,
        0.92,
        "3D goal: (3, 2, 1) µm\nInitial error: 3.742 µm\nTolerance: 0.75 µm",
        transform=axs[0, 0].transAxes,
        va="top",
    )
    axs[0, 0].legend(loc="lower right", fontsize=9)
    axs[0, 1].set(
        title="B  Actual COM motion, enlarged scale", xlabel="COM X (µm)", ylabel="COM Y (µm)", aspect="equal"
    )
    axs[0, 1].grid(alpha=0.2)
    axs[0, 1].set(xlim=(-0.02, 0.25), ylim=(-0.04, 0.12))
    axs[1, 0].set(title="C  Full 3D target error", xlabel="Time (ms)", ylabel="Target error (µm)")
    axs[1, 0].text(
        0.03,
        0.08,
        f"Final optimized error: {gate['optimized_error_um']:.4f} µm\nTarget reached: {str(gate['target_reached']).upper()}",
        transform=axs[1, 0].transAxes,
    )
    axs[1, 0].yaxis.set_major_formatter(ScalarFormatter(useOffset=False))
    a = runs["OPTIMIZED_REPLAY"][1]
    m = runs["OPTIMIZED_REPLAY"][2]
    for j, color in enumerate(["#b8544b", "#398368", "#555fa5"]):
        axs[1, 1].plot(
            a["time_s"] * 1e3,
            (m["shape_centroid_m"][:, j] - m["shape_centroid_m"][0, j]) * 1e6,
            color=color,
            label="XYZ"[j],
            lw=2,
        )
    axs[1, 1].set(
        title="D  Optimized translation components", xlabel="Time (ms)", ylabel="COM displacement (µm)"
    )
    axs[1, 1].legend()
    axs[1, 1].text(
        0.03,
        0.9,
        "Z translation is not independently demonstrated",
        transform=axs[1, 1].transAxes,
        fontsize=9,
    )
    fig.suptitle("Position-control experiment | 32 fixed IDs, zero background flow, 2 ms", fontsize=15)
    save(
        fig,
        "Fig_F10N_2_target_navigation",
        shared + [out / "runs" / n / "trajectory.h5" for n in ["OFF", "BASELINE", "OPTIMIZED_REPLAY"]],
    )

    fig = plt.figure(figsize=(15, 10), layout="constrained")
    gs = fig.add_gridspec(2, 3, height_ratios=[1.2, 1])
    x0 = runs["OFF"][1]["positions_m"][0, 0]
    for col, (name, title, x) in enumerate(
        [
            ("OFF", "A  Prescribed contact-network initial state", x0),
            (
                "OPTIMIZED_REPLAY",
                "B  Position control: unresolved layout",
                runs["OPTIMIZED_REPLAY"][1]["positions_m"][-1, 0],
            ),
            ("SHAPE_REPLAY", "C  Shape control: final state", runs["SHAPE_REPLAY"][1]["positions_m"][-1, 0]),
        ]
    ):
        ax = fig.add_subplot(gs[0, col], projection="3d")
        spheres(ax, x, p["radius_m"], title)
        ax.text2D(
            0.04,
            0.88,
            "32 actual spheres; R = 1.78 µm\nCentred, equal-scale coordinates",
            transform=ax.transAxes,
            fontsize=9,
        )
    axs = [fig.add_subplot(gs[1, i]) for i in range(3)]
    for name in ["OFF", "BASELINE", "OPTIMIZED_REPLAY", "SHAPE_REPLAY"]:
        a = runs[name][1]
        m = runs[name][2]
        t = a["time_s"] * 1e3
        axs[0].plot(t, (m["shape_rg_m"] - m["shape_rg_m"][0]) * 1e9, label=LABELS[name], color=COLORS[name])
        axs[1].plot(t, m["shape_directional_xy"] - m["shape_directional_xy"][0], color=COLORS[name])
    axs[0].set(title="D  Radius-of-gyration change", xlabel="Time (ms)", ylabel="ΔRg (nm)")
    axs[0].legend(fontsize=8)
    # Observed refinement spread: the COM can agree while individual particles do not.
    nav_variants = [runs["OPTIMIZED_REPLAY"][2]] + [
        load_run(out, n)[2] for n in ["OPTIMIZED_HALF_DT", "OPTIMIZED_FINE_CACHE"]
    ]
    tt = runs["OPTIMIZED_REPLAY"][1]["time_s"] * 1e3
    rg = np.array([(m["shape_rg_m"] - m["shape_rg_m"][0]) * 1e9 for m in nav_variants])
    ratios = np.array([m["shape_directional_xy"] - m["shape_directional_xy"][0] for m in nav_variants])
    axs[0].fill_between(tt, rg.min(0), rg.max(0), color=COLORS["OPTIMIZED_REPLAY"], alpha=0.18)
    axs[1].fill_between(tt, ratios.min(0), ratios.max(0), color=COLORS["OPTIMIZED_REPLAY"], alpha=0.18)
    inset = axs[0].inset_axes([0.48, 0.56, 0.48, 0.32])
    for name in ["OFF", "BASELINE", "SHAPE_REPLAY"]:
        a = runs[name][1]
        m = runs[name][2]
        inset.plot(a["time_s"] * 1e3, (m["shape_rg_m"] - m["shape_rg_m"][0]) * 1e9, color=COLORS[name])
    inset.tick_params(labelsize=7)
    inset.set_title("Weak controls: ΔRg (nm)", fontsize=8)
    axs[1].set(title="E  Directional shape change", xlabel="Time (ms)", ylabel="Δ√(Qxx / Qyy)")
    axs[1].ticklabel_format(axis="y", style="sci", scilimits=(-2, 2), useOffset=False)
    axs[1].text(
        0.02,
        0.04,
        "Blue band: dt/cache sensitivity\nNavigation layout is unresolved",
        transform=axs[1].transAxes,
        fontsize=9,
    )
    axs[1].text(
        0.02,
        0.94,
        f"Predefined target: 1.00 → 1.25\nFinal: {records['SHAPE_REPLAY']['directional_xy'][-1]:.6f}",
        transform=axs[1].transAxes,
        va="top",
        fontsize=9,
    )
    axs[2].axis("off")
    s = records["SHAPE_REPLAY"]
    tensors = np.array(s["tensor_initial_final_m2"]) * 1e12
    text = "F  Centred fixed-ID measurements\n\n"
    text += f"Rg: {s['rg_um'][0]:.6f} → {s['rg_um'][1]:.6f} µm\n"
    text += f"Major/minor ratio: {s['aspect'][0]:.6f} → {s['aspect'][1]:.6f}\n"
    text += f"Connected components: 1 → {s['clusters']}\n"
    text += f"Major axis identifiable: {s['axis_identifiable'][-1]}\n\n"
    text += "Q initial → final diagonal (µm²):\n"
    text += "\n".join(
        f"{label}: {tensors[0, j, j]:.6f} → {tensors[1, j, j]:.6f}"
        for j, label in enumerate(["XX", "YY", "ZZ"])
    )
    text += f"\n\nShape target reached: {str(gate['shape_target_reached']).upper()}\nConnectivity is not physical fusion."
    axs[2].text(0.02, 0.98, text, va="top", fontsize=10, linespacing=1.5)
    fig.suptitle("Morphology control | real particle geometry, translation removed", fontsize=15)
    save(
        fig,
        "Fig_F10N_3_morphology_control",
        shared
        + [out / "runs" / n / "trajectory.h5" for n in [*runs, "OPTIMIZED_HALF_DT", "OPTIMIZED_FINE_CACHE"]],
    )

    fig = plt.figure(figsize=(14, 9), layout="constrained")
    gs = fig.add_gridspec(2, 2)
    ax = fig.add_subplot(gs[0, 0])
    ax.axis("off")
    q = records["OPTIMIZED_REPLAY"]["quality"]
    r = p["radius_m"]
    cells = [
        ["Metric", "Observed", "Admission / scope"],
        [
            "Actual pressure",
            f"{q['pressure_min_Pa']:.1f}–{q['pressure_max_Pa']:.1f} Pa",
            "Table: 0.1–40000 Pa",
        ],
        ["Accepted overlap / R", f"{q['max_accepted_penetration_m'] / r:.5g}", "≤ 0.01"],
        ["Contact pass / R", f"{q['max_contact_pass_correction_m'] / r:.5g}", "≤ 0.1"],
        ["Wall correction / R", f"{q['maximum_wall_correction_m'] / r:.3g}", "≤ 0.1"],
        [
            "COM dt difference",
            f"{cmp['navigation_dt']['max_COM_difference_m'] * 1e9:.4g} nm",
            "Observed, not a proof",
        ],
        [
            "COM cache difference",
            f"{cmp['navigation_cache']['max_COM_difference_m'] * 1e9:.4g} nm",
            "3 → 1.5 µm cache",
        ],
        [
            "Same-ID dt difference",
            f"{cmp['navigation_dt']['max_fixed_id_position_difference_m'] * 1e6:.3f} µm",
            "FAIL: tolerance 0.05 µm",
        ],
        [
            "Same-ID cache difference",
            f"{cmp['navigation_cache']['max_fixed_id_position_difference_m'] * 1e6:.3f} µm",
            "Particle layout unresolved",
        ],
        [
            "Shape Δ(X/Y)",
            f"{gate['shape_directional_xy_change']:.4g}",
            f"Numerical: {gate['shape_directional_xy_numerical_envelope']:.3g}",
        ],
    ]
    table = ax.table(cellText=cells, loc="center", cellLoc="left", colWidths=[0.38, 0.28, 0.34])
    table.auto_set_font_size(False)
    table.set_fontsize(8.5)
    table.scale(1, 1.65)
    for (row, col), cell in table.get_celld().items():
        cell.set_edgecolor("#d6dce1")
        if row == 0:
            cell.set_facecolor("#e7edf1")
            cell.set_text_props(weight="bold")
    ax.set_title("A  Numerical quality and observed sensitivity", loc="left")
    ax = fig.add_subplot(gs[0, 1])
    trials = read(out / "TRIALS.json")
    failed = [x for x in trials if not x["eligible"]]
    ax.axis("off")
    t = "B  Outcomes and retained failures\n\n"
    t += f"Position target reached: {str(gate['target_reached']).upper()}\nError reduction: {gate['error_reduction_um']:.4f} µm\n"
    t += f"Shape target reached: {str(gate['shape_target_reached']).upper()}\nCompleted / attempted: {sum(x['status'] == 'COMPLETE' for x in trials)} / {len(trials)}\n\n"
    for row in failed:
        reason = read(out / "runs" / row["name"] / "SUMMARY.json")["reason"]
        t += (
            row["name"]
            + ": failed\n"
            + ("Contact budget exceeded" if "CONTACT" in str(reason) else str(reason)[:75])
            + "\nValid prefix and proposal retained.\n"
        )
    t += "\nOptimized particle layout: NOT CONVERGED.\nAbove/below-table tests: rejected without clipping."
    ax.text(0.02, 0.98, t, va="top", fontsize=10.5, linespacing=1.4)
    ax = fig.add_subplot(gs[1, 0])
    cpu = read(BASE / "analysis/CPU_CHECKS.json")
    ax.barh(
        ["CPU preflight + caches", "GPU campaign total", "All trajectory integration"],
        [cpu["wall_s"], campaign["wall_s"], campaign["integration_wall_s"]],
        color=["#899fae", "#4b83ac", "#735f99"],
    )
    ax.set(xlabel="Measured wall time (s)", title="C  Bounded computation cost")
    ax.text(
        0.98,
        0.05,
        f"{campaign['trajectory_count']} trajectories; 2 ms each if complete\nSampled device peak: {campaign['sampled_device_peak_MiB']:.0f} MiB",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=9,
    )
    ax = fig.add_subplot(gs[1, 1])
    ax.axis("off")
    ax.text(
        0.02,
        0.98,
        "D  Scientific scope\n\nHOMOGENEOUS_WATER_CONTROL_PROXY\nEXPLORATORY / PHYSICS_UNCALIBRATED\ntraining_ready = false\n\nNo PDMS reflection, free surface, CFD or streaming.\nNo experimental voltage/pressure calibration.\nNo independent signed Z translation demonstrated.\nInstantaneous graph contacts remain time-sensitive.\nDense-contact physics and fusion are not validated.",
        va="top",
        fontsize=11,
        linespacing=1.5,
    )
    fig.suptitle(
        "Control admission | measurable effects, numerical limits and failed candidates", fontsize=15
    )
    save(
        fig,
        "Fig_F10N_4_control_admission",
        shared + [out / "CAMPAIGN.json", out / "TRIALS.json", BASE / "analysis/CPU_CHECKS.json"],
    )
    write(
        folder / "FIGURE_PROVENANCE.json",
        dict(
            figures=provenance,
            source_script_sha256=sha(ROOT / "scripts/straight_tube_exp1/f10n/figures.py"),
            identity=identity(),
        ),
    )


if __name__ == "__main__":
    main()
