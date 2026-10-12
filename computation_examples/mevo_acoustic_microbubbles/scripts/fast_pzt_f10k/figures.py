"""Four publication figures from original geometry and measured arrays only."""

import re
import numpy as np
import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.text import Text
from matplotlib.ticker import NullFormatter
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from .common import ROOT, BASE, INPUT, identity, js, sha, write

BLUE = "#267aac"
RED = "#a95c60"
TEAL = "#347f77"
GOLD = "#bd9136"
GREY = "#717784"
DARK = "#25303b"


def style3d(ax, limit=None):
    ax.view_init(elev=24, azim=-65)
    for axis in [ax.xaxis, ax.yaxis, ax.zaxis]:
        axis.pane.fill = False
        axis._axinfo["grid"]["color"] = (0.7, 0.72, 0.75, 0.35)
    ax.tick_params(labelsize=7, pad=0)
    if limit:
        ax.set(xlim=(-limit, limit), ylim=(-limit, limit), zlim=(-limit, limit))
        ax.set_box_aspect((1, 1, 1))


def spheres(ax, x, r=1.78, highlight=None):
    # Actual radius and positions, no display magnification or cluster replacement.
    u, v = np.meshgrid(np.linspace(0, 2 * np.pi, 19), np.linspace(0, np.pi, 13))
    unit = np.array([np.cos(u) * np.sin(v), np.sin(u) * np.sin(v), np.cos(v)])
    for i, p in enumerate(x):
        a = p[:, None, None] + r * unit
        ax.plot_surface(
            *a,
            color=GOLD if i == highlight else BLUE,
            linewidth=0,
            antialiased=True,
            shade=True,
            rasterized=True,
        )


def mesh(ax, vertices, faces, color=RED, alpha=0.18):
    artist = Poly3DCollection(
        vertices[faces], facecolors=color, edgecolors="none", alpha=alpha, rasterized=True
    )
    ax.add_collection3d(artist)


def main():
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 11,
            "axes.labelsize": 9,
            "legend.fontsize": 8,
            "pdf.fonttype": 42,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "text.color": DARK,
            "axes.labelcolor": DARK,
            "xtick.color": DARK,
            "ytick.color": DARK,
            "savefig.facecolor": "white",
        }
    )
    out = BASE / "figures"
    out.mkdir(exist_ok=True)
    provenance = []
    geometry = np.load(INPUT / "f4/prepared_geometry.npz")
    vertices = geometry["vertices_m"]
    faces = geometry["faces"][geometry["labels"] == 1]
    admission = js(BASE / "analysis/F5_ADMISSION.json")
    state = np.load(ROOT / admission["state_file"])
    nav = js(BASE / "analysis/NAVIGATION.json")
    ns = np.load(BASE / "analysis/NAVIGATION_SERIES.npz")
    topo = js(BASE / "analysis/B2_TOPOLOGY_COMPARISON.json")
    ts = np.load(BASE / "analysis/TOPOLOGY_SERIES.npz")
    cost = js(BASE / "analysis/SCALING.json")

    def save(fig, name, sources, caption):
        fig.text(0.07, 0.024, caption, fontsize=8, va="bottom", color=DARK)
        for text_artist in fig.findobj(Text):
            text = text_artist.get_text()
            text = re.sub(
                r"\b(at|for|the|original|real|throughout|up to|dt|Actual|BLAS|CPU|GPU|radius|old)(?=\d)",
                r"\1 ",
                text,
                flags=re.I,
            )
            text = re.sub(r"([:,])(?=\S)", r"\1 ", text)
            text = text.replace(" +warmup", " + warmup")
            text_artist.set_text(text)
        for ext in ["png", "pdf"]:
            fig.savefig(out / (name + "." + ext), dpi=230, facecolor="white")
        plt.close(fig)
        provenance.append(
            dict(
                figure=name,
                caption=caption,
                source_sha256={str(p.relative_to(ROOT)): sha(p) for p in sources},
                files={ext: sha(out / (name + "." + ext)) for ext in ["png", "pdf"]},
            )
        )

    fig = plt.figure(figsize=(12, 8.8))
    gs = fig.add_gridspec(2, 2, left=0.08, right=0.96, bottom=0.13, top=0.90, hspace=0.45, wspace=0.25)
    fig.suptitle(
        "Historical split candidate: geometry is admissible, near-wall B2 is not", fontsize=15, y=0.975
    )
    ax = fig.add_subplot(gs[0, 0], projection="3d")
    mesh(ax, vertices * 1e3, faces)
    x = state["positions_m"]
    center = x.mean(0)
    ax.scatter(*center * 1e3, color=BLUE, s=35, depthshade=False)
    ax.text(*((center + np.array([0.0007, 0, 0.0002])) * 1e3), "F5 at417.5 ms", fontsize=8)
    lo = vertices.min(0) * 1e3
    hi = vertices.max(0) * 1e3
    ax.set(
        xlim=(lo[0], hi[0]),
        ylim=(lo[1], hi[1]),
        zlim=(lo[2], hi[2]),
        xlabel="X (mm)",
        ylabel="Y (mm)",
        zlabel="Z (mm)",
        title="a  Original open Y-vessel wall",
    )
    ax.set_box_aspect(hi - lo, zoom=1.2)
    style3d(ax)
    ax.view_init(elev=50, azim=-65)
    ax.set_zticks([-1, 0, 1])
    ax = fig.add_subplot(gs[0, 1], projection="3d")
    import pyvista as pv

    wall = pv.PolyData(vertices, np.c_[np.full(len(faces), 3), faces].ravel())
    bounds = np.c_[center - 18e-6, center + 18e-6].ravel()
    clipped = wall.clip_box(bounds, invert=False).extract_surface(algorithm="dataset_surface").triangulate()
    mesh(ax, (np.asarray(clipped.points) - center) * 1e6, clipped.faces.reshape(-1, 4)[:, 1:], alpha=0.35)
    spheres(ax, (x - center) * 1e6, highlight=3)
    v = state["CFD_velocity_m_s"].mean(0)
    direction = v / np.linalg.norm(v)
    ax.quiver(-12, -10, 8, *direction, length=7, color=DARK, arrow_length_ratio=0.2)
    ax.text(-12, -10, 12, "CFD direction", fontsize=8)
    style3d(ax, 18)
    ax.set(
        xlabel=r"X - $X_c$ ($\mu$m)",
        ylabel=r"Y - $Y_c$ ($\mu$m)",
        zlabel=r"Z - $Z_c$ ($\mu$m)",
        title="b  Actual32 bubbles, true radius1.78 µm",
    )
    ax = fig.add_subplot(gs[1, 0])
    gap = state["wall_gap_m"] * 1e6
    ax.bar(np.arange(32), gap, color=[GOLD if i == 3 else BLUE for i in range(32)], width=0.7)
    ax.axhline(0, color=DARK, lw=0.8)
    ax.set_ylim(0, 14)
    ax.set(
        xlabel="Fixed particle ID",
        ylabel="Bubble-surface / wall gap (µm)",
        title="c  Independent original-triangle distance",
    )
    ax.text(
        0.03,
        0.92,
        "Minimum gap ≈0\nNo resolved penetration; near-wall mobility uncalibrated",
        transform=ax.transAxes,
        va="top",
        fontsize=8,
    )
    ax = fig.add_subplot(gs[1, 1])
    with h5py.File(ROOT / admission["source"]) as h:
        t = h["time_s"][:]
        counts = h["metrics/cluster_count"][:, 0]
        mask = (t >= 0.416) & (t <= 0.422)
    ax.step(
        t[mask] * 1e3, counts[mask], where="post", color=BLUE, label="Stored M0 components (0.5 ms cadence)"
    )
    ax.axvline(418.54, color=GOLD, ls="--", label="Accepted-edge persistent split:418.54 ms")
    ax.set(
        xlabel="Historical absolute time (ms)",
        ylabel="Connected components",
        yticks=[1, 2],
        ylim=(0.7, 2.5),
        title="d  M0 history is not B2 ground truth",
    )
    ax.legend(loc="upper left")
    ax.text(
        0.02,
        0.30,
        "B2_TOPOLOGY_CASE_NOT_ADMITTED\nReason: wall contact outside interior-only B2 scope.",
        transform=ax.transAxes,
        color=RED,
        fontsize=9,
    )
    save(
        fig,
        "Fig_F10K_1_topology_candidate",
        [
            BASE / "analysis/F5_ADMISSION.json",
            ROOT / admission["state_file"],
            ROOT / admission["source"],
            INPUT / "f4/prepared_geometry.npz",
            INPUT / "result_121.vtu",
        ],
        "Original mesh and actual sphere positions; no translated cloud. Gold:ID3, later separated in M0.\nGeometric non-penetration does not validate wall hydrodynamics. EXPLORATORY / PHYSICS_UNCALIBRATED.",
    )

    fig = plt.figure(figsize=(12.8, 8.5))
    gs = fig.add_gridspec(
        2, 3, width_ratios=[1.05, 1, 1], left=0.12, right=0.97, bottom=0.16, top=0.89, wspace=0.3, hspace=0.4
    )
    fig.suptitle(
        "Scientific scope: one interior persistent merge, incomplete topology validation",
        fontsize=15,
        y=0.975,
    )
    ax = fig.add_subplot(gs[:, 0])
    labels = [
        "Fixed-state numerics",
        "Group motion",
        "Internal morphology",
        "Strict edge timing",
        "Persistent merge",
        "Persistent split",
        "Near-wall hydro.",
        "Physical calibration",
    ]
    values = np.array([[2, 2], [2, 2], [1, 1], [0, 0], [1, 2], [1, 3], [3, 3], [3, 3]])
    ax.imshow(values, aspect="auto", cmap=ListedColormap([RED, GOLD, TEAL, "#d9dce0"]), vmin=0, vmax=3)
    ax.set(
        xticks=[0, 1],
        xticklabels=["M0", "B2"],
        yticks=np.arange(8),
        yticklabels=labels,
        title="a  Evidence, not physical truth",
    )
    texts = [
        ["Checked", "Checked"],
        ["Limited\nscenarios", "Limited\nscenarios"],
        ["dt sensitive", "Edge limited"],
        ["FAIL", "FAIL"],
        ["Historical\nonly", "1 case\nchecked"],
        ["Historical\nonly", "Untested"],
        ["Uncalibrated", "Uncalibrated"],
        ["Uncalibrated", "Uncalibrated"],
    ]
    for i in range(8):
        for j in range(2):
            ax.text(
                j,
                i,
                texts[i][j],
                ha="center",
                va="center",
                fontsize=8,
                color="white" if values[i, j] in [0, 2] else DARK,
            )
    ax.tick_params(length=0)
    ax.set_ylim(7.5, -0.5)
    for col, k, title in [(1, 0, "b  B2 initial:31 +1"), (2, -1, "c  B2 at3 ms:32 connected")]:
        ax = fig.add_subplot(gs[0, col], projection="3d")
        x = ts["DT10_x"][k]
        spheres(ax, (x - x.mean(0)) * 1e6, highlight=5)
        style3d(ax, 11)
        ax.set(xlabel="X (µm)", ylabel="Y (µm)", zlabel="Z (µm)", title=title)
    ax = fig.add_subplot(gs[1, 1:])
    for tag, color, style in [("DT10", BLUE, "-"), ("DT5", GOLD, "--")]:
        labels_ = ts[tag + "_labels"]
        counts = np.array([len(np.unique(v)) for v in labels_])
        t = ts[tag + "_times"] * 1e3
        ax.step(t, counts, where="post", color=color, ls=style, lw=2, label=tag.replace("DT", "dt=") + " µs")
    ax.axvspan(0.29, 3, color=TEAL, alpha=0.08)
    ax.set(
        xlabel="Time from counterfactual initialization (ms)",
        ylabel="Connected components",
        yticks=[1, 2],
        ylim=(0.6, 2.7),
        title="d  Same merge at0.290 ms; persists for2.710 ms",
    )
    ax.legend(loc="upper right")
    m = topo["comparison"]["maxima"]
    ax.text(
        0.40,
        0.52,
        f"Max fixed-ID difference:{m['max_fixed_ID_m'] * 1e9:.2f} nm\nKabsch RMSD:{m['Kabsch_RMSD_m'] * 1e9:.2f} nm\nRaw edges:46 /47; matched timing up to60 µs",
        transform=ax.transAxes,
        fontsize=9,
    )
    save(
        fig,
        "Fig_F10K_2_scientific_scope",
        [
            BASE / "analysis/SCIENTIFIC_SCOPE.json",
            BASE / "analysis/B2_TOPOLOGY_COMPARISON.json",
            BASE / "analysis/TOPOLOGY_SERIES.npz",
        ],
        "Same B2, stiffness and controls; dt10/5 µs, original17.8/178 nm quality budgets. True-radius spheres in centroid-relative coordinates.\nA graph merge is not physical bubble fusion. Persistent splitting remains untested; old12/29 edge-timing FAIL retained; training_ready=false.",
    )

    fig = plt.figure(figsize=(12, 8.5))
    gs = fig.add_gridspec(2, 2, left=0.10, right=0.96, bottom=0.14, top=0.90, hspace=0.45, wspace=0.30)
    fig.suptitle(
        "Interior-lumen navigation: measured control time and limited transverse authority",
        fontsize=14,
        y=0.975,
    )
    c0 = np.array(nav["rows"][0]["initial_centroid_m"])
    ax = fig.add_subplot(gs[0, 0], projection="3d")
    mesh(ax, vertices * 1e3, faces, alpha=0.15)
    ax.scatter(*c0 * 1e3, color=BLUE, s=35, depthshade=False)
    ax.text(*(c0 * 1e3 + np.array([0.6, 0, 0.35])), "Interior initial group", fontsize=8)
    ax.set(
        xlim=(lo[0], hi[0]),
        ylim=(lo[1], hi[1]),
        zlim=(lo[2], hi[2]),
        xlabel="X (mm)",
        ylabel="Y (mm)",
        zlabel="Z (mm)",
        title="a  Same original Y-vessel and CFD",
    )
    ax.set_box_aspect(hi - lo, zoom=1.2)
    style3d(ax)
    ax.view_init(elev=50, azim=-65)
    ax.set_zticks([-1, 0, 1])
    inset = fig.add_axes([0.07, 0.54, 0.12, 0.15], projection="3d")
    inset.patch.set_alpha(0)
    xinit = np.asarray(nav["rows"][0]["initial_positions_m"])
    spheres(inset, (xinit - xinit.mean(0)) * 1e6)
    style3d(inset, 10)
    inset.set_axis_off()
    inset.text2D(0.0, 0.0, "Initial 32; R=1.78 µm", transform=inset.transAxes, fontsize=7)
    ax = fig.add_subplot(gs[0, 1])
    colors = {"OFF": GREY, "S_1000": RED, "O_1000": BLUE}
    for mode in colors:
        c = (ns["10_" + mode + "_c"] - c0) * 1e6
        ax.plot(c[:, 0], c[:, 2], color=colors[mode], label=mode.replace("_1000", ""))
    ax.axvline(20, color=DARK, lw=0.9, ls=":")
    ax.scatter(20, 1, marker="*", s=110, color=GOLD, label="Prospective fixed target")
    ax.annotate(
        "Original F10A +1 µm diagnostic target",
        xy=(0, 1),
        xytext=(22, 0.72),
        fontsize=7.5,
        arrowprops=dict(arrowstyle="->", color=GREY),
    )
    ax.scatter(0, 1, facecolors="none", edgecolors=GOLD, s=65)
    ax.set(
        xlim=(-4, 105),
        ylim=(-0.16, 1.15),
        xlabel="Centroid downstream displacement (µm)",
        ylabel="Centroid Z displacement (µm)",
        title="b  Actual saved centroid paths; unequal axis scales",
    )
    ax.legend(loc="center right", fontsize=7.5)
    ax = fig.add_subplot(gs[1, 0])
    for step, style in [("10", "-"), ("5", "--")]:
        t = ns[step + "_OFF_t"]
        off = ns[step + "_OFF_c"]
        for mode, color in [("O_1000", BLUE), ("S_1000", RED)]:
            delta = (ns[step + "_" + mode + "_c"] - off)[:, 2] * 1e9
            ax.plot(t * 1e3, delta, style, color=color, label=f"{mode[0]}−OFF,dt{step} µs")
    ax.axvline(2.02433, color=GREY, ls=":", label="20 µm plane:≈2.024 ms")
    ax.axhline(50, color=GOLD, ls="--", label="Original50 nm effect criterion")
    ax.set(
        xlabel="Elapsed physical time (ms)",
        ylabel="Additional target-direction motion (nm)",
        title="c  Paired effect from complete saved trajectories",
    )
    ax.legend(fontsize=7, loc="upper left", ncol=2)
    ax = fig.add_subplot(gs[1, 1])
    ax.axis("off")
    p = nav["paired"]
    txt = (
        f"At the20 µm downstream plane\n\nO−OFF: {p[0]['O_minus_OFF_plane_Z_m'] * 1e9:.2f} / {p[1]['O_minus_OFF_plane_Z_m'] * 1e9:.2f} nm (10/5 µs)\n"
        f"Paired step-size difference: {nav['dt_pair_plane_increment_difference_m'] * 1e9:.3f} nm\n"
        f"Error to fixed +1 µm target: {p[0]['O_plane_target_YZ_error_m'] * 1e6:.3f} µm\n\n"
        "O points toward +Z, but misses the50 nm effect criterion.\n"
        "Field/CFD support observed throughout10 ms;\nactual field-exit time is right-censored.\n\n"
        "Navigation to the proposed target: NOT ESTABLISHED\nNear-wall F8 target is a separate historical case."
    )
    ax.text(0, 1, "d  Control time depends on the target plane", fontweight="bold", va="top", fontsize=11)
    ax.text(0, 0.84, txt, va="top", fontsize=9, linespacing=1.4)
    save(
        fig,
        "Fig_F10K_3_navigation_feasibility",
        [
            BASE / "analysis/NAVIGATION.json",
            BASE / "analysis/NAVIGATION_SERIES.npz",
            BASE / "analysis/NAVIGATION_PREREGISTRATION.json",
            INPUT / "f4/prepared_geometry.npz",
        ],
        "1 kPa is the initial local pressure median. Curves are actual F10A trajectories; plane crossings interpolate real50 µs brackets.\nThe downstream target is a prospective protocol, not a completed historical navigation task. EXPLORATORY / training_ready=false.",
    )

    fig, aa = plt.subplots(2, 2, figsize=(12, 8.6))
    fig.subplots_adjust(left=0.09, right=0.97, bottom=0.15, top=0.88, hspace=0.45, wspace=0.3)
    fig.suptitle(
        "B2 fixed-state scaling: measured costs, numerical gates and explicit limits", fontsize=15, y=0.975
    )
    rows = cost["rows"]
    good = [r for r in rows if r["status"] == "PASS"]
    ax = aa[0, 0]
    for mode, color in [("cpu", GREY), ("gpu", BLUE)]:
        if good:
            ax.loglog(
                [r["N"] for r in good],
                [r["formal_median_s"][mode] * 1e3 for r in good],
                "o-",
                color=color,
                label=mode.upper(),
            )
    ax.set(
        xlabel="Bubbles in one connected problem,N (B=1)",
        ylabel="Complete solve cost (ms)",
        title="a  Only admitted matched timings",
    )
    ax.legend()
    ax.set_xticks([32, 128, 512], labels=["32", "128", "512"])
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax = aa[0, 1]
    if good:
        ax.bar([str(r["N"]) for r in good], [r["CPU_over_GPU_speedup"] for r in good], color=TEAL)
    ax.axhline(1, color=DARK, ls="--", lw=0.8)
    ax.set(
        xlabel="N (passed CPU and GPU)",
        ylabel="CPU / GPU cost",
        title="b  Measured speed ratio, same quality targets",
    )
    ax = aa[1, 0]
    executed = [r for r in rows if "input_contact_pairs" in r]
    if executed:
        ax.plot([r["N"] for r in executed], [r["input_contact_pairs"] for r in executed], "o-", color=GOLD)
    ax.set(
        xlabel="N actually initialized",
        ylabel="Initial active contact pairs",
        title="c  Dense single connected networks",
    )
    twin = ax.twinx()
    if executed:
        twin.plot(
            [r["N"] for r in executed],
            [
                r.get("peak_RSS_bytes", r.get("controller_sampled_peak_RSS_bytes", 0)) / 1024**2
                for r in executed
            ],
            "s--",
            color=GREY,
        )
    twin.set_ylabel("Measured process peak RSS (MiB)", color=GREY)
    twin.tick_params(axis="y", colors=GREY)
    ax = aa[1, 1]
    ax.axis("off")
    lines = []
    for r in rows:
        if r["status"] == "PASS":
            lines.append(
                f"N={r['N']}:PASS, CPU{r['formal_median_s']['cpu'] * 1e3:.2f} / GPU{r['formal_median_s']['gpu'] * 1e3:.2f} ms"
            )
        else:
            lines.append(f"N={r['N']}:{r['status']}\n  {r.get('reason', 'numerical gate failed')}")
    lines += [
        "",
        f"Solve +warmup:{cost['solve_and_warmup_s']:.3f} s /120 s",
        f"Bounded stage:{cost['stage_wall_s']:.3f} s /300 s",
        "",
        "N>32:BENCHMARK_ONLY_SYNTHETIC_STATE",
        "Regular lattice; not a realistic dense bubble trajectory.",
        "No B64 integration or AI-data admission.",
    ]
    ax.text(0, 1, "d  Cost boundary and unexecuted sizes", fontweight="bold", va="top", fontsize=11)
    ax.text(0, 0.86, "\n".join(lines), va="top", fontsize=9, linespacing=1.4)
    save(
        fig,
        "Fig_F10K_4_b2_scaling",
        [BASE / "analysis/SCALING.json", ROOT / "configs/fast_pzt/f10k/protocol.json"],
        "N=32 is real; 128/512 are regular synthetic networks with different iteration counts. BLAS 1 thread. Full timing includes preparation, residual checks and matched output.\nOne warmup + two formal solves; construction separate. Missing sizes are not measurements. Fixed-state numerical scaling is not physical validation.",
    )
    write(
        out / "FIGURE_PROVENANCE.json",
        dict(
            figures=provenance,
            identity=identity(),
            source_sha256=sha(__file__),
            no_smoothing=True,
            true_sphere_radius_m=1.78e-6,
            wall_caps_excluded=True,
            sphere_rendering="Equal-axis3D original centers and radii; centroid subtraction only in labelled local views. Original wall clipped only for local view.",
            rasterized_3D_surfaces_in_PDF=True,
            training_ready=False,
        ),
    )


if __name__ == "__main__":
    main()
