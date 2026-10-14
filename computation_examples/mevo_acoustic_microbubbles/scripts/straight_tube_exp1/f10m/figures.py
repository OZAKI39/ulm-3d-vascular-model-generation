"""Four actual English PNG/PDF groups; apparatus acoustics explicitly NOT_RUN."""

import numpy as np
import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle, Polygon
from matplotlib.colors import ListedColormap
from mevo_mb.coupling.pdms_block import apparatus
from mevo_mb.coupling.straight_tube import disc_points
from .common import BASE, ROOT, read, write, sha, identity

OUT = BASE / "figures"
DARK = "#263746"
RED = "#b85c64"
BLUE = "#267aac"
GREY = "#8197a6"
PROVENANCE = []


def save(fig, name, sources, caption):
    fig.text(0.035, 0.025, caption, fontsize=9, color=DARK, va="bottom")
    for ext in ("png", "pdf"):
        fig.savefig(OUT / (name + "." + ext), dpi=240, facecolor="white")
    plt.close(fig)
    PROVENANCE.append(
        dict(
            figure=name,
            source_sha256={str(p.relative_to(ROOT)): sha(p) for p in sources},
            file_sha256={e: sha(OUT / (name + "." + e)) for e in ("png", "pdf")},
            caption=caption,
        )
    )


def dim(ax, a, b, label, offset=(0, 0)):
    ax.annotate("", xy=b, xytext=a, arrowprops=dict(arrowstyle="<->", color=DARK, lw=0.8))
    middle = (np.asarray(a) + b) / 2 + offset
    ax.text(*middle, label, ha="center", va="bottom", fontsize=9, color=DARK)


def scene(g):
    import pyvista as pv

    lumen, elements = apparatus(g)
    p = pv.Plotter(off_screen=True, window_size=(2200, 1100))
    p.set_background("white")
    p.render_window.SetMultiSamples(0)
    p.render_window.SetAlphaBitPlanes(1)
    p.enable_depth_peeling(number_of_peels=80, occlusion_ratio=0)
    # Four outer walls plus perforated end faces: exact block minus the channel.
    for fixed in (1, 2):
        other = 3 - fixed
        for sign in (-1, 1):
            if fixed == 2 and sign == 1:
                continue
            coords = np.zeros((2, 2, 3))
            xx, vv = np.meshgrid([-16.25, 16.25], [-4.5, 4.5])
            coords[..., 0] = xx
            coords[..., other] = vv
            coords[..., fixed] = sign * 4.5
            mesh = pv.StructuredGrid(coords[..., 0], coords[..., 1], coords[..., 2])
            p.add_mesh(mesh, color="#9dacb5", opacity=0.17, specular=0.35, specular_power=50)
    theta = np.linspace(0, 2 * np.pi, 257)
    inner = 0.8 * np.ones_like(theta)
    outer = 4.5 / np.maximum(abs(np.cos(theta)), abs(np.sin(theta)))
    for sign in (-1, 1):
        rr = np.array([inner, outer])
        tt = np.broadcast_to(theta, rr.shape)
        xx = np.full(rr.shape, sign * 16.25)
        mesh = pv.StructuredGrid(xx, rr * np.cos(tt), rr * np.sin(tt))
        p.add_mesh(mesh, color="#91a0ac", opacity=0.24, specular=0.3)
    port = g["top_port"]
    px, py = np.array(port["center_xy_m"]) * 1000
    pr = port["diameter_m"] * 500
    ct, st = np.cos(theta), np.sin(theta)
    with np.errstate(divide="ignore"):
        tx = np.where(ct > 0, (16.25 - px) / ct, (-16.25 - px) / ct)
        ty = np.where(st > 0, (4.5 - py) / st, (-4.5 - py) / st)
    outer = np.minimum(abs(tx), abs(ty))
    rr = np.array([np.full_like(theta, pr), outer])
    p.add_mesh(
        pv.StructuredGrid(px + rr * ct, py + rr * st, np.full_like(rr, 4.5)),
        color="#9dacb5",
        opacity=0.17,
        specular=0.35,
    )
    xx, tt = np.meshgrid(np.unique(np.r_[np.linspace(-16.25, 16.25, 400), px - pr, px + pr]), theta)
    main = pv.StructuredGrid(xx, 0.8 * np.cos(tt), 0.8 * np.sin(tt)).extract_surface(
        algorithm="dataset_surface"
    )
    q = main.points
    main["port_outside"] = np.where(q[:, 2] >= 0, (q[:, 0] - px) ** 2 + (q[:, 1] - py) ** 2 - pr**2, 1.0)
    main = main.clip_scalar(scalars="port_outside", value=0, invert=False)
    lower = np.sqrt(0.8**2 - (py + pr * st) ** 2)
    zz = np.array([lower, np.full_like(theta, 4.5)])
    p.add_mesh(
        pv.StructuredGrid(
            np.broadcast_to(px + pr * ct, zz.shape), np.broadcast_to(py + pr * st, zz.shape), zz
        ),
        color=RED,
        opacity=0.58,
        smooth_shading=True,
        specular=0.25,
    )
    p.add_point_labels(
        [[px, py, 4.7]],
        ["Open port: 1.5 mm"],
        font_size=18,
        text_color=DARK,
        show_points=False,
        shape_opacity=0,
        always_visible=True,
    )
    p.add_mesh(
        main,
        color=RED,
        opacity=0.58,
        smooth_shading=True,
        specular=0.25,
    )
    # Block edges aid scale without altering geometry.
    p.add_mesh(
        pv.Box(bounds=(-16.25, 16.25, -4.5, 4.5, -4.5, 4.5)).extract_all_edges(),
        color="#7c8c99",
        line_width=1.1,
    )
    for e in elements:
        pos = np.array(e["position_m"]) * 1000
        n = np.array(e["normal"])
        d = pv.Disc(center=pos, inner=0, outer=3.5, normal=n, r_res=20, c_res=128)
        p.add_mesh(
            d,
            color="#9baebc",
            opacity=0.55,
            ambient=0.28,
            diffuse=0.6,
            specular=0.85,
            specular_power=70,
            smooth_shading=True,
        )
        p.add_mesh(pv.lines_from_points(disc_points(e) * 1000), color="#5e778c", line_width=1.5)
        p.add_mesh(pv.Arrow(start=pos, direction=n, scale=2.0), color=DARK)
        p.add_light(pv.Light(position=pos - n * 12 + np.array([0, 0, 12]), focal_point=pos, intensity=0.4))
    axes = p.show_bounds(
        bounds=(-17, 18, -6, 6, -6, 6),
        grid="back",
        location="outer",
        xtitle="X (mm)",
        ytitle="Y (mm)",
        ztitle="Z (mm)",
        color=DARK,
        font_size=15,
        n_xlabels=5,
        n_ylabels=3,
        n_zlabels=3,
        all_edges=False,
    )
    # The nearest X/Y corner shares two labels; retain the grid and omit one
    # redundant Y-end label so the independent coordinate values remain legible.
    axes.GetAxisLabels(1).SetValue(0, "")
    p.camera_position = [(43, -51, 36), (0, 0, 0), (0, 0, 1)]
    p.enable_parallel_projection()
    p.camera.parallel_scale = 14.5
    p.screenshot(OUT / "apparatus_render.png")
    p.close()
    return plt.imread(OUT / "apparatus_render.png")


def figure1(c):
    g = c["geometry"]
    lumen, elements = apparatus(g)
    fig = plt.figure(figsize=(13, 10))
    gs = fig.add_gridspec(
        2, 1, height_ratios=[2.2, 1.4], left=0.045, right=0.98, bottom=0.13, top=0.92, hspace=0.04
    )
    ax = fig.add_subplot(gs[0])
    ax.imshow(scene(g))
    ax.axis("off")
    ax.set_title("a  Centered PDMS channel and three 7 mm ceramic faces", loc="left", fontweight="bold")
    ax = fig.add_subplot(gs[1])
    ax.add_patch(Rectangle((-16.25, -4.5), 32.5, 9, facecolor="#d9e1e5", edgecolor=GREY))
    ax.add_patch(Rectangle((-16.25, -0.8), 32.5, 1.6, facecolor="#f5dadc", edgecolor=RED))
    for e in elements:
        q = disc_points(e) * 1000
        ax.add_patch(Polygon(q[:, :2], facecolor=GREY, edgecolor="#5e778c", alpha=0.7, lw=2))
    dim(ax, (-16.25, -5.8), (16.25, -5.8), "32.5 mm", (0, -1.0))
    dim(ax, (-6, 0.8), (-6, 4.5), "3.7 mm", (-2.0, 0))
    ax.annotate(
        "ID 1.6 mm", xy=(-11, 0), xytext=(-14, 6.4), arrowprops=dict(arrowstyle="->", color=DARK), fontsize=10
    )
    ax.annotate(
        "X+ end always closed\nTop port remains open",
        xy=(16.25, 0),
        xytext=(11, 6.4),
        arrowprops=dict(arrowstyle="->", color=DARK),
        fontsize=9,
        ha="center",
    )
    ax.text(-1.5, 7.8, "Side faces contact PDMS", ha="center", fontsize=10)
    ax.add_patch(Circle((10.85, 0), 0.75, fill=False, ec="#8c5e32", ls="--"))
    dim(ax, (10.85, -3), (16.25, -3), "5.4 mm (provisional)", (0, -2))
    ax.set(xlim=(-20, 20), ylim=(-7.3, 9), xlabel="X (mm)", ylabel="Y (mm)")
    ax.set_aspect("equal")
    ax.set_title("b  XY orthographic projection", loc="left", fontsize=11)
    fig.suptitle("F10M | Confirmed apparatus geometry", x=0.04, ha="left", fontweight="bold", fontsize=16)
    save(
        fig,
        "Fig_F10M_1_PDMS_apparatus",
        [BASE / "CONFIG_SNAPSHOT.json", BASE / "analysis/GEOMETRY.json"],
        "PDMS: 32.5 × 9 × 9 mm. Centered through-channel: ID 1.6 mm. No independent thin tube wall.\nTrue dimensional ratios; translucent display only. Normals point inward. Top port: ID 1.5 mm, open; X+ end: closed. 5.4 mm reference and lateral centering are provisional.",
    )


def material_panel(ax, axes, labels, ij, elements, grid):
    i, j = ij
    ax.pcolormesh(
        axes[i] * 1000,
        axes[j] * 1000,
        labels.T,
        cmap=ListedColormap(["#e8f3f8", "#b6bfc8", "#f0bcc2"]),
        vmin=0,
        vmax=2,
        shading="nearest",
        rasterized=True,
    )
    half = np.array([16.25, 4.5, 4.5])
    ax.add_patch(Rectangle((-half[i], -half[j]), 2 * half[i], 2 * half[j], fill=False, ec=DARK, lw=1))
    if i == 0:
        for y in [-0.8, 0.8]:
            ax.plot([-16.25, 16.25], [y, y], color=RED, lw=1)
    else:
        ax.add_patch(Circle((0, 0), 0.8, fill=False, ec=RED, lw=1))
    for e in elements:
        q = disc_points(e) * 1000
        ax.plot(q[:, i], q[:, j], color="#435e72", lw=1.6)
    h = grid["spacing_m"][0] * 1000
    end = np.array(grid["expanded_shape"]) * h / 2
    a0 = axes[i][0] * 1000
    a1 = axes[i][-1] * 1000
    b0 = axes[j][0] * 1000
    b1 = axes[j][-1] * 1000
    ax.add_patch(
        Rectangle((-end[i], -end[j]), 2 * end[i], 2 * end[j], fill=False, ec="#8b6ba2", ls="--", lw=1)
    )
    ax.add_patch(Rectangle((a0, b0), a1 - a0, b1 - b0, fill=False, ec="#8b6ba2", lw=0.8))
    ax.set(
        xlim=(-end[i] - 1, end[i] + 1),
        ylim=(-end[j] - 1, end[j] + 1),
        xlabel="XYZ"[i] + " (mm)",
        ylabel="XYZ"[j] + " (mm)",
    )
    ax.set_aspect("equal")


def figure2(c):
    grid = read(BASE / "analysis/RESOURCE_PREFLIGHT.json")
    _, elements = apparatus(c["geometry"])
    with h5py.File(BASE / "materials/material_maps.h5") as h:
        axes = [h[v + "_m"][:] for v in "xyz"]
        index = [np.argmin(abs(v)) for v in axes]
        slices = [
            h["material_id"][:, :, index[2]],
            h["material_id"][:, index[1], :],
            h["material_id"][index[0], :, :],
        ]
    fig = plt.figure(figsize=(14, 8.4))
    gs = fig.add_gridspec(
        2, 2, width_ratios=[2.6, 1], left=0.065, right=0.98, bottom=0.21, top=0.9, wspace=0.12, hspace=0.4
    )
    for pos, ij, label, data in zip(
        [(0, 0), (1, 0), (0, 1)],
        [(0, 1), (0, 2), (1, 2)],
        ["a  XY: z = 0", "b  XZ: y = 0", "c  YZ: x = 0"],
        slices,
    ):
        ax = fig.add_subplot(gs[pos])
        material_panel(ax, axes, data, ij, elements, grid)
        ax.set_title(label, loc="left", fontsize=12)
    ax = fig.add_subplot(gs[1, 1])
    ax.add_patch(Rectangle((-4.5, -4.5), 9, 9, fc="#b6bfc8", ec=DARK))
    ax.add_patch(Circle((0, 0), 3.5, fc="#879eaf", ec="#435e72", alpha=0.7))
    ax.add_patch(Circle((0, 0), 0.8, fc=RED, ec=RED, hatch="////"))
    ax.annotate(
        "100% opening\ncovered",
        xy=(0, 0),
        xytext=(0, -6.8),
        ha="center",
        fontsize=9,
        arrowprops=dict(arrowstyle="->", color=DARK),
    )
    ax.set(xlim=(-6, 6), ylim=(-8, 6), xlabel="Y (mm)", ylabel="Z (mm)")
    ax.set_aspect("equal")
    ax.set_title("d  X+ end-face interference", loc="left", fontsize=12)
    fig.suptitle(
        "F10M | Material regions and port interference", x=0.04, ha="left", fontsize=16, fontweight="bold"
    )
    from matplotlib.patches import Patch

    fig.legend(
        handles=[
            Patch(color="#b6bfc8", label="PDMS reference"),
            Patch(color="#f0bcc2", label="Channel liquid: water DEMO"),
            Patch(color="#e8f3f8", label="Exterior: water DEMO"),
            Patch(fill=False, ec="#8b6ba2", label="PML boundaries (planned)"),
        ],
        loc="lower center",
        bbox_to_anchor=(0.5, 0.082),
        ncol=4,
        frameon=False,
        fontsize=10,
    )
    save(
        fig,
        "Fig_F10M_2_material_regions",
        [
            BASE / "materials/material_maps.h5",
            BASE / "analysis/GEOMETRY.json",
            BASE / "analysis/RESOURCE_PREFLIGHT.json",
        ],
        "Slices show actual stored material labels. PZT outlines are projections, not volumetric solid-acoustic models.\nExternal medium is unconfirmed. Centered end PZT overlaps the complete 1.6 mm opening; X+ closure and open top port are user-confirmed.",
    )


def figure3():
    a = read(BASE / "adapter_validation/VALIDATION.json")
    r = read(BASE / "analysis/RESOURCE_PREFLIGHT.json")
    fig, axs = plt.subplots(2, 2, figsize=(13, 8))
    fig.subplots_adjust(left=0.15, right=0.97, top=0.88, bottom=0.15, hspace=0.44, wspace=0.3)
    ax = axs[0, 0]
    ax.barh(
        ["RTX 4090", "Apparatus estimate"],
        [24564 / 1024, r["estimated_peak_GiB"]],
        color=[GREY, "#b76b76"],
    )
    ax.axvline(20, ls="--", color=DARK, lw=1)
    ax.set_xlabel("GPU memory (GiB)")
    ax.set_title("a  Full apparatus field: NOT_RUN", loc="left", fontweight="bold")
    for i, v in enumerate([24564 / 1024, r["estimated_peak_GiB"]]):
        ax.text(v + 0.5, i, f"{v:.1f}", va="center")
    ax.set_xlim(0, 75)
    ax.text(
        0.03,
        0.05,
        "3 MHz; 700 × 270 × 270 incl. PML\nMinimum 6.19 points / wavelength",
        transform=ax.transAxes,
        fontsize=9,
    )
    ax = axs[0, 1]
    ax.axis("off")
    ax.set_title("b  Admission and model scope", loc="left", fontweight="bold")
    ax.text(
        0,
        0.95,
        "Three apparatus bases: NOT_RUN\n\n• Memory / time estimates exceed budget\n• End-solid acoustics uncalibrated\n• No 3 MHz APECSS response\n\nActual GPU adapter tests: 5 MHz\nLongitudinal equivalent fluid; no shear\nSource pressure ≠ voltage calibration",
        va="top",
        fontsize=11,
        linespacing=1.5,
        color=DARK,
    )
    ax = axs[1, 0]
    for label, color, title in [
        ("scalar_water", BLUE, "Water"),
        ("lossless_interface", GREY, "Interface, lossless"),
        ("absorbing_interface", RED, "Interface, absorbing"),
    ]:
        q = np.load(BASE / "adapter_validation" / (label + ".npz"))
        t = q["time_s"] * 1e6
        ax.plot(t, q["pressure_Pa"][2], color=color, lw=1.2, label=title)
    ax.set(xlabel="Time (µs)", ylabel="Pressure at test sensor (Pa)")
    ax.set_title("c  Actual 3D planar validation output", loc="left", fontsize=12)
    ax.legend(fontsize=8, frameon=False)
    ax.grid(alpha=0.16)
    ax = axs[1, 1]
    labels = ["|Reflection|", "Amplitude ratio"]
    xx = np.arange(2)
    expected = [abs(a["reflection_expected_pressure_coefficient"]), a["attenuation_ratio_expected"]]
    measured = [a["reflection_measured_magnitude"], a["attenuation_ratio_measured"]]
    ax.bar(xx - 0.17, expected, 0.34, color=GREY, label="Analytic reference")
    ax.bar(xx + 0.17, measured, 0.34, color=BLUE, label="Actual GPU")
    ax.set(xticks=xx, xticklabels=labels, ylabel="Dimensionless")
    ax.set_ylim(0, 1.12)
    ax.legend(frameon=False, fontsize=9)
    ax.set_title("d  Numerical interface / loss checks", loc="left", fontsize=12)
    ax.text(
        0.04,
        0.5,
        f"PDMS phase error: {a['pdms_phase_error_rad']:.4f} rad\nConstant map / scalar difference: {a['scalar_map_relative_error']:.1e}",
        transform=ax.transAxes,
        fontsize=9,
    )
    fig.suptitle(
        "F10M | Acoustic admission and measured numerical checks",
        x=0.04,
        ha="left",
        fontsize=16,
        fontweight="bold",
    )
    save(
        fig,
        "Fig_F10M_3_acoustic_admission",
        [
            BASE / "adapter_validation/VALIDATION.json",
            BASE / "analysis/RESOURCE_PREFLIGHT.json",
            BASE / "adapter_validation/absorbing_interface.npz",
            BASE / "adapter_validation/lossless_interface.npz",
            BASE / "adapter_validation/scalar_water.npz",
        ],
        "The 32.5 mm apparatus field was not computed. Panels c–d are bounded numerical tests, not apparatus pressure predictions.\nPDMS: literature reference, 3–7 MHz coverage. Top-port placement reference, actual media and source calibration remain unresolved.",
    )


def cloud(ax, x, title):
    u = np.linspace(0, 2 * np.pi, 16)
    v = np.linspace(0, np.pi, 10)
    r = 1.78
    sphere = (
        np.array(
            [
                np.outer(np.cos(u), np.sin(v)),
                np.outer(np.sin(u), np.sin(v)),
                np.outer(np.ones_like(u), np.cos(v)),
            ]
        )
        * r
    )
    for i, pt in enumerate(x * 1e6):
        q = sphere + pt[:, None, None]
        ax.plot_surface(*q, color=BLUE, linewidth=0, shade=True)
        if i in [0, 7, 19]:
            ax.text(*(pt + 3), str(i), fontsize=8, color=DARK)
    ax.set(xlim=(5, 95), ylim=(5, 95), zlim=(5, 95), xlabel="X (µm)", ylabel="Y (µm)", zlabel="Z (µm)")
    ax.set_box_aspect([1, 1, 1])
    ax.view_init(24, -58)
    ax.set_title(title, loc="left", fontsize=12)


def figure4():
    with h5py.File(BASE / "dynamics/OFF_trajectory.h5") as h:
        x = h["positions_m"][:, 0]
        t = h["time_s"][:] * 1000
        wall = h["cylinder_stats"][:]
    m = np.load(BASE / "dynamics/metrics.npz")
    s = read(BASE / "dynamics/SUMMARY.json")
    fig = plt.figure(figsize=(13, 9))
    gs = fig.add_gridspec(
        3, 3, height_ratios=[2.3, 1, 1], left=0.09, right=0.97, top=0.9, bottom=0.13, hspace=0.58, wspace=0.4
    )
    cloud(fig.add_subplot(gs[0, 0], projection="3d"), x[0], "a  Initial state")
    cloud(fig.add_subplot(gs[0, 1], projection="3d"), x[-1], "b  OFF: 1.000 ms")
    ax = fig.add_subplot(gs[0, 2])
    ax.axis("off")
    ax.text(
        0,
        0.92,
        "32 fixed particle IDs\nRadius: 1.78 µm (DEMO)\n\nZero background flow\nX+ closed; top port open\nMaximum displacement: 0\n\nThree-PZT ON: NOT_RUN\nApparatus fields unavailable\nResponse frequency unmatched",
        va="top",
        fontsize=11,
        linespacing=1.5,
        color=DARK,
    )
    eig = m["global_eigenvalues_m2"][:, 0]
    aspect = np.sqrt(eig[:, 0] / eig[:, -1])
    metrics = [
        (
            "Centroid displacement (nm)",
            np.linalg.norm(m["global_centroid_m"][:, 0] - m["global_centroid_m"][0, 0], axis=1) * 1e9,
        ),
        ("Radius of gyration (µm)", m["global_rg_m"][:, 0] * 1e6),
        ("Aspect ratio", aspect),
        ("Connected components", m["cluster_count"][:, 0]),
        ("Accepted overlap (nm)", np.full_like(t, s["quality"]["max_accepted_penetration_m"] * 1e9)),
        ("Wall correction (nm)", wall[:, :, 2].max(axis=1) * 1e9),
    ]
    for index, (title, y) in enumerate(metrics):
        ax = fig.add_subplot(gs[1 + index // 3, index % 3])
        ax.plot(t, y, color=BLUE, lw=1.7)
        ax.set(xlabel="Physical time (ms)", ylabel=title)
        ax.grid(alpha=0.17)
        if np.all(y == 0):
            ax.set_ylim(-0.05, 0.2)
    fig.suptitle(
        "F10M | No-flow baseline in the confirmed channel", x=0.04, ha="left", fontsize=16, fontweight="bold"
    )
    save(
        fig,
        "Fig_F10M_4_no_flow_bubbles",
        [BASE / "dynamics/OFF_trajectory.h5", BASE / "dynamics/metrics.npz", BASE / "dynamics/SUMMARY.json"],
        "Actual 1 ms GPU trajectory, dt = 1 µs. True sphere radii and equal coordinate scales; fixed IDs 0, 7 and 19 labeled.\nOFF positions remain bitwise unchanged. Original contact budgets retained. Main-bore check only: junction dynamics not admitted. ON: NOT_RUN.",
    )


def main():
    OUT.mkdir(exist_ok=True)
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "text.color": DARK,
            "axes.labelcolor": DARK,
            "axes.edgecolor": DARK,
            "xtick.color": DARK,
            "ytick.color": DARK,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    c = read(BASE / "CONFIG_SNAPSHOT.json")
    figure1(c)
    figure2(c)
    figure3()
    figure4()
    write(
        OUT / "FIGURE_PROVENANCE.json",
        dict(identity=identity(), source_script_sha256=sha(__file__), figures=PROVENANCE),
    )


if __name__ == "__main__":
    main()
