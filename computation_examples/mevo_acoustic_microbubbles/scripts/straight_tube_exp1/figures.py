"""English figures from resolved geometry, actual complex fields and trajectories."""

import numpy as np
import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle, Polygon
from .common import BASE, ROOT, read, write, sha, identity
from mevo_mb.coupling.straight_tube import apparatus, disc_points

RED = "#b85c64"
BLUE = "#267aac"
GREY = "#8197a6"
DARK = "#263746"
OUT = BASE / "figures"
PROVENANCE = []


def save(fig, name, sources, caption):
    fig.text(0.055, 0.024, caption, fontsize=9, color=DARK, va="bottom")
    for ext in ("png", "pdf"):
        fig.savefig(OUT / (name + "." + ext), dpi=260, facecolor="white")
    plt.close(fig)
    PROVENANCE.append(
        dict(
            figure=name,
            source_sha256={str(p.relative_to(ROOT)): sha(p) for p in sources},
            file_sha256={ext: sha(OUT / (name + "." + ext)) for ext in ("png", "pdf")},
            caption=caption,
        )
    )


def dim(ax, start, end, label, offset=(0, 0)):
    ax.annotate("", xy=end, xytext=start, arrowprops=dict(arrowstyle="<->", color=DARK, lw=0.8))
    middle = (np.array(start) + end) / 2 + offset
    ax.text(*middle, label, ha="center", va="bottom", fontsize=8, color=DARK)


def geometry3d(tube, elements, initial):
    import pyvista as pv

    p = pv.Plotter(off_screen=True, window_size=(1900, 1250))
    p.set_background("white")
    p.render_window.SetMultiSamples(0)
    p.render_window.SetAlphaBitPlanes(1)
    p.enable_depth_peeling(number_of_peels=80, occlusion_ratio=0)
    for rad, opacity in [(tube.outer_radius_m, 0.56), (tube.radius_m, 0.18)]:
        s = tube.surface(rad) * 1e3
        mesh = pv.StructuredGrid(s[..., 0], s[..., 1], s[..., 2]).extract_surface(algorithm="dataset_surface")
        p.add_mesh(mesh, color=RED, opacity=opacity, smooth_shading=True, specular=0.28, specular_power=40)
    # Annular end walls leave the lumen visibly open.
    for sign in (-1, 1):
        ring = pv.Disc(
            center=tube.world([sign * tube.length_m / 2, 0, 0]) * 1e3,
            inner=tube.radius_m * 1e3,
            outer=tube.outer_radius_m * 1e3,
            normal=tube.rotation[:, 0],
            r_res=2,
            c_res=128,
        )
        p.add_mesh(ring, color="#8a424c", opacity=0.75, specular=0.2)
    for i, e in enumerate(elements):
        center = np.array(e["position_m"]) * 1e3
        disk = pv.Disc(
            center=center, inner=0, outer=e["size_m"][0] * 500, normal=e["normal"], r_res=24, c_res=128
        )
        p.add_mesh(
            disk,
            color="#a2b0b8",
            opacity=0.5 if i != 0 else 0.33,
            ambient=0.3,
            diffuse=0.55,
            specular=0.88,
            specular_power=68,
            smooth_shading=True,
        )
        ring = disc_points(e) * 1e3
        p.add_mesh(pv.lines_from_points(ring), color="#6e8291", line_width=1.2)
        p.add_mesh(pv.Arrow(start=center, direction=e["normal"], scale=2.0), color=DARK)
        light = pv.Light(
            position=center - np.array(e["normal"]) * 10 + np.array([0, 0, 10]),
            focal_point=center,
            intensity=0.4,
        )
        p.add_light(light)
    for point in initial:
        p.add_mesh(
            pv.Sphere(radius=0.00178, center=point * 1e3, theta_resolution=12, phi_resolution=12), color=BLUE
        )
    p.show_bounds(
        bounds=(-11, 14, -5, 5, -5, 5),
        grid="back",
        location="outer",
        xtitle="X (mm)",
        ytitle="Y (mm)",
        ztitle="Z (mm)",
        color=DARK,
        font_size=14,
        n_xlabels=6,
        n_ylabels=3,
        n_zlabels=3,
        all_edges=False,
    )
    p.camera_position = [(35, -35, 27), (2, 0, 0), (0, 0, 1)]
    p.enable_parallel_projection()
    p.camera.parallel_scale = 12.4
    p.screenshot(OUT / "apparatus_3d_render.png")
    p.close()
    return plt.imread(OUT / "apparatus_3d_render.png")


def ortho(ax, tube, elements, plane, domain=False, pre=None):
    i, j = plane
    # Exact default local-coordinate projections; reject arbitrary transform here.
    if not np.allclose(tube.rotation, np.eye(3)) or np.count_nonzero(tube.center):
        raise ValueError(
            "Figure cross sections require a new registered plane choice after apparatus rotation"
        )
    L = tube.length_m * 1e3
    R = tube.radius_m * 1e3
    Ro = tube.outer_radius_m * 1e3
    if i == 0:
        ax.add_patch(Rectangle((-L / 2, -Ro), L, 2 * Ro, color=RED, alpha=0.23))
        ax.add_patch(Rectangle((-L / 2, -R), L, 2 * R, facecolor="white", edgecolor=RED, lw=0.8))
    else:
        ax.add_patch(Circle((0, 0), Ro, facecolor=RED, alpha=0.3))
        ax.add_patch(Circle((0, 0), R, facecolor="white", edgecolor=RED, lw=1))
    for e in elements:
        q = disc_points(e) * 1e3
        # Projection, not a false intersection: explicitly labelled in panels.
        ax.add_patch(Polygon(q[:, [i, j]], facecolor=GREY, edgecolor=GREY, alpha=0.28, lw=1.4))
    if domain:
        bounds = np.array(pre["grid_bounds_m"]) * 1e3
        for pad, color, style, label in [
            (0, DARK, "--", "Physical grid / PML interface"),
            (1.6, "#939393", ":", "Outer PML boundary"),
        ]:
            lo = bounds[:, 0] - pad
            hi = bounds[:, 1] + pad
            ax.add_patch(
                Rectangle(
                    (lo[i], lo[j]),
                    hi[i] - lo[i],
                    hi[j] - lo[j],
                    fill=False,
                    edgecolor=color,
                    linestyle=style,
                    lw=0.8,
                    label=label,
                )
            )
        roi = np.array(read(BASE / "RESOLVED_ACOUSTIC_CONFIG.json")["acoustics"]["roi_half_width_m"]) * 1e3
        ax.add_patch(
            Rectangle(
                (-roi[i], -roi[j]),
                2 * roi[i],
                2 * roi[j],
                fill=False,
                edgecolor="#71508a",
                lw=1,
                label="Stored acoustic ROI",
            )
        )
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("XYZ"[i] + " (mm)")
    ax.set_ylabel("XYZ"[j] + " (mm)")
    ax.grid(alpha=0.18)


def apparatus_figures():
    c = read(BASE / "CONFIG_SNAPSHOT.json")
    tube, elements = apparatus(c["demo"]["geometry"], c["acoustic"]["acoustics"]["frequency_Hz"])
    initial = np.load(BASE / "initial_state.npz")["positions_m"]
    pre = read(BASE / "analysis/PREFLIGHT.json")
    fig = plt.figure(figsize=(13, 7.8))
    fig.suptitle("Three-PZT no-flow straight-tube apparatus", fontsize=17, y=0.965)
    ax = fig.add_axes([0.01, 0.18, 0.64, 0.69])
    ax.imshow(geometry3d(tube, elements, initial))
    ax.axis("off")
    ax.set_title("a  Three-dimensional view • true dimensional ratios", loc="left", fontsize=11)
    ax = fig.add_axes([0.68, 0.50, 0.29, 0.33])
    ortho(ax, tube, elements, (0, 1))
    ax.set(xlim=(-11, 14.5), ylim=(-6, 6))
    ax.set_title("b  XY orthographic projection", fontsize=11)
    ax.text(0, 4.8, "Side Y+", ha="center", fontsize=9)
    ax.text(0, -5.1, "Side Y−", ha="center", fontsize=9)
    ax.text(13, 4.8, "End X+", ha="center", fontsize=9)
    dim(ax, (10, -2), (13, -2), "3 mm", (0, -1.2))
    dim(ax, (-2, -0.9), (-2, -3.9), "3 mm", (-2.5, -0.5))
    ax.annotate(
        "Initial 32 bubbles",
        xy=initial.mean(0)[[0, 1]] * 1e3,
        xytext=(-9, 2),
        fontsize=8,
        arrowprops=dict(arrowstyle="-", color=BLUE),
        color=BLUE,
    )
    ax = fig.add_axes([0.68, 0.19, 0.29, 0.25])
    ax.axis("off")
    ax.text(0, 0.96, "Apparatus register", fontweight="bold", fontsize=11)
    ax.text(
        0,
        0.82,
        "Inner diameter: 1.6 mm (PROVISIONAL)\nLength: 20 mm (DEMO)\nWall: 0.1 mm; OD: 1.8 mm (DEMO)\nThree active faces: Ø9 mm (PROVISIONAL)\nTwo side faces + one end face; no Z faces\nBackground flow: exactly zero\nNormals shown by dark arrows",
        va="top",
        fontsize=10,
        linespacing=1.5,
    )
    save(
        fig,
        "Fig_F10L_1_apparatus_geometry",
        [BASE / "CONFIG_SNAPSHOT.json", BASE / "initial_state.npz", BASE / "analysis/PREFLIGHT.json"],
        "DEMO geometry; open lumen ends. PZT faces have zero acoustic thickness; material reflections are not included.\nInitial spheres retain their actual radius (1.78 µm); their subpixel global size is not enlarged.",
    )
    fig = plt.figure(figsize=(13, 9))
    fig.suptitle("Geometric dimensions, acoustic domain and finite bubble radius", fontsize=16, y=0.97)
    panels = [
        ([0.06, 0.59, 0.55, 0.29], (0, 1), "a  XY tube section; source projections"),
        ([0.06, 0.20, 0.55, 0.29], (0, 2), "b  XZ tube section; source projections"),
        ([0.69, 0.57, 0.25, 0.31], (1, 2), "c  YZ tube section; source projections"),
    ]
    for box, plane, title in panels:
        ax = fig.add_axes(box)
        ortho(ax, tube, elements, plane, True, pre)
        ax.set_title(title, loc="left", fontsize=11)
        ax.set(xlim=(-19, 19) if plane[0] == 0 else (-10, 10), ylim=(-10, 10))
        if plane == (0, 1):
            dim(ax, (-10, -2), (10, -2), "20 mm (DEMO)", (0, -1.4))
            dim(ax, (10, 2.0), (13, 2.0), "3 mm", (0, 0.6))
        if plane == (1, 2):
            ax.annotate(
                "ID 1.6 mm\nOD 1.8 mm (DEMO)",
                xy=(0.8, 0),
                xytext=(1.8, -7.5),
                fontsize=8,
                arrowprops=dict(arrowstyle="-", color=DARK),
            )
            ax.text(-8, 8, "PZT Ø9 mm", fontsize=8)
    ax = fig.add_axes([0.69, 0.20, 0.25, 0.29])
    # True-scale radial-wall inset (coordinates are µm relative to inner wall).
    r = 1.78
    ax.add_patch(Rectangle((0, -15), 100, 30, color=RED, alpha=0.3))
    ax.add_patch(Circle((-r, 0), r, facecolor=BLUE, edgecolor="white", lw=0.4))
    ax.axvline(0, color=RED, lw=1)
    ax.axvline(100, color=RED, lw=1)
    dim(ax, (0, 7), (100, 7), "100 µm wall (DEMO)", (0, 1))
    ax.annotate(
        "R = 1.78 µm", xy=(-r, 0), xytext=(18, -8), fontsize=8, arrowprops=dict(arrowstyle="-", color=DARK)
    )
    ax.set(
        xlim=(-12, 112),
        ylim=(-15, 15),
        xlabel="Radial distance from lumen wall (µm)",
        ylabel="Tangential distance (µm)",
    )
    ax.set_aspect("equal")
    ax.set_title("d  Finite-radius wall constraint", loc="left", fontsize=11)
    save(
        fig,
        "Fig_F10L_2_geometry_cross_sections",
        [
            BASE / "CONFIG_SNAPSHOT.json",
            BASE / "analysis/PREFLIGHT.json",
            BASE / "RESOLVED_ACOUSTIC_CONFIG.json",
        ],
        "Equal spatial scales in every panel. Shaded grey: projected PZT faces; purple: stored ROI.\nDashed: physical grid / outside-PML interface; dotted: outer PML boundary. Wall inset is a geometric constraint illustration.",
    )


def field_figure():
    source = BASE / "acoustic/basis.h5"
    with h5py.File(source) as h:
        x, y, z = [h["axes/" + a + "_m"][:] * 1e3 for a in "xyz"]
        p = h["P_complex"][:].astype(complex)
        g = np.stack([h["gradP_complex_" + a][:] for a in "xyz"], axis=-1)
    c = read(BASE / "CONFIG_SNAPSHOT.json")
    amplitudes = np.array(c["acoustic"]["demo_control"]["amplitude_Pa"])
    phases = np.array(c["acoustic"]["demo_control"]["phase_rad"])
    controls = amplitudes * np.exp(1j * phases)
    p = np.concatenate(
        [p * controls[:, None, None, None], (controls @ p.reshape(3, -1)).reshape((1,) + p.shape[1:])]
    )
    g = np.concatenate(
        [g * controls[:, None, None, None, None], np.einsum("k,kxyzj->xyzj", controls, g)[None]], axis=0
    )
    zidx = int(np.argmin(abs(z)))
    names = ["Side Y+ only", "Side Y− only", "End X+ only", "Three PZT, in phase"]
    fig, axes = plt.subplots(4, 2, figsize=(13, 7.8))
    fig.subplots_adjust(left=0.07, right=0.88, bottom=0.16, top=0.87, hspace=0.64, wspace=0.26)
    fig.suptitle("New three-source k-Wave fields in the straight tube", fontsize=16, y=0.97)
    fig.text(
        0.07,
        0.918,
        "Amplitude and phase at Z = 0 • source input 300 Pa peak per active channel • 2.5 MHz DEMO",
        fontsize=10,
    )
    vmax = max(np.abs(p[:, :, :, zidx]).max() / 1e3, 0.01)
    for row, name in enumerate(names):
        data = p[row, :, :, zidx]
        amp = np.abs(data) / 1e3
        phase = np.angle(data)
        phase = np.where(amp >= vmax * 0.005, phase, np.nan)
        for col, image in enumerate([amp, phase]):
            ax = axes[row, col]
            artist = ax.pcolormesh(
                x,
                y,
                image.T,
                shading="nearest",
                cmap="magma" if col == 0 else "twilight",
                vmin=0 if col == 0 else -np.pi,
                vmax=vmax if col == 0 else np.pi,
                rasterized=True,
            )
            ax.set_aspect("equal")
            ax.set_ylim(-1.0, 1.0)
            ax.set_xlim(-10.2, 10.2)
            ax.axhline(0.8, color="white", lw=0.6, alpha=0.8)
            ax.axhline(-0.8, color="white", lw=0.6, alpha=0.8)
            ax.set_title(
                name + (" | pressure amplitude" if col == 0 else " | phase"), loc="left", fontsize=10
            )
            ax.set_ylabel("Y (mm)")
            ax.set_xlabel("X (mm)" if row == 3 else "")
            if col == 0:
                # Physical direction of ∇|P|, derived from the stored F6 derivative.
                grad = np.real(np.conj(data)[..., None] * g[row, :, :, zidx, :]) / np.maximum(
                    np.abs(data)[..., None], 1e-30
                )
                sx = np.arange(7, len(x) - 7, 12)
                sy = np.array([len(y) // 2 - 4, len(y) // 2, len(y) // 2 + 4])
                xx, yy = np.meshgrid(x[sx], y[sy], indexing="ij")
                vector = grad[np.ix_(sx, sy)][..., :2]
                norm = np.linalg.norm(vector, axis=-1)
                uv = vector / np.maximum(norm[..., None], 1e-30)
                ax.quiver(
                    xx,
                    yy,
                    uv[..., 0],
                    uv[..., 1],
                    color="white",
                    pivot="mid",
                    angles="xy",
                    scale_units="xy",
                    scale=2.8,
                    width=0.0028,
                )
                amp_artist = artist
            else:
                phase_artist = artist
    cax = fig.add_axes([0.91, 0.53, 0.014, 0.30])
    fig.colorbar(amp_artist, cax=cax, label=r"$|\,P\,|$ peak (kPa)")
    cax = fig.add_axes([0.91, 0.19, 0.014, 0.25])
    cb = fig.colorbar(phase_artist, cax=cax, ticks=[-np.pi, 0, np.pi], label="Phase (rad)")
    cb.ax.set_yticklabels(["−π", "0", "π"])
    save(
        fig,
        "Fig_F10L_3_acoustic_basis",
        [source, BASE / "acoustic/ACOUSTIC_GATE.json", BASE / "CONFIG_SNAPSHOT.json"],
        "Actual k-Wave outputs and their registered complex sum. Arrows: normalized ∇|P|, not force magnitude. Equal XY scales.\nWhite lines: lumen walls. Phase omitted below 0.5% of the shared amplitude maximum. Homogeneous-water demo; no experimental calibration.",
    )


def spheres(ax, x, r):
    u, v = np.meshgrid(np.linspace(0, 2 * np.pi, 17), np.linspace(0, np.pi, 11))
    unit = np.array([np.cos(u) * np.sin(v), np.sin(u) * np.sin(v), np.cos(v)])
    for point in x:
        xyz = point[:, None, None] + r * unit
        ax.plot_surface(*xyz, color=BLUE, lw=0, shade=True, rasterized=True)
    ax.set(xlim=(5, 95), ylim=(5, 95), zlim=(5, 95), xlabel="X (µm)", ylabel="Y (µm)", zlabel="Z (µm)")
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=23, azim=-58)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.fill = False
        axis._axinfo["grid"]["color"] = (0.7, 0.74, 0.78, 0.3)
    ax.tick_params(labelsize=7, pad=0)
    for index in [0, 7, 19]:
        ax.text(*(x[index] + [1, 1, 1]), str(index), fontsize=6, color=DARK)


def dynamics_figure():
    data = {}
    sources = []
    for name in ("OFF", "THREE_PZT_ON"):
        folder = BASE / "dynamics" / name
        data[name] = (
            np.load(folder / "trajectory.npz"),
            np.load(folder / "metrics.npz"),
            np.load(folder / "frame_observables.npz"),
        )
        sources.extend(folder / v for v in ("trajectory.npz", "metrics.npz", "frame_observables.npz"))
    gate = read(BASE / "dynamics/DYNAMICS_GATE.json")
    fig = plt.figure(figsize=(13, 10))
    fig.suptitle("No-flow full-particle dynamics: OFF versus three PZT ON", fontsize=16, y=0.974)
    layout = fig.add_gridspec(3, 3, left=0.07, right=0.96, bottom=0.13, top=0.91, hspace=0.55, wspace=0.36)
    off, on = data["OFF"], data["THREE_PZT_ON"]
    panels = [
        (off[0]["positions_m"][0, 0], "a  Initial 32 bubbles"),
        (off[0]["positions_m"][-1, 0], "b  OFF at 1 ms"),
        (on[0]["positions_m"][-1, 0], "c  Three PZT ON at 1 ms"),
    ]
    for j, (xx, title) in enumerate(panels):
        ax = fig.add_subplot(layout[0, j], projection="3d")
        spheres(ax, xx * 1e6, 1.78)
        ax.set_title(title, loc="left", fontsize=11, pad=4)
    axes = [fig.add_subplot(layout[i, j]) for i in (1, 2) for j in range(3)]
    for name, color, ls in [("OFF", GREY, "--"), ("THREE_PZT_ON", BLUE, "-")]:
        tr, m, obs = data[name]
        t = tr["time_s"] * 1e3
        centroid = m["global_centroid_m"][:, 0]
        eig = m["global_eigenvalues_m2"][:, 0]
        axes[0].plot(
            t,
            np.linalg.norm(centroid - centroid[0], axis=-1) * 1e9,
            color=color,
            ls=ls,
            label=name.replace("_", " "),
        )
        axes[1].plot(t, m["global_rg_m"][:, 0] * 1e6, color=color, ls=ls)
        axes[2].plot(t, np.sqrt(eig[:, 0] / eig[:, -1]), color=color, ls=ls)
        axes[3].step(t, m["cluster_count"][:, 0], where="post", color=color, ls=ls)
        quality = tr["quality_samples"][:, 0]
        axes[4].plot(t, np.max(quality[:, :, 3], axis=-1) * 1e9, color=color, ls=ls)
        wall = tr["cylinder_stats"]
        axes[5].plot(t, np.max(wall[:, :, 2], axis=-1) * 1e9, color=color, ls=ls)
    titles = [
        "d  Centroid displacement",
        "e  Radius of gyration",
        "f  Global aspect ratio",
        "g  Connected components",
        "h  Accepted overlap",
        "i  Wall correction",
    ]
    labels = [
        "Displacement (nm)",
        "Rg (µm)",
        "√(λmax / λmin)",
        "Component count",
        "Maximum overlap (nm)",
        "Maximum correction (nm)",
    ]
    for ax, title, label in zip(axes, titles, labels):
        ax.set(title=title, xlabel="Physical time (ms)", ylabel=label)
        ax.grid(alpha=0.2)
        ax.ticklabel_format(axis="y", useOffset=False, style="plain")
    axes[0].legend(fontsize=8)
    axes[4].axhline(17.8, color=RED, lw=0.8, ls=":", label="17.8 nm budget")
    axes[4].legend(fontsize=7)
    axes[5].axhline(178, color=RED, lw=0.8, ls=":", label="178 nm budget")
    axes[5].legend(fontsize=7)
    text = "Fixed IDs; actual R = 1.78 µm and equal XYZ scales. Raw saved frames; no position or time smoothing.\n"
    pressure = gate["cases"]["THREE_PZT_ON"]["pressure_peak_range_Pa"]
    text += f"ON local |P| peak: {pressure[0]:.1f}–{pressure[1]:.1f} Pa. Both trajectories: 1 ms, 1 µs steps. This short test does not establish aggregation or target control."
    save(fig, "Fig_F10L_4_no_flow_bubble_dynamics", sources + [BASE / "dynamics/DYNAMICS_GATE.json"], text)


def main():
    import argparse

    p = argparse.ArgumentParser()
    p.add_argument("--geometry-only", action="store_true")
    args = p.parse_args()
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 11,
            "pdf.fonttype": 42,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "text.color": DARK,
            "axes.labelcolor": DARK,
            "xtick.color": DARK,
            "ytick.color": DARK,
        }
    )
    OUT.mkdir(exist_ok=True)
    apparatus_figures()
    if not args.geometry_only:
        field_figure()
        dynamics_figure()
    write(
        OUT / "FIGURE_PROVENANCE.json",
        dict(
            identity=identity(),
            figures=PROVENANCE,
            geometry_style="PZT v6 muted red / grey-blue / specular faces; new true-scale analytic apparatus",
            no_smoothing=True,
            source_script_sha256=sha(__file__),
        ),
    )


if __name__ == "__main__":
    main()
