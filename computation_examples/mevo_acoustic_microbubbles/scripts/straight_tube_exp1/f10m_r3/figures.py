"""Actual geometry, stored material regions and admission evidence: four PNG/PDFs."""

import numpy as np
import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle, Patch
from matplotlib.colors import ListedColormap
from mevo_mb.coupling.pdms_cavity import apparatus
from mevo_mb.coupling.straight_tube import disc_points
from .common import BASE, ROOT, read, write, sha, identity

OUT = BASE / "figures"
DARK = "#263746"
RED = "#b85c64"
BLUE = "#8fa8bb"
GOLD = "#b59a64"
PROVENANCE = []


def save(fig, name, sources, caption):
    fig.text(0.035, 0.022, caption, fontsize=9, color=DARK, va="bottom")
    for ext in ["png", "pdf"]:
        fig.savefig(OUT / (name + "." + ext), dpi=240, facecolor="white")
    plt.close(fig)
    PROVENANCE.append(
        dict(
            figure=name,
            source_sha256={str(p.relative_to(ROOT)): sha(p) for p in sources},
            file_sha256={ext: sha(OUT / (name + "." + ext)) for ext in ["png", "pdf"]},
            caption=caption,
        )
    )


def annular_rectangle(p, center, inner, limits, plane, opacity):
    """Planar rectangular face with an actual circular hole; display mesh in mm."""
    import pyvista as pv

    lo, hi, low, high = limits
    corners = np.array([[x, y] for x in [lo, hi] for y in [low, high]]) - center
    theta = np.unique(
        np.r_[np.linspace(0, 2 * np.pi, 257), np.mod(np.arctan2(corners[:, 1], corners[:, 0]), 2 * np.pi)]
    )
    ct, st = np.cos(theta), np.sin(theta)
    with np.errstate(divide="ignore"):
        tx = np.where(ct > 0, (hi - center[0]) / ct, (lo - center[0]) / ct)
        ty = np.where(st > 0, (high - center[1]) / st, (low - center[1]) / st)
    rr = np.array([np.full(len(theta), inner), np.minimum(abs(tx), abs(ty))])
    u = center[0] + rr * ct
    v = center[1] + rr * st
    if plane[0] == "z":
        xyz = (u, v, np.full_like(u, plane[1]))
    else:
        xyz = (np.full_like(u, plane[1]), u, v)
    p.add_mesh(pv.StructuredGrid(*xyz), color="#98a7b1", opacity=opacity, specular=0.35, specular_power=50)


def union_mesh(cavity):
    import pyvista as pv

    geom = BASE / "geometry"
    geom.mkdir(exist_ok=True)
    h = 0.05
    # A visualization discretization only; all particle walls are queried from
    # exposed analytic patches and curves, independently of this mesh.
    axes = [np.arange(-16.35, 16.351, h), np.arange(-0.9, 0.901, h), np.arange(-0.9, 4.601, h)]
    grid = pv.ImageData(dimensions=[len(a) for a in axes], spacing=(h, h, h), origin=[a[0] for a in axes])
    q = np.stack(np.meshgrid(*axes, indexing="ij"), -1)
    grid["membership"] = cavity.implicit(q * 1e-3).ravel(order="F") * 1e3
    mesh = grid.contour([0], scalars="membership")
    cc = mesh.cell_centers().points
    # Top plane is geometrically open; keep the X cap surfaces sealed.
    mesh = (
        mesh.extract_cells(cc[:, 2] < 4.5 - 1e-6).extract_surface(algorithm="dataset_surface").triangulate()
    )
    mesh.save(geom / "cavity_open_top_sealed_x.vtp")
    write(
        geom / "MESH_SCOPE.json",
        dict(
            scope="visualization only; no particle wall discretization",
            spacing_m=5e-5,
            top_caps_removed=True,
            X_caps_retained=True,
            input_geometry_sha256=sha(ROOT / "configs/straight_tube_exp1/f10m_r3/geometry.yaml"),
            mesh_sha256=sha(geom / "cavity_open_top_sealed_x.vtp"),
        ),
    )
    return mesh


def scene(cavity, elements, mesh, show_block=True, show_pzt=True):
    import pyvista as pv

    p = pv.Plotter(off_screen=True, window_size=(2400, 1100))
    p.set_background("white")
    p.render_window.SetMultiSamples(0)
    p.render_window.SetAlphaBitPlanes(1)
    p.enable_depth_peeling(number_of_peels=80, occlusion_ratio=0)
    if show_block:
        for sign in [-1, 1]:
            x, z = np.meshgrid([-16.25, 16.25], [-4.5, 4.5])
            p.add_mesh(pv.StructuredGrid(x, np.full_like(x, sign * 4.5), z), color="#99aab5", opacity=0.14)
            annular_rectangle(p, np.array([0.0, 0.0]), 0.8, (-4.5, 4.5, -4.5, 4.5), ("x", sign * 16.25), 0.22)
        x, y = np.meshgrid([-16.25, 16.25], [-4.5, 4.5])
        p.add_mesh(pv.StructuredGrid(x, y, np.full_like(x, -4.5)), color="#99aab5", opacity=0.14)
        for cx, limits in [(-10.85, (-16.25, 0, -4.5, 4.5)), (10.85, (0, 16.25, -4.5, 4.5))]:
            annular_rectangle(p, np.array([cx, 0.0]), 0.75, limits, ("z", 4.5), 0.19)
        p.add_mesh(
            pv.Box(bounds=(-16.25, 16.25, -4.5, 4.5, -4.5, 4.5)).extract_all_edges(),
            color="#7c8c99",
            line_width=1,
        )
    p.add_mesh(
        mesh,
        color=RED,
        opacity=0.63 if show_block else 0.9,
        smooth_shading=True,
        specular=0.3,
        specular_power=35,
    )
    # Rings mark opening rims; no cap or fixed meniscus is drawn.
    for cx in cavity.ports * 1000:
        angle = np.linspace(0, 2 * np.pi, 193)
        ring = np.c_[cx + 0.75 * np.cos(angle), 0.75 * np.sin(angle), np.full_like(angle, 4.5)]
        p.add_mesh(pv.lines_from_points(ring), color=DARK, line_width=1.5)
    if show_pzt:
        for e in elements:
            pos = np.array(e["position_m"]) * 1000
            n = np.array(e["normal"])
            p.add_mesh(
                pv.Disc(center=pos, inner=0, outer=3.5, normal=n, r_res=20, c_res=128),
                color=BLUE,
                opacity=0.55,
                ambient=0.28,
                diffuse=0.6,
                specular=0.85,
                specular_power=70,
                smooth_shading=True,
            )
            p.add_mesh(pv.lines_from_points(disc_points(e) * 1000), color="#597387", line_width=1.6)
            p.add_mesh(pv.Arrow(start=pos, direction=n, scale=2), color=DARK)
            p.add_light(pv.Light(position=pos - n * 12 + [0, 0, 12], focal_point=pos, intensity=0.4))
    else:
        for sign in [-1, 1]:
            p.add_mesh(
                pv.Disc(center=(sign * 16.25, 0, 0), inner=0, outer=0.8, normal=(1, 0, 0), c_res=128),
                color=GOLD,
                opacity=1,
            )
    ax = p.show_bounds(
        bounds=(-18, 18, -6, 6, -5, 6),
        grid="back",
        location="outer",
        xtitle="X (mm)",
        ytitle="Y (mm)",
        ztitle="Z (mm)",
        color=DARK,
        font_size=16,
        n_xlabels=5,
        n_ylabels=3,
        n_zlabels=3,
        all_edges=False,
    )
    ax.GetAxisLabels(1).SetValue(0, "")
    p.camera_position = [(43, -53, 39), (0, 0, 0), (0, 0, 1)]
    p.enable_parallel_projection()
    p.camera.parallel_scale = 14.5
    path = OUT / ("apparatus_render.png" if show_pzt else "cavity_render.png")
    p.screenshot(path)
    p.close()
    return plt.imread(path)


def dim(ax, a, b, text, dy=0.15):
    ax.annotate("", xy=b, xytext=a, arrowprops=dict(arrowstyle="<->", color=DARK, lw=0.8))
    q = (np.array(a) + b) / 2
    ax.text(q[0], q[1] + dy, text, ha="center", va="bottom", fontsize=10)


def figure1(cavity, elements, mesh):
    fig = plt.figure(figsize=(14, 9))
    ax = fig.add_axes([0.02, 0.24, 0.96, 0.64])
    ax.imshow(scene(cavity, elements, mesh))
    ax.axis("off")
    fig.suptitle("F10M-R3 | Four-PZT PDMS apparatus", x=0.04, ha="left", fontsize=17, fontweight="bold")
    fig.text(0.06, 0.2, "Confirmed dimensions", fontsize=12, fontweight="bold")
    fig.text(
        0.06,
        0.18,
        "PDMS: 32.5 × 9 × 9 mm\nMain channel ID: 1.6 mm\nFour ceramic faces: 7 mm",
        fontsize=11,
        linespacing=1.6,
        va="top",
    )
    fig.text(0.39, 0.2, "Mechanical boundaries", fontsize=12, fontweight="bold")
    fig.text(
        0.39,
        0.18,
        "Two top holes: ID 1.5 mm; open\nBoth X ends: always sealed\nSide PDMS path: 3.7 mm",
        fontsize=11,
        linespacing=1.6,
        va="top",
    )
    fig.text(0.72, 0.2, "Provisional hole reference", fontsize=12, fontweight="bold")
    fig.text(
        0.72,
        0.18,
        "End face to hole center: 5.4 mm\nHole axes: X = ±10.85 mm\nCenter-versus-edge reference unconfirmed",
        fontsize=10,
        linespacing=1.6,
        va="top",
    )
    save(
        fig,
        "Fig_F10M_R3_1_apparatus",
        [
            ROOT / "configs/straight_tube_exp1/f10m_r3/geometry.yaml",
            BASE / "geometry/cavity_open_top_sealed_x.vtp",
        ],
        "True dimensional ratios. Red surfaces bound the geometric cavity; top holes have no caps or assumed menisci.\nArrows point inward. PZT faces are numerical source locations; end-seal materials and acoustic coupling remain unknown.",
    )


def figure2(cavity, elements, mesh):
    fig = plt.figure(figsize=(14, 9))
    ax = fig.add_axes([0.01, 0.46, 0.98, 0.44])
    ax.imshow(scene(cavity, elements, mesh, False, False))
    ax.axis("off")
    fig.suptitle(
        "F10M-R3 | Connected cavity and exposed walls", x=0.04, ha="left", fontsize=17, fontweight="bold"
    )
    ax = fig.add_axes([0.075, 0.13, 0.56, 0.3])
    x = np.linspace(-16.5, 16.5, 1101)
    z = np.linspace(-1.0, 4.7, 191)
    xx, zz = np.meshgrid(x, z, indexing="ij")
    mask = cavity.implicit(np.stack([xx, np.zeros_like(xx), zz], -1) * 1e-3) >= 0
    ax.contourf(x, z, mask.T, levels=[0.5, 1.5], colors=["#e6b7bd"])
    for sign in [-1, 1]:
        ax.plot([sign * 16.25] * 2, [-0.8, 0.8], color=GOLD, lw=3)
    for cx in cavity.ports * 1000:
        ax.annotate(
            "Open",
            xy=(cx, 4.5),
            xytext=(cx, 5.1),
            ha="center",
            fontsize=10,
            arrowprops=dict(arrowstyle="->", color=DARK),
        )
    ax.text(0, 2.5, "One connected void\nTwo open top mouths; sealed X faces", ha="center", fontsize=11)
    ax.set(xlim=(-17.5, 17.5), ylim=(-1.5, 5.6), xlabel="X (mm)", ylabel="Z (mm)")
    ax.set_aspect("equal")
    ax.set_title("a  Actual union section: Y = 0", loc="left", fontsize=12)
    ax = fig.add_axes([0.73, 0.17, 0.24, 0.25])
    path = np.array([[cavity.ports[1], 0, z] for z in np.linspace(0.00045, 0.0013, 65)])
    gap = np.array([cavity.query(q, 1.78e-6)["finite_radius_gap_m"] for q in path])
    np.savetxt(
        BASE / "analysis/JUNCTION_CLEARANCE.csv",
        np.c_[path[:, 2], gap],
        delimiter=",",
        header="z_m,finite_sphere_union_clearance_m",
        comments="",
    )
    ax.plot(path[:, 2] * 1000, gap * 1e6, color="#2c779e", lw=2, label="Exposed union wall")
    ax.axvline(0.8, color=RED, ls="--", lw=1, label="Old cylindrical roof")
    ax.set(xlabel="Z along hole axis (mm)", ylabel="Sphere clearance (µm)", ylim=(0, 900))
    ax.ticklabel_format(axis="y", style="plain", useOffset=False)
    ax.set_title("b  Finite-radius junction query", loc="left", fontsize=11)
    ax.grid(alpha=0.15)
    ax.legend(fontsize=8, frameon=False)
    save(
        fig,
        "Fig_F10M_R3_2_cavity_connectivity",
        [
            BASE / "analysis/CONNECTIVITY.json",
            BASE / "analysis/JUNCTION_CLEARANCE.csv",
            ROOT / "src/mevo_mb/coupling/pdms_cavity.py",
        ],
        "Six-neighbor flood fill: one connected cavity component. The wall query includes exposed surfaces and intersection curves.\nClearance is geometric (bubble radius 1.78 µm); interior wetting is a DEMO assumption. Free-surface dynamics are not validated.",
    )


def figure3(cavity, elements):
    fig = plt.figure(figsize=(14, 8.5))
    fig.suptitle(
        "F10M-R3 | Four sources, mirror ports and sealed ends",
        x=0.04,
        ha="left",
        fontsize=17,
        fontweight="bold",
    )
    ax = fig.add_axes([0.075, 0.49, 0.64, 0.38])
    ax.add_patch(Rectangle((-16.25, -4.5), 32.5, 9, fc="#e2e7eb", ec="#6b7f8e"))
    ax.add_patch(Rectangle((-16.25, -0.8), 32.5, 1.6, fc="#edc6ca", ec=RED))
    for cx in cavity.ports * 1000:
        ax.add_patch(Circle((cx, 0), 0.75, fc="white", ec=DARK, lw=1.3))
    for e in elements:
        p = np.array(e["position_m"]) * 1000
        n = np.array(e["normal"])
        if n[0]:
            ax.plot([p[0]] * 2, [-3.5, 3.5], color="#55778d", lw=4)
        else:
            ax.plot([-3.5, 3.5], [p[1]] * 2, color="#55778d", lw=4)
        ax.annotate("", xy=(p[:2] + n[:2] * 2), xytext=p[:2], arrowprops=dict(arrowstyle="->", color=DARK))
    for sign in [-1, 1]:
        dim(ax, (sign * 10.85, -2.1), (sign * 16.25, -2.1), "5.4 mm*", dy=0.2)
        dim(ax, (sign * (10.85 + 0.75), -3.9), (sign * 16.25, -3.9), "4.65 mm", dy=0.05)
    ax.set(xlim=(-18, 18), ylim=(-6, 6), xlabel="X (mm)", ylabel="Y (mm)")
    ax.set_aspect("equal")
    ax.set_title("a  XY projection: four source faces and two top mouths", loc="left", fontsize=12)
    ax = fig.add_axes([0.78, 0.48, 0.2, 0.39])
    ax.add_patch(Rectangle((-4.5, -4.5), 9, 9, fc="#d7dfe5", ec=DARK))
    ax.add_patch(Circle((0, 0), 3.5, fc=BLUE, ec="#55778d", alpha=0.65))
    ax.add_patch(Circle((0, 0), 0.8, fc=GOLD, ec="#7f693e", hatch="///"))
    ax.set(xlim=(-5.5, 5.5), ylim=(-5.5, 5.5), xlabel="Y (mm)", ylabel="Z (mm)")
    ax.set_aspect("equal")
    ax.set_title("b  Both X faces", loc="left", fontsize=12)
    ax = fig.add_axes([0.08, 0.15, 0.56, 0.25])
    ax.axis("off")
    rows = [
        [e["id"], str(np.round(np.array(e["position_m"]) * 1000, 2).tolist()), str(e["normal"])]
        for e in elements
    ]
    table = ax.table(
        cellText=rows,
        colLabels=["PZT face", "Center (mm)", "Inward unit normal"],
        loc="center",
        cellLoc="center",
        colWidths=[0.29, 0.38, 0.33],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 1.7)
    for (r, c), cell in table.get_celld().items():
        cell.set_edgecolor("#bdc7d0")
        cell.set_text_props(color=DARK)
        if r == 0:
            cell.set_facecolor("#edf1f4")
    fig.text(0.71, 0.32, "Both X ends: mechanically sealed", fontsize=11, fontweight="bold")
    fig.text(
        0.71,
        0.23,
        "Gold: seal interface, material unknown\nNo assumed PDMS plug thickness\nNo automatic rigid acoustic reflection",
        fontsize=10,
        linespacing=1.7,
    )
    save(
        fig,
        "Fig_F10M_R3_3_source_port_geometry",
        [BASE / "analysis/GEOMETRY.json"],
        "* Provisional reference: end-PZT face plane to adjacent hole CENTER; center-versus-edge measurement remains unconfirmed.\n4.65 mm is the geometric face-to-hole-volume clearance under that reference. PZT backing and seal-stack thickness are unknown.",
    )


def figure4():
    path = BASE / "materials/regions_and_unresolved_media.h5"
    colors = ["#e8edf3", "#b6c1ca", "#edbdc4", "#d6c9e6", "#d6c9e6", "#bfab7d"]
    with h5py.File(path) as h:
        axes = [h[k + "_m"][:] * 1000 for k in "xyz"]
        i = [np.argmin(abs(a)) for a in axes]
        slices = [h["region_id"][:, :, i[2]], h["region_id"][:, i[1], :], h["region_id"][i[0], :, :]]
        bounds = [h["boundary_id"][:, :, i[2]], h["boundary_id"][:, i[1], :], h["boundary_id"][i[0], :, :]]
    fig = plt.figure(figsize=(14, 9))
    fig.suptitle(
        "F10M-R3 | Material regions and acoustic admission", x=0.04, ha="left", fontsize=17, fontweight="bold"
    )
    rects = [[0.075, 0.57, 0.55, 0.29], [0.075, 0.20, 0.55, 0.29], [0.73, 0.59, 0.22, 0.28]]
    for rect, ij, data, bound, title in zip(
        rects, [(0, 1), (0, 2), (1, 2)], slices, bounds, ["a  XY: Z = 0", "b  XZ: Y = 0", "c  YZ: X ≈ 0"]
    ):
        ax = fig.add_axes(rect)
        marked = np.where(np.isin(bound, [1, 2]), 5, data)
        ax.pcolormesh(
            axes[ij[0]],
            axes[ij[1]],
            marked.T,
            cmap=ListedColormap(colors),
            vmin=-0.5,
            vmax=5.5,
            shading="nearest",
            rasterized=True,
        )
        if ij[0] == 0:
            for sign in [-1, 1]:
                ax.plot([sign * 16.25] * 2, [-0.8, 0.8], color=GOLD, lw=3)
        ax.set_aspect("equal")
        ax.set_xlabel("XYZ"[ij[0]] + " (mm)")
        ax.set_ylabel("XYZ"[ij[1]] + " (mm)")
        ax.set_title(title, loc="left", fontsize=12)
    fig.text(0.70, 0.48, "Four apparatus fields: NOT RUN", fontsize=13, fontweight="bold", color="#985864")
    fig.text(
        0.70,
        0.38,
        "Top fill / exterior medium: UNKNOWN\nEnd seals / acoustic transfer: UNKNOWN\nDrive frequency / pressure calibration: UNKNOWN",
        fontsize=10,
        linespacing=1.8,
    )
    fig.text(
        0.70,
        0.25,
        "Incomplete maps intentionally contain NaN\nSolver input rejection: PASS\nNo four-field reuse from old three-source data",
        fontsize=10,
        linespacing=1.8,
    )
    fig.legend(
        handles=[
            Patch(color=colors[1], label="PDMS reference"),
            Patch(color=colors[2], label="Main liquid DEMO"),
            Patch(color=colors[3], label="Top void: phase unknown"),
            Patch(color=colors[0], label="Exterior unknown"),
            Patch(color=colors[5], label="Seal interface unknown"),
        ],
        loc="lower center",
        bbox_to_anchor=(0.5, 0.09),
        ncol=3,
        frameon=False,
        fontsize=10,
    )
    save(
        fig,
        "Fig_F10M_R3_4_materials_and_admission",
        [path, BASE / "analysis/ACOUSTIC_ADMISSION.json", BASE / "analysis/CPU_CHECKS.json"],
        "Geometry grid: 100 µm, for material bookkeeping only. Port colors label VOID regions, not an assumed water fill.\nMechanical closure does not specify seal material or acoustic reflection. Zero active flow does not model streaming, menisci or residual injection flow.",
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
    cavity, elements = apparatus(c["geometry"])
    mesh = union_mesh(cavity)
    figure1(cavity, elements, mesh)
    figure2(cavity, elements, mesh)
    figure3(cavity, elements)
    figure4()
    write(
        OUT / "FIGURE_PROVENANCE.json",
        dict(identity=identity(), source_script_sha256=sha(__file__), figures=PROVENANCE),
    )


if __name__ == "__main__":
    main()
