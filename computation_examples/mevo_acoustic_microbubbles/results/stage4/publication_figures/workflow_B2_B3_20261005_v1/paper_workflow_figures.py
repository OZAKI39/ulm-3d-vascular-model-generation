"""Publication figures from sealed Stage 4I B2/B3; no simulation or trajectory edits."""

import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager, colors
from matplotlib.collections import LineCollection
from matplotlib.patches import Circle, FancyArrowPatch
from matplotlib.cm import ScalarMappable
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection, Text3D
import numpy as np
from PIL import Image
import pyvista as pv
import yaml

from mevo_mb.acoustics.basis import ComplexAcousticField
from mevo_mb.bubbles.adaptive_table import AdaptiveResponseTableBubbleModel
from mevo_mb.bubbles.forces import ClassicalFarFieldBjerknes, primary_bjerknes

ROOT = Path(__file__).resolve().parents[2]
GOLD = ROOT / "results/stage4/performance_refactor/gold_baseline/GOLD_BASELINE_MANIFEST.json"
BLUE, ORANGE, INK, GREY = "#0072B2", "#D55E00", "#20262D", "#73808D"
CMAP = "magma"
FONT_PATH = Path("/mnt/c/Windows/Fonts/arial.ttf")


def sha(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def style():
    if FONT_PATH.is_file():
        for font_file in FONT_PATH.parent.glob("arial*.ttf"):
            font_manager.fontManager.addfont(str(font_file))
    font = "Arial" if FONT_PATH.is_file() else "DejaVu Sans"
    plt.rcParams.update(
        {
            "font.family": font,
            "font.size": 8,
            "axes.labelsize": 8,
            "axes.titlesize": 9,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 8,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "text.color": INK,
            "axes.labelcolor": INK,
            "axes.edgecolor": INK,
            "xtick.color": INK,
            "ytick.color": INK,
            "axes.linewidth": 0.65,
            "lines.linewidth": 0.8,
            "xtick.major.width": 0.65,
            "ytick.major.width": 0.65,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "mathtext.fontset": "stixsans",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "legend.frameon": False,
        }
    )
    return font


def mesh_lines(mesh, normal="z", origin=(0, 0, 0)):
    part = mesh.slice(normal=normal, origin=origin)
    cells, offset, result = part.lines, 0, []
    while offset < len(cells):
        count = cells[offset]
        points = part.points[cells[offset + 1 : offset + 1 + count]] * 1000
        result.extend(np.stack([points[:-1], points[1:]], axis=1))
        offset += count + 1
    return np.asarray(result)


def load_data():
    gold = json.loads(GOLD.read_text())
    assert gold["status"] == "SEALED"
    protected = {str(GOLD.relative_to(ROOT)): sha(GOLD)}
    basis_path = ROOT / "results/stage4/acoustic/20261002T143722411057Z/acoustic/basis.h5"
    protected[str(basis_path.relative_to(ROOT))] = sha(basis_path)
    with h5py.File(basis_path) as f:
        basis = f["P_complex"][:]
    data = {}
    for name in ["B2", "B3"]:
        run = ROOT / gold["cases"][name]["result_id"]
        for rel in [
            "particles/trajectories.h5",
            "particles/SUMMARY.json",
            "config.yaml",
            "BRANCH_METRICS.json",
            "PHYSICS_BANNER.json",
            "logs/source_snapshot.zip",
        ]:
            p = run / rel
            assert sha(p) == gold["files"][f"cases/{name}/{rel}"]["sha256"], p
            protected[str(p.relative_to(ROOT))] = sha(p)
        config = yaml.safe_load((run / "config.yaml").read_text())
        with h5py.File(run / "particles/trajectories.h5") as h:
            assert h.attrs["complete"]
            d = {key: h[key][:] for key in ["time_s", "positions_m", "particle_id", "radii_m", "active"]}
            d.update(
                {
                    key: h[key][0]
                    for key in [
                        "primary_force_N",
                        "secondary_force_N",
                        "local_pressure_complex_Pa",
                        "volume_response_complex_m3",
                        "flow_velocity_m_s",
                        "free_total_velocity_m_s",
                        "total_velocity_m_s",
                    ]
                }
            )
        assert d["positions_m"].shape == (301, 256, 3)
        assert np.array_equal(d["time_s"], np.arange(301) * 0.02)
        field_path = ROOT / config["effective_acoustic_field"]["path"]
        assert sha(field_path) == config["pzt_control"]["field_sha256"]
        protected[str(field_path.relative_to(ROOT))] = sha(field_path)
        field = ComplexAcousticField.load(field_path)
        control = np.array(
            [e["amplitude_Pa"] * np.exp(1j * e["phase_rad"]) * e["enabled"] for e in config["pzt"]]
        )
        reconstructed = np.einsum("i,ixyz->xyz", control, basis)
        basis_error = np.linalg.norm(reconstructed - field.pressure_Pa) / np.linalg.norm(field.pressure_Pa)
        assert basis_error < 1e-12
        surface = ROOT / Path(config["geometry"]["mesh_surface"]).relative_to(
            "/workspace/mevo_acoustic_microbubbles"
        )
        protected[str(surface.relative_to(ROOT))] = sha(surface)
        mesh = pv.read(surface)
        outline = mesh_lines(mesh)
        wire = [outline]
        for x_mm in [-3, -2, 0, 2, 4, 6, 8, 10]:
            wire.append(mesh_lines(mesh, "x", (x_mm * 0.001, 0, 0)))
        wire = np.concatenate(wire)
        keep = np.all((wire[..., 0] >= -3.4) & (wire[..., 0] <= 10.8), axis=1)
        wire = wire[keep]
        table_path = ROOT / config["bubble_response"]["path"]
        protected[str(table_path.relative_to(ROOT))] = sha(table_path)
        table = AdaptiveResponseTableBubbleModel(table_path)
        response = table.harmonic_response(d["radii_m"], field.frequency_Hz, d["local_pressure_complex_Pa"])[
            "harmonics_m3"
        ]
        np.testing.assert_allclose(response[:, 0], d["volume_response_complex_m3"], rtol=1e-12, atol=0)
        grad = field.interpolate_complex_gradient(d["positions_m"][0])
        fp = primary_bjerknes(response[:, 0], grad)
        fs = ClassicalFarFieldBjerknes(enforce_far_field=False).forces(
            d["positions_m"][0], d["radii_m"], response, field.frequency_Hz, config["medium"]["density_kg_m3"]
        )
        fp_error = np.linalg.norm(fp - d["primary_force_N"]) / np.linalg.norm(fp)
        fs_error = np.linalg.norm(fs - d["secondary_force_N"]) / np.linalg.norm(fs)
        assert fp_error < 1e-12 and fs_error < 1e-11, (fp_error, fs_error)
        expected_free = d["flow_velocity_m_s"] + (fp + fs) / (
            6 * np.pi * config["medium"]["viscosity_Pa_s"] * d["radii_m"][:, None]
        )
        np.testing.assert_allclose(expected_free, d["free_total_velocity_m_s"], rtol=1e-11, atol=1e-16)
        distances = np.linalg.norm(d["positions_m"][0, :, None] - d["positions_m"][0, None, :], axis=-1)
        ratios = distances / (d["radii_m"][:, None] + d["radii_m"][None, :])
        scores = np.where(np.triu(ratios >= 5, 1), distances, np.inf)
        i, j = np.unravel_index(np.argmin(scores), scores.shape)
        h = np.arange(1, response.shape[-1] + 1)
        corr = (
            (2 * np.pi * field.frequency_Hz) ** 2
            / 2
            * np.real(np.sum(h * h * response[i] * np.conj(response[j])))
        )
        pair_force = config["medium"]["density_kg_m3"] * corr / (4 * np.pi * distances[i, j] ** 2)
        branches = json.loads((run / "BRANCH_METRICS.json").read_text())
        outcomes = np.zeros(256, int)
        exit_time = np.full(256, np.inf)
        mapping = {int(pid): k for k, pid in enumerate(d["particle_id"])}
        for event in branches["crossings"]:
            k = mapping[event["particle_id"]]
            outcomes[k] = 1 if event["branch"] == "BRANCH_A" else 2
            exit_time[k] = event["crossing_time_s"]
        assert np.array_equal(d["active"], d["time_s"][:, None] < exit_time)
        d.update(
            name=name,
            config=config,
            field=field,
            control=control,
            surface=surface,
            outline=outline,
            wire=wire,
            response_harmonics=response,
            gradient=grad,
            pair=dict(
                indices=[int(i), int(j)],
                ids=d["particle_id"][[i, j]].tolist(),
                distance_m=float(distances[i, j]),
                separation_ratio=float(ratios[i, j]),
                force_N=float(pair_force),
                harmonics=response.shape[-1],
            ),
            branches=branches,
            outcome=outcomes,
            exit_time=exit_time,
            run=str(run.relative_to(ROOT)),
            validation=dict(
                basis_relative_L2=float(basis_error),
                primary_relative_L2=float(fp_error),
                secondary_relative_L2=float(fs_error),
                stored_frames=301,
                initial_particles=256,
            ),
        )
        data[name] = d
    np.testing.assert_array_equal(data["B2"]["positions_m"][0], data["B3"]["positions_m"][0])
    return data, protected


def region(fig, box, relative, projection=None):
    x, y, w, h = box
    rx, ry, rw, rh = relative
    kwargs = dict(projection=projection) if projection else {}
    return fig.add_axes([x + w * rx, y + h * ry, w * rw, h * rh], **kwargs)


def label(fig, box, x, y, text, **kwargs):
    return fig.text(box[0] + box[2] * x, box[1] + box[3] * y, text, **kwargs)


def axes3d(ax, bounds, elev=24, azim=-64):
    ax.set_proj_type("ortho")
    ax.view_init(elev=elev, azim=azim)
    ax.set(xlim=bounds[0], ylim=bounds[1], zlim=bounds[2])
    ax.set_box_aspect([b[1] - b[0] for b in bounds], zoom=1.04)
    for axis, name in zip([ax.xaxis, ax.yaxis, ax.zaxis], ["x (mm)", "y (mm)", "z (mm)"]):
        axis.set_label_text(name, fontsize=8)
        axis.labelpad = 0
        axis.set_pane_color((1, 1, 1, 0))
        axis.line.set_linewidth(0.55)
        axis._axinfo["grid"]["color"] = "#E0E5E9"
        axis._axinfo["grid"]["linewidth"] = 0.35
        axis._axinfo["tick"]["inward_factor"] = 0
        axis._axinfo["tick"]["outward_factor"] = 0.15
    ax.tick_params(pad=0, labelsize=8)


def pzt_geometry(ax, d):
    bounds = [(-3.5, 6), (-7, 7), (-7, 7)]
    axes3d(ax, bounds, elev=22, azim=-60)
    wire = d["wire"]
    take = np.all(wire[..., 0] < 6, axis=1)
    ax.add_collection3d(Line3DCollection(wire[take], colors="#7B8994", linewidths=0.6, alpha=0.6))
    t = np.linspace(0, 2 * np.pi, 97)
    for e in d["config"]["pzt"]:
        c = np.array(e["position_m"]) * 1000
        n = np.array(e["normal"])
        u = np.array([1.0, 0, 0])
        v = np.cross(n, u)
        pts = c + e["size_m"][0] * 500 * (np.cos(t)[:, None] * u + np.sin(t)[:, None] * v)
        color = BLUE if e["phase_rad"] < 0.1 else ORANGE
        ax.add_collection3d(
            Poly3DCollection([pts], facecolors=colors.to_rgba(color, 0.12), edgecolors=color, linewidths=0.8)
        )
        pos = c * 1.25
        ax.text(
            *pos,
            e["id"].replace("_plus", "+").replace("_minus", "-").upper(),
            fontsize=8,
            color=color,
            ha="center",
        )
        ax.quiver(*c, *(n * 1.5), arrow_length_ratio=0.22, color=color, linewidth=0.7)
    ax.set_xticks([0, 5])
    ax.set_yticks([-5, 0, 5])
    ax.set_zticks([-5, 0, 5])


def draw_setup(fig, box, d, letter="a"):
    label(fig, box, 0, 1, f"{letter}  PZT control setup", fontsize=10, fontweight="bold", va="top")
    label(
        fig,
        box,
        0.99,
        1,
        f"{d['name']} | Ideal acoustic sources",
        fontsize=8,
        color=GREY,
        ha="right",
        va="top",
    )
    ax = region(fig, box, [0.015, 0.20, 0.49, 0.69], "3d")
    pzt_geometry(ax, d)
    ax = region(fig, box, [0.555, 0.16, 0.43, 0.62])
    ax.axis("off")
    elements = d["config"]["pzt"]
    cell = [
        [
            e["id"].replace("_plus", "+").replace("_minus", "-").upper(),
            f"{e['amplitude_Pa'] / 1000:.2f}",
            f"{np.rad2deg(e['phase_rad']):.0f}",
            "ON" if e["enabled"] else "OFF",
        ]
        for e in elements
    ]
    table = ax.table(
        cellText=cell,
        colLabels=["Element", "A (kPa)", r"$\phi$ (deg)", "State"],
        cellLoc="center",
        colLoc="center",
        bbox=[0, 0, 1, 1],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    for (row, col), cell in table.get_celld().items():
        cell.set_edgecolor("#D5DDE3")
        cell.set_linewidth(0.45)
        cell.set_facecolor("#EEF2F5" if row == 0 else "white")
        if row == 0:
            cell.get_text().set_weight("bold")
    label(fig, box, 0.555, 0.86, r"$u_i=s_i A_i e^{\mathrm{i}\phi_i}$", fontsize=11)
    label(fig, box, 0.555, 0.015, "f = 2.5 MHz   |   Disc diameter = 6 mm", fontsize=8)


def pressure3d(ax, d):
    field = d["field"]
    axes = [v * 1000 for v in field.axes_m]
    p = abs(field.pressure_Pa) / 1000
    ix = np.argmin(abs(axes[0] + 2))
    iy = np.argmin(abs(axes[1]))
    iz = np.argmin(abs(axes[2]))
    mx = (axes[0] >= -3) & (axes[0] <= 3)
    my = (axes[1] >= -2) & (axes[1] <= 2)
    norm = colors.Normalize(0, 39)
    cmap = plt.get_cmap(CMAP)
    X, Y = np.meshgrid(axes[0][mx], axes[1][my], indexing="ij")
    ax.plot_surface(
        X,
        Y,
        np.zeros_like(X),
        facecolors=cmap(norm(p[mx][:, my, iz])),
        rstride=1,
        cstride=1,
        shade=False,
        linewidth=0,
        antialiased=False,
    )
    X, Z = np.meshgrid(axes[0][mx], axes[2], indexing="ij")
    ax.plot_surface(
        X,
        np.zeros_like(X),
        Z,
        facecolors=cmap(norm(p[mx, iy, :])),
        rstride=1,
        cstride=1,
        shade=False,
        linewidth=0,
        antialiased=False,
        alpha=1.0,
    )
    Y, Z = np.meshgrid(axes[1][my], axes[2], indexing="ij")
    ax.plot_surface(
        np.full_like(Y, axes[0][ix]),
        Y,
        Z,
        facecolors=cmap(norm(p[ix, my, :])),
        rstride=1,
        cstride=1,
        shade=False,
        linewidth=0,
        antialiased=False,
        alpha=1.0,
    )
    axes3d(ax, [(-3, 3), (-2, 2), (-1.6, 1.6)], elev=28, azim=-62)
    ax.set_xticks([-2, 0, 2])
    ax.set_yticks([-2, 0, 2])
    ax.set_zticks([-1, 0, 1])


def pressure2d(ax, d, extent=(-3.2, 10.2, -6.2, 6.2)):
    field = d["field"]
    x, y, z = [a * 1000 for a in field.axes_m]
    iz = np.argmin(abs(z))
    dx, dy = x[1] - x[0], y[1] - y[0]
    image = ax.imshow(
        (abs(field.pressure_Pa[:, :, iz]) / 1000).T,
        origin="lower",
        extent=[x[0] - dx / 2, x[-1] + dx / 2, y[0] - dy / 2, y[-1] + dy / 2],
        vmin=0,
        vmax=39,
        cmap=CMAP,
        interpolation="nearest",
        aspect="equal",
    )
    ax.add_collection(LineCollection(d["outline"][..., :2], colors="white", linewidths=0.65))
    ax.set(xlim=extent[:2], ylim=extent[2:], xlabel="x (mm)", ylabel="y (mm)")
    ax.set_xticks([0, 5, 10])
    ax.set_yticks([-5, 0, 5])
    return image


def draw_acoustics(fig, box, d, letter="b"):
    label(
        fig,
        box,
        0,
        1,
        f"{letter}  Three-dimensional acoustic field",
        fontsize=10,
        fontweight="bold",
        va="top",
    )
    label(fig, box, 0.99, 1, "k-Wave basis superposition", fontsize=8, color=GREY, ha="right", va="top")
    ax = region(fig, box, [0.005, 0.23, 0.46, 0.61], "3d")
    pressure3d(ax, d)
    label(fig, box, 0.065, 0.855, "Orthogonal slices", fontsize=8)
    ax = region(fig, box, [0.54, 0.35, 0.35, 0.49])
    pressure2d(ax, d)
    ax.text(0.035, 0.94, "z = 0 mm", transform=ax.transAxes, va="top", color="white", fontsize=8)
    cax = region(fig, box, [0.93, 0.37, 0.016, 0.45])
    cbar = fig.colorbar(ScalarMappable(norm=colors.Normalize(0, 39), cmap=CMAP), cax=cax)
    cbar.set_ticks([0, 20, 39])
    cbar.set_label(r"$|P|$ (kPa)", labelpad=3)
    label(fig, box, 0.58, 0.015, r"$P(\mathbf{x})=\sum_i u_i H_i(\mathbf{x})$", fontsize=11)


def force_quiver(ax, d):
    x = d["positions_m"][0] * 1000
    force = d["primary_force_N"] * 1e12
    ax.scatter(x[:, 1], x[:, 2], s=1.6, color="#AEBCC6", alpha=0.75, zorder=1)
    indices = np.arange(0, 256, 4)
    q = ax.quiver(
        x[indices, 1],
        x[indices, 2],
        force[indices, 1],
        force[indices, 2],
        angles="xy",
        scale_units="xy",
        scale=400,
        color=BLUE,
        width=0.006,
        headwidth=3.5,
        headlength=4,
    )
    ax.quiverkey(q, 0.59, 1.04, 100, "100 pN", labelpos="E", coordinates="axes", fontproperties={"size": 8})
    ax.set(xlim=(-1.35, 1.35), ylim=(-1.25, 1.25), xlabel="y (mm)", ylabel="z (mm)")
    ax.set_xticks([-1, 0, 1])
    ax.set_yticks([-1, 0, 1])
    ax.set_aspect("equal")


def pair_diagram(ax, d):
    pair = d["pair"]
    dist = pair["distance_m"] * 1e6
    radius = d["radii_m"][pair["indices"][0]] * 1e6
    ax.set(xlim=(-dist * 0.75, dist * 0.75), ylim=(-8, 8))
    ax.set_aspect("equal")
    ax.axis("off")
    for k, c in enumerate([-dist / 2, dist / 2]):
        ax.add_patch(Circle((c, 0), radius, facecolor="#D7E9F3", edgecolor=BLUE, linewidth=0.8))
        ax.plot([c], [0], ".", ms=2, color=BLUE)
        ax.text(c, 3.1, f"ID {pair['ids'][k]}", ha="center", fontsize=8)
        direction = (1 if k == 0 else -1) * np.sign(pair["force_N"])
        start = c + direction * (radius + 0.35)
        end = start + direction * 3.4
        ax.add_patch(
            FancyArrowPatch(
                (start, 0), (end, 0), arrowstyle="-|>", mutation_scale=9, color=ORANGE, linewidth=0.9
            )
        )
    ax.annotate(
        "",
        xy=(dist / 2, -3.4),
        xytext=(-dist / 2, -3.4),
        arrowprops=dict(arrowstyle="<->", lw=0.65, color=INK),
    )
    ax.text(0, -5.35, f"d = {dist:.2f} µm", ha="center", fontsize=8)
    ax.text(
        0, 6.0, f"Pair force = {abs(pair['force_N']) * 1e12:.2f} pN", ha="center", fontsize=8, color=ORANGE
    )


def draw_forces(fig, box, d, letter="c"):
    label(fig, box, 0, 1, f"{letter}  Microbubble acoustic forces", fontsize=10, fontweight="bold", va="top")
    label(fig, box, 0.99, 1, "Cycle-averaged; t = 0", fontsize=8, color=GREY, ha="right", va="top")
    label(fig, box, 0.13, 0.86, "Primary Bjerknes force", fontsize=9, color=BLUE, fontweight="bold")
    label(fig, box, 0.59, 0.86, "Secondary Bjerknes force", fontsize=9, color=ORANGE, fontweight="bold")
    ax = region(fig, box, [0.08, 0.29, 0.35, 0.48])
    force_quiver(ax, d)
    ax = region(fig, box, [0.54, 0.22, 0.43, 0.49])
    pair_diagram(ax, d)
    label(fig, box, 0.53, 0.165, r"$R_0=1.78\ \mathrm{\mu m}$" + "; true pair geometry; t = 0", fontsize=8)
    label(
        fig,
        box,
        0.06,
        0.015,
        r"$\mathbf{F}_{P,i}=-\frac{1}{2}\operatorname{Re}[\hat V_{i,1}\nabla P_i^*]$",
        fontsize=10,
    )
    label(
        fig,
        box,
        0.53,
        0.015,
        r"$\mathbf{F}_{S,ij}=\frac{\rho\langle\dot V_i\dot V_j\rangle}{4\pi d_{ij}^2}\,\mathbf{e}_{ij}$",
        fontsize=10,
    )


def track_segments(d):
    points = d["positions_m"] * 1000
    t = d["time_s"]
    segments = []
    times = []
    for k in range(256):
        stop = np.searchsorted(t, d["exit_time"][k], side="left")
        stop = min(stop, len(t) - 1)
        p = points[: stop + 1, k]
        pt = t[: stop + 1].copy()
        if np.isfinite(d["exit_time"][k]):
            pt[-1] = d["exit_time"][k]
        segments.extend(np.stack([p[:-1], p[1:]], axis=1))
        times.extend((pt[:-1] + pt[1:]) / 2)
    return np.asarray(segments), np.asarray(times)


def trajectory3d(ax, d):
    segments, times = track_segments(d)
    ax.add_collection3d(Line3DCollection(d["wire"], colors="#B4BEC6", linewidths=0.5, alpha=0.65))
    lc = Line3DCollection(segments, cmap="viridis", norm=colors.Normalize(0, 6), linewidths=0.5, alpha=0.8)
    lc.set_array(times)
    ax.add_collection3d(lc)
    x = d["positions_m"][0] * 1000
    ax.scatter(*x.T, s=1.1, color=BLUE, alpha=0.7, depthshade=False)
    axes3d(ax, [(-3, 10.8), (-6.5, 6.5), (-1.5, 1.5)], elev=44, azim=-70)
    ax.set_xticks([0, 5, 10])
    ax.set_yticks([-5, 0, 5])
    ax.set_zticks([0])
    ax.text(8.8, 5.9, 0.4, "A", fontsize=9, fontweight="bold")
    ax.text(8.8, -5.9, 0.4, "B", fontsize=9, fontweight="bold")
    return lc


def draw_tracking(fig, box, d, letter="d"):
    label(fig, box, 0, 1, f"{letter}  Lagrangian particle tracking", fontsize=10, fontweight="bold", va="top")
    label(fig, box, 0.99, 1, "256 bubbles; 0-6 s", fontsize=8, color=GREY, ha="right", va="top")
    ax = region(fig, box, [0, 0.32, 0.55, 0.59], "3d")
    trajectory3d(ax, d)
    cax = region(fig, box, [0.09, 0.215, 0.32, 0.02])
    cbar = fig.colorbar(
        ScalarMappable(norm=colors.Normalize(0, 6), cmap="viridis"), cax=cax, orientation="horizontal"
    )
    cbar.set_ticks([0, 3, 6])
    cbar.set_label("")
    label(fig, box, 0.43, 0.21, "Time (s)", fontsize=8)
    ax = region(fig, box, [0.65, 0.32, 0.32, 0.32])
    t = d["time_s"]
    for outcome, title, color, ls in [
        (0, "Unassigned", GREY, "-"),
        (1, "Branch A", BLUE, "-"),
        (2, "Branch B", ORANGE, "--"),
    ]:
        values = (
            np.sum(t[:, None] < d["exit_time"][None, :], axis=1)
            if outcome == 0
            else np.sum((t[:, None] >= d["exit_time"][None, :]) & (d["outcome"][None, :] == outcome), axis=1)
        )
        ax.step(t, values, where="post", label=title, color=color, ls=ls, lw=0.9)
    ax.set(xlim=(0, 6), ylim=(-5, 270), xlabel="Time (s)", ylabel="Bubble count")
    ax.set_xticks([0, 3, 6])
    ax.set_yticks([0, 128, 256])
    ax.legend(
        loc="upper left", bbox_to_anchor=(0, 1.8), ncol=1, borderaxespad=0, handlelength=1.4, labelspacing=0.1
    )
    label(
        fig,
        box,
        0.10,
        -0.055,
        r"$\mathbf{v}_i^{\mathrm{free}}=\mathbf{u}_{\mathrm{flow}}+\frac{\mathbf{F}_{P,i}+\sum_{j\ne i}\mathbf{F}_{S,ij}}{6\pi\mu R_i}$",
        fontsize=10,
    )
    label(fig, box, 0.66, -0.055, "RK2 + geometric constraints", fontsize=8)


def save(fig, out, name, preview=False):
    fig.canvas.draw()
    texts = []
    # Check figure-level annotations; 3D axis labels are additionally inspected in previews.
    fw, fh = fig.canvas.get_width_height()
    renderer = fig.canvas.get_renderer()
    for t in fig.texts:
        b = t.get_window_extent(renderer)
        assert b.x0 >= 0 and b.y0 >= 0 and b.x1 <= fw and b.y1 <= fh, (name, t.get_text(), b.bounds, (fw, fh))
        texts.append(
            dict(text=t.get_text(), font_size_pt=float(t.get_fontsize()), bounds_pixels=list(b.bounds))
        )
    clipped = []
    for text_artist in fig.findobj(matplotlib.text.Text):
        if isinstance(text_artist, Text3D):
            continue  # Text3D extents use raw data coordinates; inspect projected labels in previews.
        if not text_artist.get_visible() or not text_artist.get_text().strip():
            continue
        bb = text_artist.get_window_extent(renderer)
        if bb.x0 < -1 or bb.y0 < -1 or bb.x1 > fw + 1 or bb.y1 > fh + 1:
            clipped.append(dict(text=text_artist.get_text(), bounds=list(bb.bounds)))
    assert not clipped, (name, clipped)
    if preview:
        fig.savefig(out / (name + "_preview.png"), dpi=180)
    else:
        fig.savefig(
            out / (name + ".pdf"),
            dpi=600,
            metadata={"Creator": "MEVO sealed-data figure renderer", "Title": name},
        )
        fig.savefig(out / (name + ".svg"), dpi=600)
        fig.savefig(out / (name + ".png"), dpi=600)
        with Image.open(out / (name + ".png")) as im:
            im.convert("RGB").save(out / (name + ".tif"), compression="tiff_lzw", dpi=(600, 600))
    result = dict(
        name=name,
        size_mm=(fig.get_size_inches() * 25.4).tolist(),
        texts=texts,
        visible_2d_text_within_page=True,
        projected_3d_text_requires_visual_review=True,
    )
    plt.close(fig)
    return result


def comparison(data, out, preview):
    fig = plt.figure(figsize=(180 / 25.4, 200 / 25.4))
    for j, name in enumerate(["B2", "B3"]):
        d = data[name]
        x = 0.10 + j * 0.48
        fig.text(
            x,
            0.965,
            ("a" if j == 0 else "b") + "  " + name + " acoustic field",
            fontsize=10,
            fontweight="bold",
        )
        ax = fig.add_axes([x, 0.61, 0.365, 0.30])
        pressure2d(ax, d)
        ax.text(0.03, 0.94, "z = 0 mm", color="white", transform=ax.transAxes, va="top", fontsize=8)
        fig.text(
            x,
            0.47,
            ("c" if j == 0 else "d") + "  " + name + " bubble trajectories",
            fontsize=10,
            fontweight="bold",
        )
        ax = fig.add_axes([x, 0.135, 0.365, 0.30])
        segments, times = track_segments(d)
        ax.add_collection(LineCollection(d["outline"][..., :2], colors="#A5B0B9", linewidths=0.65))
        lc = LineCollection(
            segments[..., :2], cmap="viridis", norm=colors.Normalize(0, 6), linewidths=0.45, alpha=0.8
        )
        lc.set_array(times)
        ax.add_collection(lc)
        ax.set(xlim=(-3.2, 10.5), ylim=(-6.3, 6.3), xlabel="x (mm)", ylabel="y (mm)", aspect="equal")
        ax.set_xticks([0, 5, 10])
        ax.set_yticks([-5, 0, 5])
        ax.text(0.03, 0.95, "x-y projection", transform=ax.transAxes, va="top", fontsize=8)
        fig.text(
            x,
            0.067,
            f"A = {d['branches']['N_A']}; B = {d['branches']['N_B']}; unassigned = {d['branches']['N_unassigned']}",
            fontsize=8,
        )
    cax = fig.add_axes([0.28, 0.545, 0.46, 0.011])
    cb = fig.colorbar(
        ScalarMappable(norm=colors.Normalize(0, 39), cmap=CMAP), cax=cax, orientation="horizontal"
    )
    cb.set_ticks([0, 20, 39])
    cb.set_label(r"$|P|$ (kPa)", labelpad=0)
    cax = fig.add_axes([0.28, 0.03, 0.46, 0.011])
    cb = fig.colorbar(
        ScalarMappable(norm=colors.Normalize(0, 6), cmap="viridis"), cax=cax, orientation="horizontal"
    )
    cb.set_ticks([0, 3, 6])
    fig.text(0.77, 0.028, "Time (s)", fontsize=8)
    return save(fig, out, "FigS1_B2_B3_comparison", preview)


def captions(data):
    pair = data["B2"]["pair"]
    return f"""Figure 1. Data-grounded PZT-to-microbubble workflow illustrated with the sealed Stage 4I B2 case. (a) Actual positions, normals and 6 mm diameters of four ideal disc sources, with source-pressure amplitudes, phases and enabled states. Blue and orange distinguish 0 and 180 degree source phases. The phasor uses p(t) = Re[P exp(i omega t)]; s_i is the source switch. This is an acoustic boundary-source representation, not an electromechanical PZT model. (b) Native k-Wave complex-basis superposition at 2.5 MHz in homogeneous water. Orthogonal slices in the indicated local region and the z = 0 mm amplitude section use the same linear 0-39 kPa scale. Vessel outlines provide geometric context and are not acoustic material interfaces. No spatial smoothing or intensity rescaling is applied. (c) Saved primary-force y-z components at t = 0: all 256 bubble centers are shown and arrows use every fourth global particle index (64 arrows) to reduce overlap. Arrow scale is 100 pN. The secondary diagram shows actual radii and separation of particle IDs {pair["ids"][0]} and {pair["ids"][1]}, the closest initial pair with d/(R_i+R_j) >= 5. A rigid pair-axis coordinate transformation preserves geometry. Pair-force magnitude ({abs(pair["force_N"]) * 1e12:.2f} pN) is reconstructed from the same response table, retaining all {pair["harmonics"]} volume harmonics and their h-squared derivative weights; arrow directions encode the computed attraction/repulsion and their graphical lengths are schematic. This selected separation does not establish complete far-field or experimental validity. (d) All 256 recorded trajectories over 0-6 s, with color denoting time and a cropped vessel view of the tracking region. Lines join saved states (20 ms spacing) and terminate at each recorded first-crossing event; they are not recovered integration substeps. The nominal integration timestep is 0.5 ms with adaptive displacement constraints. The equation describes unconstrained Stokes migration; RK2 with wall/contact constraints yields the stored motion. Counts distinguish branch arrival from observation-censored, unassigned particles. B2 ends with A = 76, B = 0 and unassigned = 180. Microbubble position markers are enlarged for visibility; trajectory coordinates are not magnified.

Figures 1a-1d. Separately exported components of Figure 1, with the same data, force conventions, scales and display qualifications.

Figure S1. B2/B3 comparison using identical initial positions, particle identities and radii. Top: native z = 0 mm pressure amplitudes on a shared scale. Bottom: x-y projections of all recorded trajectories with a shared 0-6 s color scale. The orthographic projections can overlap tracks separated in z. B2 ends with 76 branch-A arrivals, 0 branch-B arrivals and 180 unassigned particles; B3 ends with 0, 88 and 168, respectively. These are single matched-seed simulations, not independent experimental replicates or an estimate of population uncertainty.

Model scope. The source runs are Stage 4I, not the earlier Stage-1 software demo. Existing limits on near-contact forces, near-wall hydrodynamics, quasi-steady drag and finite slip Reynolds number remain. Contact enforces geometric nonpenetration only and does not model adhesion or fusion. No figure establishes local physical accuracy. No new acoustic, CFD or trajectory simulation was run for these figures.
"""


def export_data(data, out):
    with h5py.File(out / "source_data.h5", "w") as h:
        h.attrs["scope"] = "Plot data copied/derived from sealed Stage 4I B2/B3; original inputs unchanged"
        for name, d in data.items():
            g = h.create_group(name)
            for key in [
                "time_s",
                "positions_m",
                "particle_id",
                "radii_m",
                "active",
                "primary_force_N",
                "secondary_force_N",
                "flow_velocity_m_s",
                "free_total_velocity_m_s",
                "response_harmonics",
                "control",
                "outcome",
                "exit_time",
            ]:
                g.create_dataset(key, data=d[key], compression="gzip")
            for key, arr in zip("xyz", d["field"].axes_m):
                g.create_dataset("axes/" + key + "_m", data=arr)
            iz = np.argmin(abs(d["field"].axes_m[2]))
            g.create_dataset("pressure_z0_Pa", data=d["field"].pressure_Pa[:, :, iz], compression="gzip")
            ix = np.argmin(abs(d["field"].axes_m[0] + 0.002))
            iy = np.argmin(abs(d["field"].axes_m[1]))
            g.create_dataset("pressure_xminus2mm_Pa", data=d["field"].pressure_Pa[ix], compression="gzip")
            g.create_dataset("pressure_y0_Pa", data=d["field"].pressure_Pa[:, iy], compression="gzip")
            g.create_dataset("vessel_wire_mm", data=d["wire"], compression="gzip")
            g.create_dataset("vessel_z0_outline_mm", data=d["outline"], compression="gzip")
            g.attrs["pair_example_json"] = json.dumps(d["pair"])
            g.attrs["pzt_configuration_json"] = json.dumps(d["config"]["pzt"])
            g.attrs["frequency_Hz"] = d["field"].frequency_Hz
            g.attrs["force_and_velocity_time_s"] = 0.0
            g.attrs["response_harmonics_units"] = "m^3"
            g.attrs["outcome_codes"] = "0: unassigned at 6 s; 1: branch A; 2: branch B"
            g.attrs["exit_time_units"] = "s; infinity means no recorded crossing"
            g.attrs["control_units"] = "Pa; Re[P exp(+i omega t)] convention"
    with (out / "pzt_controls.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "case",
                "element",
                "x_m",
                "y_m",
                "z_m",
                "diameter_m",
                "amplitude_Pa",
                "phase_rad",
                "frequency_Hz",
                "enabled",
            ]
        )
        for name, d in data.items():
            for e in d["config"]["pzt"]:
                writer.writerow(
                    [
                        name,
                        e["id"],
                        *e["position_m"],
                        e["size_m"][0],
                        e["amplitude_Pa"],
                        e["phase_rad"],
                        d["field"].frequency_Hz,
                        e["enabled"],
                    ]
                )
    with (out / "initial_forces.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "case",
                "particle_id",
                "x_m",
                "y_m",
                "z_m",
                "Fp_x_N",
                "Fp_y_N",
                "Fp_z_N",
                "Fs_x_N",
                "Fs_y_N",
                "Fs_z_N",
            ]
        )
        for name, d in data.items():
            for k, pid in enumerate(d["particle_id"]):
                writer.writerow(
                    [name, pid, *d["positions_m"][0, k], *d["primary_force_N"][k], *d["secondary_force_N"][k]]
                )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--preview", action="store_true")
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    font = style()
    data, protected = load_data()
    d = data["B2"]
    plots = []
    fig = plt.figure(figsize=(180 / 25.4, 218 / 25.4))
    funcs = [draw_setup, draw_acoustics, draw_forces, draw_tracking]
    boxes = [
        [0.055, 0.755, 0.89, 0.215],
        [0.055, 0.517, 0.89, 0.215],
        [0.055, 0.279, 0.89, 0.215],
        [0.055, 0.041, 0.89, 0.215],
    ]
    for fn, box, letter in zip(funcs, boxes, "abcd"):
        fn(fig, box, d, letter)
    for y in [0.742, 0.504, 0.266]:
        fig.add_artist(
            plt.Line2D([0.055, 0.945], [y, y], transform=fig.transFigure, color="#D6DEE4", linewidth=0.55)
        )
        fig.add_artist(
            FancyArrowPatch(
                (0.033, y + 0.009),
                (0.033, y - 0.009),
                transform=fig.transFigure,
                arrowstyle="-|>",
                mutation_scale=8,
                color=BLUE,
                linewidth=0.8,
            )
        )
    plots.append(save(fig, out, "Fig1_B2_workflow", args.preview))
    for fn, letter, name in zip(
        funcs, "abcd", ["PZT_control", "acoustic_field", "microbubble_forces", "particle_tracking"]
    ):
        fig = plt.figure(figsize=(180 / 25.4, 70 / 25.4))
        fn(fig, [0.055, 0.10, 0.89, 0.82], d, letter)
        plots.append(save(fig, out, f"Fig1{letter}_{name}", args.preview))
    plots.append(comparison(data, out, args.preview))
    assert all(sha(ROOT / p) == digest for p, digest in protected.items())
    report = dict(
        status="PREVIEW" if args.preview else "RENDERED",
        cases={n: dict(run=d["run"], checks=d["validation"], pair=d["pair"]) for n, d in data.items()},
        protected_inputs=protected,
        renderer_sha256=sha(__file__),
        font=font,
        formats=["PDF", "SVG", "PNG", "TIFF"],
        raster_dpi=600,
        scope="Existing Stage 4I; figures are derived displays, no new simulation or physical validation",
        typography="English; 8-11 pt; true-type fonts embedded in PDF, editable SVG text",
        style_references=[
            "https://journals.plos.org/plosone/s/figures",
            "https://www.nature.com/nature/for-authors/final-submission",
        ],
        plots=plots,
    )
    (out / "FIGURE_PROVENANCE.json").write_text(json.dumps(report, indent=2) + "\n")
    if not args.preview:
        (out / "CAPTIONS.txt").write_text(captions(data))
        export_data(data, out)
        shutil.copy2(__file__, out / "paper_workflow_figures.py")
    print(json.dumps(dict(status=report["status"], output=str(out), checks=report["cases"])), flush=True)


if __name__ == "__main__":
    main()
