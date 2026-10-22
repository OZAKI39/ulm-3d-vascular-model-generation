"""Two data-grounded mechanism stories with actual PZT layout and true-radius bubbles."""

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
from time import perf_counter

import h5py
import numpy as np
import pyvista as pv
from PIL import Image, ImageDraw, ImageFont
import imageio_ffmpeg
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/stage4/visualization/aggregation_translation_PZT_v1_20261008"
REF = ROOT / "results/stage4/publication_figures/PZT_control_reference_v6"
F4 = ROOT / "results/fast_pzt/f4"
F3 = ROOT / "results/fast_pzt/f3"
R = 1.78e-6
BLUE = "#3288d0"
ORANGE = "#c78352"
BG = "#080d15"
FONT = Path("/mnt/c/Windows/Fonts/arial.ttf")
if FONT.exists():
    font_manager.fontManager.addfont(str(FONT))
plt.rcParams.update(
    {
        "font.family": "Arial" if FONT.exists() else "DejaVu Sans",
        "font.size": 11,
        "text.color": "#263746",
        "axes.labelcolor": "#263746",
        "xtick.color": "#263746",
        "ytick.color": "#263746",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
    }
)


def write(p, d):
    p.write_text(
        json.dumps(d, indent=2, default=lambda a: a.tolist() if hasattr(a, "tolist") else str(a)) + "\n"
    )


def sha(p):
    with Path(p).open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def js(p):
    return json.loads(Path(p).read_text())


def read(path):
    with h5py.File(path / "trajectory.h5") as h:
        return dict(
            time=h["time_s"][:],
            x=h["positions_m"][:, 0],
            id=h["particle_id"][:],
            r=h["radii_m"][:],
            nmax=h["metrics/largest_size"][:, 0],
            nc=h["metrics/cluster_count"][:, 0],
            rg=h["metrics/global_rg_m"][:, 0],
            eigen=h["metrics/global_eigenvalues_m2"][:, 0],
        )


def load():
    paths = {
        "aggregation_on": F4 / "aggregation_v1/parent_wall_on",
        "aggregation_off": OUT / "data/secondary_off",
        "translation_baseline": F3 / "repeats_v1/stable_reduction_0",
        "translation_controlled": F3 / "controls_v1/amplitude",
    }
    data = {k: read(p) for k, p in paths.items()}
    for a, b in [("aggregation_on", "aggregation_off"), ("translation_baseline", "translation_controlled")]:
        np.testing.assert_array_equal(data[a]["x"][0], data[b]["x"][0])
        np.testing.assert_array_equal(data[a]["time"], data[b]["time"])
    for key, d in data.items():
        assert np.isfinite(d["x"]).all()
        np.testing.assert_array_equal(d["id"], np.arange(32))
        assert (d["r"] == R).all()
        assert js(paths[key] / "SUMMARY.json")["numerical_quality"] == "PASS"
    return paths, data


def reference():
    info = js(REF / "PROVENANCE.json")
    file = REF / "paper_pzt_reference.py"
    assert sha(file) == info["renderer_sha256"]
    spec = importlib.util.spec_from_file_location("v6_style", file)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return info, mod


def vessel():
    mesh = pv.read(F4 / "inputs/labelled_lumen_m.vtp")
    return mesh.extract_cells(mesh["ModelFaceID"] == 1).extract_surface(algorithm="dataset_surface")


def setup(p, white):
    p.set_background("white" if white else BG)
    p.render_window.SetMultiSamples(0)
    p.render_window.SetAlphaBitPlanes(1)
    p.enable_depth_peeling(number_of_peels=80, occlusion_ratio=0)


def font(size):
    return ImageFont.truetype(str(FONT) if FONT.exists() else "DejaVuSans.ttf", size)


def global_image(event, region, white=True, size=(1400, 1150)):
    info, ref = reference()
    p = pv.Plotter(off_screen=True, window_size=size, lighting="three lights")
    setup(p, white)
    wall = vessel()
    wall.points = wall.points.astype(float) * 1000
    p.add_mesh(
        wall,
        color="#b85c64",
        opacity=0.88,
        smooth_shading=True,
        specular=0.24,
        specular_power=40,
        ambient=0.2,
        diffuse=0.72,
        backface_params=dict(
            color="#44212b", opacity=0.96, ambient=0.12, diffuse=0.45, specular=0.03, specular_power=30
        ),
    )
    rims = wall.extract_feature_edges(
        boundary_edges=True, non_manifold_edges=False, feature_edges=False, manifold_edges=False
    )
    p.add_mesh(rims, color="#c57981", lighting=False, line_width=1)
    elements = info["source_elements"]
    labels = []
    for i, e in enumerate(elements):
        disc = ref.disc_surface(e)
        n = np.array(e["normal"])
        back = disc.copy(deep=True)
        back.points -= n * 0.35
        ff = back.faces.reshape(-1, 4).copy()
        ff[:, [2, 3]] = ff[:, [3, 2]]
        back.faces = ff
        rim = disc.points[1:]
        rear = rim - n * 0.35
        count = len(rim)
        side = pv.PolyData(
            np.vstack([rim, rear]),
            np.array([[4, j, j + count, (j + 1) % count + count, (j + 1) % count] for j in range(count)]),
        )
        prominent = e["id"] in ["y_minus", "z_plus"]
        factor = 1 if event["enabled"][i] else 0.20
        p.add_mesh(
            disc,
            color="#8aa7b9" if i == 1 else "#a2b0b8",
            opacity=(0.42 if prominent else 0.27) * factor,
            ambient=0.55,
            diffuse=0.45,
            specular=0.30,
            specular_power=50,
            smooth_shading=True,
        )
        p.add_mesh(
            back,
            color="#b4c2cb",
            opacity=(0.55 if prominent else 0.30) * factor,
            ambient=0.32,
            diffuse=0.50,
            specular=0.88,
            specular_power=68,
            smooth_shading=True,
        )
        p.add_mesh(
            side,
            color="#8499a5",
            opacity=(0.65 if prominent else 0.45) * factor,
            specular=0.30,
            specular_power=30,
            smooth_shading=True,
        )
        for ring in [rim, rear]:
            curve = pv.lines_from_points(np.vstack([ring, ring[:1]]), close=False)
            p.add_mesh(
                curve.tube(radius=0.018, n_sides=8),
                color="#7c909e",
                opacity=0.75 * factor,
                smooth_shading=True,
                specular=0.55,
                specular_power=50,
            )
        labels.append(
            (
                np.array(e["position_m"]) * 1000,
                ["Y+", "Y−", "Z+", "Z−"][i] + (" ON" if event["enabled"][i] else " OFF"),
            )
        )
    # Same reference scene lighting, including the two broad outward reflections.
    lookup = {e["id"]: e for e in elements}
    for li in info["render"]["outward_reflection_lights"]:
        e = lookup[li["element"]]
        center = np.array(e["position_m"]) * 1000 - np.array(e["normal"]) * 0.35
        p.add_light(
            pv.Light(
                position=li["position_mm"],
                focal_point=center,
                color="#eef4ff",
                intensity=0.65,
                positional=True,
                cone_angle=180,
                attenuation_values=(1, 0, 0),
            )
        )
    xmin, xmax, ymin, ymax, zmin, zmax = -14.5, 11.5, -7.5, 7.5, -5, 5
    segments = []
    for x in [-14.5, -10, -5, 0, 5, 10, 11.5]:
        if x != xmax:
            segments.extend([[[x, ymin, zmin], [x, ymax, zmin]], [[x, ymax, zmin], [x, ymax, zmax]]])
    for y in [-7.5, -4, 0, 4, 7.5]:
        if y != ymin:
            segments.append([[xmin, y, zmin], [xmax, y, zmin]])
        segments.append([[xmin, y, zmin], [xmin, y, zmax]])
    for z in [-5, -2, 0, 2, 5]:
        segments.extend([[[xmin, ymax, z], [xmax, ymax, z]], [[xmin, ymin, z], [xmin, ymax, z]]])
    grid = pv.PolyData(
        np.array(segments).reshape(-1, 3),
        lines=np.array([[2, 2 * i, 2 * i + 1] for i in range(len(segments))]),
    )
    p.add_mesh(grid, color="#bec9d1" if white else "#475565", line_width=0.65)
    cam = info["render"]
    p.camera_position = [cam["camera_eye_mm"], cam["camera_focal_mm"], cam["camera_up"]]
    p.enable_parallel_projection()
    pts = np.array(segments).reshape(-1, 3)
    right = np.array(cam["camera_right"])
    up = np.array(cam["camera_up"])
    scale = max(np.ptp(pts @ up) / 1.40, np.ptp(pts @ right) / (1.45 * size[0] / size[1]))
    p.camera.parallel_scale = scale
    p.show(auto_close=False)
    im = Image.fromarray(p.screenshot(return_img=True))
    p.close()
    draw = ImageDraw.Draw(im)

    def xy(point):
        delta = np.array(point) - np.array(cam["camera_focal_mm"])
        return np.array(
            [
                size[0] / 2 + delta @ right / (2 * scale) * size[1],
                size[1] / 2 - delta @ up / (2 * scale) * size[1],
            ]
        )

    color = "#263746" if white else "#d4dfeb"
    fs = max(16, round(size[0] / 48))
    f = font(fs)
    offsets = [(48, -33), (-108, 2), (-58, -72), (-52, 73)]
    for (point, text), offset in zip(labels, offsets):
        pos = xy(point)
        off = np.array(offset) * size[0] / 1000
        end = pos + off
        draw.line([tuple(pos), tuple(end)], fill="#637e93", width=max(1, size[0] // 700))
        draw.text(tuple(end + [-15, -10]), text, font=f, fill=color, anchor="lm" if off[0] > 0 else "rm")
    # Screen-space locator at the actual projected origin; its size is not a
    # bubble diameter or an outline of the microscopic observation volume.
    locator = xy(np.asarray(region) * 1000)
    reach = max(4, size[0] / 150)
    draw.line([tuple(locator - [reach, 0]), tuple(locator + [reach, 0])], fill=color, width=2)
    draw.line([tuple(locator - [0, reach]), tuple(locator + [0, reach])], fill=color, width=2)
    for point, text, offset in [
        ([-10, -7.5, -5], "-10", (0, 14)),
        ([0, -7.5, -5], "0", (0, 14)),
        ([10, -7.5, -5], "10", (0, 14)),
        ([11.5, -4, -5], "-4", (16, 8)),
        ([11.5, 0, -5], "0", (16, 8)),
        ([11.5, 4, -5], "4", (16, 8)),
        ([11.5, 7.5, -4], "-4", (18, 0)),
        ([11.5, 7.5, 0], "0", (18, 0)),
        ([11.5, 7.5, 4], "4", (18, 0)),
        ([-1.5, -7.5, -5], "X (mm)", (0, 42)),
        ([11.5, 0, -5], "Y (mm)", (52, 32)),
        ([11.5, 7.5, 0], "Z (mm)", (72, 0)),
    ]:
        pp = xy(point) + np.array(offset) * size[0] / 1000
        draw.text(tuple(pp), text, font=font(max(14, round(fs * 0.8))), fill=color, anchor="mm")
    return im


def local_scene(states, wall_enabled, white, size=(1100, 900)):
    # One fixed origin and camera per comparison; no recentering per video frame.
    origin = states[0].mean(axis=0)
    coords = (states - origin) * 1e6
    lo = coords.min(axis=(0, 1)) - R * 1e6 - 3
    hi = coords.max(axis=(0, 1)) + R * 1e6 + 3
    extent = max(20, float(np.max(hi - lo)))
    center = (lo + hi) / 2
    bounds = np.array([center - extent / 2, center + extent / 2]).T.ravel()
    p = pv.Plotter(off_screen=True, window_size=size, lighting="three lights", border=False)
    setup(p, white)
    if wall_enabled:
        wall = vessel()
        wall.points = (wall.points.astype(float) - origin) * 1e6
        wall = wall.clip_box(bounds, invert=False)
        if wall.n_cells:
            p.add_mesh(wall, color="#b85c64", opacity=0.30, smooth_shading=False, ambient=0.25, specular=0.24)
    sphere = pv.Sphere(radius=R * 1e6, theta_resolution=30, phi_resolution=24)
    actors = []
    for pos in coords[0]:
        a = p.add_mesh(
            sphere,
            color="#2875b7" if white else BLUE,
            smooth_shading=True,
            ambient=0.23,
            diffuse=0.72,
            specular=0.50,
            specular_power=36,
            reset_camera=False,
        )
        a.position = pos
        actors.append(a)
    p.camera_position = [center + extent * np.array([0.85, -1.1, 1.5]), center, [0, 0, 1]]
    p.enable_parallel_projection()
    p.camera.parallel_scale = extent * 0.98
    xmin, xmax, ymin, ymax, zmin, zmax = bounds
    ticks = [np.linspace(lo, hi, 3) for lo, hi in bounds.reshape(3, 2)]
    segments = []
    for x in ticks[0]:
        segments.extend([[[x, ymin, zmin], [x, ymax, zmin]], [[x, ymax, zmin], [x, ymax, zmax]]])
    for y in ticks[1]:
        segments.extend([[[xmin, y, zmin], [xmax, y, zmin]], [[xmin, y, zmin], [xmin, y, zmax]]])
    for z in ticks[2]:
        segments.extend([[[xmin, ymin, z], [xmin, ymax, z]], [[xmin, ymax, z], [xmax, ymax, z]]])
    grid = pv.PolyData(
        np.array(segments).reshape(-1, 3),
        lines=np.array([[2, 2 * i, 2 * i + 1] for i in range(len(segments))]),
    )
    p.add_mesh(grid, color="#b3bdc8" if white else "#566473", line_width=0.7, lighting=False)
    p.show(auto_close=False)
    return (
        p,
        actors,
        origin,
        dict(
            origin_m=origin,
            bounds_um=bounds,
            camera=list(p.camera_position),
            parallel_scale_um=p.camera.parallel_scale,
            radius_m=R,
            background="white" if white else BG,
        ),
    )


def image_state(scene, state, mark_centroid=False):
    p, actors, origin, _ = scene
    for actor, pos in zip(actors, (state - origin) * 1e6):
        actor.position = pos
    p.render()
    im = Image.fromarray(p.screenshot(return_img=True))

    def screen(point):
        p.renderer.SetWorldPoint(*point, 1)
        p.renderer.WorldToDisplay()
        xx, yy, _ = p.renderer.GetDisplayPoint()
        return np.array([xx, im.height - 1 - yy])

    # Horizontal labels stay readable and inside each panel. The three open
    # grid planes and tick locations remain in actual local micrometres.
    d = ImageDraw.Draw(im)
    xmin, xmax, ymin, ymax, zmin, zmax = scene[3]["bounds_um"]
    color = "#526170" if scene[3]["background"] == "white" else "#a6b6c6"
    fs = max(16, round(im.width / 37))

    def label(point, text, offset):
        pt = screen(point) + offset
        bb = d.textbbox((0, 0), text, font=font(fs))
        half = (bb[2] - bb[0]) / 2 + 3
        pt[0] = np.clip(pt[0], half, im.width - half)
        pt[1] = np.clip(pt[1], fs, im.height - fs)
        d.text(tuple(pt), text, font=font(fs), fill=color, anchor="mm")

    for x in np.linspace(xmin, xmax, 3):
        label([x, ymin, zmin], f"{x:.1f}", [0, 14])
    for y in np.linspace(ymin, ymax, 3)[1:]:
        label([xmax, y, zmin], f"{y:.1f}", [18, 6])
    for z in np.linspace(zmin, zmax, 3)[1:]:
        label([xmax, ymax, z], f"{z:.1f}", [21, 0])
    label([(xmin + xmax) / 2, ymin, zmin], "x (µm)", [-15, 42])
    label([xmax, (ymin + ymax) / 2, zmin], "y (µm)", [36, 30])
    label([xmax, ymax, (zmin + zmax) / 2], "z (µm)", [48, -26])
    if mark_centroid:
        # Markers use the same camera projection as the spheres, with no shift
        # magnification; origin is the common initial centroid.
        start = screen([0, 0, 0])
        end = screen((state.mean(0) - origin) * 1e6)
        d.line([tuple(start), tuple(end)], fill="#dea75a", width=2)
        color = "#493e58" if scene[3]["background"] == "white" else "#f6dd9a"
        for point, s, col in [(start, 8, "#91949b"), (end, 5, color)]:
            for direction in [np.array([s, 0]), np.array([0, s])]:
                d.line([tuple(point - direction), tuple(point + direction)], fill=col, width=2)
    return im


def selected(kind, data):
    return (
        (data["aggregation_on"], data["aggregation_off"])
        if kind == "aggregation"
        else (data["translation_baseline"], data["translation_controlled"])
    )


def evidence(paths, data):
    OUT.mkdir(exist_ok=True, parents=True)
    for name, d in data.items():
        np.savez_compressed(OUT / "data" / (name + ".npz"), **d)
    a, b = selected("aggregation", data)
    c, d = selected("translation", data)
    info = js(F3 / "MORPHODYNAMICS_COMPARISON.json")["effects"]["amplitude"]
    report = dict(
        sources={
            k: dict(path=str(p), sha256=sha(p / "trajectory.h5"), config=js(p / "CONFIG.json"))
            for k, p in paths.items()
        },
        reference_png=dict(
            path=str(REF / "Fig1a_PZT_control.png"), sha256=sha(REF / "Fig1a_PZT_control.png")
        ),
        aggregation=dict(
            duration_s=float(a["time"][-1]),
            same_initial_positions=True,
            same_controls=True,
            initial_cluster_count=int(a["nc"][0]),
            on_final_largest=int(a["nmax"][-1]),
            off_final_largest=int(b["nmax"][-1]),
            initial_rg_um=float(a["rg"][0] * 1e6),
            on_final_rg_um=float(a["rg"][-1] * 1e6),
            off_final_rg_um=float(b["rg"][-1] * 1e6),
            meaning="In this matched finite-duration real-PZT F4 experiment, Secondary ON aggregates; OFF does not. Primary, mobility, contact and WALL stay enabled.",
        ),
        translation=dict(
            duration_s=float(c["time"][-1]),
            same_initial_positions=True,
            endpoint_centroid_contrast_um=float(
                np.linalg.norm(d["x"][-1].mean(0) - c["x"][-1].mean(0)) * 1e6
            ),
            endpoint_rg_difference_percent=float(abs(d["rg"][-1] / c["rg"][-1] - 1) * 100),
            inherited_shape_assessment=info,
            meaning="Whole-cluster centroid translation changes under PZT control; robust shape control is not established. Actual rotation and internal rearrangement are retained, not suppressed.",
        ),
        scope="EXPLORATORY. First story: F4 geometric WALL only. Second: archived F3, zero flow and no wall forces; global vessel is geometric context.",
        visualization="Actual saved positions, all 32 IDs, true 1.78 um radii, fixed cameras, no trajectory interpolation or manufactured displacement",
        PZT="Actual four source centres/normals/diameters; original v6 material style. Thin display backing/rims are illustrative. OFF sources are faded, not removed.",
        mesh=dict(
            path=str(F4 / "inputs/labelled_lumen_m.vtp"), sha256=sha(F4 / "inputs/labelled_lumen_m.vtp")
        ),
    )
    assert report["aggregation"]["on_final_largest"] == 32 and report["aggregation"]["off_final_largest"] == 1
    write(OUT / "EVIDENCE.json", report)


def figures(paths, data):
    render_meta = {}
    for kind in ["aggregation", "translation"]:
        a, b = selected(kind, data)
        states = np.stack([a["x"][0], a["x"][-1], b["x"][-1]])
        scene = local_scene(states, kind == "aggregation", True, size=(800, 900))
        imgs = [image_state(scene, x, kind == "translation") for x in states]
        scene[0].close()
        render_meta[kind] = scene[3]
        event = js(
            paths["aggregation_on" if kind == "aggregation" else "translation_controlled"] / "CONFIG.json"
        )["events"][0]
        global_img = global_image(event, a["x"][0].mean(0))
        fig = plt.figure(figsize=(16, 8.5), layout="constrained")
        fig.get_layout_engine().set(rect=(0.018, 0.01, 0.965, 0.975))
        gs = fig.add_gridspec(3, 4, height_ratios=[0.13, 1, 0.64], width_ratios=[1.35, 1, 1, 1])
        top = fig.add_subplot(gs[0, :])
        top.axis("off")
        title = (
            "Bubble–bubble attraction enables aggregation"
            if kind == "aggregation"
            else "PZT shifts the cluster; robust shape control is not established"
        )
        top.text(0, 1, title, fontsize=21, va="top", transform=top.transAxes)
        subtitle = (
            "Matched initial positions and PZT drive; only Secondary Bjerknes force is switched OFF"
            if kind == "aggregation"
            else "Same mature cluster, two PZT settings; actual internal rearrangements are retained"
        )
        top.text(0, 0, subtitle, fontsize=11, va="bottom", transform=top.transAxes, color="#526170")
        ax = fig.add_subplot(gs[1, 0])
        ax.imshow(global_img)
        ax.axis("off")
        ax.set_title("PZT layout and vessel context", loc="left", fontsize=13)
        ax = fig.add_subplot(gs[2, 0])
        ax.axis("off")
        ax.text(
            0.02,
            0.92,
            "Four discs: 6 mm diameter\nCentres: Y = ±4 mm, Z = ±4 mm\nCross: observation location\n\nLocal panels share a fixed camera\nand physical scale.",
            transform=ax.transAxes,
            fontsize=10,
            va="top",
        )
        ax.text(
            0.02,
            0.32,
            "All four PZT elements stay ON."
            if kind == "aggregation"
            else "Shown: controlled PZT setting.\nUnchanged: all four elements ON.\n\nLocal cross: cluster centroid.\nGrey cross: initial centroid.",
            transform=ax.transAxes,
            fontsize=10,
            va="top",
        )
        titles = (
            ["Same dispersed initial state", "Attraction ON · 600 ms", "Attraction OFF · 600 ms"]
            if kind == "aggregation"
            else ["Same mature initial state", "Unchanged PZT · 120 ms", "Controlled PZT · 120 ms"]
        )
        for j, (img, title) in enumerate(zip(imgs, titles), 1):
            ax = fig.add_subplot(gs[1, j])
            ax.imshow(img)
            ax.axis("off")
            ax.set_title(title, fontsize=12)
            n = [a["nmax"][0], a["nmax"][-1], b["nmax"][-1]][j - 1]
            ax.text(0.5, 0.01, f"Largest cluster: {n}/32", ha="center", transform=ax.transAxes, fontsize=10)
        axes = [fig.add_subplot(gs[2, j]) for j in [1, 2, 3]]
        t = a["time"] * 1000
        if kind == "aggregation":
            for d, label, color, ls in [(a, "Attraction ON", BLUE, "-"), (b, "Attraction OFF", ORANGE, "--")]:
                axes[0].plot(t, d["nmax"], ls, color=color, label=label, lw=1.5)
                axes[1].plot(t, d["rg"] * 1e6, ls, color=color, lw=1.5)
            axes[0].set(ylabel="Largest cluster size", ylim=(0, 34))
            axes[0].legend(frameon=False, fontsize=9)
            axes[1].set(ylabel="Radius of gyration (µm)")
            axes[2].axis("off")
            axes[2].text(0, 0.87, "ON: 32-bubble cluster\n\nOFF: 32 separate bubbles", fontsize=14, va="top")
            axes[2].text(
                0,
                0.28,
                "PZT remains ON in both cases.\nObservation window: 600 ms.\nNo flow; geometric WALL only.",
                fontsize=10,
                va="top",
            )
            footer = "All 32 actual bubbles • R = 1.78 µm • Local views retain physical relative sizes • EXPLORATORY"
        else:
            centa = a["x"].mean(1)
            centb = b["x"].mean(1)
            delta = np.linalg.norm(centb - centa, axis=1) * 1e6
            axes[0].plot(t, delta, color="#6d568e", lw=1.5)
            axes[0].set(ylabel="Control centroid contrast (µm)")
            for d, label, color, ls in [(a, "Unchanged", BLUE, "-"), (b, "Controlled", ORANGE, "--")]:
                axes[1].plot(t, d["rg"] * 1e6, ls, color=color, label=label, lw=1.4)
            axes[1].set(ylabel="Radius of gyration (µm)")
            axes[1].legend(frameon=False, fontsize=9)
            axes[2].axis("off")
            axes[2].text(0, 0.92, "Translation contrast: 1.34 µm\nRg contrast: 1.59%", fontsize=14, va="top")
            axes[2].text(
                0,
                0.57,
                "No sustained shape response\nbeyond the existing sensitivity\ncomparison was established.\n\nTransient changes and internal\nrearrangements remain visible.",
                fontsize=10,
                va="top",
            )
            footer = "F3 actual trajectories • R = 1.78 µm • Zero flow; no wall forces (global vessel is context) • EXPLORATORY"
        for ax in axes[:2]:
            ax.set_xlabel("Time (ms)")
            ax.grid(alpha=0.15)
            ax.ticklabel_format(axis="y", useOffset=False)
        fig.supxlabel(footer, fontsize=10)
        stem = "Fig_1_attraction_on_off" if kind == "aggregation" else "Fig_2_PZT_translation_shape_limit"
        for ext in ["png", "pdf", "svg"]:
            fig.savefig(OUT / (stem + "." + ext), dpi=300, facecolor="white")
        plt.close(fig)
        with Image.open(OUT / (stem + ".png")) as im:
            im.thumbnail((1600, 1000))
            im.save(OUT / (stem + "_preview.png"))
    write(OUT / "FIGURE_RENDER.json", render_meta)
    print("Figures complete", flush=True)


def movie(kind, paths, data, preview=False):
    a, b = selected(kind, data)
    stride = 5 if kind == "aggregation" else 1
    indices = np.arange(0, len(a["time"]), stride)
    assert len(indices) == 601
    # Both local panels use the exact same camera/bounds built from both histories.
    all_states = np.concatenate([a["x"], b["x"]])
    s1 = local_scene(all_states, kind == "aggregation", False, size=(600, 720))
    s2 = local_scene(all_states, kind == "aggregation", False, size=(600, 720))
    event = js(
        paths["aggregation_on" if kind == "aggregation" else "translation_controlled"] / "CONFIG.json"
    )["events"][0]
    global_img = global_image(event, a["x"][0].mean(0), False, size=(640, 720))
    title = (
        "Attraction ON: aggregation   /   Attraction OFF: no aggregation"
        if kind == "aggregation"
        else "PZT moves the cluster; robust shape control remains unproven"
    )
    headings = (
        ["Attraction ON", "Attraction OFF"] if kind == "aggregation" else ["Unchanged PZT", "Controlled PZT"]
    )
    name = (
        "Video_1_attraction_on_off_60fps"
        if kind == "aggregation"
        else "Video_2_PZT_translation_shape_limit_60fps"
    )

    def frame(index):
        canvas = Image.new("RGB", (1920, 1080), BG)
        d = ImageDraw.Draw(canvas)
        d.text((30, 26), title, font=font(39), fill="#edf3fa")
        d.text(
            (30, 88),
            f"t = {a['time'][index] * 1000:.1f} ms   |   N = 32   |   R = 1.78 µm",
            font=font(27),
            fill="#bac9da",
        )
        d.text((30, 144), "PZT layout", font=font(29), fill="#dbe6f3")
        for x, heading in zip([675, 1300], headings):
            d.text((x, 144), heading, font=font(29), fill="#dbe6f3")
        canvas.paste(global_img, (10, 190))
        canvas.paste(image_state(s1, a["x"][index], kind == "translation"), (655, 190))
        canvas.paste(image_state(s2, b["x"][index], kind == "translation"), (1280, 190))
        d = ImageDraw.Draw(canvas)
        d.text(
            (30, 925),
            (
                "Grey-blue discs: actual PZT layout\nCross: location of the local views"
                if kind == "aggregation"
                else "PZT layout: controlled setting\nF3: no wall forces; vessel is context"
            ),
            font=font(21),
            fill="#aebed0",
        )
        for x, dd in [(675, a), (1300, b)]:
            d.text(
                (x, 925),
                f"Largest cluster: {dd['nmax'][index]}/32\nRg = {dd['rg'][index] * 1e6:.2f} µm",
                font=font(24),
                fill="#dbe6f3",
            )
        if kind == "aggregation":
            text = "Same PZT drive; Primary and WALL remain ON. Only bubble–bubble attraction is disabled."
        else:
            shift = np.linalg.norm(b["x"][index].mean(0) - a["x"][index].mean(0)) * 1e6
            text = f"Control centroid contrast: {shift:.3f} µm. Internal rearrangement remains; stable shape control is not established."
        d.text((30, 1018), text, font=font(23), fill="#c9d7e7")
        d.text((1880, 88), "EXPLORATORY", font=font(22), fill="#aebed0", anchor="ra")
        d.text(
            (990, 88),
            "Local cross: centroid" if kind == "translation" else "Local views: actual relative sizes",
            font=font(23),
            fill="#bac9da",
        )
        return canvas

    if preview:
        for i in [0, len(a["time"]) // 2, len(a["time"]) - 1]:
            frame(i).save(OUT / (name + f"_preview_{i:04d}.png"))
        s1[0].close()
        s2[0].close()
        return
    begin = perf_counter()
    writer = imageio_ffmpeg.write_frames(
        str(OUT / (name + ".mp4")),
        (1920, 1080),
        fps=60,
        codec="libx264",
        pix_fmt_in="rgb24",
        pix_fmt_out="yuv420p",
        macro_block_size=1,
        output_params=["-crf", "19", "-preset", "medium", "-movflags", "+faststart"],
    )
    writer.send(None)
    for f, i in enumerate(indices):
        im = frame(i)
        writer.send(np.asarray(im))
        if f in [0, 300, 600]:
            im.save(OUT / (name + f"_frame_{f:04d}.jpg"), quality=92)
        if f % 120 == 0:
            print(f"{kind}: {f}/601", flush=True)
    writer.close()
    s1[0].close()
    s2[0].close()
    np.savez_compressed(
        OUT / (name + "_FRAME_MAP.npz"),
        video_frame=np.arange(len(indices)),
        source_frame_index=indices,
        physical_time_s=a["time"][indices],
        particle_id=a["id"],
    )
    write(
        OUT / (name + "_IDENTITY.json"),
        dict(
            frames=len(indices),
            fps=60,
            physical_time_s=[float(a["time"][0]), float(a["time"][-1])],
            physical_interval_s=float(np.diff(a["time"][indices])[0]),
            source_names=["aggregation_on", "aggregation_off"]
            if kind == "aggregation"
            else ["translation_baseline", "translation_controlled"],
            local_scene=s1[3],
            position_interpolation=False,
            bubble_radius_scale=1,
            spatial_rescaling=False,
            centroid_markers=kind == "translation",
            locator="Screen-space cross at the projected initial centroid; not a bubble or region boundary",
            source_mesh_sha256=sha(F4 / "inputs/labelled_lumen_m.vtp"),
            render_wall_s=perf_counter() - begin,
            PZT_enabled=event["enabled"],
            PZT_amplitudes_Pa=event["amplitude_Pa"],
            PZT_phases_rad=event["phase_rad"],
        ),
    )
    print(kind + " video complete", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--video", choices=["aggregation", "translation"])
    ap.add_argument("--preview", action="store_true")
    ap.add_argument(
        "--reuse-evidence", action="store_true", help="Reuse prepared evidence for independent render workers"
    )
    args = ap.parse_args()
    paths, data = load()
    if args.reuse_evidence:
        assert (OUT / "EVIDENCE.json").exists()
    else:
        evidence(paths, data)
    if args.video:
        movie(args.video, paths, data, args.preview)
    else:
        figures(paths, data)
    write(
        OUT / ("RENDER_SOURCE_" + (args.video or "figures") + ".json"),
        dict(file=str(Path(__file__)), sha256=sha(__file__), pyvista=pv.__version__, vtk=pv.vtk_version_info),
    )


if __name__ == "__main__":
    main()
