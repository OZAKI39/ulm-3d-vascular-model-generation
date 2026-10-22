"""Render the actual B2 vessel surface and source apertures with opaque depth testing."""

import argparse
import hashlib
import itertools
import json
from pathlib import Path
import shutil

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import FancyArrowPatch
import numpy as np
from PIL import Image
import pyvista as pv
import yaml

ROOT = Path(__file__).resolve().parents[2]
PRIOR = ROOT / "results/stage4/publication_figures/workflow_B2_B3_20261005_v1"
COLORS = ["#d55e00", "#0072b2", "#d55e00", "#d55e00"]
NAMES = ["Y+", "Y−", "Z+", "Z−"]


def sha(path):
    with Path(path).open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def disc_surface(element):
    center = np.array(element["position_m"]) * 1000
    n = np.array(element["normal"], dtype=float)
    n /= np.linalg.norm(n)
    u = np.cross(n, np.eye(3)[np.argmin(abs(n))])
    u /= np.linalg.norm(u)
    v = np.cross(n, u)
    count = 256
    t = np.arange(count) * 2 * np.pi / count
    radius = element["size_m"][0] * 500
    ring = center + radius * (np.cos(t)[:, None] * u + np.sin(t)[:, None] * v)
    pts = np.vstack([center, ring])
    faces = np.array([[3, 0, i + 1, (i + 1) % count + 1] for i in range(count)])
    mesh = pv.PolyData(pts, faces)
    np.testing.assert_allclose((ring - center) @ n, 0, atol=1e-12)
    np.testing.assert_allclose(np.linalg.norm(ring - center, axis=1), radius, atol=1e-12)
    return mesh


def render(config, vessel, dpi):
    box = np.array([0.025, 0.125, 0.95, 0.81])
    figsize = np.array([180, 118]) / 25.4
    size = np.rint(figsize * dpi * box[2:]).astype(int)
    plotter = pv.Plotter(off_screen=True, window_size=size.tolist(), lighting="three lights")
    plotter.set_background("white")
    plotter.enable_anti_aliasing("ssaa")
    vessel_mm = vessel.copy(deep=True)
    vessel_mm.points = vessel.points.astype(float) * 1000
    elements = config["pzt"]
    discs = [disc_surface(e) for e in elements]
    checks = []
    for i, j in itertools.combinations(range(4), 2):
        _, count = discs[i].collision(discs[j], cell_tolerance=0.0, box_tolerance=0.0)
        assert count == 0, (i, j, count)
        delta = abs(np.array(elements[i]["position_m"]) - np.array(elements[j]["position_m"])) * 1000
        ri, rj = elements[i]["size_m"][0] * 500, elements[j]["size_m"][0] * 500
        ni, nj = np.abs(elements[i]["normal"]), np.abs(elements[j]["normal"])
        if np.dot(ni, nj) > 0.99:
            gap = float(delta @ ni)
        else:
            gap = float(np.hypot(delta @ ni - rj, delta @ nj - ri))
        checks.append(dict(pair=[NAMES[i], NAMES[j]], intersection_cells=count, analytic_aperture_gap_mm=gap))
    for i, disc in enumerate(discs):
        _, count = disc.collision(vessel_mm, cell_tolerance=0.0, box_tolerance=0.0)
        assert count == 0, (NAMES[i], "vessel", count)
    actors = [
        plotter.add_mesh(
            vessel_mm,
            color="#aab9c4",
            opacity=1,
            smooth_shading=True,
            specular=0.13,
            specular_power=24,
            ambient=0.18,
            diffuse=0.78,
        )
    ]
    for disc, color in zip(discs, COLORS):
        actors.append(
            plotter.add_mesh(
                disc, color=color, opacity=1, smooth_shading=False, ambient=0.45, diffuse=0.55, specular=0.08
            )
        )
    # Equal physical scale, orthographic view, camera-only fitting of the full geometry.
    eye = np.array([12.0, -55.0, 43.0])
    focal = np.array([-1.0, 0.0, 0.0])
    direction = focal - eye
    direction /= np.linalg.norm(direction)
    right = np.cross(direction, [0, 0, 1])
    right /= np.linalg.norm(right)
    up = np.cross(right, direction)
    points = np.concatenate([vessel_mm.points, *[disc.points for disc in discs]])
    horiz, vert = (points - focal) @ right, (points - focal) @ up
    shift = (horiz.min() + horiz.max()) / 2 * right + (vert.min() + vert.max()) / 2 * up
    eye += shift
    focal += shift
    scale = max(np.ptp(vert) / (2 * 0.76), np.ptp(horiz) / (2 * 0.91 * size[0] / size[1]))
    plotter.camera_position = [eye.tolist(), focal.tolist(), up.tolist()]
    plotter.enable_parallel_projection()
    plotter.camera.parallel_scale = float(scale)
    plotter.show(auto_close=False)
    # Visibility is determined by the same per-pixel depth buffer. Opaque ID pass
    # must be invariant to reversed actor insertion order, unlike painter sorting.
    plotter.disable_anti_aliasing()
    id_colors = [(128, 128, 128), (230, 30, 30), (30, 230, 30), (30, 30, 230), (230, 30, 230)]
    props = []
    for actor, color in zip(actors, id_colors):
        prop = actor.GetProperty()
        props.append((prop.GetColor(), prop.GetLighting()))
        prop.SetColor(*(np.array(color) / 255))
        prop.LightingOff()
    plotter.render()
    mask_image = plotter.screenshot(return_img=True)
    for actor in actors:
        plotter.renderer.RemoveActor(actor)
    for actor in reversed(actors):
        plotter.renderer.AddActor(actor)
    plotter.render()
    reversed_image = plotter.screenshot(return_img=True)
    np.testing.assert_array_equal(mask_image, reversed_image)
    anchors = []
    visible_pixels = {}
    for name, color in zip(NAMES, id_colors[1:]):
        yy, xx = np.where(np.all(mask_image[:, :, :3] == color, axis=-1))
        assert len(xx) > 100
        visible_pixels[name] = int(len(xx))
        # A point within the visible mask, not a possibly hidden world-space center.
        k = np.argmin((xx - xx.mean()) ** 2 + (yy - yy.mean()) ** 2)
        anchors.append([box[0] + box[2] * xx[k] / size[0], box[1] + box[3] * (1 - yy[k] / size[1])])
    for actor, (color, lighting) in zip(actors, props):
        actor.GetProperty().SetColor(color)
        actor.GetProperty().SetLighting(lighting)
    plotter.enable_anti_aliasing("ssaa")
    plotter.render()
    image = plotter.screenshot(return_img=True)
    plotter.close()
    return image, dict(
        box=box.tolist(),
        figsize=figsize.tolist(),
        size=size.tolist(),
        anchors=anchors,
        camera_eye_mm=eye.tolist(),
        camera_focal_mm=focal.tolist(),
        camera_up=up.tolist(),
        camera_right=right.tolist(),
        parallel_scale_mm=float(scale),
        source_pair_checks=checks,
        source_vessel_intersections=0,
        insertion_order_invariant=True,
        visible_source_pixels=visible_pixels,
        full_vessel_cells=vessel.n_cells,
        full_vessel_points=vessel.n_points,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--preview", action="store_true")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    provenance = json.loads((PRIOR / "FIGURE_PROVENANCE.json").read_text())
    config_path = ROOT / provenance["cases"]["B2"]["run"] / "config.yaml"
    surface_path = ROOT / next(p for p in provenance["protected_inputs"] if p.endswith(".vtp"))
    sources = {str(p.relative_to(ROOT)): sha(p) for p in [config_path, surface_path]}
    assert all(digest == provenance["protected_inputs"][rel] for rel, digest in sources.items())
    config = yaml.safe_load(config_path.read_text())
    vessel = pv.read(surface_path)
    dpi = 180 if args.preview else 600
    image, audit = render(config, vessel, dpi)
    fonts = Path("/mnt/c/Windows/Fonts")
    for path in fonts.glob("arial*.ttf"):
        font_manager.fontManager.addfont(str(path))
    plt.rcParams.update(
        {
            "font.family": "Arial" if (fonts / "arial.ttf").is_file() else "DejaVu Sans",
            "font.size": 10,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
            "text.color": "#20262d",
        }
    )
    fig = plt.figure(figsize=audit["figsize"], dpi=dpi, facecolor="white")
    ax = fig.add_axes(audit["box"])
    ax.imshow(image, interpolation="none")
    ax.axis("off")
    fig.text(0.04, 0.963, "a  PZT control setup", fontsize=12, fontweight="bold", va="top")
    fig.text(0.96, 0.963, "B2", fontsize=10, ha="right", va="top", color="#667581")
    # Leader labels stay outside source faces; control values are intentionally omitted.
    offsets = [(42, 24), (-54, -7), (-42, 38), (35, -32)]
    for name, xy, offset, color in zip(NAMES, audit["anchors"], offsets, COLORS):
        ax.annotate(
            name,
            xy=xy,
            xycoords=fig.transFigure,
            xytext=offset,
            textcoords="offset points",
            ha="left" if offset[0] > 0 else "right",
            va="center",
            color=color,
            fontsize=11,
            fontweight="bold",
            arrowprops=dict(arrowstyle="-", color=color, lw=0.75),
            annotation_clip=False,
        )
    # A compact orientation triad replaces a surrounding frame. Its screen-space
    # directions are the actual orthographic camera projections of the SI axes.
    origin = np.array([0.10, 0.12])
    right, up = np.array(audit["camera_right"]), np.array(audit["camera_up"])
    aspect = audit["figsize"][0] / audit["figsize"][1]
    for i, name in enumerate("xyz"):
        delta = 0.045 * np.array([right[i], up[i] * aspect])
        end = origin + delta
        fig.add_artist(
            FancyArrowPatch(
                origin,
                end,
                transform=fig.transFigure,
                arrowstyle="-|>",
                mutation_scale=8,
                lw=0.8,
                color="#536675",
            )
        )
        textpos = origin + 1.3 * delta
        fig.text(*textpos, name, fontsize=9, ha="center", va="center", color="#536675")
    box = audit["box"]
    bar_width = box[2] * 2 / (2 * audit["parallel_scale_mm"] * audit["size"][0] / audit["size"][1])
    start = 0.90 - bar_width
    fig.add_artist(
        plt.Line2D([start, 0.90], [0.09, 0.09], transform=fig.transFigure, color="#20262d", lw=1.5)
    )
    for x in [start, 0.90]:
        fig.add_artist(plt.Line2D([x, x], [0.084, 0.096], transform=fig.transFigure, color="#20262d", lw=0.8))
    fig.text((start + 0.90) / 2, 0.107, "2 mm", ha="center", fontsize=9)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    w, h = fig.canvas.get_width_height()
    for artist in fig.findobj(matplotlib.text.Text):
        if artist.get_visible() and artist.get_text():
            bb = artist.get_window_extent(renderer)
            assert bb.x0 >= 0 and bb.x1 <= w + 1 and bb.y0 >= 0 and bb.y1 <= h + 1, artist.get_text()
    stem = "Fig1a_PZT_control"
    png = args.output / (stem + ("_preview.png" if args.preview else ".png"))
    fig.savefig(png, dpi=dpi)
    if not args.preview:
        fig.savefig(args.output / (stem + ".pdf"), dpi=600, metadata={"Title": "B2 PZT control geometry"})
        fig.savefig(args.output / (stem + ".svg"), dpi=600)
        with Image.open(png) as im:
            im.convert("RGB").save(args.output / (stem + ".tif"), compression="tiff_lzw", dpi=(600, 600))
            im.convert("RGB").resize(tuple(round(v * 0.3) for v in im.size), Image.Resampling.LANCZOS).save(
                args.output / (stem + "_preview.png"), dpi=(180, 180)
            )
    plt.close(fig)
    assert all(sha(ROOT / rel) == digest for rel, digest in sources.items())
    report = dict(
        status="PREVIEW" if args.preview else "PASS",
        protected_inputs=sources,
        renderer_sha256=sha(__file__),
        render=audit,
        source_elements=config["pzt"],
        geometry="Full original vessel surface; opaque zero-thickness acoustic apertures",
        export="600 dpi raster 3D scene; editable vector labels in PDF/SVG",
        constraints="No position changes, mesh decimation, geometry smoothing, transparency, or new simulation",
    )
    (args.output / "PROVENANCE.json").write_text(json.dumps(report, indent=2) + "\n")
    if not args.preview:
        shutil.copy2(__file__, args.output / "paper_pzt_surface.py")
        (args.output / "CAPTION.txt").write_text(
            "B2 PZT control geometry. The full, original vessel-domain surface is shown with four opaque ideal "
            "acoustic source apertures. Source positions, normals and diameters match the sealed B2 configuration. "
            "The aperture surfaces have zero thickness; no unmeasured transducer housing is added. "
            "Blue denotes the 0-degree source and orange the 180-degree sources. Per-pixel depth testing "
            "determines all surface occlusion. Opaque object visibility was verified to be independent of "
            "actor insertion order; all four apertures remain partly visible. No source pairs or source-vessel "
            "surfaces intersect. Neighbouring perpendicular source apertures have an analytic minimum gap "
            "of sqrt(2) mm. Orthographic projection preserves physical proportions; the scale bar is 2 mm. "
            "Smooth normals affect illumination only; the original vessel vertices and connectivity are unchanged. "
            "Control amplitudes, phases, frequency and enabled states are provided separately, outside the figure. "
            "The PNG/TIFF masters are 600 dpi. PDF/SVG contain the rendered 3D raster scene and vector annotations.\n"
        )
    print(
        json.dumps(
            dict(
                status=report["status"],
                output=str(png),
                geometry_checks=audit["source_pair_checks"],
                visible_pixels=audit["visible_source_pixels"],
            )
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
