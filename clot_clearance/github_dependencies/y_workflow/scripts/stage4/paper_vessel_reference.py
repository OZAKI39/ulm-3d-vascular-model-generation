"""Render the exact v6 vessel/camera/materials without any PZT display geometry."""

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
from PIL import Image
import pyvista as pv

ROOT = Path(__file__).resolve().parents[2]
REFERENCE = ROOT / "results/stage4/publication_figures/PZT_control_reference_v6"


def sha(path):
    with Path(path).open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def load_reference():
    info = json.loads((REFERENCE / "PROVENANCE.json").read_text())
    script = REFERENCE / "paper_pzt_reference.py"
    assert sha(script) == info["renderer_sha256"]
    spec = importlib.util.spec_from_file_location("archived_pzt_v6", script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return info, module


def render(vessel, info, dpi):
    audit = dict(info["render"])
    audit["size"] = (
        np.rint(np.array(audit["figsize"]) * dpi * np.array(audit["box"][2:])).astype(int).tolist()
    )
    p = pv.Plotter(off_screen=True, window_size=audit["size"], lighting="three lights")
    p.set_background("white")
    p.render_window.SetMultiSamples(0)
    p.render_window.SetAlphaBitPlanes(1)
    p.enable_depth_peeling(number_of_peels=100, occlusion_ratio=0)
    wall = vessel.copy(deep=True)
    wall.points = vessel.points.astype(float) * 1000
    actors = [
        p.add_mesh(
            wall,
            color="#b85c64",
            opacity=0.88,
            smooth_shading=True,
            specular=0.24,
            specular_power=40,
            ambient=0.20,
            diffuse=0.72,
            backface_params=dict(
                color="#44212b", opacity=0.96, ambient=0.12, diffuse=0.45, specular=0.03, specular_power=30
            ),
        )
    ]
    rims = wall.extract_feature_edges(
        boundary_edges=True, non_manifold_edges=False, feature_edges=False, manifold_edges=False
    )
    actors.append(
        p.add_mesh(rims, color="#c57981", opacity=1, lighting=False, line_width=max(0.9, dpi / 180 * 0.9))
    )
    xmin, xmax, ymin, ymax, zmin, zmax = -14.5, 11.5, -7.5, 7.5, -5.0, 5.0
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
        np.asarray(segments).reshape(-1, 3),
        lines=np.array([[2, 2 * i, 2 * i + 1] for i in range(len(segments))]),
    )
    actors.append(p.add_mesh(grid, color="#bec9d1", line_width=max(0.6, dpi / 180 * 0.55)))
    edges = [
        [[xmin, ymin, zmin], [xmax, ymin, zmin]],
        [[xmax, ymin, zmin], [xmax, ymax, zmin]],
        [[xmax, ymax, zmin], [xmax, ymax, zmax]],
    ]
    for edge in edges:
        actors.append(p.add_mesh(pv.Line(*edge), color="#647887", line_width=max(0.8, dpi / 180 * 0.7)))
    p.camera_position = [audit["camera_eye_mm"], audit["camera_focal_mm"], audit["camera_up"]]
    p.enable_parallel_projection()
    p.camera.parallel_scale = audit["parallel_scale_mm"]
    # The existing two broad fill/reflection lights also illuminate the vessel.
    # Preserve these exact scene lights, without creating any PZT actor.
    elements = {e["id"]: e for e in info["source_elements"]}
    for light in audit["outward_reflection_lights"]:
        e = elements[light["element"]]
        center = np.array(e["position_m"]) * 1000 - np.array(e["normal"]) * 0.35
        p.add_light(
            pv.Light(
                position=light["position_mm"],
                focal_point=center,
                color="#eef4ff",
                intensity=0.65,
                positional=True,
                cone_angle=180,
                attenuation_values=(1, 0, 0),
            )
        )
    p.show(auto_close=False)
    image = p.screenshot(return_img=True)
    assert p.renderer.GetLastRenderingUsedDepthPeeling()
    for a in actors:
        p.renderer.RemoveActor(a)
    for a in reversed(actors):
        p.renderer.AddActor(a)
    p.render()
    reverse = p.screenshot(return_img=True)
    delta = abs(image.astype(np.int16) - reverse.astype(np.int16))
    assert delta.max() <= 2
    p.close()
    return (
        image,
        audit,
        dict(
            rendered_PZT_actor_count=0,
            total_mesh_actors=len(actors),
            reverse_actor_max_rgb_delta=int(delta.max()),
            depth_peeling_used=True,
        ),
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--preview", action="store_true")
    args = ap.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    info, reference = load_reference()
    surface = ROOT / next(
        path for path in info["protected_inputs"] if path.endswith("mesh-complete.exterior.vtp")
    )
    wall, paths, wall_audit = reference.load_open_wall(surface)
    protected = {str(p.relative_to(ROOT)): sha(p) for p in [surface, *paths]}
    assert all(info["protected_inputs"][key] == digest for key, digest in protected.items())
    original_png_sha = sha(REFERENCE / "Fig1a_PZT_control.png")
    dpi = 180 if args.preview else 600
    image, audit, render_check = render(wall, info, dpi)
    fonts = Path("/mnt/c/Windows/Fonts")
    for font in fonts.glob("arial*.ttf"):
        font_manager.fontManager.addfont(str(font))
    plt.rcParams.update(
        {
            "font.family": "Arial" if (fonts / "arial.ttf").exists() else "DejaVu Sans",
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
    # Reuse the exact world-registered v6 axis/tick layout. Empty sources create no leaders.
    reference.annotations(fig, ax, {"pzt": []}, audit)
    fig.texts[0].set_text("Y-shaped vessel")
    fig.texts[1].remove()
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    w, h = fig.canvas.get_width_height()
    for artist in fig.findobj(matplotlib.text.Text):
        if artist.get_visible() and artist.get_text():
            bb = artist.get_window_extent(renderer)
            assert bb.x0 >= 0 and bb.x1 <= w + 1 and bb.y0 >= 0 and bb.y1 <= h + 1, artist.get_text()
    stem = "Fig1a_vessel_only"
    png = args.output / (stem + ("_preview.png" if args.preview else ".png"))
    fig.savefig(png, dpi=dpi)
    if not args.preview:
        fig.savefig(
            args.output / (stem + ".pdf"),
            dpi=dpi,
            metadata={"Title": "Y-shaped vessel, original v6 style without PZT"},
        )
        fig.savefig(args.output / (stem + ".svg"), dpi=dpi)
        with Image.open(png) as im:
            im.convert("RGB").save(args.output / (stem + ".tif"), compression="tiff_lzw", dpi=(dpi, dpi))
            im.convert("RGB").resize(tuple(round(v * 0.3) for v in im.size), Image.Resampling.LANCZOS).save(
                args.output / (stem + "_preview.png"), dpi=(180, 180)
            )
    plt.close(fig)
    assert all(sha(ROOT / key) == value for key, value in protected.items())
    assert sha(REFERENCE / "Fig1a_PZT_control.png") == original_png_sha
    report = dict(
        status="PREVIEW" if args.preview else "PASS",
        reference=str(REFERENCE),
        reference_png_sha256=original_png_sha,
        reference_renderer_sha256=info["renderer_sha256"],
        renderer_sha256=sha(__file__),
        protected_inputs=protected,
        camera_and_style_reference=audit,
        render_check=render_check,
        open_wall_verification=wall_audit,
        changes=[
            "All PZT discs, backing, rims, labels and leaders omitted",
            "Title: Y-shaped vessel; B2 tag omitted",
        ],
        preserved="Exact v6 WALL coordinates/connectivity, three open ends, red material/opacity, inner-wall shading, scene lights, camera, physical scale, grid and axis typography",
        dpi=dpi,
        simulation_run=False,
        geometry_source="Reference v6 SimVascular WALL patch; not substituted with a different vessel mesh",
    )
    (args.output / "PROVENANCE.json").write_text(json.dumps(report, indent=2) + "\n")
    if not args.preview:
        shutil.copy2(__file__, args.output / Path(__file__).name)
        (args.output / "README.txt").write_text(
            "Y-shaped vessel without PZT, preserving the PZT_control_reference_v6 camera and visual style.\n"
            "Original WALL.vtp only; all three ends remain open. No new simulation and no geometric changes.\n"
            "PNG/TIFF: 600 dpi; PDF/SVG: embedded 3D scene and editable English labels.\n"
            "Reproduce from project root with the shared scientific Python environment:\n"
            "PYTHONPATH=src:. python -m scripts.stage4.paper_vessel_reference --output results/stage4/publication_figures/vessel_only_reference_v1\n"
        )
        files = {
            p.name: dict(bytes=p.stat().st_size, sha256=sha(p))
            for p in sorted(args.output.iterdir())
            if p.is_file() and p.name not in ["DELIVERY_SHA256.json", "render.log"]
        }
        (args.output / "DELIVERY_SHA256.json").write_text(json.dumps(files, indent=2) + "\n")
    print(json.dumps(dict(status=report["status"], output=str(png), **render_check)), flush=True)


if __name__ == "__main__":
    main()
