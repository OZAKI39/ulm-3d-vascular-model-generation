"""White-background, reference-style B2 PZT scene with depth-correct transparency."""

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
import numpy as np
from PIL import Image
import pyvista as pv
import yaml

ROOT = Path(__file__).resolve().parents[2]
PRIOR = ROOT / "results/stage4/publication_figures/workflow_B2_B3_20261005_v1"
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
    box = np.array([0.025, 0.115, 0.94, 0.80])
    figsize = np.array([180, 132]) / 25.4
    size = np.rint(figsize * dpi * box[2:]).astype(int)
    plotter = pv.Plotter(off_screen=True, window_size=size.tolist(), lighting="three lights")
    plotter.set_background("white")
    # Multisampling must be disabled for reliable transparent depth peeling.
    plotter.render_window.SetMultiSamples(0)
    plotter.render_window.SetAlphaBitPlanes(1)
    plotter.enable_depth_peeling(number_of_peels=100, occlusion_ratio=0.0)
    vessel_mm = vessel.copy(deep=True)
    vessel_mm.points = vessel.points.astype(float) * 1000
    elements = config["pzt"]
    discs = [disc_surface(e) for e in elements]
    geometry_checks = []
    for i, j in itertools.combinations(range(4), 2):
        _, count = discs[i].collision(discs[j], cell_tolerance=0.0, box_tolerance=0.0)
        assert count == 0
        geometry_checks.append(dict(pair=[NAMES[i], NAMES[j]], intersection_cells=count))
    actors = [
        plotter.add_mesh(
            vessel_mm,
            color="#b85c64",
            opacity=0.88,
            smooth_shading=True,
            specular=0.24,
            specular_power=40,
            ambient=0.20,
            diffuse=0.72,
        )
    ]
    display_bodies = []
    source_materials = []
    for index, (element, disc) in enumerate(zip(elements, discs)):
        n = np.asarray(element["normal"], dtype=float)
        prominent = element["id"] in ["z_plus", "y_minus"]
        front_opacity = 0.42 if prominent else 0.27
        back_opacity = 0.55 if prominent else 0.30
        side_opacity = 0.65 if prominent else 0.45
        source_materials.append(
            dict(
                element=element["id"],
                front_opacity=front_opacity,
                back_opacity=back_opacity,
                side_opacity=side_opacity,
            )
        )
        color = "#8aa7b9" if index == 1 else "#a2b0b8"
        # The acoustic front face remains exactly on its configured plane. A
        # 0.35 mm back extrusion is display-only, not an identified PZT dimension.
        back = disc.copy(deep=True)
        back.points = back.points - n * 0.35
        faces = back.faces.reshape(-1, 4).copy()
        faces[:, [2, 3]] = faces[:, [3, 2]]
        back.faces = faces
        rim = disc.points[1:]
        back_rim = rim - n * 0.35
        count = len(rim)
        sides = np.array([[4, i, i + count, (i + 1) % count + count, (i + 1) % count] for i in range(count)])
        side = pv.PolyData(np.vstack([rim, back_rim]), sides)
        body = (
            disc.merge([back, side], merge_points=True)
            .extract_surface(algorithm="dataset_surface")
            .triangulate()
        )
        display_bodies.append(body)
        actors.append(
            plotter.add_mesh(
                disc,
                color=color,
                opacity=front_opacity,
                smooth_shading=True,
                ambient=0.55,
                diffuse=0.45,
                specular=0.30,
                specular_power=50,
            )
        )
        actors.append(
            plotter.add_mesh(
                back,
                color="#b4c2cb",
                opacity=back_opacity,
                smooth_shading=True,
                ambient=0.32,
                diffuse=0.50,
                specular=0.88,
                specular_power=68,
            )
        )
        actors.append(
            plotter.add_mesh(
                side,
                color="#8499a5",
                opacity=side_opacity,
                smooth_shading=True,
                specular=0.30,
                specular_power=30,
            )
        )
        for ring in [rim, back_rim]:
            curve = pv.lines_from_points(np.vstack([ring, ring[:1]]), close=False)
            actors.append(
                plotter.add_mesh(
                    curve.tube(radius=0.018, n_sides=8),
                    color="#7c909e",
                    opacity=0.75,
                    smooth_shading=True,
                    specular=0.55,
                    specular_power=50,
                )
            )
    for i, body in enumerate(display_bodies):
        _, count = body.collision(vessel_mm, cell_tolerance=0.0, box_tolerance=0.0)
        assert count == 0, (NAMES[i], "vessel", count)
        for j in range(i):
            _, count = body.collision(display_bodies[j], cell_tolerance=0.0, box_tolerance=0.0)
            assert count == 0, (NAMES[i], NAMES[j], count)
    # Three faint coordinate planes echo the reference video on a white page.
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
    points = np.asarray(segments).reshape(-1, 3)
    grid = pv.PolyData(points, lines=np.array([[2, 2 * i, 2 * i + 1] for i in range(len(segments))]))
    actors.append(plotter.add_mesh(grid, color="#e0e6eb", line_width=max(0.6, dpi / 180 * 0.55)))
    edges = np.array(
        [
            [[xmin, ymin, zmin], [xmax, ymin, zmin]],
            [[xmax, ymin, zmin], [xmax, ymax, zmin]],
            [[xmax, ymax, zmin], [xmax, ymax, zmax]],
        ]
    )
    for edge in edges:
        actors.append(plotter.add_mesh(pv.Line(*edge), color="#a9b7c2", line_width=max(0.8, dpi / 180 * 0.7)))
    azim, elev = -48.0, 26.0
    az, el = np.deg2rad([azim, elev])
    direction = -np.array([np.cos(el) * np.cos(az), np.cos(el) * np.sin(az), np.sin(el)])
    right = np.cross(direction, [0, 0, 1])
    right /= np.linalg.norm(right)
    up = np.cross(right, direction)
    focal = np.zeros(3)
    allpoints = np.concatenate(
        [vessel_mm.points, *[d.points for d in display_bodies], points, edges.reshape(-1, 3)]
    )
    horizontal, vertical = allpoints @ right, allpoints @ up
    focal += (horizontal.max() + horizontal.min()) / 2 * right + (vertical.max() + vertical.min()) / 2 * up
    scale = max(np.ptp(vertical) / (2 * 0.83), np.ptp(horizontal) / (2 * 0.88 * size[0] / size[1]))
    scale /= 1.17
    eye = focal - direction * 60
    plotter.camera_position = [eye.tolist(), focal.tolist(), up.tolist()]
    plotter.enable_parallel_projection()
    plotter.camera.parallel_scale = float(scale)
    outward_lights = []
    for element in elements:
        outside_normal = -np.asarray(element["normal"], dtype=float)
        center = np.asarray(element["position_m"]) * 1000 + outside_normal * 0.35
        to_eye = eye - center
        to_eye /= np.linalg.norm(to_eye)
        if np.dot(outside_normal, to_eye) > 0:
            reflection_direction = 2 * np.dot(outside_normal, to_eye) * outside_normal - to_eye
            light_position = center + reflection_direction * 18
            light = pv.Light(
                position=light_position,
                focal_point=center,
                color="#eef4ff",
                intensity=0.65,
                positional=True,
                cone_angle=180,
                attenuation_values=(1, 0, 0),
            )
            plotter.add_light(light)
            outward_lights.append(dict(element=element["id"], position_mm=light_position.tolist()))
    plotter.show(auto_close=False)
    image = plotter.screenshot(return_img=True)
    assert plotter.renderer.GetLastRenderingUsedDepthPeeling()
    # Read the actual composite, not only a renderer capability flag: transparent
    # layers must remain unchanged when actor insertion order is reversed.
    assert np.count_nonzero(np.min(image[:, :, :3], axis=-1) < 180) > 1000
    for actor in actors:
        plotter.renderer.RemoveActor(actor)
    for actor in reversed(actors):
        plotter.renderer.AddActor(actor)
    plotter.render()
    reversed_image = plotter.screenshot(return_img=True)
    delta = abs(image.astype(np.int16) - reversed_image.astype(np.int16))
    if delta.max() > 2:
        Image.fromarray(image).save("/tmp/pzt_v3_order_first.png")
        Image.fromarray(reversed_image).save("/tmp/pzt_v3_order_reverse.png")
    assert delta.max() <= 2, (delta.max(), delta.mean())
    plotter.close()
    return image, dict(
        box=box.tolist(),
        figsize=figsize.tolist(),
        size=size.tolist(),
        camera_eye_mm=eye.tolist(),
        camera_focal_mm=focal.tolist(),
        camera_up=up.tolist(),
        camera_right=right.tolist(),
        parallel_scale_mm=float(scale),
        camera_azimuth_deg=azim,
        camera_elevation_deg=elev,
        source_pair_checks=geometry_checks,
        display_body_intersections=0,
        depth_peeling_used=True,
        reverse_actor_max_rgb_delta=int(delta.max()),
        reverse_actor_mean_rgb_delta=float(delta.mean()),
        full_vessel_cells=vessel.n_cells,
        full_vessel_points=vessel.n_points,
        source_materials=source_materials,
        vessel_color="#b85c64",
        vessel_opacity=0.88,
        outward_reflection_lights=outward_lights,
        back_specular=0.88,
        back_specular_power=68,
        display_only_backing_thickness_mm=0.35,
        display_only_rim_radius_mm=0.018,
    )


def project(point, audit):
    delta = np.asarray(point) - np.array(audit["camera_focal_mm"])
    scale = audit["parallel_scale_mm"]
    aspect = audit["size"][0] / audit["size"][1]
    u = 0.5 + delta @ np.array(audit["camera_right"]) / (2 * scale * aspect)
    v = 0.5 + delta @ np.array(audit["camera_up"]) / (2 * scale)
    x, y, w, h = audit["box"]
    return np.array([x + w * u, y + h * v])


def annotations(fig, ax, config, audit):
    fig.text(0.045, 0.965, "PZT control setup", fontsize=12, va="top", fontweight="normal")
    fig.text(0.955, 0.965, "B2", fontsize=10, ha="right", va="top", color="#667581")
    offsets = [(40, 23), (-44, -7), (-28, 37), (-12, -38)]
    for e, name, offset in zip(config["pzt"], NAMES, offsets):
        xy = project(np.array(e["position_m"]) * 1000, audit)
        ax.annotate(
            name,
            xy=xy,
            xycoords=fig.transFigure,
            xytext=offset,
            textcoords="offset points",
            fontsize=10,
            color="#4f6574",
            fontweight="normal",
            ha="left" if offset[0] > 0 else "right",
            va="center",
            arrowprops=dict(arrowstyle="-", color="#738997", lw=0.65),
            annotation_clip=False,
        )
    # World-registered, editable English axis labels and sparse ticks.
    entries = []
    for x in [-10, 0, 10]:
        entries.append(([x, -7.5, -5], str(x), (0, -12), "center"))
    for y in [-4, 0, 4]:
        entries.append(([11.5, y, -5], str(y), (10, -6), "left"))
    for z in [-4, 0, 4]:
        entries.append(([11.5, 7.5, z], str(z), (9, 0), "left"))
    for p, text, off, ha in entries:
        ax.annotate(
            text,
            xy=project(p, audit),
            xycoords=fig.transFigure,
            xytext=off,
            textcoords="offset points",
            fontsize=8,
            color="#637582",
            ha=ha,
            va="center",
            annotation_clip=False,
        )
    for p, text, off in [
        ([-1.5, -7.5, -5], "X (mm)", (0, -27)),
        ([11.5, 0, -5], "Y (mm)", (31, -19)),
        ([11.5, 7.5, 0], "Z (mm)", (35, 0)),
    ]:
        ax.annotate(
            text,
            xy=project(p, audit),
            xycoords=fig.transFigure,
            xytext=off,
            textcoords="offset points",
            fontsize=9,
            color="#4f6574",
            ha="center",
            va="center",
            annotation_clip=False,
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
    annotations(fig, ax, config, audit)
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
        geometry="Full original vessel surface and configured acoustic front faces; display-only thin backing and rims",
        export="600 dpi raster 3D scene; editable vector labels in PDF/SVG",
        constraints="No changes to acoustic planes, aperture diameters, vessel vertices or connectivity; no new simulation",
        reference_video_sha256=sha(
            "/home/lzy/projects/clot_clearance/visualization/streaming_fragmentation_demo/damage_source_v3/particle_damage_load_source_60fps.mp4"
        ),
        reference_style="Arial regular; orthographic scene; fine three-dimensional coordinate grid, adapted to white background",
    )
    (args.output / "PROVENANCE.json").write_text(json.dumps(report, indent=2) + "\n")
    if not args.preview:
        shutil.copy2(__file__, args.output / "paper_pzt_reference.py")
        (args.output / "CAPTION.txt").write_text(
            "B2 PZT control geometry. Full original vessel-domain surface and four configured acoustic apertures, "
            "shown in a white-background orthographic view with a fine coordinate grid and regular Arial labels. "
            "Low-saturation grey-blue source materials and a muted red vessel are display choices. Vessel opacity "
            "is 0.88. Y-minus and Z-plus front/back/side opacities are 0.42/0.55/0.65; Y-plus and Z-minus "
            "retain 0.27/0.30/0.45. Overlapping surfaces compound opacity. Colors do not encode field values. "
            "A 0.35 mm backing and 0.018 mm edge radius are illustrative display geometry only; no measured "
            "PZT thickness, housing design or electromechanical model is implied. The authoritative acoustic "
            "front planes, centers, normals and 6 mm diameters are unchanged. No acoustic apertures or displayed "
            "thin bodies intersect one another or the vessel. Transparency uses depth peeling; reverse actor "
            "insertion checks the actual composited image. Vessel lighting uses smooth normals without moving "
            "mesh vertices or changing connectivity. Camera azimuth is -48 degrees and elevation 26 degrees. "
            "Outward-facing backs use a soft specular material with two broad point-light highlights. "
            "Control values remain outside the figure. PDF/SVG contain a 600 dpi rendered 3D scene with "
            "editable vector text. The stylistic reference is particle_damage_load_source_60fps.mp4; its "
            "simulation data are not used in this figure. No simulation was run for this display.\n"
        )
    print(
        json.dumps(
            dict(
                status=report["status"],
                output=str(png),
                geometry_checks=audit["source_pair_checks"],
                reverse_actor_max_rgb_delta=audit["reverse_actor_max_rgb_delta"],
            )
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
