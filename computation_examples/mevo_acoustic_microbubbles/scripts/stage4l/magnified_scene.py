"""Global movie context with the publication v6 PZT geometry and materials.

PZT recipe transcribed from the frozen paper_pzt_reference.py, SHA256
207c650fcc2a7d2714d42d6f7cb6aab9d4d26bcb0e29c66d88f91803fcd101ab.
The backing thickness is illustrative; acoustic front faces remain exact.
"""

from __future__ import annotations

import itertools
import numpy as np
import pyvista as pv
from PIL import Image, ImageDraw, ImageFont
from matplotlib import colormaps

from scripts.stage4l.render_motion import style

WIDTH, HEIGHT = 3840, 2160
BOX = np.array([65, 140, 1790, 835]) * 2
BOUNDS = np.array([[-14.5, 11.5], [-7.5, 7.5], [-2.2, 2.2]])


def font(size):
    return ImageFont.truetype(str(style.FONT), size * 2)


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
    return pv.PolyData(pts, faces)


def add_sources(plotter, elements, vessel):
    actors, bodies, materials = [], [], []
    discs = [disc_surface(e) for e in elements]
    for i, j in itertools.combinations(range(4), 2):
        _, count = discs[i].collision(discs[j], cell_tolerance=0.0, box_tolerance=0.0)
        assert count == 0
    for index, (element, disc) in enumerate(zip(elements, discs, strict=True)):
        n = np.asarray(element["normal"], dtype=float)
        prominent = element["id"] in ["z_plus", "y_minus"]
        front_opacity, back_opacity, side_opacity = (0.42, 0.55, 0.65) if prominent else (0.27, 0.30, 0.45)
        materials.append(
            dict(
                element=element["id"],
                front_opacity=front_opacity,
                back_opacity=back_opacity,
                side_opacity=side_opacity,
            )
        )
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
        bodies.append(body)
        actors.append(
            plotter.add_mesh(
                disc,
                color="#8aa7b9" if index == 1 else "#a2b0b8",
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
    for i, body in enumerate(bodies):
        _, count = body.collision(vessel, cell_tolerance=0.0, box_tolerance=0.0)
        assert count == 0
        for other in bodies[:i]:
            _, count = body.collision(other, cell_tolerance=0.0, box_tolerance=0.0)
            assert count == 0
    return actors, bodies, materials


def lines(segments):
    points = np.asarray(segments).reshape(-1, 3)
    return pv.PolyData(points, lines=np.array([[2, 2 * i, 2 * i + 1] for i in range(len(segments))]))


def make_projector(audit):
    right, up, eye = (np.asarray(audit[k]) for k in ["right", "up", "eye_direction"])
    focal = np.asarray(audit["focal_mm"])
    box = np.asarray(audit["viewport_px"])
    factor = box[3] / (2 * audit["parallel_scale_mm"])

    def project(positions_m):
        delta = np.asarray(positions_m) * 1000 - focal
        centers = np.column_stack(
            [box[0] + box[2] / 2 + delta @ right * factor, box[1] + box[3] / 2 - delta @ up * factor]
        )
        return centers, factor * 1000, right, -up, eye

    return project


def scene(d, wall, scale):
    plotter = pv.Plotter(off_screen=True, window_size=BOX[2:].tolist(), lighting="three lights")
    plotter.set_background("black")
    plotter.render_window.SetMultiSamples(0)
    plotter.render_window.SetAlphaBitPlanes(1)
    plotter.enable_depth_peeling(number_of_peels=100, occlusion_ratio=0.0)
    vessel = wall.copy(deep=True)
    vessel.points = vessel.points.astype(float) * 1000
    actors = [
        plotter.add_mesh(
            vessel,
            color="#b2c2d5",
            opacity=0.08,
            smooth_shading=True,
            specular=0.24,
            specular_power=40,
            ambient=0.40,
            diffuse=0.60,
        )
    ]
    ports = vessel.extract_feature_edges(
        boundary_edges=True, non_manifold_edges=False, feature_edges=False, manifold_edges=False
    )
    actors.append(plotter.add_mesh(ports, color="#d8e2ed", opacity=0.8, lighting=False, line_width=1.5))
    sections = style.base.wall_sections(d["surface"], bounds=BOUNDS)
    actors.append(
        plotter.add_mesh(lines(sections), color="#d8e2ed", opacity=0.52, lighting=False, line_width=1.2)
    )
    source_actors, bodies, materials = add_sources(plotter, d["config"]["pzt"], vessel)
    actors.extend(source_actors)
    x, y, z = d["axes"]
    iz, ix = np.argmin(abs(z)), np.argmin(abs(x))
    X, Y = np.meshgrid(x, y, indexing="ij")
    Y2, Z = np.meshgrid(y, z, indexing="ij")
    for coordinates, values, opacity in [
        ((X, Y, np.zeros_like(X)), abs(d["P"][:, :, iz]) / 1000, 0.67),
        ((np.zeros_like(Y2), Y2, Z), abs(d["P"][ix]) / 1000, 0.46),
    ]:
        field = pv.StructuredGrid(*coordinates)
        field["pressure_kPa"] = values.ravel(order="F")
        actors.append(
            plotter.add_mesh(
                field,
                scalars="pressure_kPa",
                cmap="magma",
                clim=[0, 39],
                opacity=opacity,
                lighting=False,
                show_scalar_bar=False,
                interpolate_before_map=True,
            )
        )
    xmin, xmax, ymin, ymax, zmin, zmax = BOUNDS.ravel()
    segments = []
    for x in [-14.5, -10, -5, 0, 5, 10, 11.5]:
        segments.extend([[[x, ymin, zmin], [x, ymax, zmin]], [[x, ymax, zmin], [x, ymax, zmax]]])
    for y in [-7.5, -5, 0, 5, 7.5]:
        segments.extend([[[xmin, y, zmin], [xmax, y, zmin]], [[xmin, y, zmin], [xmin, y, zmax]]])
    for z in [-2, 0, 2]:
        segments.extend([[[xmin, ymax, z], [xmax, ymax, z]], [[xmin, ymin, z], [xmin, ymax, z]]])
    grid = lines(segments)
    actors.append(plotter.add_mesh(grid, color="#344557", opacity=0.68, lighting=False, line_width=1.0))
    edges = np.array(
        [
            [[xmin, ymin, zmin], [xmax, ymin, zmin]],
            [[xmax, ymin, zmin], [xmax, ymax, zmin]],
            [[xmax, ymax, zmin], [xmax, ymax, zmax]],
        ]
    )
    actors.append(plotter.add_mesh(lines(edges), color="#75899d", lighting=False, line_width=1.5))
    az, el = np.deg2rad([-68.0, 24.0])
    eye_direction = np.array([np.cos(el) * np.cos(az), np.cos(el) * np.sin(az), np.sin(el)])
    right = np.cross(-eye_direction, [0, 0, 1])
    right /= np.linalg.norm(right)
    up = np.cross(right, -eye_direction)
    allpoints = np.concatenate([vessel.points, *[b.points for b in bodies], grid.points])
    horizontal, vertical = allpoints @ right, allpoints @ up
    focal = (horizontal.max() + horizontal.min()) / 2 * right + (vertical.max() + vertical.min()) / 2 * up
    parallel_scale = max(np.ptp(vertical) / 2, np.ptp(horizontal) / 2 * BOX[3] / BOX[2]) * 1.06
    eye = focal + eye_direction * 60
    plotter.camera_position = [eye.tolist(), focal.tolist(), up.tolist()]
    plotter.enable_parallel_projection()
    plotter.camera.parallel_scale = float(parallel_scale)
    lights = []
    for e in d["config"]["pzt"]:
        n = -np.asarray(e["normal"], float)
        center = np.asarray(e["position_m"]) * 1000 + n * 0.35
        view = eye - center
        view /= np.linalg.norm(view)
        if n @ view > 0:
            position = center + 18 * (2 * (n @ view) * n - view)
            plotter.add_light(
                pv.Light(
                    position=position,
                    focal_point=center,
                    color="#eef4ff",
                    intensity=0.65,
                    positional=True,
                    cone_angle=180,
                    attenuation_values=(1, 0, 0),
                )
            )
            lights.append(dict(element=e["id"], position_mm=position.tolist()))
    plotter.show(auto_close=False)
    rendered = plotter.screenshot(return_img=True)[:, :, :3]
    assert plotter.renderer.GetLastRenderingUsedDepthPeeling()
    for actor in actors:
        plotter.renderer.RemoveActor(actor)
    for actor in reversed(actors):
        plotter.renderer.AddActor(actor)
    plotter.render()
    reverse = plotter.screenshot(return_img=True)[:, :, :3]
    delta = abs(rendered.astype(np.int16) - reverse.astype(np.int16))
    assert delta.max() <= 2, (delta.max(), delta.mean())
    plotter.close()
    audit = dict(
        right=right.tolist(),
        up=up.tolist(),
        eye_direction=eye_direction.tolist(),
        focal_mm=focal.tolist(),
        parallel_scale_mm=float(parallel_scale),
        viewport_px=BOX.tolist(),
        bounds_mm=BOUNDS.tolist(),
        camera_moves=False,
        camera_azimuth_deg=-68,
        camera_elevation_deg=24,
        full_wall_points=wall.n_points,
        full_wall_triangles=wall.n_cells,
        wall_decimation=False,
        open_boundary_edges=ports.n_cells,
        depth_peeling_used=True,
        reverse_actor_max_rgb_delta=int(delta.max()),
        pzt_body_intersections=0,
        pzt_backing_display_thickness_mm=0.35,
        source_materials=materials,
        outward_reflection_lights=lights,
        back_specular=0.88,
        back_specular_power=68,
        glyph_layer="Analytic spheres composited over translucent context; sphere-sphere depth tested",
    )
    project = make_projector(audit)
    uv = project(wall.points)[0]
    assert np.all(uv.min(axis=0) > [100, BOX[1]]) and np.all(uv.max(axis=0) < [WIDTH - 100, BOX[1] + BOX[3]])
    audit["wall_projected_bounds_px"] = [uv.min(axis=0).tolist(), uv.max(axis=0).tolist()]
    image = Image.new("RGB", (WIDTH, HEIGHT), "black")
    image.paste(Image.fromarray(rendered), (int(BOX[0]), int(BOX[1])))
    draw = ImageDraw.Draw(image)

    def text(x, y, label, size=16, color="#c5d3e1", anchor="la"):
        draw.text((x * 2, y * 2), label, font=font(size), fill=color, anchor=anchor)

    text(67, 35, "Acoustic field & bubble transport", 30, "white")
    text(67, 82, "2.5 MHz  |  Flow + primary force  |  256 bubbles", 19)
    for row, pair in enumerate([d["config"]["pzt"][:2], d["config"]["pzt"][2:]]):
        label = "     ".join(
            f"{e['id'].replace('_plus', '+').replace('_minus', '-').upper()}: "
            f"{e['amplitude_Pa'] / 1000:.2f} kPa / {np.rad2deg(e['phase_rad']):.0f}°"
            for e in pair
        )
        text(67, 125 + row * 30, label, 16, "#b0c1ce")
    text(67, 244, f"Morphology ×{scale:g} (display)  |  0.30× playback", 17, "#bed0e0")
    for values, func, offset in [
        ([-10, -5, 0, 5, 10], lambda v: [v, ymin, zmin], [0, 13]),
        ([-5, 0, 5], lambda v: [xmax, v, zmin], [20, 8]),
        ([-2, 0, 2], lambda v: [xmax, ymax, v], [22, 0]),
    ]:
        for value in values:
            p = project(np.array([func(value)]) / 1000)[0][0] / 2 + offset
            text(*p, f"{value:g}", 13, anchor="mm")
    for label, p, offset in [
        ("X (mm)", [0, ymin, zmin], [0, 47]),
        ("Y (mm)", [xmax, 0, zmin], [52, 34]),
        ("Z (mm)", [xmax, ymax, 0], [74, 0]),
    ]:
        xy = project(np.array([p]) / 1000)[0][0] / 2 + offset
        text(*xy, label, 17, anchor="mm")
    for e in d["config"]["pzt"]:
        n = np.asarray(e["normal"])
        point = np.asarray(e["position_m"]) - n * 0.0006
        xy = project(point[None])[0][0] / 2
        text(*xy, e["id"].replace("_plus", "+").replace("_minus", "-").upper(), 16, "#d0dbe2", anchor="mm")
    grad = (colormaps["magma"](np.linspace(0, 1, 1600))[:, :3] * 255).astype("uint8")
    image.paste(Image.fromarray(np.tile(grad[None], (22, 1, 1))), (1080, 2000))
    for val in [0, 10, 20, 30, 39]:
        text(540 + 800 * val / 39, 1028, str(val), 13, anchor="mm")
    text(1370, 998, "| P | (kPa)", 17)
    return image, project, audit
