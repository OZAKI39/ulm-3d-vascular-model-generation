"""F5 native scene styles adapted from the existing v6 PZT and F4 renderers."""

from pathlib import Path
from types import SimpleNamespace
import json
import numpy as np
import pyvista as pv
from PIL import Image, ImageDraw, ImageFont
from .common import INPUT

F4 = INPUT / "f4"
R = 1.78e-6
BLUE = "#3288d0"
BG = "#000000"
FONT = Path("/mnt/c/Windows/Fonts/arial.ttf")


def reference():
    data = json.loads((Path(__file__).parent / "assets/PZT_v6_style.json").read_text())
    return data, SimpleNamespace(disc_surface=disc_surface)


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


def vessel():
    mesh = pv.read(F4 / "labelled_lumen_m.vtp")
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
