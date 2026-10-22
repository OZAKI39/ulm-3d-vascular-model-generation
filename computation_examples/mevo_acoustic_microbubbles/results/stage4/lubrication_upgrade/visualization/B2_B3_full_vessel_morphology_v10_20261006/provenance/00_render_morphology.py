"""Full-vessel overview of every stored microbubble at its true center and radius.

Scientific geometry is never inflated, interpolated, grouped into a surrogate,
filtered or projected back to the wall. A fixed common camera preserves scale.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
from pathlib import Path
import shutil
import time
import xml.etree.ElementTree as ET

import h5py
import imageio.v2 as imageio
from matplotlib import colors, font_manager
from matplotlib.cm import ScalarMappable
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
from PIL import Image, ImageDraw, ImageFont
import pyvista as pv
from scipy.spatial.distance import pdist

from mevo_mb.coupling.bifurcation import BifurcationGeometry
from scripts.stage4l import render_motion as previous

ROOT, style = previous.ROOT, previous.style
write, sha = previous.write, previous.sha
WIDTH, HEIGHT, FPS, SS = 3840, 2160, 60, 32
SVG = "http://www.w3.org/2000/svg"
ET.register_namespace("", SVG)


def projection(ax, positions_m):
    centers = style.base.screen_points(ax, positions_m * 1000)
    centers[:, 1] = HEIGHT - centers[:, 1]
    origin = positions_m.mean(axis=0) * 1000
    probe = style.base.screen_points(ax, np.vstack([origin, origin + np.eye(3)]))
    matrix = (probe[1:] - probe[0]).T
    matrix[1] *= -1
    singular = np.linalg.svd(matrix, compute_uv=False)
    np.testing.assert_allclose(singular, singular.mean(), rtol=1e-10, atol=1e-9)
    scale = float(singular.mean()) * 1000
    u, v = matrix / singular.mean()
    elev, azim = np.deg2rad([ax.elev, ax.azim])
    eye = np.array([np.cos(elev) * np.cos(azim), np.cos(elev) * np.sin(azim), np.sin(elev)])
    return centers, scale, u, v, eye


def analytic_spheres(background, ax, positions_m, radii_m, ids, supersampling=SS):
    """Analytic orthographic spheres, sparse per-pixel depth and 32x32 coverage.

    No minimum glyph diameter. Fully hidden spheres remain submitted. Sparse
    buffers avoid allocating a supersampled image of the entire vessel.
    """
    centers, scale, u, v, eye = projection(ax, positions_m)
    radius_px = radii_m * scale
    affine = np.stack([u, v]) * scale
    offset = centers[0] - affine @ positions_m[0]
    light = np.array([-0.4, -0.6, 1.0])
    light /= np.linalg.norm(light)
    light_camera = np.array([u @ light, v @ light, eye @ light])
    blue = np.array(colors.to_rgb(style.BLUE)) * 255
    offsets = (np.arange(supersampling) + 0.5) / supersampling
    buffers = {}
    for i, (center, radius) in enumerate(zip(centers, radius_px, strict=True)):
        lower = np.floor(center - radius).astype(int)
        upper = np.ceil(center + radius).astype(int)
        assert np.all(lower >= 0) and np.all(upper <= [WIDTH, HEIGHT]), "Global view crops a bubble"
        for y in range(lower[1], upper[1]):
            for x in range(lower[0], upper[0]):
                nx = (x + offsets[None] - center[0]) / radius
                ny = (y + offsets[:, None] - center[1]) / radius
                disc = nx * nx + ny * ny
                hit = disc <= 1
                if not hit.any():
                    continue
                key = (y, x)
                if key not in buffers:
                    buffers[key] = (
                        np.full(disc.shape, -np.inf),
                        np.full(disc.shape, -1, dtype=np.int32),
                        np.zeros(disc.shape + (3,), dtype=np.float64),
                    )
                depth, owner, rgb = buffers[key]
                nz = np.sqrt(np.maximum(0, 1 - disc))
                z = positions_m[i] @ eye + radii_m[i] * nz
                update = hit & (z > depth)
                shade = 0.40 + 0.60 * np.maximum(
                    0, nx * light_camera[0] + ny * light_camera[1] + nz * light_camera[2]
                )
                depth[update] = z[update]
                owner[update] = int(ids[i])
                rgb[update] = (shade[..., None] * blue)[update]
    result = np.asarray(background).copy()
    covered = set()
    coverage_area = 0.0
    for (y, x), (_, owner, rgb) in buffers.items():
        hit = owner >= 0
        covered.update(np.unique(owner[hit]).tolist())
        coverage = float(hit.mean())
        coverage_area += coverage
        value = rgb.mean(axis=(0, 1)) + result[y, x] * (1 - coverage)
        result[y, x] = np.clip(np.rint(value), 0, 255).astype(np.uint8)
    return Image.fromarray(result), dict(
        submitted_ids=ids.tolist(),
        covered_ids=sorted(covered),
        centers_px=centers.tolist(),
        radii_px=radius_px.tolist(),
        pixels_per_m=scale,
        covered_pixel_area=coverage_area,
        minimum_glyph_size=False,
        projection_matrix_px_per_m=affine.tolist(),
        projection_offset_px=offset.tolist(),
    )


def scene(d, wall):
    mesh_mm = np.asarray(wall.points) * 1000
    bounds = np.array([[-14.7, 11.7], [-7.4, 7.4], [-5.0, 5.0]])
    assert np.all(mesh_mm.min(axis=0) > bounds[:, 0]) and np.all(mesh_mm.max(axis=0) < bounds[:, 1])
    d["view_bounds_mm"] = bounds
    fig = style.base.plt.figure(figsize=(19.2, 10.8), dpi=200, facecolor="black")
    fig.text(0.035, 0.946, "Acoustic field & bubble transport", color="white", fontsize=28)
    fig.text(0.035, 0.902, "2.5 MHz  |  Flow + primary force  |  256 bubbles", color="#bed0e0", fontsize=17)
    fig.text(
        0.035, 0.737, "Normal lubrication  |  Primary only  |  0.30× playback", color="#bed0e0", fontsize=14
    )
    for row, pair in enumerate([d["config"]["pzt"][:2], d["config"]["pzt"][2:]]):
        label = "     ".join(
            f"{e['id'].replace('_plus', '+').replace('_minus', '-').upper()}: "
            f"{e['amplitude_Pa'] / 1000:.2f} kPa / {np.rad2deg(e['phase_rad']):.0f}°"
            for e in pair
        )
        fig.text(0.035, 0.848 - row * 0.030, label, color="#dec09a", fontsize=14)
    ax = fig.add_axes([0.13, 0.15, 0.87, 0.80], projection="3d", facecolor="black", computed_zorder=False)
    style.original.axes_style(ax, bounds, ["X (mm)", "Y (mm)", "Z (mm)"], zoom=1.18)
    ax.set_xticks([-10, -5, 0, 5, 10])
    ax.set_yticks([-5, 0, 5])
    ax.set_zticks([-4, 0, 4])
    for axis in [ax.xaxis, ax.yaxis, ax.zaxis]:
        axis.label.set_fontsize(18)
        axis.labelpad = 9
    ax.tick_params(labelsize=14)
    norm, cmap = style.original.source_and_field(ax, d)
    triangles = wall.faces.reshape(-1, 4)[:, 1:]
    ax.add_collection3d(
        Poly3DCollection(
            mesh_mm[triangles],
            facecolors="#b2c2d5",
            edgecolors="none",
            linewidths=0,
            alpha=0.045,
            antialiased=False,
            zorder=2.5,
        )
    )
    bar_ax = fig.add_axes([0.24, 0.090, 0.48, 0.014])
    bar = fig.colorbar(
        ScalarMappable(norm=norm, cmap=cmap), cax=bar_ax, orientation="horizontal", ticks=[0, 10, 20, 30, 39]
    )
    bar.ax.tick_params(colors="white", labelsize=12, pad=3)
    bar.outline.set_visible(False)
    fig.text(0.74, 0.088, "|P| (kPa)", color="#c2d2e1", fontsize=14)
    # Matplotlib's square 3-D axes patch can crop a projected vessel endpoint
    # despite the point lying inside the 3-D limits. Keep the complete scene;
    # the explicit full-canvas bounds check below guards the actual output.
    for artist in ax.collections:
        artist.set_clip_on(False)
    fig.canvas.draw()
    uv = style.base.screen_points(ax, mesh_mm)
    assert np.all(uv.min(axis=0) > [60, 250]) and np.all(uv.max(axis=0) < [3780, 1870]), uv.min(axis=0)
    background = Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[:, :, :3].copy())
    vector = io.StringIO()
    fig.savefig(vector, format="svg", facecolor="black")
    return (
        fig,
        ax,
        background,
        vector.getvalue(),
        dict(
            bounds_mm=bounds.tolist(),
            camera="Fixed orthographic; elevation 24 deg, azimuth -68 deg",
            camera_moves=False,
            full_wall_vertices=wall.n_points,
            full_wall_triangles=wall.n_cells,
            mesh_screen_bounds_bottom_origin=[uv.min(axis=0).tolist(), uv.max(axis=0).tolist()],
            open_ports=True,
            wall_decimation=False,
        ),
    )


def label_geometry(d, k, ax):
    group = d["clusters"][k][0]
    ix = np.array([d["id_index"][pid] for pid in group["particle_ids"]])
    uv, scale, _, _, _ = projection(ax, d["x"][k, ix])
    center = uv.mean(axis=0)
    radius = float(np.max(np.linalg.norm(uv - center, axis=1) + d["radii"][ix] * scale))
    tx = float(np.clip(center[0] + radius + 48, 100, WIDTH - 650))
    ty = float(np.clip(center[1] - radius - 48, 650, HEIGHT - 480))
    direction = np.array([tx - 12, ty + 16]) - center
    direction /= np.linalg.norm(direction)
    # Preserve the tiny physical cluster silhouette: its label leader starts
    # outside every projected member sphere, with six pixels of clearance.
    leader_start = center + direction * (radius + 6)
    return group, leader_start, (tx, ty)


def annotate(image, d, k, ax):
    draw = ImageDraw.Draw(image)
    font = ImageFont.truetype(str(style.FONT), 66)
    count_font = ImageFont.truetype(str(style.FONT), 48)
    draw.text((134, 448), f"t = {d['t'][k] * 1000:.1f} ms", font=font, fill=style.MINT, anchor="la")
    group, center, (tx, ty) = label_geometry(d, k, ax)
    if group["count"] > 1:
        label = f"Largest: {group['count']} bubbles"
        draw.line([tuple(center), (tx - 12, ty + 16)], fill=style.BLUE, width=2)
        box = draw.textbbox((tx, ty), label, font=count_font, anchor="lm")
        draw.rounded_rectangle((box[0] - 12, box[1] - 8, box[2] + 12, box[3] + 8), radius=8, fill="#080e15")
        draw.text((tx, ty), label, font=count_font, fill="#c1e3ff", anchor="lm")
    return image


def save_vector(path, background_svg, d, k, ax):
    """Resolution-independent full view, with one exact projected circle per ID."""
    root = ET.fromstring(background_svg)
    defs = root.find(f"{{{SVG}}}defs")
    gradient = ET.SubElement(defs, f"{{{SVG}}}radialGradient", id="bubble_light", cx="36%", cy="28%", r="72%")
    for offset, color in [("0%", "#339dff"), ("60%", "#287ccb"), ("100%", "#143f66")]:
        ET.SubElement(gradient, f"{{{SVG}}}stop", offset=offset, attrib={"stop-color": color})
    width = float(root.attrib["viewBox"].split()[2])
    factor = width / WIDTH
    layer = ET.SubElement(root, f"{{{SVG}}}g", id="recorded_bubbles", transform=f"scale({factor:.17g})")
    active = d["states"][k]["active_mask"]
    x, radii, ids = d["x"][k, active], d["radii"][active], d["ids"][active]
    uv, scale, _, _, eye = projection(ax, x)
    for j in np.argsort(x @ eye, kind="stable"):
        ET.SubElement(
            layer,
            f"{{{SVG}}}circle",
            cx=f"{uv[j, 0]:.17g}",
            cy=f"{uv[j, 1]:.17g}",
            r=f"{radii[j] * scale:.17g}",
            fill="url(#bubble_light)",
            attrib={
                "data-particle-id": str(int(ids[j])),
                "data-radius-m": repr(float(radii[j])),
                "data-position-m": json.dumps(x[j].tolist()),
            },
        )
    label = ET.SubElement(
        layer,
        f"{{{SVG}}}text",
        x="134",
        y="506",
        fill=style.MINT,
        attrib={"font-family": "Arial, sans-serif", "font-size": "66"},
    )
    label.text = f"t = {d['t'][k] * 1000:.1f} ms"
    group, center, (tx, ty) = label_geometry(d, k, ax)
    if group["count"] > 1:
        ET.SubElement(
            layer,
            f"{{{SVG}}}line",
            x1=str(center[0]),
            y1=str(center[1]),
            x2=str(tx - 12),
            y2=str(ty + 16),
            stroke=style.BLUE,
            attrib={"stroke-width": "2"},
        )
        ET.SubElement(
            layer,
            f"{{{SVG}}}rect",
            x=str(tx - 12),
            y=str(ty - 36),
            width="590",
            height="70",
            rx="8",
            fill="#080e15",
        )
        label = ET.SubElement(
            layer,
            f"{{{SVG}}}text",
            x=str(tx),
            y=str(ty + 18),
            fill="#c1e3ff",
            attrib={"font-family": "Arial, sans-serif", "font-size": "48"},
        )
        label.text = f"Largest: {group['count']} bubbles"
    description = ET.SubElement(root, f"{{{SVG}}}desc")
    description.text = "Full vessel, fixed global view. Every active recorded bubble, exact center and radius. No cluster glyph or synthetic packing. Vector shading is illustrative; geometry is exact."
    ET.ElementTree(root).write(path, encoding="utf-8", xml_declaration=True)


def audit_native(d, info):
    geometry = BifurcationGeometry(d["config"]["geometry"])
    geometry.bind_mesh_surface(pv.read(d["surface"]))
    active = np.stack([s["active_mask"] for s in d["states"]])
    radii = np.broadcast_to(d["radii"], active.shape)
    gap, _, valid = geometry.wall_gap_and_normal(d["x"][active], radii[active])
    assert valid.all() and gap.min() >= -1.001e-12
    pair_gap = np.array([pdist(x[a]).min() - 2 * d["radii"][0] for x, a in zip(d["x"], active, strict=True)])
    assert pair_gap.min() > -1.01e-11
    with h5py.File(info["trajectory"]) as h:
        for value, name in [(d["x"], "positions_m"), (d["radii"], "radii_m"), (active, "active")]:
            assert np.array_equal(value, h[name][:])
    return dict(
        saved_positions_bitwise_exact=True,
        radii_bitwise_exact=True,
        radius_scale=1.0,
        position_scale=1.0,
        stored_states=len(d["t"]),
        minimum_wall_gap_m=float(gap.min()),
        minimum_pair_gap_m=float(pair_gap.min()),
        all_active_centers_inside=True,
        particle_interpolation=False,
        cluster_reconstruction=False,
        particle_filtering=False,
        wall_symbol_clipping=False,
        camera_moves=False,
        supersampling=SS,
    )


def viewer(out, videos):
    body = """<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Full vessel · Recorded microbubbles</title><style>body{background:#070a0f;color:#eef4ff;font:17px/1.7 Arial,sans-serif;max-width:1500px;margin:30px auto;padding:0 24px}h1,h2{font-weight:400}p{color:#b7cee2}video{width:100%;background:black}a{color:#66e0c2}</style><h1>完整血管全景 · 真实比例微泡聚集</h1><p>固定单一全局视角，保留黑底、蓝色球体、PZT、紫色声场、英文坐标和泡数标注。完整血管包括延伸至 X=−14 mm 的入口段和两个出口。壁面采用原始三角网格，不做减面；端口保持开放。画面不显示 B2/B3。</p><p>每个活动微泡直接来自保存数据，R=1.78 µm；泡心、半径、间距和血管共用同一物理比例。没有合并微泡、重构团簇、筛除粒子、放大球径、移动群体、靠壁裁切或添加聚集圆环。原始微泡在毫米级全景中小于一个像素，因此视频不能分辨全部微米级间隙。另附同一完整全景的 SVG 矢量图，可在图形软件中放大检查逐泡几何，不改变视频取景。</p><p>3840×2160，60 fps，完整 0–6 s 主力场，0.30× 匀速，约 20 s。原始状态每 20 ms 保存一次，每个状态显示四帧，时间文字显示保存时刻。没有插值微泡坐标或混合两种聚集形态；60 fps 编码不代表有 60 fps 的独立计算状态。</p><p>原始数据未改动，没有启动新模拟。原求解器约 10 pm 的球对几何容差保留，视觉平稳与完整解码不证明局部物理精度。视频采用解析球面、逐采样点深度遮挡和 32×32 像素覆盖采样；SVG 的圆形投影半径精确，渐变着色仅为球体光照示意。</p>"""
    for key, value in videos.items():
        body += f'<h2>{key}</h2><video controls playsinline preload="metadata" poster="keyframes/{key}_000600.png" src="{value["path"]}"></video><p>'
        body += " · ".join(
            f'<a href="vectors/{key}_{k:03d}.svg">完整全景矢量图 t={k * 0.02:g} s</a>' for k in [0, 150, 300]
        )
        body += "</p>"
    body += '<p><a href="MORPHOLOGY_CHECKS.json">几何检查</a> · <a href="RENDER_MANIFEST.json">数据来源</a> · <a href="DELIVERY_VERIFICATION.json">交付核验</a></p></html>'
    (out / "OPEN_RESULTS.html").write_text(body)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--case", choices=["B2", "B3", "all"], default="all")
    parser.add_argument("--preview-only", action="store_true")
    args = parser.parse_args()
    out = args.output.resolve()
    assert not out.exists(), "Preserve earlier deliveries"
    for folder in ["keyframes", "vectors", "provenance"]:
        (out / folder).mkdir(parents=True)
    runs, protected = previous.selected_runs()
    sources = [
        Path(__file__),
        Path(previous.__file__),
        Path(style.__file__),
        Path(style.original.__file__),
        Path(style.base.__file__),
        Path(style.playback.__file__),
        ROOT / "scripts/stage3/render_motion.py",
        ROOT / "src/mevo_mb/coupling/bifurcation.py",
        ROOT / "src/mevo_mb/bubbles/metrics.py",
        ROOT / "scripts/stage4l/check_morphology.py",
        ROOT / "tests/test_stage4l_morphology_rendering.py",
    ]
    source_sha = {str(p): sha(p) for p in sources}
    for i, p in enumerate(sources):
        shutil.copy2(p, out / "provenance" / f"{i:02d}_{p.name}")
    protected[str(style.FONT)] = sha(style.FONT)
    font_manager.fontManager.addfont(str(style.FONT))
    style.base.plt.rcParams.update({"font.family": "Arial", "font.weight": "normal", "font.size": 14})
    videos, checks, rows = {}, {}, []
    started = time.perf_counter()
    for case in ["B2", "B3"] if args.case == "all" else [args.case]:
        d, _, info = previous.load(case, "primary", runs["primary"][case], protected)
        wall_path = d["surface"].parent / "mesh-surfaces/WALL.vtp"
        protected[str(wall_path)] = sha(wall_path)
        wall = pv.read(wall_path).triangulate()
        exterior = pv.read(d["surface"])
        assert set(wall["GlobalNodeID"]).issubset(set(exterior["GlobalNodeID"]))
        nodes = {int(pid): i for i, pid in enumerate(exterior["GlobalNodeID"])}
        source_points = exterior.points[[nodes[int(pid)] for pid in wall["GlobalNodeID"]]]
        assert np.array_equal(wall.points, source_points), "Wall differs from the authoritative mesh"
        checks[case] = audit_native(d, info)
        fig, ax, background, vector, screen = scene(d, wall)
        checks[case]["scene"] = screen
        background.save(out / f"{case}_BACKGROUND.png")
        active_all = np.stack([s["active_mask"] for s in d["states"]])
        np.savez_compressed(
            out / f"{case}_RECORDED_SPHERES.npz",
            time_s=d["t"],
            positions_m=d["x"],
            radii_m=d["radii"],
            active=active_all,
            particle_id=d["ids"],
        )
        selected = [0, 40, 384, 385, 400, 600, 608, 609, 800, 1000, 1200]
        frames = selected if args.preview_only else range(1201)
        evidence_by_state = {}
        last_state = None
        path = out / f"stage4L_{case}_full_vessel_true_morphology_4k60.mp4"
        encoder = (
            None
            if args.preview_only
            else imageio.get_writer(
                path,
                fps=FPS,
                codec="libx264",
                pixelformat="yuv420p",
                macro_block_size=None,
                ffmpeg_log_level="error",
                output_params=[
                    "-preset",
                    "medium",
                    "-crf",
                    "12",
                    "-threads",
                    "4",
                    "-profile:v",
                    "high",
                    "-level:v",
                    "5.2",
                    "-g",
                    "120",
                    "-movflags",
                    "+faststart",
                ],
            )
        )
        try:
            for frame in frames:
                k = min(frame // 4, 300)
                if k != last_state:
                    active = active_all[k]
                    im, evidence = analytic_spheres(
                        background, ax, d["x"][k, active], d["radii"][active], d["ids"][active]
                    )
                    assert evidence["submitted_ids"] == d["ids"][active].tolist()
                    group = d["clusters"][k][0]
                    evidence_by_state[str(k)] = dict(
                        source=k, largest_cluster_ids=group["particle_ids"], **evidence
                    )
                    im = annotate(im, d, k, ax)
                    if k in [0, 150, 300]:
                        save_vector(out / "vectors" / f"{case}_{k:03d}.svg", vector, d, k, ax)
                    last_state = k
                rows.append([case, frame, frame / FPS, frame * 0.005, k, float(d["t"][k])])
                if frame in selected:
                    im.save(out / "keyframes" / f"{case}_{frame:06d}.png")
                if encoder is not None:
                    encoder.append_data(np.asarray(im))
                if frame % 240 == 0:
                    print(
                        json.dumps(
                            dict(case=case, frame=frame, total=1201, elapsed_s=time.perf_counter() - started)
                        ),
                        flush=True,
                    )
        finally:
            if encoder is not None:
                encoder.close()
            style.base.plt.close(fig)
        write(out / f"{case}_SPHERE_PROJECTIONS.json", evidence_by_state)
        info["sampling"] = (
            "All 301 actual saved states; four-frame hold; no particle interpolation or image blending"
        )
        videos[case] = dict(
            path=path.name,
            source=info,
            camera="fixed_full_vessel",
            fps=FPS,
            frames=1201,
            resolution=[WIDTH, HEIGHT],
            playback_speed=0.3,
            sha256=sha(path) if path.exists() else None,
        )
    assert all(sha(p) == h for p, h in (protected | source_sha).items())
    write(out / "MORPHOLOGY_CHECKS.json", checks)
    with (out / "frame_map.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            ["case", "frame", "video_time_s", "requested_time_s", "source_index", "source_time_s"]
        )
        writer.writerows(rows)
    write(
        out / "RENDER_MANIFEST.json",
        dict(
            status="PREVIEW" if args.preview_only else "RENDERED_PENDING_DELIVERY_CHECK",
            videos=videos,
            protected_inputs=protected,
            original_source_sha256=source_sha,
            scientific_data_modified=False,
            simulations_launched=False,
            particle_interpolation=False,
            particle_radius_scale=1.0,
            cluster_geometry_transform=False,
            case_identifiers_in_video=False,
            camera_moves=False,
            physical_state_cadence_s=0.02,
            playback_speed=0.3,
            elapsed_s=time.perf_counter() - started,
        ),
    )
    viewer(out, videos)


if __name__ == "__main__":
    main()
