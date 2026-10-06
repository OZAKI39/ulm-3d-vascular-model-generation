"""One enlarged PZT/field scene, with explicit ring or spherical cluster displays."""

import argparse
import csv
import json
from pathlib import Path
import shutil
import time

from matplotlib import colors, font_manager
from matplotlib.cm import ScalarMappable
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
from PIL import Image, ImageDraw, ImageFont

from scripts.stage4 import render_aggregation as original
from scripts.stage4 import render_linked_zoom as playback
from mevo_mb.bubbles.metrics import aggregation_metrics

base, ROOT = original.base, original.ROOT
BLUE, MINT = original.BLUE, original.MINT
FPS, FONT = playback.FPS, playback.FONT
AREA_PER_MEMBER_PT2 = 18.0
UNIFORM_STATE_FRAMES = np.arange(301) * 4


def uniform_timeline():
    """Constant 0.30x playback: each video frame advances physical time by 5 ms."""
    for frame in range(1201):
        left = frame // 4
        yield left, min(left + 1, 300), (frame % 4) / 4, frame * 0.005


def prepare_clusters(d):
    """Connected components from the frozen scientific threshold, active particles only."""
    d["clusters"] = []
    for k, state in enumerate(d["states"]):
        indices = np.flatnonzero(state["active_mask"])
        stats, labels, _, _ = aggregation_metrics(
            d["x"][k, indices], d["radii"][indices], d["config"]["metrics"]
        )
        assert stats["cluster_count"] == state["stats"]["cluster_count"]
        groups = []
        for label in np.unique(labels):
            members = indices[labels == label]
            groups.append(
                dict(
                    particle_ids=d["ids"][members].tolist(),
                    count=len(members),
                    center_m=d["x"][k, members].mean(axis=0).tolist(),
                    marker_area_pt2=AREA_PER_MEMBER_PT2 * len(members),
                )
            )
        groups.sort(key=lambda g: (-g["count"], min(g["particle_ids"])))
        assert sum(g["count"] for g in groups) == state["active_particles"]
        assert groups[0]["count"] == stats["largest_cluster_size"]
        d["clusters"].append(groups)


def scene(name, d, sphere_mode=False, magnification=1.0):
    fig = base.plt.figure(figsize=(19.2, 10.8), dpi=100, facecolor="black")
    fig.text(0.035, 0.946, f"{name}  |  Acoustic field & bubble aggregation", color="white", fontsize=28)
    fig.text(
        0.035,
        0.902,
        "2.5 MHz  |  Flow + primary + secondary force  |  256 bubbles",
        color="#bed0e0",
        fontsize=17,
    )
    # Source controls occupy the free upper-left margin, without a second scene.
    elements = d["config"]["pzt"]
    for row, pair in enumerate([elements[:2], elements[2:]]):
        label = "     ".join(
            f"{e['id'].replace('_plus', '+').replace('_minus', '-').upper()}: "
            f"{e['amplitude_Pa'] / 1000:.2f} kPa / {np.rad2deg(e['phase_rad']):.0f}°"
            for e in pair
        )
        fig.text(0.035, 0.848 - row * 0.030, label, color="#dec09a", fontsize=14)
    ax = fig.add_axes([0.13, 0.15, 0.87, 0.80], projection="3d", facecolor="black", computed_zorder=False)
    original.axes_style(ax, original.FIELD_BOUNDS, ["X (mm)", "Y (mm)", "Z (mm)"], zoom=1.31)
    ax.set_xticks([0, 5, 10])
    ax.set_yticks([-5, 0, 5])
    ax.set_zticks([-4, 0, 4])
    for axis in [ax.xaxis, ax.yaxis, ax.zaxis]:
        axis.label.set_fontsize(18)
        axis.labelpad = 9
    ax.tick_params(labelsize=14)
    norm, cmap = original.source_and_field(ax, d)
    rings = (
        None
        if sphere_mode
        else ax.scatter(
            [],
            [],
            [],
            facecolors=(0.08, 0.32, 0.65, 0.18),
            edgecolors=(0.24, 0.66, 1.0, 0.85),
            linewidths=1.0,
            depthshade=False,
            zorder=8,
        )
    )
    exits = ax.scatter(
        [],
        [],
        [],
        facecolors="#64756f" if sphere_mode else "none",
        edgecolors="#74948d",
        linewidths=0.5,
        s=5 if sphere_mode else 12,
        alpha=0.5,
        depthshade=False,
        zorder=9,
    )
    if sphere_mode:
        points = Poly3DCollection([], edgecolors="none", antialiased=False, zsort="average", zorder=10)
        ax.add_collection3d(points)
        template, intensity = original.sphere_faces()
        face_rgb = intensity[:, None] * np.array(colors.to_rgb(BLUE))
    else:
        points = ax.scatter(
            [], [], [], c=BLUE, edgecolors="#c1e3ff", linewidths=0.23, s=22, depthshade=False, zorder=10
        )
    bar_ax = fig.add_axes([0.24, 0.090, 0.48, 0.014])
    bar = fig.colorbar(
        ScalarMappable(norm=norm, cmap=cmap), cax=bar_ax, orientation="horizontal", ticks=[0, 10, 20, 30, 39]
    )
    bar.ax.tick_params(colors="white", labelsize=12, pad=3)
    bar.outline.set_visible(False)
    fig.text(0.74, 0.088, "|P| (kPa)", color="#c2d2e1", fontsize=14)
    dynamic = [exits, points] if sphere_mode else [rings, exits, points]
    for artist in dynamic:
        artist.set_visible(False)
    fig.canvas.draw()
    checks = dict(three_dimensional_axes=1, zoom_panel=False, roi_connectors=False)
    checks["blue_rings"] = not sphere_mode
    checks["right_explanatory_text"] = False
    checks["clock_position"] = "upper_left"
    boxes = {}
    for axis in [ax.xaxis, ax.yaxis, ax.zaxis]:
        box = axis.label.get_window_extent(fig.canvas.get_renderer())
        assert box.x0 > 8 and box.x1 < 1912 and box.y0 > 100 and box.y1 < 850, box
        boxes[axis.label.get_text()] = list(box.bounds)
    checks["axis_label_bounds"] = boxes
    uv = base.screen_points(ax, np.vstack([np.zeros(3), np.eye(3)]))
    singular = np.linalg.svd((uv[1:] - uv[0]).T, compute_uv=False)
    assert np.isclose(singular[0], singular[1], rtol=1e-10)
    checks["pixels_per_mm"] = float(singular.mean())
    checks["scale_relative_to_v6"] = float(singular.mean() / 45.620670734403305)
    assert checks["scale_relative_to_v6"] > 1.1
    if sphere_mode:
        checks["sphere_diameter_pixels"] = float(0.00356 * magnification * singular.mean())
    checks["scale_relative_to_previous_left_view"] = float(singular.mean() / 35.67770403587953)
    assert checks["scale_relative_to_previous_left_view"] > 1.2
    # Check all source rims as well as every original bubble coordinate against safe margins.
    source_points = []
    for e in elements:
        normal = np.array(e["normal"], dtype=float)
        u = np.cross(normal, np.eye(3)[int(np.argmin(abs(normal)))])
        u /= np.linalg.norm(u)
        v = np.cross(normal, u)
        angles = np.linspace(0, 2 * np.pi, 97)
        rim = np.array(e["position_m"]) * 1000 + e["size_m"][0] * 500 * (
            np.cos(angles)[:, None] * u + np.sin(angles)[:, None] * v
        )
        source_points.extend(rim)
    geometry = np.vstack([d["x"].reshape(-1, 3) * 1000, source_points])
    projected = base.screen_points(ax, geometry)
    lo, hi = projected.min(axis=0), projected.max(axis=0)
    assert np.all(lo > [45, 160]) and np.all(hi < [1875, 920]), (lo, hi)
    checks["source_and_particle_screen_bounds_bottom_origin"] = [lo.tolist(), hi.tolist()]
    background = fig.canvas.copy_from_bbox(fig.bbox)
    for artist in dynamic:
        artist.set_visible(True)
    for k, groups in enumerate(d["clusters"]):
        for group in groups:
            center = base.screen_points(ax, np.asarray(group["center_m"])[None] * 1000)[0]
            group["screen_center_top_origin"] = [float(center[0]), float(1080 - center[1])]
            if sphere_mode:
                indices = [d["id_index"][pid] for pid in group["particle_ids"]]
                center_m = np.asarray(group["center_m"])
                real = d["x"][k, indices]
                display = center_m + magnification * (real - center_m)
                radii = magnification * d["radii"][indices]
                # Uniform magnification preserves packing, distance/radius ratios, and centroid.
                np.testing.assert_allclose(display.mean(axis=0), center_m, atol=1e-15, rtol=0)
                np.testing.assert_allclose(
                    (display[:, None] - display[None, :]) / magnification,
                    real[:, None] - real[None, :],
                    atol=1e-15,
                    rtol=0,
                )
                group["display_positions_m"] = display.tolist()
                group["display_radii_m"] = radii.tolist()
                uv = base.screen_points(ax, display * 1000)
                radius_px = radii * 1000 * singular.mean()
                lower = (uv - radius_px[:, None]).min(axis=0)
                upper = (uv + radius_px[:, None]).max(axis=0)
                assert np.all(lower > [40, 160]) and np.all(upper < [1880, 920]), (lower, upper)
                group["label_radius_px"] = float(np.max(np.linalg.norm(uv - center, axis=1) + radius_px))
                group.pop("marker_area_pt2")
        if groups[0]["count"] > 1:
            center = np.array(groups[0]["screen_center_top_origin"])
            radius_px = (
                groups[0]["label_radius_px"]
                if sphere_mode
                else np.sqrt(groups[0]["marker_area_pt2"]) * 100 / 72 / 2
            )
            assert center[0] - radius_px > 40 and center[0] + radius_px < 1880

    def render(k):
        fig.canvas.restore_region(background)
        state = d["states"][k]
        x = d["x"][k] * 1000
        active = state["active_mask"]
        exits._offsets3d = tuple(x[~active].T)
        if sphere_mode:
            meshes = []
            for group in d["clusters"][k]:
                for center, radius in zip(
                    group["display_positions_m"], group["display_radii_m"], strict=True
                ):
                    meshes.append((template * radius + np.asarray(center)) * 1000)
            assert len(meshes) == int(active.sum())
            points.set_verts(np.concatenate(meshes))
            points.set_facecolor(np.tile(np.c_[face_rgb, np.ones(len(face_rgb))], (len(meshes), 1)))
        else:
            points._offsets3d = tuple(x[active].T)
            groups = [g for g in d["clusters"][k] if g["count"] > 1]
            centers = np.array([g["center_m"] for g in groups]).reshape(-1, 3) * 1000
            rings._offsets3d = tuple(centers.T)
            rings.set_sizes(np.array([g["marker_area_pt2"] for g in groups]))
        for artist in dynamic:
            artist.do_3d_projection()
            ax.draw_artist(artist)
        return Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[:, :, :3].copy())

    return fig, render, checks


def overlay(im, d, k, physical_time):
    draw = ImageDraw.Draw(im)
    font = ImageFont.truetype(str(FONT), 33)
    count_font = ImageFont.truetype(str(FONT), 24)
    draw.text((67, 224), f"t = {physical_time * 1000:.1f} ms", font=font, fill=MINT, anchor="la")
    # The label refers to the most recent saved state, never an interpolated group count.
    group = d["clusters"][k][0]
    if group["count"] > 1:
        x, y = group["screen_center_top_origin"]
        radius = (
            group["label_radius_px"]
            if "label_radius_px" in group
            else np.sqrt(group["marker_area_pt2"]) * 100 / 72 / 2
        )
        tx, ty = x + radius + 24, y - radius - 24
        draw.line([(x + radius * 0.6, y - radius * 0.6), (tx - 6, ty + 12)], fill=BLUE, width=2)
        label = (
            f"Largest: {group['count']} bubbles"
            if "label_radius_px" in group
            else f"{group['count']} bubbles"
        )
        box = draw.textbbox((tx, ty), label, font=count_font, anchor="lm")
        draw.rounded_rectangle((box[0] - 6, box[1] - 4, box[2] + 6, box[3] + 4), radius=4, fill="#080e15")
        draw.text((tx, ty), label, font=count_font, fill="#c1e3ff", anchor="lm")
    return im


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--preview-only", action="store_true")
    parser.add_argument("--sphere-clusters", action="store_true")
    parser.add_argument("--structure-magnification", type=float, default=1.0)
    parser.add_argument("--constant-speed", action="store_true")
    args = parser.parse_args()
    if args.sphere_clusters:
        args.constant_speed = True
    assert np.isfinite(args.structure_magnification) and args.structure_magnification >= 1
    sample_timeline = list(uniform_timeline()) if args.constant_speed else playback.timeline()
    exact_frames = UNIFORM_STATE_FRAMES if args.constant_speed else playback.STATE_FRAMES
    frame_count = len(sample_timeline)
    out = args.output.resolve()
    assert not list(out.glob("*.mp4")), "Use a new directory to preserve earlier exports"
    (out / "keyframes").mkdir(parents=True, exist_ok=True)
    (out / "provenance").mkdir(exist_ok=True)
    font_manager.fontManager.addfont(str(FONT))
    base.plt.rcParams.update({"font.family": "Arial", "font.weight": "normal", "font.size": 14})
    data, protected = original.load(range(301), half_width_um=25)
    protected[str(FONT)] = base.sha(FONT)
    sources = [
        Path(__file__),
        Path(original.__file__),
        Path(playback.__file__),
        Path(base.__file__),
        ROOT / "scripts/stage3/render_motion.py",
        ROOT / "src/mevo_mb/bubbles/metrics.py",
    ]
    source_hashes = {str(p): base.sha(p) for p in sources}
    for path in sources:
        shutil.copy2(path, out / "provenance" / f"{path.parent.name}_{path.name}")
    videos, checks, rows = {}, {}, []
    started = time.perf_counter()
    snapshots = [0, 2, 10, 70, 72, 100, 300]
    for name, d in data.items():
        prepare_clusters(d)
        fig, render, checks[name] = scene(name, d, args.sphere_clusters, args.structure_magnification)
        if args.preview_only:
            for k in snapshots:
                overlay(render(k), d, k, d["t"][k]).save(out / "keyframes" / f"{name}_{k:03d}.png")
            base.plt.close(fig)
            continue
        variant = "spheres_uniform" if args.sphere_clusters and args.constant_speed else "single_overview"
        video = out / f"stage4_{name}_{variant}_60fps.mp4"
        left_image, right_image, previous = render(0), render(1), 0
        with playback.writer(video) as encoder:
            for frame, (left, right, alpha, t) in enumerate(sample_timeline):
                if left != previous:
                    left_image, right_image = right_image, render(right) if right > left else right_image
                    previous = left
                im = Image.blend(left_image, right_image, alpha) if alpha else left_image.copy()
                im = overlay(im, d, left, t)
                encoder.append_data(np.asarray(im))
                if alpha == 0 and left in snapshots:
                    im.save(out / "keyframes" / f"{name}_{left:03d}.png")
                if frame in ([1, 2, 3, 7, 9, 41, 489] if args.constant_speed else [12, 36, 252, 489]):
                    im.save(out / "keyframes" / f"{name}_transition_{frame:04d}.png")
                rows.append([name, frame, frame / FPS, t, left, right, alpha])
                if frame % 240 == 0:
                    print(
                        json.dumps(
                            dict(
                                case=name,
                                frame=frame,
                                total=frame_count,
                                elapsed_s=time.perf_counter() - started,
                            )
                        ),
                        flush=True,
                    )
        base.plt.close(fig)
        videos[name] = dict(
            path=video.name,
            fps=FPS,
            frames=frame_count,
            duration_s=frame_count / FPS,
            resolution=[1920, 1080],
            sha256=base.sha(video),
        )
    (out / "SCREEN_CHECKS.json").write_text(json.dumps(checks, indent=2) + "\n")
    if args.preview_only:
        return
    with (out / "frame_map.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "case",
                "frame",
                "video_time_s",
                "display_time_s",
                "source_left",
                "source_right",
                "image_blend_fraction",
            ]
        )
        writer.writerows(rows)
    assert all(base.sha(path) == digest for path, digest in (protected | source_hashes).items())
    manifest = dict(
        status="RENDERED_PENDING_DELIVERY_CHECK",
        layout="single_overview_spheres" if args.sphere_clusters else "single_overview",
        videos=videos,
        protected_inputs=protected,
        renderer_sources=source_hashes,
        source_runs={name: d["result_id"] for name, d in data.items()},
        temporal_display=dict(
            original_states=301,
            original_spacing_s=0.02,
            exact_state_video_frames=exact_frames.tolist(),
            method="Cross-dissolve saved raster states, no particle-coordinate interpolation",
            speed="Constant 0.30x across 0–6 s" if args.constant_speed else "0–0.2 s: 0.05x; 0.2–6 s: 0.30x",
            constant_speed=args.constant_speed,
            physical_step_per_video_frame_s=0.005 if args.constant_speed else None,
            count_label="Largest group of the most recent saved state; not an interpolated count",
        ),
        aggregation_display=dict(
            center="Arithmetic mean of actual active group member positions",
            threshold_m=1e-5,
            area_pt2_per_member=AREA_PER_MEMBER_PT2,
            scope="All active groups shown at true centers; ring area encodes count, not physical diameter",
            particle_symbols="Enlarged for visibility; every active particle remains at its actual coordinate",
            radii="Original R0=1.78 um unchanged; no fusion or radial dynamics",
        ),
        field=dict(
            quantity="Steady |P|",
            range_kPa=[0, 39],
            frequency_Hz=2.5e6,
            slices="Native XY z=0 and YZ x=0",
            source_geometry="Four actual 6 mm PZT discs",
            inputs={n: str(d["field_path"].relative_to(ROOT)) for n, d in data.items()},
        ),
        states={
            name: [
                {key: val for key, val in s.items() if key not in ["local", "roi", "group", "active_mask"]}
                | dict(
                    visible_particle_ids=d["ids"][s["active_mask"]].tolist(),
                    cluster_particle_ids=d["ids"][s["group"]].tolist(),
                    overview_groups=d["clusters"][k],
                )
                for k, s in enumerate(d["states"])
            ]
            for name, d in data.items()
        },
        science_scope="Sealed Stage 4I B2/B3 only. Existing near-contact, near-wall and finite-Re limits remain.",
        screen_checks=checks,
        elapsed_s=time.perf_counter() - started,
    )
    if args.sphere_clusters:
        manifest["aggregation_display"] = dict(
            representation="One blue shaded sphere per active physical particle, no ring glyphs or synthetic packing",
            structure_magnification=args.structure_magnification,
            transform="X_display = cluster_centroid + scale * (X_real - cluster_centroid); R_display = scale * R_real",
            true_geometry="Actual saved group membership, centroid, orientation and every distance/radius ratio preserved",
            scope="In-place magnified cluster geometry is a display glyph, not millimetre-scale physical extent. Overlap between distinct magnified groups is display occlusion, not simulated contact.",
            threshold_m=1e-5,
            radius_m=1.78e-6,
        )
    (out / "RENDER_MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n")
    html = """<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>B2/B3 单视图微泡聚集</title><style>body{background:#070a0f;color:#eef4ff;font:17px/1.7 Arial,sans-serif;max-width:1500px;margin:30px auto;padding:0 24px}h1,h2{font-weight:400}p{color:#b7cee2}video{width:100%;background:black}a{color:#66e0c2}</style><h1>B2/B3：声场中的微泡聚集</h1><p>每段只保留一张三维全景，删除局部放大窗口、定位框和连接线。原全景按相同相机方向等比例放大。四个 PZT 声源、声压幅值截面和微泡位置均取自正式结果。</p><p>蓝点表示实际位置的微泡；蓝色聚集环位于泡群真实质心，面积正比于群内泡数。数字标注最大群体的微泡数量。粒子符号和聚集环是可见性标记，不表示物理半径，也不代表融合；未移动粒子或夸大群内距离。聚集依据仍为原中心距离 ≤10 µm 的连通分组，出流微泡不再计入活动泡群。</p><p>使用全部 301 个原始保存状态，60 fps、23.35 秒。帧间采用画面渐变；没有用位置插值重建未保存的近接触轨迹。数字表示最近的原始保存状态。原始 R0=1.78 µm、轨迹、近壁/近接触及有限雷诺数限制保持。</p>"""
    if args.sphere_clusters:
        html = (
            html[: html.index("<h1>")]
            + f"""<h1>B2/B3：球形微泡聚集</h1><p>保留单一三维场景、PZT 声源和稳态声场；蓝色圆环已删除。每个活动微泡对应一个有立体明暗的圆球，最大泡群保留泡数标注。</p><p>真实直径 3.56 µm 在当前全景仅约 0.16 像素。因此球半径和群内相对坐标共同放大 {args.structure_magnification:g} 倍，泡群质心仍在真实位置；逐泡身份、群内排布、方向和距离/半径比均来自保存数据。不同放大泡群的画面重叠不是物理接触，显示尺寸不用于读取实际毫米距离。</p><p>全程统一 0.30×，60 fps，每帧物理时间推进 5 ms。模拟 0–6 s 对应视频约 20 秒，不再对开头单独慢放。原始保存间隔仍为 20 ms；采用相邻原始球体画面的短渐变，不捏造中间接触轨迹。全帧解码和统一倍率分别检查，帧数正确不再被等同于观感流畅。</p>"""
        )
    for name in ["B2", "B3"]:
        html += f'<h2>{name}</h2><video controls playsinline preload="metadata" poster="keyframes/{name}_010.png" src="{videos[name]["path"]}"></video>'
    html += """<p>最终 B2：A=76、B=0、截止未分配=180；B3：A=0、B=88、截止未分配=168。总数各为 256。</p><p><a href="RENDER_MANIFEST.json">来源与逐状态聚集数据</a> · <a href="frame_map.csv">逐帧时间映射</a> · <a href="DELIVERY_VERIFICATION.json">视频检查</a></p><script>const videos=[...document.querySelectorAll('video')];videos.forEach(v=>v.addEventListener('play',()=>videos.forEach(o=>{if(o!==v)o.pause()})));</script></html>"""
    (out / "OPEN_RESULTS.html").write_text(html, encoding="utf-8")
    (out / "REPRODUCE.txt").write_text(
        "PYTHONPATH=src:. /home/lzy/projects/computation_examples/_shared/python_environment/bin/python "
        "-m scripts.stage4.render_single_overview --output <new-directory>"
        + (
            f" --sphere-clusters --structure-magnification {args.structure_magnification:g}"
            if args.sphere_clusters
            else ""
        )
        + (" --constant-speed" if args.constant_speed else "")
        + "\n"
        "Then run scripts.stage4.check_linked_zoom on the output directory.\n"
        "Visualization only; no development controller resumed.\n"
    )
    print(json.dumps(dict(status="RENDERED", output=str(out))), flush=True)


if __name__ == "__main__":
    main()
