"""Linked overview/zoom using all sealed states and explicit image cross-dissolves.

Intermediate video frames blend images, never particle coordinates: twenty-ms
trajectory samples cannot reconstruct near-contact motion. The overview ROI
locator alone follows a display chord; it is not a new particle trajectory.
"""

import argparse
import csv
import itertools
import json
from pathlib import Path
import shutil
import subprocess
import time

import imageio.v2 as imageio
import imageio_ffmpeg
from matplotlib import colors, font_manager
from matplotlib.cm import ScalarMappable
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
from PIL import Image, ImageDraw, ImageFont

from scripts.stage4 import render_aggregation as original

base = original.base
ROOT, FPS = base.ROOT, 60
FONT = Path("/mnt/c/Windows/Fonts/arial.ttf")
BLUE, MINT = original.BLUE, original.MINT
INDICES = list(range(301))
HALF_UM = 25.0
MICRO_BOUNDS = np.array([[-HALF_UM, HALF_UM]] * 3)
# Every original sample lands on an exact video frame. No held keyframes.
INTERVAL_FRAMES = np.array([24] * 10 + [4] * 290)
STATE_FRAMES = np.r_[0, np.cumsum(INTERVAL_FRAMES)]
FRAME_COUNT = int(STATE_FRAMES[-1]) + 1


def timeline():
    rows = []
    for left, n in enumerate(INTERVAL_FRAMES):
        for j in range(int(n)):
            alpha = j / n
            rows.append((left, left + 1, float(alpha), (left + alpha) * 0.02))
    rows.append((300, 300, 0.0, 6.0))
    return rows


def scene(name, d):
    fig = base.plt.figure(figsize=(19.2, 10.8), dpi=100, facecolor="black")
    fig.text(0.035, 0.946, f"{name}  |  Acoustic field & bubble aggregation", color="white", fontsize=28)
    fig.text(
        0.035,
        0.902,
        "2.5 MHz  |  Flow + primary + secondary force  |  256 bubbles",
        color="#bed0e0",
        fontsize=17,
    )
    fig.text(0.035, 0.848, "PZT sources + pressure amplitude", color="white", fontsize=21)
    fig.text(0.610, 0.848, "Zoom of highlighted region", color=MINT, fontsize=21)
    elements = d["config"]["pzt"]
    for row, pair in enumerate([elements[:2], elements[2:]]):
        label = "     ".join(
            f"{e['id'].replace('_plus', '+').replace('_minus', '-').upper()}: "
            f"{e['amplitude_Pa'] / 1000:.2f} kPa / {np.rad2deg(e['phase_rad']):.0f}°"
            for e in pair
        )
        fig.text(0.035, 0.812 - row * 0.027, label, color="#dec09a", fontsize=13)
    global_ax = fig.add_axes(
        [0.002, 0.16, 0.566, 0.60], projection="3d", facecolor="black", computed_zorder=False
    )
    original.axes_style(global_ax, original.FIELD_BOUNDS, ["X (mm)", "Y (mm)", "Z (mm)"], zoom=1.22)
    global_ax.set_xticks([0, 5, 10])
    global_ax.set_yticks([-5, 0, 5])
    global_ax.set_zticks([-4, 0, 4])
    norm, cmap = original.source_and_field(global_ax, d)
    pop = global_ax.scatter(
        [], [], [], c=BLUE, edgecolors="#c1e3ff", linewidths=0.22, s=13, depthshade=False, zorder=10
    )
    exits = global_ax.scatter(
        [],
        [],
        [],
        facecolors="none",
        edgecolors=MINT,
        linewidths=0.5,
        s=13,
        alpha=0.6,
        depthshade=False,
        zorder=9,
    )
    micro_ax = fig.add_axes(
        [0.584, 0.14, 0.394, 0.64], projection="3d", facecolor="black", computed_zorder=False
    )
    original.axes_style(micro_ax, MICRO_BOUNDS, ["ΔX (µm)", "ΔY (µm)", "ΔZ (µm)"], zoom=1.00)
    for axis in [micro_ax.xaxis, micro_ax.yaxis, micro_ax.zaxis]:
        axis.set_ticks([-20, 0, 20])
    particles = Poly3DCollection([], edgecolors="none", antialiased=False, zsort="average", zorder=8)
    micro_ax.add_collection3d(particles)
    bar_ax = fig.add_axes([0.045, 0.104, 0.435, 0.013])
    bar = fig.colorbar(
        ScalarMappable(norm=norm, cmap=cmap), cax=bar_ax, orientation="horizontal", ticks=[0, 10, 20, 30, 39]
    )
    bar.ax.tick_params(colors="white", labelsize=12, pad=3)
    bar.outline.set_visible(False)
    fig.text(0.495, 0.102, "|P| (kPa)", color="#c2d2e1", fontsize=14)
    # Explanation and population statistics live in the viewer/manifest, not a footer.
    dynamic = [(global_ax, exits), (global_ax, pop), (micro_ax, particles)]
    for _, artist in dynamic:
        artist.set_visible(False)
    fig.canvas.draw()
    checks = {}
    for label, ax in [("overview", global_ax), ("zoom", micro_ax)]:
        boxes = {}
        for axis in [ax.xaxis, ax.yaxis, ax.zaxis]:
            box = axis.label.get_window_extent(fig.canvas.get_renderer())
            assert box.x0 > 3 and box.x1 < 1917 and box.y0 > 3 and box.y1 < 1077, (label, box)
            boxes[axis.label.get_text()] = list(box.bounds)
        uv = base.screen_points(ax, np.vstack([np.zeros(3), np.eye(3)]))
        scale = np.linalg.svd((uv[1:] - uv[0]).T, compute_uv=False)
        assert np.isclose(scale[0], scale[1], rtol=1e-10)
        checks[label] = dict(axis_label_bounds=boxes, pixels_per_axis_unit=float(scale.mean()))
    checks["zoom"]["equilibrium_diameter_pixels"] = 3.56 * checks["zoom"]["pixels_per_axis_unit"]
    checks["screen_magnification"] = (
        1000 * checks["zoom"]["pixels_per_axis_unit"] / checks["overview"]["pixels_per_axis_unit"]
    )
    fig.text(
        0.610,
        0.805,
        f"50 µm cube  |  ×{checks['screen_magnification']:.0f}  |  Bubble #0",
        color="#bed0e0",
        fontsize=15,
    )
    fig.canvas.draw()
    background = fig.canvas.copy_from_bbox(fig.bbox)
    for _, artist in dynamic:
        artist.set_visible(True)
    template, intensity = original.sphere_faces()
    vertices = np.array(list(itertools.product([-HALF_UM, HALF_UM], repeat=3)))
    corners = base.screen_points(micro_ax, vertices)
    corners[:, 1] = 1080 - corners[:, 1]
    checks["zoom_cube_bounds_top_origin"] = [corners.min(axis=0).tolist(), corners.max(axis=0).tolist()]
    anchor = base.screen_points(global_ax, d["x"][:, d["anchor"]] * 1000)
    anchor[:, 1] = 1080 - anchor[:, 1]
    checks["anchor_saved_screen_positions"] = anchor.tolist()
    # True-size ROI cube projection, separate from the deliberately enlarged locator.
    offsets = base.screen_points(global_ax, vertices * 0.001) - base.screen_points(
        global_ax, np.zeros((1, 3))
    )
    offsets[:, 1] *= -1
    checks["actual_roi_screen_offsets"] = offsets.tolist()

    def render(k):
        state = d["states"][k]
        fig.canvas.restore_region(background)
        x = d["x"][k] * 1000
        active = state["active_mask"]
        pop._offsets3d, exits._offsets3d = tuple(x[active].T), tuple(x[~active].T)
        meshes, facecolors = [], []
        for i in np.flatnonzero(state["roi"]):
            meshes.append(template * d["radii"][i] * 1e6 + state["local"][i])
            rgb = np.array(colors.to_rgb(MINT if i == d["anchor"] else BLUE))
            facecolors.append(np.c_[intensity[:, None] * rgb, np.ones(len(intensity))])
        particles.set_visible(bool(meshes))
        if meshes:
            verts = np.concatenate(meshes)
            particles.set_verts(verts)
            particles.set_facecolor(np.concatenate(facecolors))
            screen = base.screen_points(micro_ax, verts.reshape(-1, 3))
            assert np.all(screen.min(axis=0) > [1120, 170]), screen.min(axis=0)
            assert np.all(screen.max(axis=0) < [1900, 840]), screen.max(axis=0)
        for ax, artist in dynamic:
            if artist.get_visible():
                artist.do_3d_projection()
                ax.draw_artist(artist)
        return Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[:, :, :3].copy())

    return fig, render, checks


def overlay(im, d, checks, left, right, alpha, physical_time):
    """Crisp clock and a common ROI callout; no scientific glyph repositioning."""
    draw = ImageDraw.Draw(im)
    clock_font = ImageFont.truetype(str(FONT), 33)
    small_font = ImageFont.truetype(str(FONT), 21)
    draw.text((1852, 28), f"t = {physical_time * 1000:.1f} ms", font=clock_font, fill=MINT, anchor="ra")
    end = d["exit_time"][d["anchor"]]
    live = physical_time < end
    positions = np.asarray(checks["anchor_saved_screen_positions"])
    loc = (1 - alpha) * positions[left] + alpha * positions[right]
    if np.isfinite(end) and d["t"][left] < end <= d["t"][right]:
        fraction = min(1.0, (physical_time - d["t"][left]) / (end - d["t"][left]))
        loc = (1 - fraction) * positions[left] + fraction * positions[right]
    x, y = loc
    (xmin, ymin), (xmax, ymax) = checks["zoom_cube_bounds_top_origin"]
    xmin, ymin, xmax, ymax = xmin - 20, ymin - 16, xmax + 20, ymax + 16
    color, faint = ("#66e0c2", "#35685e") if live else ("#758780", "#263c37")
    # Open corners and two leaders form a single enlarged view without an outer frame.
    for cy in [ymin, ymax]:
        sy = -9 if cy == ymin else 9
        draw.line([(x + 10, y + sy), (xmin, cy)], fill=faint, width=2)
    for cx, sign_x in [(xmin, 1), (xmax, -1)]:
        for cy, sign_y in [(ymin, 1), (ymax, -1)]:
            draw.line([(cx + sign_x * 22, cy), (cx, cy), (cx, cy + sign_y * 22)], fill=color, width=2)
    draw.rectangle([x - 9, y - 9, x + 9, y + 9], outline=color, width=2)
    # A subpixel cube is genuinely 50 um across; the 18-px bracket is only a locator.
    offsets = np.asarray(checks["actual_roi_screen_offsets"])
    for a, b in itertools.combinations(range(8), 2):
        if (a ^ b) in (1, 2, 4):
            draw.line([tuple(loc + offsets[a]), tuple(loc + offsets[b])], fill=color, width=1)
    draw.text((x - 13, y - 37), "ROI", font=small_font, fill=color)
    if not live:
        text = f"Bubble #0 crossed B at {end:.3f} s"
        draw.text(((xmin + xmax) / 2, (ymin + ymax) / 2), text, font=small_font, fill="#bed0e0", anchor="mm")
    return im


def blend(left_image, right_image, alpha, d, left, right, physical_time):
    """Cross-dissolve actual raster states; ensure no ghost after an absorbing exit."""
    if alpha == 0:
        return left_image.copy(), 0.0
    frame = Image.blend(left_image, right_image, alpha)
    micro_alpha = alpha
    end = d["exit_time"][d["anchor"]]
    if d["t"][left] < end <= d["t"][right]:
        micro_alpha = min(1.0, (physical_time - d["t"][left]) / (end - d["t"][left]))
        box = (1120, 240, 1920, 950)
        frame.paste(Image.blend(left_image.crop(box), right_image.crop(box), micro_alpha), box)
    return frame, micro_alpha


def writer(path):
    return imageio.get_writer(
        path,
        fps=FPS,
        codec="libx264",
        quality=None,
        macro_block_size=1,
        ffmpeg_log_level="error",
        pixelformat="yuv420p",
        output_params=[
            "-preset",
            "medium",
            "-crf",
            "18",
            "-threads",
            "4",
            "-profile:v",
            "high",
            "-level:v",
            "4.2",
            "-g",
            "120",
            "-maxrate",
            "12M",
            "-bufsize",
            "24M",
            "-movflags",
            "+faststart",
        ],
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--preview-only", action="store_true")
    args = parser.parse_args()
    out = args.output.resolve()
    assert not list(out.glob("*.mp4")), "Use a new directory; preserve previous exports"
    (out / "keyframes").mkdir(parents=True, exist_ok=True)
    (out / "provenance").mkdir(exist_ok=True)
    font_manager.fontManager.addfont(str(FONT))
    base.plt.rcParams.update(
        {"font.family": "Arial", "font.weight": "normal", "font.size": 14, "axes.labelweight": "normal"}
    )
    data, protected = original.load(INDICES, half_width_um=HALF_UM)
    protected[str(FONT)] = base.sha(FONT)
    sources = [
        Path(__file__),
        Path(original.__file__),
        Path(base.__file__),
        ROOT / "scripts/stage3/render_motion.py",
        ROOT / "src/mevo_mb/bubbles/metrics.py",
    ]
    for path in sources:
        shutil.copy2(path, out / "provenance" / f"{path.parent.name}_{path.name}")
    source_hashes = {str(path): base.sha(path) for path in sources}
    checks, videos, frame_rows = {}, {}, []
    start = time.perf_counter()
    for name, d in data.items():
        fig, render, checks[name] = scene(name, d)
        snapshots = [0, 2, 10, 70, 72, 100, 300]
        if args.preview_only:
            for k in snapshots:
                im = overlay(render(k), d, checks[name], k, k, 0.0, d["t"][k])
                im.save(out / "keyframes" / f"{name}_{k:03d}.png")
            base.plt.close(fig)
            continue
        video = out / f"stage4_{name}_linked_zoom_60fps.mp4"
        left_image, right_image = render(0), render(1)
        previous = 0
        with writer(video) as encoder:
            for frame_index, (left, right, alpha, physical_time) in enumerate(timeline()):
                if left != previous:
                    left_image = right_image
                    right_image = render(right) if right > left else left_image
                    previous = left
                im, micro_alpha = blend(left_image, right_image, alpha, d, left, right, physical_time)
                im = overlay(im, d, checks[name], left, right, alpha, physical_time)
                encoder.append_data(np.asarray(im))
                if alpha == 0 and left in snapshots:
                    im.save(out / "keyframes" / f"{name}_{left:03d}.png")
                if frame_index in [12, 36, 252, 489]:
                    im.save(out / "keyframes" / f"{name}_transition_{frame_index:04d}.png")
                frame_rows.append(
                    [
                        name,
                        frame_index,
                        frame_index / FPS,
                        physical_time,
                        left,
                        right,
                        alpha,
                        micro_alpha,
                        physical_time < d["exit_time"][d["anchor"]],
                    ]
                )
                if frame_index % 120 == 0:
                    print(
                        json.dumps(
                            dict(
                                case=name,
                                frame=frame_index,
                                total=FRAME_COUNT,
                                elapsed_s=time.perf_counter() - start,
                            )
                        ),
                        flush=True,
                    )
        base.plt.close(fig)
        videos[name] = dict(
            path=video.name,
            frames=FRAME_COUNT,
            fps=FPS,
            resolution=[1920, 1080],
            duration_s=FRAME_COUNT / FPS,
            sha256=base.sha(video),
        )
        print(json.dumps(dict(case=name, status="ENCODED")), flush=True)
    (out / "SCREEN_CHECKS.json").write_text(json.dumps(checks, indent=2) + "\n")
    if args.preview_only:
        return
    # Deliver a 1080p comparison to reduce decoder load; singles retain full detail.
    comparison = out / "stage4_B2_B3_linked_zoom_60fps.mp4"
    subprocess.run(
        [
            imageio_ffmpeg.get_ffmpeg_exe(),
            "-v",
            "error",
            "-n",
            "-threads",
            "2",
            "-i",
            str(out / videos["B2"]["path"]),
            "-threads",
            "2",
            "-i",
            str(out / videos["B3"]["path"]),
            "-filter_complex_threads",
            "2",
            "-filter_complex",
            "[0:v]scale=960:540[a];[1:v]scale=960:540[b];[a][b]hstack[v]",
            "-map",
            "[v]",
            "-c:v",
            "libx264",
            "-crf",
            "18",
            "-preset",
            "medium",
            "-threads",
            "4",
            "-pix_fmt",
            "yuv420p",
            "-g",
            "120",
            "-maxrate",
            "12M",
            "-bufsize",
            "24M",
            "-movflags",
            "+faststart",
            str(comparison),
        ],
        check=True,
    )
    videos["B2_B3"] = dict(
        path=comparison.name,
        frames=FRAME_COUNT,
        fps=FPS,
        resolution=[1920, 540],
        duration_s=FRAME_COUNT / FPS,
        sha256=base.sha(comparison),
    )
    with (out / "frame_map.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "case",
                "frame",
                "video_time_s",
                "display_time_s",
                "source_left",
                "source_right",
                "image_blend_fraction",
                "zoom_image_blend_fraction",
                "anchor_active",
            ]
        )
        w.writerows(frame_rows)
    assert all(base.sha(path) == digest for path, digest in (protected | source_hashes).items())
    manifest = dict(
        status="RENDERED_PENDING_DELIVERY_CHECK",
        videos=videos,
        protected_inputs=protected,
        renderer_sources=source_hashes,
        source_runs={n: d["result_id"] for n, d in data.items()},
        temporal_display=dict(
            original_states=301,
            original_spacing_s=0.02,
            interval_frames=INTERVAL_FRAMES.tolist(),
            exact_state_video_frames=STATE_FRAMES.tolist(),
            method="Cross-dissolve adjacent original raster states. No particle-position interpolation or new contact dynamics.",
            clock="Display time between the two bracketing saved times; not a newly integrated state.",
            speed="0–0.2 s: 0.05x; 0.2–6 s: 0.30x. No long still-image holds.",
            exit="Zoom dissolve ends exactly at the recorded absorbing crossing; no post-exit ghost.",
        ),
        zoom=dict(
            anchor_particle_id=0,
            half_width_um=HALF_UM,
            equilibrium_radius_um=1.78,
            transform="Translation only x_i-x_0, identical camera orientation, isotropic scale",
            connection="True 50 um cube plus enlarged locator at its center; two leaders to zoom corners",
            locator_motion="Display chord only, clipped at the recorded crossing; no effect on particle geometry",
        ),
        field=dict(
            quantity="Steady pressure amplitude |P|",
            range_kPa=[0, 39],
            frequency_Hz=2.5e6,
            slices="XY z=0 and YZ x=0 at native resolution",
            source_geometry="Four actual 6 mm PZT discs",
            inputs={n: str(d["field_path"].relative_to(ROOT)) for n, d in data.items()},
        ),
        states={
            name: [
                {key: val for key, val in s.items() if key not in ["local", "roi", "group", "active_mask"]}
                | dict(
                    visible_particle_ids=d["ids"][s["roi"]].tolist(),
                    cluster_particle_ids=d["ids"][s["group"]].tolist(),
                )
                for s in d["states"]
            ]
            for name, d in data.items()
        },
        science_scope="Original sealed Stage 4I B2/B3. Spatial aggregation is not fusion. Existing near-contact, near-wall and finite-Re limitations remain.",
        screen_checks=checks,
        elapsed_s=time.perf_counter() - start,
    )
    (out / "RENDER_MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n")
    html = """<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>B2/B3 局部放大动画</title><style>body{background:#070a0f;color:#eef4ff;font:17px/1.7 Arial,sans-serif;max-width:1500px;margin:30px auto;padding:0 24px}h1,h2{font-weight:400}p{color:#b7cee2}video{width:100%;background:black}a{color:#66e0c2}</style><h1>B2/B3：声场与局部放大</h1><p>右侧与左侧绿色 ROI 指示位置相连，显示同一位置的 50 µm 立方区域，相机方向相同，球半径保持 1.78 µm。左侧定位框为了可见性放大；实际 ROI 在毫米尺度下不足几个像素。声源几何和稳态声压幅值来自正式计算。</p><p>每例使用完整 301 个原始保存状态，输出恒定 60 fps；前 200 ms 慢放，此后播放至 6 s。帧间使用原始画面的交叉渐变，可能短暂呈现双重影像；它仅改善播放连续性，不代表求解器输出了中间轨迹。没有用会造成微泡穿透的位置插值重建近接触动力学。</p>"""
    for name in ["B2", "B3"]:
        html += f'<h2>{name}</h2><video controls playsinline preload="metadata" poster="keyframes/{name}_010.png" src="{videos[name]["path"]}"></video>'
    html += f'<h2>B2 / B3 同步对照</h2><video controls playsinline preload="metadata" src="{comparison.name}"></video>'
    html += """<p>右侧始终跟随 #0。B3 的 #0 在 1.434823 s 穿过 B 判定平面后，局部泡群停止显示，全景继续。最终 B2：A=76、B=0、截止未分配=180；B3：A=0、B=88、截止未分配=168。聚集不代表融合，近壁、近接触和有限雷诺数的模型限制保持。</p><p><a href="RENDER_MANIFEST.json">来源与逐状态校验</a> · <a href="frame_map.csv">逐帧时间及渐变映射</a></p><script>const videos=[...document.querySelectorAll('video')];videos.forEach(v=>v.addEventListener('play',()=>videos.forEach(o=>{if(o!==v)o.pause()})));</script></html>"""
    (out / "OPEN_RESULTS.html").write_text(html, encoding="utf-8")
    (out / "REPRODUCE.txt").write_text(
        "PYTHONPATH=src:. /home/lzy/projects/computation_examples/_shared/python_environment/bin/python "
        "-m scripts.stage4.render_linked_zoom --output <new-directory>\n"
        "Development remains paused; this is visualization of sealed Stage 4I only.\n"
    )
    print(json.dumps(dict(status="RENDERED", output=str(out))), flush=True)


if __name__ == "__main__":
    main()
