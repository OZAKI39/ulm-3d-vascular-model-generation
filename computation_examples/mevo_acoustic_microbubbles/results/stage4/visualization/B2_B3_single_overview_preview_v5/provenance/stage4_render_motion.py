"""Read-only animation of sealed Stage 4I B2/B3 trajectories."""

import argparse
import csv
import json
from pathlib import Path
import shutil
import time

import h5py
import imageio.v2 as imageio
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import numpy as np
from PIL import Image
import pyvista as pv
import yaml

from scripts.stage3.render_motion import sample_positions, screen_points, sha


ROOT = Path(__file__).resolve().parents[2]
REFERENCE = Path("/home/lzy/projects/clot_clearance/visualization/streaming_fragmentation_demo/"
                 "damage_source_v3/particle_damage_load_source_60fps.mp4")
BOUNDS = np.array([[-3.5, 11.5], [-7.2, 7.2], [-2., 2.]])
BLUE, MINT, GRID = "#339dff", "#66e0c2", "#24364c"
FPS, TRAIL_S = 60, .12


def sample_case(d, requested):
    """Display chords stop at recorded event time, not the next saved frame."""
    requested = np.asarray(requested, dtype=float)
    if requested.ndim != 1 or not np.isfinite(requested).all():
        raise ValueError("Display times must be a finite one-dimensional array")
    x, left, right, alpha = sample_positions(d["t"], d["x"], requested)
    for event in d["events"]:
        index = d["id_index"][event["particle_id"]]
        te = event["crossing_time_s"]
        before = np.searchsorted(d["t"], te, side="left") - 1
        if before < 0:
            raise ValueError("Crossing at/before the initial saved state is unsupported")
        mask = (requested > d["t"][before]) & (requested <= d["t"][before + 1])
        fraction = np.clip((requested[mask] - d["t"][before]) / (te - d["t"][before]), 0, 1)
        origin, target = d["x"][before, index], np.asarray(event["position_m"])
        x[mask, index] = origin + fraction[:, None] * (target - origin)
        # Preserve the recorded terminal coordinate without cancellation roundoff.
        x[requested >= te, index] = target
    return x, left, right, alpha


def load_data(display_times):
    manifest_path = ROOT / "results/stage4/performance_refactor/gold_baseline/GOLD_BASELINE_MANIFEST.json"
    gold = json.loads(manifest_path.read_text())
    assert gold["status"] == "SEALED"
    protected = {str(manifest_path): sha(manifest_path), str(REFERENCE): sha(REFERENCE)}
    data = {}
    for name in ["B2", "B3"]:
        case = gold["cases"][name]
        run = ROOT / case["result_id"]
        for relative in ["particles/trajectories.h5", "BRANCH_METRICS.json", "config.yaml",
                         "SOURCE_IDENTITY.json", "PHYSICS_BANNER.json", "INITIAL_POPULATION.npz"]:
            p = run / relative
            digest = sha(p)
            assert digest == gold["files"][f"cases/{name}/{relative}"]["sha256"], p
            protected[str(p)] = digest
        with h5py.File(run / "particles/trajectories.h5") as h:
            assert h.attrs["complete"]
            t, x, ids, radii, active = (h[k][:] for k in
                ["time_s", "positions_m", "particle_id", "radii_m", "active"])
        assert np.all(np.diff(t) > 0) and np.isfinite(x).all()
        assert t[0] == 0 and abs(t[-1] - 6.) < 1e-12
        branches = json.loads((run / "BRANCH_METRICS.json").read_text())
        assert branches["N_wall_loss"] == 0
        mapping = {int(pid): i for i, pid in enumerate(ids)}
        exit_time, outcome = np.full(len(ids), np.inf), np.zeros(len(ids), dtype=int)
        for event in branches["crossings"]:
            i = mapping[event["particle_id"]]
            assert not np.isfinite(exit_time[i])
            exit_time[i] = event["crossing_time_s"]
            outcome[i] = 1 if event["branch"] == "BRANCH_A" else 2
        assert np.array_equal(active, t[:, None] < exit_time[None, :])
        d = dict(t=t, x=x, ids=ids, radii=radii, events=branches["crossings"], id_index=mapping,
                 exit_time=exit_time, outcome=outcome, branches=branches, result_id=case["result_id"])
        recovered = sample_case(d, t)[0]
        assert np.array_equal(recovered, x), "Saved coordinates must be reconstructed exactly"
        d["pos"], d["left"], d["right"], d["alpha"] = sample_case(d, display_times)
        config = yaml.safe_load((run / "config.yaml").read_text())
        surface = ROOT / Path(config["geometry"]["mesh_surface"]).relative_to("/workspace/mevo_acoustic_microbubbles")
        d["surface"] = surface
        protected[str(surface)] = sha(surface)
        data[name] = d
    assert np.array_equal(data["B2"]["x"][0], data["B3"]["x"][0])
    assert data["B2"]["surface"] == data["B3"]["surface"]
    return data, protected


def wall_sections(path):
    """Sections of the actual particle-wall mesh, clipped only for the camera."""
    mesh = pv.read(path)
    lines = []
    planes = [("z", [0, 0, z]) for z in [-.0007, 0, .0007]]
    planes += [("x", [x / 1000, 0, 0]) for x in [-3, -1, 1, 3, 5, 7, 9, 10]]
    for normal, origin in planes:
        section = mesh.slice(normal=normal, origin=origin)
        cells, offset = section.lines, 0
        while offset < len(cells):
            n = cells[offset]
            pts = section.points[cells[offset + 1:offset + n + 1]] * 1000
            for a, b in zip(pts[:-1], pts[1:]):
                if np.all(np.minimum(a, b) >= BOUNDS[:, 0]) and np.all(np.maximum(a, b) <= BOUNDS[:, 1]):
                    lines.append([a, b])
            offset += n + 1
    return np.asarray(lines)


def trail_segments(d, now):
    start = max(0, now - TRAIL_S)
    saved = d["t"][(d["t"] > start) & (d["t"] < now)]
    times = np.concatenate([[start], saved, [now]])
    x = sample_case(d, times)[0]
    keep = d["exit_time"] >= start
    return np.swapaxes(x[:, keep] * 1000, 0, 1)


def make_scene(data, names, wall):
    comparison = len(names) == 2
    fig = plt.figure(figsize=(19.2, 10.8), dpi=100, facecolor="black")
    slots = [(.012, .19, .48, .62), (.512, .19, .48, .62)] if comparison else [(.015, .225, .93, .66)]
    title = "B2 / B3  |  Microbubble branch steering" if comparison else f"{names[0]}  |  Microbubble branch steering"
    fig.text(.042, .946, title, color="white", fontsize=29, weight="normal")
    fig.text(.043, .899, "Flow + primary + secondary force  |  256 bubbles per case  |  2.5 MHz",
             color="#b5c7dc", fontsize=17)
    timer = fig.text(.954, .946, "t = 0.000 s", ha="right", color=MINT, fontsize=24)
    dynamic, checks = [], {}
    for rect, name in zip(slots, names, strict=True):
        d = data[name]
        ax = fig.add_axes(rect, projection="3d", facecolor="black", computed_zorder=False)
        ax.set_proj_type("ortho")
        ax.view_init(elev=24, azim=-68)
        ax.set(xlim=BOUNDS[0], ylim=BOUNDS[1], zlim=BOUNDS[2])
        ax.set_box_aspect(BOUNDS[:, 1] - BOUNDS[:, 0], zoom=1.23 if comparison else 1.60)
        ax.set_xticks([0, 5, 10])
        ax.set_yticks([-6, -3, 0, 3, 6])
        ax.set_zticks([-1, 0, 1])
        for axis, label in zip([ax.xaxis, ax.yaxis, ax.zaxis], ["X (mm)", "Y (mm)", "Z (mm)"], strict=True):
            axis.set_label_text(label, color="white", fontsize=15 if comparison else 19)
            axis.labelpad = 8 if comparison else 15
            if comparison and label == "Z (mm)":
                axis.labelpad = 0
            axis.line.set_color("#42556c")
            axis.set_pane_color((0, 0, 0, 1))
            axis._axinfo["grid"]["color"] = GRID
            axis._axinfo["grid"]["linewidth"] = .65
        ax.tick_params(colors="white", labelsize=12 if comparison else 16, pad=1)
        ax.add_collection3d(Line3DCollection(wall, colors="#58758c", linewidths=.65, alpha=.35, zorder=1))
        ax.text(10.25, 5.9, .1, "A", color="#c4d6e4", fontsize=18, zorder=2)
        ax.text(10.25, -6.4, .1, "B", color="#c4d6e4", fontsize=18, zorder=2)
        if comparison:
            fig.text(rect[0] + .035, .828, f"{name}  |  Target: branch {'A' if name == 'B2' else 'B'}",
                     color="white", fontsize=21)
        else:
            fig.text(.043, .850, f"Target: branch {'A' if name == 'B2' else 'B'}", color="#c4d6e4", fontsize=18)
        traces = Line3DCollection(np.zeros((1, 2, 3)), colors=MINT,
                                 linewidths=.7 if comparison else 1., alpha=.48, zorder=3)
        ax.add_collection3d(traces)
        bubbles = ax.scatter(*(d["x"][0] * 1000).T, s=19 if comparison else 29,
                             c=BLUE, edgecolors="#badfff", linewidths=.35, depthshade=False, zorder=5)
        exits = ax.scatter([], [], [], s=17 if comparison else 24, facecolors="none",
                           edgecolors=MINT, linewidths=.6, depthshade=False, alpha=.65, zorder=4)
        for item in [traces, bubbles, exits]:
            item.set_clip_on(False)
        counter = fig.text(rect[0] + .035, .164, "", color="#d6e4f1", fontsize=15 if comparison else 19)
        dynamic.append((ax, name, traces, bubbles, exits, counter))
    handles = [Line2D([], [], marker="o", color="none", markerfacecolor=BLUE, markeredgecolor="#badfff",
                      markersize=7, label="Active bubbles"),
               Line2D([], [], color=MINT, lw=1.3, label="Recent path (0.12 s)"),
               Line2D([], [], marker="o", color="none", markeredgecolor=MINT,
                      markerfacecolor="none", markersize=7, label="Recorded branch crossing")]
    fig.legend(handles=handles, loc="center", bbox_to_anchor=(.5, .094), ncol=3, frameon=False,
               labelcolor="white", fontsize=16, columnspacing=2.4)
    fig.text(.042, .032, "Distal vessel view  |  Playback: 0.30×  |  Bubble symbols enlarged",
             color="#91a4b8", fontsize=13)
    fig.text(.958, .032, "Saved trajectories  |  Near-contact model approximated", ha="right",
             color="#91a4b8", fontsize=13)
    fig.canvas.draw()
    for ax, name, traces, bubbles, exits, counter in dynamic:
        projected = screen_points(ax, data[name]["x"].reshape(-1, 3) * 1000)
        lo, hi = projected.min(axis=0), projected.max(axis=0)
        assert np.all(lo > [15, 220]) and np.all(hi < [1905, 900]), (name, lo, hi)
        # Verify identical millimetre scale on all three axes.
        origin = BOUNDS.mean(axis=1)
        uv = screen_points(ax, np.vstack([origin, origin + np.eye(3)]))
        singular = np.linalg.svd((uv[1:] - uv[0]).T, compute_uv=False)
        assert np.isclose(singular[0], singular[1], rtol=1e-10)
        checks[name] = dict(screen_bounds_pixels_bottom_origin=[lo.tolist(), hi.tolist()],
                            pixels_per_mm=float(singular.mean()),
                            true_equilibrium_diameter_pixels=(data[name]["radii"] * 2000 * singular.mean()).tolist())
        labels = {}
        for axis in [ax.xaxis, ax.yaxis, ax.zaxis]:
            box = axis.label.get_window_extent(fig.canvas.get_renderer())
            assert box.x0 > 5 and box.x1 < 1915 and box.y0 > 5 and box.y1 < 1075, box
            labels[axis.label.get_text()] = list(box.bounds)
        checks[name]["axis_label_bounds_pixels_bottom_origin"] = labels
        for item in [traces, bubbles, exits, counter]:
            item.set_visible(False)
    timer.set_visible(False)
    fig.canvas.draw()
    background = fig.canvas.copy_from_bbox(fig.bbox)
    for _, _, traces, bubbles, exits, counter in dynamic:
        for item in [traces, bubbles, exits, counter]:
            item.set_visible(True)
    timer.set_visible(True)

    def render(frame, now):
        fig.canvas.restore_region(background)
        timer.set_text(f"t = {now:.3f} s")
        for ax, name, traces, bubbles, exits, counter in dynamic:
            d = data[name]
            current = d["pos"][frame] * 1000
            active = now < d["exit_time"]
            a = np.count_nonzero((~active) & (d["outcome"] == 1))
            b = np.count_nonzero((~active) & (d["outcome"] == 2))
            assert active.sum() + a + b == len(d["ids"])
            counter.set_text(f"A: {a:3d}    B: {b:3d}    {'Unassigned at cutoff' if now == 6 else 'Not yet crossed'}: {active.sum():3d}")
            bubbles._offsets3d = tuple(current[active].T)
            exits._offsets3d = tuple(current[~active].T)
            traces.set_segments(trail_segments(d, now))
            for artist in [traces, exits, bubbles]:
                artist.do_3d_projection()
                ax.draw_artist(artist)
            fig.draw_artist(counter)
        fig.draw_artist(timer)
        return np.asarray(fig.canvas.buffer_rgba())[:, :, :3].copy()

    return fig, render, checks


def write_viewer(out):
    sections = []
    for key, title in [("B2_B3", "B2 / B3 同步对照"), ("B2", "B2：目标 A 支"), ("B3", "B3：目标 B 支")]:
        stem = f"stage4_{key}_motion_3d_60fps"
        sections.append(f'<h2>{title}</h2><video controls playsinline preload="metadata" poster="keyframes/{stem}_0630.png" src="{stem}.mp4"></video><p><a href="{stem}.mp4" download>下载 MP4</a></p>')
    html = '''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>B2 / B3 微泡三维运动</title><style>body{background:#070a0f;color:#e8f2ff;font:17px/1.7 Arial,sans-serif;max-width:1440px;margin:32px auto;padding:0 24px}h1,h2{font-weight:400}video{display:block;width:100%;background:black}a{color:#66e0c2}p{color:#b5c7dc}</style>
<h1>B2 / B3 · 微泡三维运动</h1><p>使用已封存的 Stage 4I 正式结果，流动 + 主辐射力 + 次辐射力。每组 256 个相同初始微泡，物理时间 0–6 s；1920 × 1080，60 fps，0.30 倍速播放。</p>
<p>蓝色圆点表示仍在计算域内的微泡；青绿色空心圆表示已经记录的首次分支穿越位置，并非停留的微泡。最终 B2：A 76、B 0、未出流 180；B3：A 0、B 88、未出流 168。两组均无壁面损失。</p>'''
    html += "\n".join(sections)
    html += '''<p>为使微米级微泡在毫米级血管视野中可见，圆点符号统一放大；位置、血管几何和三轴空间比例未放大或变形。血管线条来自实际壁面网格的切面，画面只展示远端母管和分叉。显示帧按保存轨迹插值，并遵守记录的分支穿越时刻；没有重新计算、平滑轨迹或过滤粒子。</p>
<p>群聚中的微泡可能在投影中重叠，全部 256 个微泡始终纳入显示和计数。保存状态间隔 20 ms，动画不能分辨求解器的微秒级接触过程。近接触、近壁和有限雷诺数适用性警告仍然保留；动画不是新的物理验证，也不属于 J 敏感性算例。</p>
<p><a href="RENDER_MANIFEST.json">数据身份与渲染记录</a> · <a href="frame_map.csv">逐帧时间映射</a></p></html>'''
    (out / "OPEN_RESULTS.html").write_text(html, encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--preview-only", action="store_true")
    args = parser.parse_args()
    out = args.output.resolve()
    if list(out.glob("*.mp4")):
        raise FileExistsError("Preserve existing videos; use a new output directory")
    out.mkdir(parents=True, exist_ok=True)
    (out / "keyframes").mkdir(exist_ok=True)
    (out / "provenance").mkdir(exist_ok=True)
    font = Path("/mnt/c/Windows/Fonts/arial.ttf")
    font_manager.fontManager.addfont(str(font))
    plt.rcParams.update({"font.family": "Arial", "font.weight": "normal", "font.size": 14,
                         "axes.labelweight": "normal", "axes.titleweight": "normal"})
    assert sha(font_manager.findfont("Arial", fallback_to_default=False)) == sha(font)
    times = np.r_[np.zeros(30), np.arange(1201) * .005, np.full(90, 6.)]
    data, protected = load_data(times)
    protected[str(font)] = sha(font)
    wall = wall_sections(data["B2"]["surface"])
    (out / "provenance/INPUT_FILES_SHA256.json").write_text(json.dumps(protected, indent=2) + "\n")
    for script in [Path(__file__), ROOT / "scripts/stage3/render_motion.py"]:
        shutil.copy2(script, out / "provenance" / (script.parent.name + "_" + script.name))
    videos, checks = {}, {}
    started = time.perf_counter()
    for names in [["B2"], ["B3"], ["B2", "B3"]]:
        stem = "stage4_" + "_".join(names) + "_motion_3d_60fps"
        fig, render, bounds = make_scene(data, names, wall)
        checks[stem] = bounds
        for frame in [30, 230, 630, 1230]:
            imageio.imwrite(out / "keyframes" / f"{stem}_{frame:04d}.png", render(frame, times[frame]))
        if args.preview_only:
            plt.close(fig)
            continue
        with imageio.get_writer(out / f"{stem}.mp4", fps=FPS, codec="libx264", quality=8,
                                macro_block_size=1, pixelformat="yuv420p",
                                ffmpeg_params=["-movflags", "+faststart"]) as writer:
            for frame, now in enumerate(times):
                writer.append_data(render(frame, now))
                if frame % 180 == 0 or frame == len(times) - 1:
                    print(json.dumps(dict(video=stem, frame=frame, total=len(times),
                                          elapsed_s=round(time.perf_counter() - started, 2))), flush=True)
        plt.close(fig)
        reader = imageio.get_reader(out / f"{stem}.mp4")
        count, meta = reader.count_frames(), reader.get_meta_data()
        assert count == len(times) and meta["fps"] == FPS and meta["size"] == (1920, 1080)
        for frame in [0, len(times) // 2, len(times) - 1]:
            assert reader.get_data(frame).shape == (1080, 1920, 3)
        reader.close()
        videos[stem] = dict(frames=count, fps=FPS, duration_s=count / FPS, resolution=[1920, 1080],
                            sha256=sha(out / f"{stem}.mp4"))
    if args.preview_only:
        (out / "PREVIEW_CHECKS.json").write_text(json.dumps(checks, indent=2) + "\n")
        return
    with (out / "frame_map.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["case", "frame", "video_time_s", "physical_time_s", "left_saved_state", "right_saved_state", "linear_alpha", "A", "B", "not_crossed"])
        for name, d in data.items():
            for i, now in enumerate(times):
                crossed = now >= d["exit_time"]
                writer.writerow([name, i, i / FPS, now, d["left"][i], d["right"][i], d["alpha"][i],
                    np.count_nonzero(crossed & (d["outcome"] == 1)), np.count_nonzero(crossed & (d["outcome"] == 2)), np.count_nonzero(~crossed)])
    np.savez_compressed(out / "display_samples.npz", physical_time_s=times,
                        **{f"{name}_positions_m": d["pos"] for name, d in data.items()},
                        **{f"{name}_particle_id": d["ids"] for name, d in data.items()})
    changed = [path for path, digest in protected.items() if sha(path) != digest]
    assert not changed
    manifest = dict(status="PASS", reference_video=str(REFERENCE), reference_sha256=sha(REFERENCE),
        renderer_sha256=sha(__file__), videos=videos, physical_interval_s=[0, 6], playback_speed=.30,
        initial_hold_s=.5, final_hold_s=1.5, font=dict(family="Arial", weight="normal", sha256=sha(font)),
        camera=dict(projection="orthographic", elevation_deg=24, azimuth_deg=-68, limits_mm=BOUNDS.tolist(),
                    spatial_axes_scale="Equal millimetre scale; distal vessel crop contains every recorded particle coordinate"),
        bubble_glyph="Enlarged screen symbols, not physical radii; coordinates and displacement unchanged",
        colors=dict(background="#000000", bubbles=BLUE, trails=MINT, grid=GRID),
        interpolation="Piecewise linear display sampling with exact event-time terminal coordinates; all saved states reconstructed bitwise",
        exit_semantics="Active glyph changes to recorded hollow crossing marker at exact saved event time; counters use the same event time",
        trajectory_scale=1, smoothing=False, particle_filtering=False, new_simulation=False,
        all_ids_accounted_for_every_frame=True, wall_mesh=str(data["B2"]["surface"].relative_to(ROOT)),
        trail_s=TRAIL_S, screen_checks=checks,
        sources={name: dict(result_id=d["result_id"], saved_states=len(d["t"]), particle_count=len(d["ids"]),
                            A=d["branches"]["N_A"], B=d["branches"]["N_B"], unassigned=d["branches"]["N_unassigned"],
                            wall_loss=d["branches"]["N_wall_loss"]) for name, d in data.items()},
        scope="Sealed Stage 4I results, not incomplete J. 20 ms saved states do not resolve microsecond contacts; display chords may cross curved walls or overlap. No geometry/contact correction or new physical validation. Finite-Re, near-wall and near-contact limitations remain.",
        changed_inputs=changed, elapsed_s=time.perf_counter() - started)
    (out / "RENDER_MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n")
    write_viewer(out)
    sheet = Image.new("RGB", (1440, 540), "black")
    for row, name in enumerate(["B2", "B3"]):
        for col, frame in enumerate([30, 630, 1230]):
            with Image.open(out / "keyframes" / f"stage4_{name}_motion_3d_60fps_{frame:04d}.png") as im:
                sheet.paste(im.resize((480, 270), Image.Resampling.LANCZOS), (col * 480, row * 270))
    sheet.save(out / "contact_sheet.png")
    delivery = {str(p.relative_to(out)): dict(bytes=p.stat().st_size, sha256=sha(p))
                for p in sorted(out.rglob("*")) if p.is_file() and p.name != "render.log"}
    (out / "DELIVERY_SHA256.json").write_text(json.dumps(delivery, indent=2) + "\n")
    print(json.dumps(dict(status="PASS", output=str(out), elapsed_s=manifest["elapsed_s"])), flush=True)


if __name__ == "__main__":
    main()
