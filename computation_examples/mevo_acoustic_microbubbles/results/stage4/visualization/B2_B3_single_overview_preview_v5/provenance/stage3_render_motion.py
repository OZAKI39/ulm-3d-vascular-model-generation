"""Render saved Stage 3 trajectories; no solver imports or physics modifications."""

import argparse
import csv
import hashlib
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
from mpl_toolkits.mplot3d import proj3d
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import numpy as np
from PIL import Image


ROOT = Path(__file__).resolve().parents[2]
REFERENCE = Path(
    "/home/lzy/projects/clot_clearance/visualization/streaming_fragmentation_demo/"
    "streaming_traction_v4/particle_damage_streaming_traction_60fps.mp4"
)
BLUE = "#339dff"
MINT = "#66e0c2"
GRID = "#24364c"
LABELS = {
    "C1": "C1  |  Flow only",
    "C2": "C2  |  Flow + primary force",
    "C3": "C3  |  Flow + primary + secondary",
    "C4": "C4  |  Same forces, PZT phase shifted",
}
# A cloud close-up, not the 20 mm vessel boundary. Equal physical scale on all axes.
BOUNDS_UM = np.array([[-90, 90], [-10, 170], [-50, 130]])


def sha(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def sample_positions(times, positions, requested):
    """Piecewise linear display sampling, exact endpoints, no extrapolation."""
    requested = np.asarray(requested)
    if np.any(requested < times[0]) or np.any(requested > times[-1]):
        raise ValueError("Display sampling outside the saved time interval")
    left = np.clip(np.searchsorted(times, requested, side="right") - 1, 0, len(times) - 1)
    right = np.minimum(left + 1, len(times) - 1)
    dt = times[right] - times[left]
    alpha = np.divide(requested - times[left], dt, out=np.zeros_like(requested), where=dt > 0)
    # This form reproduces exact saved values when alpha == 0.
    sampled = positions[left] + alpha[..., None, None] * (positions[right] - positions[left])
    return sampled, left, right, alpha


def screen_points(ax, xyz):
    u, v, _ = proj3d.proj_transform(*np.asarray(xyz).T, ax.get_proj())
    return ax.transData.transform(np.column_stack([u, v]))


def write_viewer(out):
    (out / "OPEN_RESULTS.html").write_text('''<!doctype html>
<html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Stage 3 微泡三维运动</title>
<style>
body{background:#070a0f;color:#e8f2ff;font:17px/1.7 Arial,sans-serif;margin:32px auto;padding:0 24px;max-width:1440px}
h1{font-size:30px;font-weight:400}h2{font-size:23px;font-weight:400;margin-top:36px}
video{display:block;width:100%;background:black;border-radius:5px}a{color:#66e0c2}
p{color:#b5c7dc}small{color:#91a4b8}.links{display:flex;gap:24px;flex-wrap:wrap}
</style><h1>Stage 3 · 微泡三维运动</h1>
<p>正式圆管算例的 64 个微泡，物理时间 0–2 ms。黑底、正交视角、蓝色微泡、青绿色短尾迹；
真实空间比例，10,000 倍慢放。每段约 22 秒，1920 × 1080，60 fps。</p>
<h2>C3：流动 + 主辐射力 + 次辐射力</h2>
<video controls playsinline preload="metadata" poster="preview.png" src="stage3_C3_motion_3d_60fps.mp4"></video>
<h2>C1–C4 同步对照</h2>
<p>C1 仅流动；C2 加主辐射力；C3 再加次辐射力；C4 在 C3 基础上将第二个 PZT 相位改为 90°。
四组初始泡群、时钟、镜头和坐标尺度相同。</p>
<video controls playsinline preload="metadata" poster="keyframes/stage3_C1_C4_comparison_3d_60fps_1230.png"
src="stage3_C1_C4_comparison_3d_60fps.mp4"></video>
<p class="links"><a href="stage3_C3_motion_3d_60fps.mp4" download>下载主动画</a>
<a href="stage3_C1_C4_comparison_3d_60fps.mp4" download>下载对照动画</a>
<a href="RENDER_MANIFEST.json">数据来源与渲染记录</a><a href="frame_map.csv">逐帧时间映射</a></p>
<small>画面是泡群局部特写，坐标单位 µm；尾迹保留最近 0.25 ms。只对显示帧线性插值，没有重新运行计算、
放大位移、过滤微泡或模拟径向脉动。近接触动力学仍沿用 Stage 3 的近似适用范围。</small></html>
''', encoding="utf-8")


def load_data(display_times):
    state_path = ROOT / "results/stage3/cylinder/STATE.json"
    state = json.loads(state_path.read_text())
    if state["status"] != "PASS":
        raise ValueError("Formal Stage 3 cylinder acceptance is not PASS")
    protected = [state_path, ROOT / "PHYSICS_ACCEPTANCE.json", REFERENCE]
    data = {}
    for row in state["rows"]:
        run = ROOT / row["run"]
        path = run / "particles/trajectories.h5"
        protected.extend([path, run / "config.yaml", run / "SOURCE_IDENTITY.json",
                          run / "PHYSICS_BANNER.json", run / "acoustic/CONTROL.json"])
        with h5py.File(path) as h:
            if not h.attrs["complete"]:
                raise ValueError(f"Incomplete trajectory: {path}")
            t, x, radii, ids = (h[k][:] for k in ["time_s", "positions_m", "radii_m", "particle_id"])
        assert np.all(np.diff(t) > 0) and np.isfinite(x).all()
        pos, left, right, alpha = sample_positions(t, x, display_times)
        recovered, _, _, _ = sample_positions(t, x, t)
        assert np.array_equal(recovered, x)
        assert np.array_equal(pos[0], x[0]) and np.array_equal(pos[-1], x[-1])
        ii, jj = np.triu_indices(len(radii), 1)
        gaps = np.linalg.norm(pos[:, ii] - pos[:, jj], axis=2) - (radii[ii] + radii[jj])
        # Includes saved vertices: trails cannot cut corners across accepted solver states.
        data[row["scenario"]] = dict(
            t=t, x=x, radius=radii, ids=ids, pos=pos, left=left, right=right, alpha=alpha,
            row=row, source=str(path.relative_to(ROOT)),
            minimum_display_gap_m=float(gaps.min()),
        )
    assert all(np.array_equal(d["x"][0], data["C3"]["x"][0]) for d in data.values())
    return data, {str(p): sha(p) for p in protected}


def trail_segments(d, physical_time, current):
    start = max(d["t"][0], physical_time - 0.00025)
    a = np.searchsorted(d["t"], start, side="right")
    b = np.searchsorted(d["t"], physical_time, side="left")
    first = sample_positions(d["t"], d["x"], np.array([start]))[0][0]
    points = np.concatenate([first[None], d["x"][a:b], current[None]], axis=0) * 1e6
    return np.swapaxes(points, 0, 1)


def make_scene(data, comparison):
    fig = plt.figure(figsize=(19.2, 10.8), dpi=100, facecolor="black")
    if comparison:
        slots = [(0.035, 0.54, 0.43, 0.31), (0.535, 0.54, 0.43, 0.31),
                 (0.035, 0.13, 0.43, 0.31), (0.535, 0.13, 0.43, 0.31)]
        names = list(LABELS)
        title = "Stage 3  |  Microbubble motion: four controls"
        subtitle = "Identical initial cloud, camera and scale  |  64 bubbles per case  |  Flow along +X"
    else:
        slots = [(0.035, 0.15, 0.93, 0.68)]
        names = ["C3"]
        title = "Stage 3  |  Microbubble migration & aggregation"
        subtitle = "C3: flow + primary + secondary force  |  64 bubbles  |  2.5 MHz  |  Flow along +X"
    fig.text(0.042, 0.947, title, color="white", fontsize=28, weight="normal")
    fig.text(0.043, 0.900, subtitle, color="#b5c7dc", fontsize=16)
    timer = fig.text(0.954, 0.947, "t = 0.000 ms", ha="right", color=MINT, fontsize=23)
    axes, dynamic = [], []
    for rect, name in zip(slots, names):
        d = data[name]
        ax = fig.add_axes(rect, projection="3d", facecolor="black", computed_zorder=False)
        ax.set_proj_type("ortho")
        ax.view_init(elev=24, azim=-68)
        ax.set(xlim=BOUNDS_UM[0], ylim=BOUNDS_UM[1], zlim=BOUNDS_UM[2])
        ax.set_box_aspect((1, 1, 1), zoom=1.05 if comparison else 1.18)
        ax.set_xticks([-80, 0, 80])
        ax.set_yticks([0, 80, 160])
        ax.set_zticks([-40, 40, 120])
        fontsize = 13 if comparison else 19
        for axis, label in zip([ax.xaxis, ax.yaxis, ax.zaxis], ["X (µm)", "Y (µm)", "Z (µm)"]):
            axis.set_label_text(label, color="white", fontsize=fontsize)
            axis.line.set_color("#42556c")
            axis.set_pane_color((0, 0, 0, 1))
            axis._axinfo["grid"]["color"] = GRID
            axis._axinfo["grid"]["linewidth"] = 0.65
            axis.labelpad = 6 if comparison else 16
        ax.tick_params(colors="white", labelsize=11 if comparison else 16, pad=1)
        if comparison:
            fig.text(rect[0] + 0.025, rect[1] + rect[3] + 0.012,
                     LABELS[name], color="white", fontsize=18)
            phase = "[0°, 90°, 0°, 0°]" if name == "C4" else "[0°, 0°, 0°, 0°]"
            note = "PZT off" if name == "C1" else f"PZT phases: {phase}"
            fig.text(rect[0] + 0.025, rect[1] + rect[3] - 0.018, note,
                     color="#90a9c2", fontsize=12, ha="left")
        initial_segments = np.repeat((d["pos"][0] * 1e6)[:, None, :], 2, axis=1)
        traces = Line3DCollection(initial_segments, colors=MINT, linewidths=0.85 if comparison else 1.25,
                                  alpha=0.55, zorder=3)
        ax.add_collection3d(traces)
        bubbles = ax.scatter(*(d["pos"][0] * 1e6).T, s=30, color=BLUE, edgecolors="#b5dcff",
                             linewidths=0.4 if comparison else 0.65, depthshade=False, zorder=5)
        bubbles.set_clip_on(False)
        traces.set_clip_on(False)
        axes.append((ax, name))
        dynamic.append((ax, name, traces, bubbles))
    handles = [Line2D([], [], marker="o", color="none", markerfacecolor=BLUE,
                      markeredgecolor="#b5dcff", markersize=8, label="Microbubbles (R0 = 1.78 µm)"),
               Line2D([], [], color=MINT, lw=1.7, label="Recent path (0.25 ms)")]
    fig.legend(handles=handles, loc="center", bbox_to_anchor=(0.5, 0.065), ncol=2,
               frameon=False, labelcolor="white", fontsize=15, columnspacing=3)
    fig.text(0.042, 0.02, "Cloud close-up  |  10,000× slow motion  |  Display interpolation only",
             fontsize=12, color="#91a4b8")
    fig.text(0.958, 0.02, "Scoped model: near-contact motion approximated", ha="right",
             fontsize=12, color="#91a4b8")
    fig.canvas.draw()
    screen_checks = {}
    for ax, name, traces, bubbles in dynamic:
        # Orthographic equal-axis projection: screen discs have the actual equilibrium diameter.
        origin = np.mean(BOUNDS_UM, axis=1)
        projected = screen_points(ax, np.vstack([origin, origin + np.eye(3)]))
        singular = np.linalg.svd((projected[1:] - projected[0]).T, compute_uv=False)
        assert np.isclose(singular[0], singular[1], rtol=1e-10)
        px_per_um = singular.mean()
        bubbles.set_sizes((2 * data[name]["radius"] * 1e6 * px_per_um * 72 / fig.dpi) ** 2)
        xy = screen_points(ax, data[name]["x"].reshape(-1, 3) * 1e6)
        lo, hi = xy.min(axis=0), xy.max(axis=0)
        assert np.all(lo > [20, 100]) and np.all(hi < [1900, 935]), (name, lo, hi)
        screen_checks[name] = dict(pixel_bounds_bottom_origin=[lo.tolist(), hi.tolist()],
                                   pixels_per_um=float(px_per_um),
                                   bubble_diameters_px=(2 * data[name]["radius"] * 1e6 * px_per_um).tolist())
        bubbles.set_visible(False)
        traces.set_visible(False)
    timer.set_visible(False)
    fig.canvas.draw()
    background = fig.canvas.copy_from_bbox(fig.bbox)
    for _, _, traces, bubbles in dynamic:
        traces.set_visible(True)
        bubbles.set_visible(True)
    timer.set_visible(True)

    def render(frame, physical_time):
        fig.canvas.restore_region(background)
        timer.set_text(f"t = {physical_time * 1e3:.3f} ms")
        for ax, name, traces, bubbles in dynamic:
            current = data[name]["pos"][frame]
            traces.set_segments(trail_segments(data[name], physical_time, current))
            bubbles._offsets3d = tuple((current * 1e6).T)
            for artist in [traces, bubbles]:
                artist.do_3d_projection()
                ax.draw_artist(artist)
        fig.draw_artist(timer)
        return np.asarray(fig.canvas.buffer_rgba())[:, :, :3].copy()

    return fig, render, screen_checks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "results/stage3/visualization/reference_style_v1")
    parser.add_argument("--font", type=Path, default=Path("/mnt/c/Windows/Fonts/arial.ttf"))
    parser.add_argument("--preview-only", action="store_true")
    args = parser.parse_args()
    out = args.output.resolve()
    if not args.preview_only and list(out.glob("*.mp4")):
        raise FileExistsError("Choose a new output directory to preserve existing videos")
    out.mkdir(parents=True, exist_ok=True)
    (out / "keyframes").mkdir(exist_ok=True)
    (out / "provenance").mkdir(exist_ok=True)
    font_manager.fontManager.addfont(str(args.font))
    plt.rcParams.update({"font.family": "Arial", "font.weight": "normal", "font.size": 14,
                         "axes.labelweight": "normal", "axes.titleweight": "normal"})
    actual_font = Path(font_manager.findfont("Arial", fallback_to_default=False))
    assert sha(actual_font) == sha(args.font)
    fps, moving_frames = 60, 1201
    physical_times = np.concatenate([np.zeros(30), np.linspace(0, 0.002, moving_frames), np.full(90, 0.002)])
    data, before = load_data(physical_times)
    before[str(args.font)] = sha(args.font)
    (out / "provenance/INPUT_FILES_SHA256.json").write_text(json.dumps(before, indent=2) + "\n")
    shutil.copy2(__file__, out / "provenance/render_motion.py")
    selections = [30, 330, 630, 930, 1230]
    videos = {}
    started = time.perf_counter()
    for kind, comparison in [("stage3_C3_motion_3d_60fps", False), ("stage3_C1_C4_comparison_3d_60fps", True)]:
        fig, render, bounds = make_scene(data, comparison)
        for frame in selections:
            imageio.imwrite(out / "keyframes" / f"{kind}_{frame:04d}.png", render(frame, physical_times[frame]))
        if args.preview_only:
            plt.close(fig)
            print(json.dumps(dict(preview=kind, screen_checks=bounds)), flush=True)
            continue
        with imageio.get_writer(out / f"{kind}.mp4", fps=fps, codec="libx264", quality=8,
                                macro_block_size=1, pixelformat="yuv420p",
                                ffmpeg_params=["-movflags", "+faststart"]) as writer:
            for frame, physical_time in enumerate(physical_times):
                writer.append_data(render(frame, physical_time))
                if frame % 180 == 0 or frame == len(physical_times) - 1:
                    print(json.dumps(dict(video=kind, frame=frame, total=len(physical_times),
                                          elapsed_s=round(time.perf_counter() - started, 2))), flush=True)
        plt.close(fig)
        reader = imageio.get_reader(out / f"{kind}.mp4")
        count = reader.count_frames()
        meta = reader.get_meta_data()
        assert count == len(physical_times) and meta["size"] == (1920, 1080) and meta["fps"] == fps
        for frame in [0, len(physical_times) // 2, len(physical_times) - 1]:
            assert reader.get_data(frame).shape == (1080, 1920, 3)
        reader.close()
        videos[kind] = dict(frame_count=count, fps=fps, duration_s=count / fps, resolution=[1920, 1080],
                            sha256=sha(out / f"{kind}.mp4"), screen_checks=bounds)
    if args.preview_only:
        return
    arrays = dict(physical_time_s=physical_times, video_time_s=np.arange(len(physical_times)) / fps)
    with (out / "frame_map.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["scenario", "frame", "video_time_s", "physical_time_s", "left_state", "right_state", "alpha"])
        for name, d in data.items():
            arrays[f"{name}_positions_m"] = d["pos"]
            arrays[f"{name}_particle_id"] = d["ids"]
            arrays[f"{name}_radii_m"] = d["radius"]
            for i, physical_time in enumerate(physical_times):
                writer.writerow([name, i, i / fps, physical_time, d["left"][i], d["right"][i], d["alpha"][i]])
    np.savez_compressed(out / "display_samples.npz", **arrays)
    shutil.copy2(out / "keyframes/stage3_C3_motion_3d_60fps_0630.png", out / "preview.png")
    images = [Image.open(out / "keyframes" / f"stage3_C3_motion_3d_60fps_{f:04d}.png") for f in [30, 630, 1230]]
    sheet = Image.new("RGB", (1440, 270), "black")
    for i, im in enumerate(images):
        sheet.paste(im.resize((480, 270), Image.Resampling.LANCZOS), (i * 480, 0))
        im.close()
    sheet.save(out / "contact_sheet.png")
    changed = [p for p, digest in before.items() if sha(p) != digest]
    assert not changed, changed
    manifest = dict(
        status="PASS", reference_video=str(REFERENCE), reference_video_sha256=sha(REFERENCE),
        renderer_sha256=sha(__file__), font=dict(family="Arial", sha256=sha(actual_font)), videos=videos,
        physical_interval_s=[0, 0.002], playback_slowdown=10000, initial_hold_s=0.5, final_hold_s=1.5,
        camera=dict(elevation_deg=24, azimuth_deg=-68, projection="orthographic", limits_um=BOUNDS_UM.tolist()),
        colors=dict(background="#000000", bubbles=BLUE, trails=MINT, grid=GRID),
        trajectory_scale=1, glyph_radius="Actual equilibrium radius; orthographic physical scale; no radial pulsation",
        particle_filtering=False, smoothing=False, extrapolation=False, new_simulation=False,
        interpolation="Piecewise linear at uniform physical time; exact saved-time reconstruction verified",
        trails="Previous 0.25 ms, all saved vertices plus interpolated boundary points",
        overlap_note="Linear display chords can introduce nanometre-scale overlap; never corrected or used as solver data",
        sources={name: dict(trajectory=d["source"], sha256=sha(ROOT / d["source"]),
                            saved_states=len(d["t"]), particles=len(d["radius"]),
                            min_display_surface_gap_m=d["minimum_display_gap_m"],
                            exact_saved_state_reconstruction=True) for name, d in data.items()},
        physics_scope="Saved Stage 3 engineering trajectories; near-contact approximation; no new physical validation",
        protected_inputs=len(before), changed_preexisting_files=changed, elapsed_s=time.perf_counter() - started,
    )
    (out / "RENDER_MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n")
    write_viewer(out)
    delivery = {str(p.relative_to(out)): dict(bytes=p.stat().st_size, sha256=sha(p))
                for p in sorted(out.rglob("*")) if p.is_file() and p.name != "render.log"}
    (out / "DELIVERY_SHA256.json").write_text(json.dumps(delivery, indent=2) + "\n")
    print(json.dumps(dict(status=manifest["status"], output=str(out), elapsed_s=manifest["elapsed_s"])), flush=True)


if __name__ == "__main__":
    main()
