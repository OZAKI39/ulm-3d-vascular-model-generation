"""Clean black-background rendering of one continuous Secondary ON/OFF run."""

import argparse
from pathlib import Path
import shutil
from time import perf_counter

import h5py
import imageio_ffmpeg
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
import numpy as np
from PIL import Image, ImageDraw

from . import visual_mechanism_stories as style

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/stage4/visualization/aggregation_sequence_v2_20261008"
RUN = OUT / "data/run"
STEM = "Video_attraction_initial_on_off_60fps"
style.BG = "#000000"


def inputs():
    data = style.read(RUN)
    cfg = style.js(RUN / "CONFIG.json")
    summary = style.js(RUN / "SUMMARY.json")
    assert summary["numerical_quality"] == "PASS"
    assert len(data["time"]) == 1201
    assert cfg["events"][1]["time_s"] == 0.3
    np.testing.assert_array_equal(data["id"], np.arange(32))
    assert np.isfinite(data["x"]).all() and (data["r"] == style.R).all()
    with h5py.File(RUN / "trajectory.h5") as h:
        np.testing.assert_array_equal(h["secondary_force_enabled"][:], np.arange(1201) < 600)
        assert h["active"][:].all()
    return data, cfg


def figure(data, cfg):
    scene = style.local_scene(data["x"], True, False, size=(800, 900))
    indices = [0, 599, 1200]
    images = [style.image_state(scene, data["x"][i]) for i in indices]
    scene[0].close()
    global_img = style.global_image(cfg["events"][0], data["x"][0].mean(0), False)
    fig = plt.figure(figsize=(16, 7), facecolor="black")
    fig.text(
        0.025, 0.96, "Initial state → Attraction ON → Attraction OFF", color="#edf3fa", fontsize=23, va="top"
    )
    fig.text(0.975, 0.905, "Attraction switched OFF at 300 ms", color="#b7c9dc", fontsize=11, ha="right")
    ax = fig.add_axes([0.012, 0.10, 0.28, 0.72])
    ax.imshow(global_img)
    ax.axis("off")
    ax.set_title("PZT layout", color="#dbe6f3", fontsize=14, pad=12)
    for j, (idx, im, label) in enumerate(
        zip(indices, images, ["Initial state", "Attraction ON", "Attraction OFF"])
    ):
        x = 0.305 + j * 0.227
        ax = fig.add_axes([x, 0.09, 0.215, 0.73])
        ax.imshow(im)
        ax.axis("off")
        ax.set_title(label, color="#dbe6f3", fontsize=14, pad=12)
        fig.text(
            x + 0.1075,
            0.076,
            f"t = {data['time'][idx] * 1000:.1f} ms",
            color="#dbe6f3",
            ha="center",
            fontsize=11,
        )
        fig.text(
            x + 0.1075,
            0.036,
            f"Largest cluster: {data['nmax'][idx]}/32",
            color="#dbe6f3",
            ha="center",
            fontsize=11,
        )
    for x in [0.521, 0.748]:
        fig.add_artist(
            FancyArrowPatch(
                (x, 0.48),
                (x + 0.012, 0.48),
                transform=fig.transFigure,
                arrowstyle="->",
                mutation_scale=15,
                lw=1.2,
                color="#b7c9dc",
            )
        )
    for ext in ["png", "pdf", "svg"]:
        fig.savefig(OUT / ("Fig_attraction_initial_on_off." + ext), dpi=300, facecolor="black")
    plt.close(fig)
    with Image.open(OUT / "Fig_attraction_initial_on_off.png") as im:
        im.thumbnail((1600, 900))
        im.save(OUT / "Fig_attraction_initial_on_off_preview.png")
    style.write(
        OUT / "FIGURE_IDENTITY.json",
        dict(indices=indices, local_scene=scene[3], trajectory_sha256=style.sha(RUN / "trajectory.h5")),
    )


def video(data, cfg, preview=False):
    scene = style.local_scene(data["x"], True, False, size=(1170, 870))
    global_img = style.global_image(cfg["events"][0], data["x"][0].mean(0), False, size=(660, 810))

    def frame(index):
        im = Image.new("RGB", (1920, 1080), "black")
        d = ImageDraw.Draw(im)
        current = 0 if index == 0 else 1 if index < 600 else 2
        labels = ["Initial state", "Attraction ON", "Attraction OFF"]
        positions = [30, 425, 900]
        for n, (label, x) in enumerate(zip(labels, positions)):
            d.text((x, 28), label, fill="#ecf4ff" if n == current else "#8593a2", font=style.font(38))
            if n == current:
                length = d.textlength(label, font=style.font(38))
                d.line([(x, 80), (x + length, 80)], fill="#7eadd7", width=3)
        for x in [340, 800]:
            d.text((x, 28), "→", font=style.font(38), fill="#a5b4c5")
        d.text((30, 108), f"t = {data['time'][index] * 1000:.1f} ms", font=style.font(30), fill="#dbe6f3")
        d.text((1890, 32), "OFF at 300 ms", font=style.font(25), fill="#aebfd2", anchor="ra")
        d.text((30, 164), "PZT layout", font=style.font(28), fill="#dbe6f3")
        d.text(
            (1890, 105),
            f"N = 32    |    Largest cluster: {data['nmax'][index]}/32",
            font=style.font(27),
            fill="#dbe6f3",
            anchor="ra",
        )
        im.paste(global_img, (15, 210))
        im.paste(style.image_state(scene, data["x"][index]), (720, 180))
        return im

    selected = [0, 200, 599, 600, 1200]
    if preview:
        for idx in selected:
            frame(idx).save(OUT / (STEM + f"_preview_{idx:04d}.png"))
        scene[0].close()
        return
    started = perf_counter()
    writer = imageio_ffmpeg.write_frames(
        str(OUT / (STEM + ".mp4")),
        (1920, 1080),
        fps=60,
        codec="libx264",
        pix_fmt_in="rgb24",
        pix_fmt_out="yuv420p",
        macro_block_size=1,
        output_params=["-crf", "18", "-preset", "medium", "-movflags", "+faststart"],
    )
    writer.send(None)
    for idx in range(1201):
        im = frame(idx)
        writer.send(np.asarray(im))
        if idx in selected:
            im.save(OUT / (STEM + f"_frame_{idx:04d}.png"))
        if idx % 200 == 0:
            print(f"video {idx}/1201", flush=True)
    writer.close()
    scene[0].close()
    np.savez_compressed(
        OUT / "FRAME_MAP.npz",
        video_frame=np.arange(1201),
        source_frame=np.arange(1201),
        time_s=data["time"],
        secondary_enabled=np.arange(1201) < 600,
        particle_id=data["id"],
    )
    style.write(
        OUT / "VIDEO_IDENTITY.json",
        dict(
            frames=1201,
            fps=60,
            duration_s=1201 / 60,
            switch_frame=600,
            switch_time_s=0.3,
            physical_interval_s=0.0005,
            trajectory_sha256=style.sha(RUN / "trajectory.h5"),
            local_scene=scene[3],
            background_rgb=[0, 0, 0],
            position_interpolation=False,
            bubble_radius_scale=1,
            cluster_proxy=False,
            camera_tracking=False,
            render_wall_s=perf_counter() - started,
        ),
    )
    print("Video complete", flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--preview", action="store_true")
    parser.add_argument("--video-only", action="store_true")
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    data, cfg = inputs()
    if not args.video_only:
        figure(data, cfg)
    video(data, cfg, args.preview)
    (OUT / "source").mkdir(exist_ok=True)
    for path in [Path(__file__), Path(style.__file__), ROOT / "scripts/stage4/visual_sequential_ablation.py"]:
        shutil.copy2(path, OUT / "source" / path.name)
    style.write(
        OUT / "RENDER_SOURCE.json",
        dict(
            renderer=style.sha(__file__),
            style_dependency=style.sha(style.__file__),
            reference_png=style.sha(style.REF / "Fig1a_PZT_control.png"),
        ),
    )


if __name__ == "__main__":
    main()
