"""True-radius same-camera F3 shapes, force-vector panels, and actual-frame comparison MP4."""

import argparse
import hashlib
from time import perf_counter
import numpy as np
import h5py
import pyvista as pv
from PIL import Image
import imageio_ffmpeg
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scripts.fast_pzt_f2.visualize import plotter, snapshots
from scripts.fast_pzt_f2.analyze import read
from .common import BASE, write, js
from .analyze import BRANCHES

FIG = BASE / "figures"
VIS = BASE / "visualization"


def geometry_panels():
    FIG.mkdir(exist_ok=True)
    VIS.mkdir(exist_ok=True)
    data = {n: read(p) for n, p in BRANCHES.items()}
    x = np.stack(
        [data["unchanged"]["positions_m"][0, 0]]
        + [data[n]["positions_m"][-1, 0] for n in ["unchanged", "phase", "amplitude"]]
    )
    images = snapshots(x * 1e6)
    fig, axes = plt.subplots(2, 2, figsize=(10, 9), layout="constrained")
    for ax, img, title in zip(
        axes.ravel(),
        images,
        [
            "a  Mature state · 102.2 ms",
            "b  Unchanged · 222.2 ms",
            "c  Phase only · 222.2 ms",
            "d  Amplitude / enable · 222.2 ms",
        ],
    ):
        ax.imshow(img)
        ax.axis("off")
        ax.set_title(title, fontsize=12)
    fig.suptitle(
        "F3-4 | Actual cluster morphology\n32 independent bubbles · R = 1.78 µm · identical camera and scale · EXPLORATORY",
        fontsize=13,
    )
    for ext in ["png", "pdf"]:
        fig.savefig(FIG / f"Fig_F3_4_actual_morphology.{ext}", dpi=230)
    plt.close(fig)
    np.savez_compressed(
        FIG / "Fig_F3_4_data.npz", positions_m=x, particle_id=np.arange(32), radii_m=np.full(32, 1.78e-6)
    )
    driving = js(BASE / "driving_v3/RESULT.json")
    names = ["unchanged", driving["selected"]["phase"], driving["selected"]["amplitude"]]
    scales = {}
    with h5py.File(BASE / "driving_v3/FORCE_FIELD_DATA.h5") as h:
        points = h["positions_m"][:] * 1e6
        for name in names:
            force = (
                h[name + "/primary_N"][:]
                if name == "unchanged"
                else h[name + "/control_delta_centered_primary_N"][:]
            )
            factor = 4.0 / np.max(np.linalg.norm(force, axis=1))  # output micrometers per newton, arrows only
            p, _ = plotter(np.stack([points, points]), size=(950, 850), white=True)
            p.camera.parallel_scale *= 1.15
            for pos, f in zip(points, force):
                norm = np.linalg.norm(f)
                if norm == 0:
                    continue
                direction = f / norm
                start = pos + 1.78 * direction
                arrow = pv.Arrow(
                    start=start, direction=direction, scale=norm * factor, tip_radius=0.08, shaft_radius=0.025
                )
                p.add_mesh(arrow, color="#ac513b", smooth_shading=True, reset_camera=False)
            p.render()
            p.screenshot(FIG / f"force_{name}.png")
            p.close()
            scales[name] = dict(
                arrow_micrometers_per_pN=factor * 1e-12,
                arrow_origin="Bubble surface in force direction",
                sphere_radius_um=1.78,
                coordinate_scale="Same geometry/camera in all three panels",
                kind="Primary" if name == "unchanged" else "Centered control-induced Primary difference",
            )
    write(FIG / "VECTOR_SCALES.json", scales)


def video(preview=False):
    data = {n: read(p) for n, p in BRANCHES.items()}
    names = ["unchanged", "phase", "amplitude"]
    t = data["unchanged"]["time_s"]
    x = np.stack([data[n]["positions_m"][:, 0] for n in names]) * 1e6
    for n in names:
        np.testing.assert_array_equal(data[n]["time_s"], t)
    low = x.min(axis=(0, 1, 2)) - 3.56
    high = x.max(axis=(0, 1, 2)) + 3.56
    center = (low + high) / 2
    extent = float(np.max(high - low))
    bounds = np.array([center - extent / 2, center + extent / 2]).T.ravel()
    size = (1800, 800)
    p = pv.Plotter(off_screen=True, shape=(1, 3), window_size=size, lighting="three lights", border=False)
    actors = []
    texts = []
    sphere = pv.Sphere(radius=1.78, theta_resolution=28, phi_resolution=24)
    selection = js(BASE / "driving_v3/SELECTED.json")
    for k, (name, title) in enumerate(
        zip(names, ["Unchanged PZT", "Phase-only control", "Amplitude / enable control"])
    ):
        p.subplot(0, k)
        p.set_background("#070b12")
        aa = []
        for point in x[k, 0]:
            a = p.add_mesh(
                sphere,
                color="#3288e8",
                smooth_shading=True,
                ambient=0.22,
                diffuse=0.78,
                specular=0.5,
                specular_power=28,
            )
            a.position = point
            aa.append(a)
        actors.append(aa)
        p.camera_position = [center + extent * np.array([1.4, -1.8, 1.05]), center, (0, 0, 1)]
        p.enable_parallel_projection()
        p.camera.parallel_scale = extent * 0.92
        p.show_bounds(
            bounds=bounds,
            grid="back",
            location="outer",
            ticks="outside",
            n_xlabels=3,
            n_ylabels=3,
            n_zlabels=3,
            font_size=12,
            color="#aab6c6",
            all_edges=False,
            xtitle="x (um)",
            ytitle="y (um)",
            ztitle="z (um)",
            use_2d=False,
            bold=False,
            padding=0,
        )
        texts.append(p.add_text("", position="upper_left", font_size=12, color="#eef3fb"))
        subtitle = [
            "All four ON, 151.34 Pa each",
            "Phases: 30 / 45.67 / 0 / 0 deg",
            "Y- and Z- ON, 246.26 Pa each",
        ][k]
        p.add_text(
            subtitle + "\nN = 32 | R = 1.78 um | Actual scale\nStable GPU reduction | EXPLORATORY",
            position="lower_left",
            font_size=10,
            color="#b1bdce",
        )
    target = VIS / "F3_mature_PZT_comparison_actual_radius_60fps.mp4"
    writer = None
    if not preview:
        writer = imageio_ffmpeg.write_frames(
            str(target),
            size,
            fps=60,
            codec="libx264",
            pix_fmt_in="rgb24",
            pix_fmt_out="yuv420p",
            macro_block_size=2,
            output_params=["-preset", "medium", "-crf", "19", "-movflags", "+faststart"],
        )
        writer.send(None)
    started = perf_counter()
    chosen = [0, 100, 300, 600]
    for i in chosen if preview else range(len(t)):
        for k, name in enumerate(names):
            for a, point in zip(actors[k], x[k, i]):
                a.position = point
            title = ["Unchanged PZT", "Phase-only control", "Amplitude / enable control"][k]
            texts[k].SetText(
                2,
                f"{title}\nt = {102.2 + t[i] * 1e3:.1f} ms | After switch: {t[i] * 1e3:.1f} ms\nComponents: {data[name]['cluster_count'][i, 0]} | Largest: {data[name]['largest_size'][i, 0]}",
            )
        p.render()
        frame = p.screenshot(return_img=True)
        if writer is not None:
            writer.send(np.ascontiguousarray(frame[:, :, :3]))
        if i in chosen:
            Image.fromarray(frame).save(VIS / f"frame_{i:04d}.png")
        if i % 150 == 0:
            print(f"F3 render {i}/{len(t)} frames; {perf_counter() - started:.1f} s", flush=True)
    if writer is not None:
        writer.close()
    p.close()
    if not preview:
        stream = imageio_ffmpeg.read_frames(str(target))
        meta = next(stream)
        count = sum(1 for _ in stream)
        assert count == len(t) and meta["fps"] == 60
        np.testing.assert_allclose(np.diff(t), 0.0002, rtol=0, atol=1e-14)
        write(
            VIS / "VIDEO_VALIDATION.json",
            dict(
                status="PASS",
                file=target.name,
                sha256=hashlib.sha256(target.read_bytes()).hexdigest(),
                decoded_frames=count,
                fps=60,
                pixel_dimensions=size,
                physical_duration_s=0.12,
                absolute_time_offset_s=0.1022,
                encoded_duration_s=count / 60,
                physical_interval_s=0.0002,
                interpolation=False,
                particle_subset=False,
                cluster_rescaling=False,
                actual_radius_m=1.78e-6,
                camera="Same fixed parallel projection, absolute bounds and isotropic coordinate scaling for all views",
                sources={n: str((path / "trajectory.h5").relative_to(BASE)) for n, path in BRANCHES.items()},
                controls=selection,
                render_wall_s=perf_counter() - started,
            ),
        )
        np.save(VIS / "frame_absolute_time_s.npy", 0.1022 + t)
        print("F3 VIDEO COMPLETE", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--preview", action="store_true")
    a = ap.parse_args()
    geometry_panels()
    video(a.preview)
