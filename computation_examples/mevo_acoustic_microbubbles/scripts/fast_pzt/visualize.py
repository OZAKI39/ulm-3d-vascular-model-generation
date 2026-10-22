"""Actual-radius 3D snapshots and a uniformly timed video of one continuous run."""

import argparse
import hashlib
from time import perf_counter
import numpy as np
import pyvista as pv
import imageio_ffmpeg
import matplotlib.pyplot as plt
from .common import BASE, write
from .analyze import FIG, ON, C, read, style, finish


def plotter(x, radius, size, white=False):
    # All geometry in micrometers, using one common isotropic unit conversion.
    p = pv.Plotter(off_screen=True, window_size=size, lighting="three lights")
    p.set_background("white" if white else "#070b12")
    low, high = x.min(axis=(0, 1)) - 2 * radius, x.max(axis=(0, 1)) + 2 * radius
    center = (low + high) / 2
    extent = float(np.max(high - low))
    bounds = [v for pair in zip(center - extent / 2, center + extent / 2) for v in pair]
    actors = []
    sphere = pv.Sphere(radius=radius, theta_resolution=28, phi_resolution=24)
    for point in x[0]:
        actor = p.add_mesh(
            sphere.copy(),
            color="#246dc9" if white else "#3288e8",
            smooth_shading=True,
            ambient=0.22,
            diffuse=0.78,
            specular=0.5,
            specular_power=28,
            show_scalar_bar=False,
        )
        actor.position = point
        actors.append(actor)
    p.camera_position = [center + extent * np.array([1.4, -1.8, 1.05]), center, (0, 0, 1)]
    p.enable_parallel_projection()
    p.camera.parallel_scale = extent * (0.88 if white else 1.0)
    p.show_bounds(
        bounds=bounds,
        grid="back",
        location="outer",
        ticks="outside",
        n_xlabels=3,
        n_ylabels=3,
        n_zlabels=3,
        font_size=18 if white else 16,
        color="#353e48" if white else "#9aa9be",
        all_edges=False,
        xtitle="x (um)",
        ytitle="y (um)",
        ztitle="z (um)",
        use_2d=False,
        bold=False,
        padding=0,
    )
    p.render()
    return p, actors


def formation():
    d = read(ON)
    x = d["positions_m"][:, 0] * 1e6
    radius = d["radii_m"][0] * 1e6
    idx = [int(np.argmin(abs(d["time_s"] - t))) for t in [0, 0.01, 0.05]]
    p, actors = plotter(x, radius, (900, 780), white=True)
    panels = []
    for i in idx:
        for a, point in zip(actors, x[i]):
            a.position = point
        p.render()
        panels.append(p.screenshot(return_img=True))
    p.close()
    style()
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.5), layout="constrained")
    for ax, panel, i, name in zip(axes, panels, idx, ["Dispersed", "Encounters", "Aggregated"]):
        ax.imshow(panel)
        ax.axis("off")
        ax.set_title(
            f"{name} · {d['time_s'][i] * 1e3:.0f} ms\n"
            f"Components: {d['cluster_count'][i, 0]} | Largest: {d['largest_size'][i, 0]}",
            fontsize=11,
        )
    fig.suptitle(
        "Dispersed-to-cluster formation | N = 32 | R = 1.78 µm | actual scale\n"
        "Uniform 1 kPa, 2.5 MHz · EXPLORATORY",
        fontsize=13,
    )
    finish(fig, "Fig1_aggregation_formation")
    np.savez_compressed(
        FIG / "Fig1_data.npz",
        positions_m=d["positions_m"][idx],
        time_s=d["time_s"][idx],
        radii_m=d["radii_m"],
        particle_id=np.arange(32),
        cluster_labels=d["labels"][idx],
    )


def video(preview=False):
    d = read(C)
    times = d["time_s"]
    x = d["positions_m"][:, 0] * 1e6
    np.testing.assert_allclose(np.diff(times), 0.0002, rtol=0, atol=1e-14)
    out = BASE / "visualization"
    out.mkdir(exist_ok=True)
    target = out / "F1_PZT_aggregation_actual_radius_60fps.mp4"
    p, actors = plotter(x, d["radii_m"][0] * 1e6, (1600, 1000))
    title = p.add_text("", position="upper_left", font_size=15, color="#e4eaf5", name="title")
    p.add_text(
        "N = 32  |  R = 1.78 um  |  Actual sphere scale\n"
        "EXPLORATORY  |  60 fps  |  0.2 ms physical time / frame",
        position="lower_left",
        font_size=12,
        color="#a4b1c4",
    )
    writer = None
    if not preview:
        writer = imageio_ffmpeg.write_frames(
            str(target),
            (1600, 1000),
            fps=60,
            codec="libx264",
            pix_fmt_in="rgb24",
            pix_fmt_out="yuv420p",
            quality=8,
            macro_block_size=2,
            output_params=["-preset", "medium", "-crf", "18", "-movflags", "+faststart"],
        )
        writer.send(None)
    selected = [0, 395, 400, 650, 1000]
    start = perf_counter()
    for i in selected if preview else range(len(times)):
        for a, point in zip(actors, x[i]):
            a.position = point
        event = d["controls"][int(d["control_index"][i])]
        control_name = "Initial PZT controls" if times[i] < 0.08 else "Changed amplitudes / phases / Z- OFF"
        amp = [event["amplitude_Pa"][k] if event["enabled"][k] else 0 for k in range(4)]
        phase = np.array(event["phase_rad"]) * 180 / np.pi
        title.SetText(
            2,
            "PZT-controlled microbubble aggregation\n"
            f"t = {times[i] * 1e3:6.1f} ms  |  Components: {d['cluster_count'][i, 0]}  |  Largest: {d['largest_size'][i, 0]}\n"
            f"{control_name}\n"
            "Elements: Y+ / Y- / Z+ / Z-\n"
            + "Amplitude (Pa): "
            + " / ".join(f"{a:.0f}" for a in amp)
            + "\n"
            + "Phase (deg): "
            + " / ".join(f"{a:.0f}" for a in phase),
        )
        p.render()
        frame = p.screenshot(return_img=True)
        if writer is not None:
            writer.send(np.ascontiguousarray(frame[:, :, :3]))
        if i in selected:
            from PIL import Image

            Image.fromarray(frame).save(out / f"frame_{i:04d}.png")
        if not preview and i % 200 == 0:
            print(f"Rendered {i}/{len(times)} actual frames; {perf_counter() - start:.1f} s", flush=True)
    if writer is not None:
        writer.close()
    p.close()
    if not preview:
        reader = imageio_ffmpeg.read_frames(str(target), pix_fmt="rgb24")
        metadata = next(reader)
        decoded = sum(1 for _ in reader)
        assert decoded == len(times)
        assert abs(metadata["fps"] - 60) < 1e-8
        np.save(out / "frame_physical_times_s.npy", times)
        write(
            out / "VIDEO_VALIDATION.json",
            dict(
                status="PASS",
                file=target.name,
                sha256=hashlib.sha256(target.read_bytes()).hexdigest(),
                source=str((C / "trajectory.h5").relative_to(BASE)),
                decoded_frames=decoded,
                fps=60,
                pixel_dimensions=[1600, 1000],
                physical_interval_s=0.0002,
                physical_duration_s=float(times[-1]),
                encoded_duration_s=decoded / 60,
                actual_radius_m=float(d["radii_m"][0]),
                interpolation=False,
                cluster_rescaling=False,
                particle_subset=False,
                camera="Fixed parallel projection; identical isotropic conversion for coordinates and radii",
                position_rule="Exactly one saved continuous-run frame per encoded video frame",
                render_wall_s=perf_counter() - start,
            ),
        )


if __name__ == "__main__":
    a = argparse.ArgumentParser()
    a.add_argument("--preview", action="store_true")
    q = a.parse_args()
    FIG.mkdir(exist_ok=True)
    formation()
    video(q.preview)
