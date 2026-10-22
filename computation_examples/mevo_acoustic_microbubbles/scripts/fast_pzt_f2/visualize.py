"""Actual-radius fixed-camera F2 panels and video; one observed state per frame."""

import argparse
import hashlib
from time import perf_counter
import numpy as np
import h5py
import pyvista as pv
import imageio_ffmpeg
from PIL import Image
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from .common import BASE, write
from .analyze import read, js, BRANCHES

FIG = BASE / "figures"
VIS = BASE / "visualization"


def plotter(x, size=(1000, 900), white=False):
    p = pv.Plotter(off_screen=True, window_size=size, lighting="three lights")
    p.set_background("white" if white else "#070b12")
    low = x.min(axis=(0, 1)) - 3.56
    high = x.max(axis=(0, 1)) + 3.56
    center = (low + high) / 2
    extent = float(np.max(high - low))
    bounds = np.array([center - extent / 2, center + extent / 2]).T.ravel()
    sphere = pv.Sphere(radius=1.78, theta_resolution=28, phi_resolution=24)
    actors = []
    for pos in x[0]:
        a = p.add_mesh(
            sphere,
            color="#286db8" if white else "#3288e8",
            smooth_shading=True,
            ambient=0.22,
            diffuse=0.78,
            specular=0.5,
            specular_power=28,
            show_scalar_bar=False,
        )
        a.position = pos
        actors.append(a)
    p.camera_position = [center + extent * np.array([1.4, -1.8, 1.05]), center, (0, 0, 1)]
    p.enable_parallel_projection()
    p.camera.parallel_scale = extent * 0.86
    p.show_bounds(
        bounds=bounds,
        grid="back",
        location="outer",
        ticks="outside",
        n_xlabels=3,
        n_ylabels=3,
        n_zlabels=3,
        font_size=17 if white else 14,
        color="#353e48" if white else "#b0bbcb",
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


def snapshots(x):
    p, actors = plotter(x, white=True)
    panels = []
    for state in x:
        for a, pos in zip(actors, state):
            a.position = pos
        p.render()
        panels.append(p.screenshot(return_img=True))
    p.close()
    return panels


def comparison():
    FIG.mkdir(exist_ok=True)
    d = read(BASE / "mature_v1/formation")
    idx = js(BASE / "mature_v1/MATURITY.json")["frame_index"]
    x = np.stack(
        [
            d["positions_m"][0, 0],
            d["positions_m"][idx, 0],
            read(BRANCHES["amplitude"])["positions_m"][-1, 0],
            read(BRANCHES["unchanged"])["positions_m"][-1, 0],
        ]
    )
    panels = snapshots(x * 1e6)
    fig, axes = plt.subplots(2, 2, figsize=(10, 9), layout="constrained")
    for ax, img, title in zip(
        axes.ravel(),
        panels,
        [
            "a  Dispersed · 0 ms",
            "b  Mature · 102.2 ms",
            "c  Amplitude / enable change · 222.2 ms",
            "d  Unchanged control · 222.2 ms",
        ],
    ):
        ax.imshow(img)
        ax.axis("off")
        ax.set_title(title, fontsize=12)
    fig.suptitle(
        "F2-4 | Actual 3D cluster shapes\n32 bubbles · R = 1.78 µm · identical camera and coordinate scale · EXPLORATORY",
        fontsize=13,
    )
    for ext in ["png", "pdf"]:
        fig.savefig(FIG / f"Fig_F2_4_3D_cluster_shapes.{ext}", dpi=220)
    plt.close(fig)
    np.savez_compressed(
        FIG / "Fig_F2_4_data.npz",
        positions_m=x,
        particle_id=np.arange(32),
        radii_m=np.full(32, 1.78e-6),
        time_s=[0, 0.1022, 0.2222, 0.2222],
    )
    contact = np.stack(
        [read(BASE / "sensitivity_v1" / n)["positions_m"][-1, 0] for n in ["near_5us", "contact_5us"]]
    )
    # Both contact panels share absolute coordinate bounds, projection and sphere radius.
    for name, img in zip(["contact_F1_PBD", "contact_F2_candidate"], snapshots(contact * 1e6)):
        Image.fromarray(img).save(FIG / f"{name}.png")


def video(preview=False):
    with h5py.File(VIS / "continuous_actual_states.h5") as h:
        x = h["positions_m"][:, 0] * 1e6
        t = h["time_s"][:]
        nc = h["cluster_count"][:, 0]
        nm = h["largest_size"][:, 0]
    idx = js(BASE / "mature_v1/MATURITY.json")["frame_index"]
    size = (1200, 900)
    p, actors = plotter(x, size=size)
    title = p.add_text("", position="upper_left", font_size=14, color="#eef3fb")
    p.add_text(
        "N = 32  |  R = 1.78 um  |  Actual sphere scale  |  EXPLORATORY\n60 fps  |  0.2 ms physical time per frame  |  Fixed camera",
        position="lower_left",
        font_size=11,
        color="#b0bbcb",
    )
    target = VIS / "F2_PZT_mature_cluster_actual_radius_60fps.mp4"
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
    start = perf_counter()
    chosen = [0, idx, len(t) - 1]
    for i in chosen if preview else range(len(t)):
        for a, pos in zip(actors, x[i]):
            a.position = pos
        stage = "Fixed PZT: cluster formation" if i < idx else "Mature cluster: amplitude / enable switched"
        control = (
            "Y+ / Y- / Z+ / Z-: 151 / 151 / 151 / 151 Pa"
            if i < idx
            else "Y+ / Y- / Z+ / Z-: OFF / 246 / OFF / 246 Pa"
        )
        title.SetText(
            2,
            f"PZT-controlled microbubble aggregation\nt = {t[i] * 1e3:5.1f} ms  |  Components: {nc[i]}  |  Largest: {nm[i]}\n{stage}\n{control}  |  All phases: 0 deg",
        )
        p.render()
        frame = p.screenshot(return_img=True)
        if writer is not None:
            writer.send(np.ascontiguousarray(frame[:, :, :3]))
        if i in chosen:
            Image.fromarray(frame).save(VIS / f"frame_{i:04d}.png")
        if not preview and i % 200 == 0:
            print(f"Rendered {i}/{len(t)} actual states; {perf_counter() - start:.1f} s", flush=True)
    if writer is not None:
        writer.close()
    p.close()
    if not preview:
        reader = imageio_ffmpeg.read_frames(str(target), pix_fmt="rgb24")
        meta = next(reader)
        count = sum(1 for _ in reader)
        assert count == len(t) and meta["fps"] == 60
        np.testing.assert_allclose(np.diff(t), 0.0002, rtol=0, atol=1e-14)
        np.save(VIS / "frame_physical_times_s.npy", t)
        write(
            VIS / "VIDEO_VALIDATION.json",
            dict(
                status="PASS",
                file=target.name,
                sha256=hashlib.sha256(target.read_bytes()).hexdigest(),
                decoded_frames=count,
                fps=60,
                pixel_dimensions=size,
                physical_interval_s=0.0002,
                physical_duration_s=float(t[-1]),
                encoded_duration_s=count / 60,
                actual_radius_m=1.78e-6,
                interpolation=False,
                cluster_rescaling=False,
                particle_subset=False,
                camera="Fixed parallel projection; identical isotropic coordinate/radius conversion",
                source="continuous_actual_states.h5",
                position_rule="Actual common prefix and exact-state continuation; restart equality verified; all 32 IDs",
                render_wall_s=perf_counter() - start,
            ),
        )
        print("VIDEO COMPLETE", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--preview", action="store_true")
    a = ap.parse_args()
    comparison()
    video(a.preview)
