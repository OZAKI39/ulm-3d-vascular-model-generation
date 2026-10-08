"""Actual-frame F5 figures and constant-rate video; no trajectory interpolation."""

import argparse
from time import perf_counter
import numpy as np
import pyvista as pv
from PIL import Image, ImageDraw
import imageio_ffmpeg
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from . import style
from .common import BASE, INPUT, js, write, sha
from .analyze import read

if style.FONT.exists():
    font_manager.fontManager.addfont(str(style.FONT))
plt.rcParams.update(
    {
        "font.family": "Arial" if style.FONT.exists() else "DejaVu Sans",
        "font.size": 11,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.labelcolor": "#263746",
        "text.color": "#263746",
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
    }
)


def save(fig, out, name):
    for ext in ["png", "pdf", "svg"]:
        fig.savefig(out / (name + "." + ext), dpi=300, facecolor=fig.get_facecolor())
    plt.close(fig)
    with Image.open(out / (name + ".png")) as im:
        im.thumbnail((1700, 1000))
        im.save(out / (name + "_preview.png"))


def flow_scene(x):
    p = pv.Plotter(off_screen=True, window_size=(1550, 1050), lighting="three lights")
    style.setup(p, True)
    wall = style.vessel()
    wall.points = wall.points.astype(float) * 1000
    p.add_mesh(wall, color="#b85c64", opacity=0.20, smooth_shading=True, specular=0.25, ambient=0.2)
    rims = wall.extract_feature_edges(
        boundary_edges=True, non_manifold_edges=False, feature_edges=False, manifold_edges=False
    )
    p.add_mesh(rims, color="#b85c64", lighting=False, line_width=1.5)
    grid = pv.read(INPUT / "result_121.vtu")
    grid["Speed (mm/s)"] = np.linalg.norm(grid["Velocity"], axis=1) * 1000
    section = grid.slice(normal=[0, 0, 1], origin=[0, 0, 0])
    section.points = section.points.astype(float) * 1000
    p.add_mesh(
        section,
        scalars="Speed (mm/s)",
        cmap="magma",
        clim=(0, 10.074017175615313),
        show_scalar_bar=False,
        lighting=False,
    )
    points = section.points[::90]
    u = section["Velocity"][::90]
    valid = np.linalg.norm(u, axis=1) > 0.0008
    cloud = pv.PolyData(points[valid])
    cloud["direction"] = u[valid]
    arrows = cloud.glyph(orient="direction", scale=False, factor=0.70)
    p.add_mesh(arrows, color="#d6ecf5", lighting=False)
    # Original stored acoustic support, with no extrapolated pressure shown.
    box = pv.Box(bounds=(-10.2, 10.2, -1.8, 1.8, -1.8, 1.8))
    p.add_mesh(box, style="wireframe", color="#497aab", line_width=1.6, opacity=0.9)
    dots = pv.PolyData(x * 1000).glyph(geom=pv.Sphere(radius=style.R * 1000), scale=False, orient=False)
    p.add_mesh(dots, color=style.BLUE, smooth_shading=True)
    p.add_point_labels(
        [x.mean(0) * 1000],
        ["Initial cloud"],
        font_size=22,
        point_size=5,
        point_color=style.BLUE,
        text_color="#263746",
        shape_opacity=0,
        always_visible=True,
    )
    info, _ = style.reference()
    cam = info["render"]
    p.camera_position = [cam["camera_eye_mm"], cam["camera_focal_mm"], cam["camera_up"]]
    p.enable_parallel_projection()
    p.camera.parallel_scale = 10.7
    segments = []
    for xx in [-14.5, -10, -5, 0, 5, 10, 11.5]:
        segments.extend([[[xx, -7.5, -2], [xx, 7.5, -2]], [[xx, 7.5, -2], [xx, 7.5, 2]]])
    for yy in [-7.5, -4, 0, 4, 7.5]:
        segments.extend([[[-14.5, yy, -2], [11.5, yy, -2]], [[-14.5, yy, -2], [-14.5, yy, 2]]])
    for zz in [-2, 0, 2]:
        segments.extend([[[-14.5, -7.5, zz], [-14.5, 7.5, zz]], [[-14.5, 7.5, zz], [11.5, 7.5, zz]]])
    mesh = pv.PolyData(
        np.array(segments).reshape(-1, 3),
        lines=np.array([[2, 2 * i, 2 * i + 1] for i in range(len(segments))]),
    )
    p.add_mesh(mesh, color="#a8b3be", lighting=False, line_width=0.7)
    p.show(auto_close=False)
    im = Image.fromarray(p.screenshot(return_img=True))
    draw = ImageDraw.Draw(im)
    labels = []
    for xx in [-10, 0, 10]:
        labels.append(([xx, -7.5, -2], str(xx), (0, 20)))
    for yy in [-4, 0, 4]:
        labels.append(([11.5, yy, -2], str(yy), (24, 12)))
    for zz in [-2, 0, 2]:
        labels.append(([-14.5, -7.5, zz], str(zz), (-28, 0)))
    labels.extend(
        [
            ([-1.5, -7.5, -2], "X (mm)", (0, 60)),
            ([11.5, 0, -2], "Y (mm)", (76, 43)),
            ([-14.5, -7.5, 2], "Z (mm)", (16, -45)),
        ]
    )
    for point, label, offset in labels:
        p.renderer.SetWorldPoint(*point, 1)
        p.renderer.WorldToDisplay()
        xx, yy, _ = p.renderer.GetDisplayPoint()
        draw.text(
            (xx + offset[0], im.height - 1 - yy + offset[1]),
            label,
            font=style.font(24),
            fill="#364b5d",
            anchor="mm",
        )
    p.close()
    return im


def figures(out, d, cfg, analysis):
    x = d["x"]
    t = d["time_s"]
    scene = style.local_scene(x, True, True, size=(1000, 900))
    # One common local camera for initial, switch and final positions.
    panels = [style.image_state(scene, x[i]) for i in [0, 600, 1200]]
    write(out / "LOCAL_CAMERA.json", scene[3])
    scene[0].close()
    initial = style.local_scene(x[:1], True, True, size=(850, 900))
    initial_image = style.image_state(initial, x[0])
    initial[0].close()
    fig = plt.figure(figsize=(14, 7.2), facecolor="white")
    fig.text(0.025, 0.963, "F5-1  |  Prescribed flow in the original Y vessel", fontsize=20, va="top")
    ax = fig.add_axes([0.01, 0.12, 0.72, 0.76])
    ax.imshow(flow_scene(x[0]))
    ax.axis("off")
    ax = fig.add_axes([0.735, 0.19, 0.255, 0.65])
    ax.imshow(initial_image)
    ax.axis("off")
    ax.set_title("Initial cloud: 32 bubbles", fontsize=12)
    cax = fig.add_axes([0.13, 0.105, 0.40, 0.024])
    fig.colorbar(
        plt.cm.ScalarMappable(norm=plt.Normalize(0, 10.074017175615313), cmap="magma"),
        cax=cax,
        orientation="horizontal",
        label="CFD speed on Z = 0 plane (mm/s)",
    )
    fig.text(0.71, 0.095, "Blue box: stored acoustic domain", fontsize=10)
    fig.text(0.71, 0.060, "Spheres: R = 1.78 µm; arrows: flow direction", fontsize=10)
    save(fig, out, "Fig_F5_1_flow_original_Y")

    fig = plt.figure(figsize=(15, 6.3), facecolor="white")
    fig.text(
        0.025,
        0.96,
        "F5-2  |  Initial state → Acoustic attraction → Flow-driven release",
        fontsize=19,
        va="top",
    )
    fig.text(0.975, 0.87, "Ultrasound switched OFF at 300 ms", ha="right", fontsize=11)
    names = ["Initial state", "Acoustic ON + Flow ON", "Ultrasound OFF + Flow ON"]
    for k, (im, index, name) in enumerate(zip(panels, [0, 600, 1200], names)):
        left = 0.018 + k * 0.33
        ax = fig.add_axes([left, 0.18, 0.32, 0.63])
        ax.imshow(im)
        ax.axis("off")
        ax.set_title(name, fontsize=12)
        rg = d["metrics"]["global_rg_m"][index, 0] * 1e6
        nc = d["metrics"]["cluster_count"][index, 0]
        largest = d["metrics"]["largest_size"][index, 0]
        delta = np.linalg.norm(x[index].mean(0) - x[0].mean(0)) * 1e6
        fig.text(
            left + 0.16,
            0.145,
            f"t = {t[index] * 1000:.0f} ms" + (" (switch state)" if index == 600 else ""),
            ha="center",
            fontsize=11,
        )
        fig.text(left + 0.16, 0.100, f"Clusters: {nc}     Largest: {largest}/32", ha="center", fontsize=11)
        fig.text(left + 0.16, 0.052, f"Rg = {rg:.2f} µm     |Δxc| = {delta:.2f} µm", ha="center", fontsize=11)
    save(fig, out, "Fig_F5_2_attraction_release")

    fig, axes = plt.subplots(2, 3, figsize=(14.5, 8.1), layout="constrained")
    fig.suptitle("F5-3  |  Release from the same 300 ms state", fontsize=19)
    for name, color, label in [
        ("A", style.BLUE, "Ultrasound OFF + Flow ON"),
        ("B", "#b06c46", "Ultrasound OFF + Flow OFF"),
    ]:
        a = analysis[name]
        ms = a["time_s"] * 1000
        axes[0, 0].plot(ms, a["streamwise_displacement_m"] * 1e6, color=color, label=label)
        axes[0, 0].plot(ms, a["transverse_displacement_m"] * 1e6, color=color, ls=":", alpha=0.9)
        axes[0, 1].plot(ms, a["fixed_rg_m"] * 1e6, color=color)
        for k, ls in enumerate(["-", "--", ":"]):
            axes[0, 2].plot(ms, a["fixed_eigenvalues_m2"][:, k] * 1e12, color=color, ls=ls)
        axes[1, 0].plot(ms, a["aligned_internal_rmsd_m"] * 1e6, color=color)
        axes[1, 1].plot(ms, a["nearest_neighbor_m"].mean(1) * 1e6, color=color)
        axes[1, 2].step(ms, a["largest_size"], where="post", color=color)
        axes[1, 2].step(ms, a["cluster_count"], where="post", color=color, ls=":")
    titles = [
        "Centroid displacement",
        "Fixed-ID radius of gyration",
        "Shape eigenvalues: λ1, λ2, λ3",
        "Internal RMSD after rigid alignment",
        "Mean nearest-neighbor distance",
        "Cluster connectivity",
    ]
    labels = [
        "Displacement (µm)",
        "Rg (µm)",
        "Eigenvalue (µm²)",
        "Internal RMSD (µm)",
        "Distance (µm)",
        "Number of bubbles / clusters",
    ]
    for ax, title, label in zip(axes.flat, titles, labels):
        ax.set(title=title, xlabel="Time after ultrasound OFF (ms)", ylabel=label, xlim=(0, 300))
        ax.grid(alpha=0.17)
    axes[0, 0].legend(fontsize=9, loc="upper left")
    axes[0, 0].text(0.04, 0.62, "Solid: X; dotted: transverse", transform=axes[0, 0].transAxes, fontsize=9)
    axes[0, 2].text(
        0.04, 0.91, "λ1: solid   λ2: dashed   λ3: dotted", transform=axes[0, 2].transAxes, fontsize=9
    )
    axes[1, 2].text(
        0.04, 0.90, "Solid: largest group; dotted: clusters", transform=axes[1, 2].transAxes, fontsize=9
    )
    A = js(BASE / "ON_OFF_FLOW_COMPARISON.json")["A"]["topology_counts"]
    axes[1, 2].text(
        0.04,
        0.79,
        f"Persistent split events: {A['split_persistent']} (≥ 2 ms)",
        transform=axes[1, 2].transAxes,
        fontsize=9,
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=2),
    )
    save(fig, out, "Fig_F5_3_release_comparison")


def video(out, d, cfg, a, preview):
    x, t = d["x"], d["time_s"]
    scene = style.local_scene(x, True, False, size=(920, 820))
    global_images = [
        style.global_image(event, x[0].mean(0), False, size=(650, 790)) for event in cfg["events"]
    ]
    W, H = 1920, 1080

    def frame(index):
        im = Image.new("RGB", (W, H), "black")
        draw = ImageDraw.Draw(im)
        stage = 0 if index == 0 else 1 if index < 600 else 2
        labels = ["Initial state", "Acoustic ON + Flow ON", "Ultrasound OFF + Flow ON"]
        for j, (label, left) in enumerate(zip(labels, [32, 390, 1030])):
            color = "#edf4ff" if j == stage else "#8292a1"
            draw.text((left, 28), label, font=style.font(33), fill=color)
            if j == stage:
                width = draw.textlength(label, font=style.font(33))
                draw.line([(left, 78), (left + width, 78)], fill="#75abd8", width=3)
        for xx in [300, 950]:
            draw.text((xx, 25), "→", font=style.font(35), fill="#91a7ba")
        draw.text((32, 107), f"t = {t[index] * 1000:.1f} ms", font=style.font(30), fill="#dce8f5")
        draw.text((1885, 108), "Ultrasound OFF at 300 ms", font=style.font(25), fill="#c1d0df", anchor="ra")
        im.paste(global_images[int(index >= 600)], (5, 235))
        im.paste(style.image_state(scene, x[index]), (650, 200))
        draw = ImageDraw.Draw(im)
        draw.text((32, 181), "Y vessel and PZT layout", font=style.font(25), fill="#dce8f5")
        draw.text((720, 181), "Local view  ·  R = 1.78 µm", font=style.font(25), fill="#dce8f5")
        left = 1585
        largest = int(a["largest_size"][index])
        nc = int(a["cluster_count"][index])
        lines = [
            "N = 32",
            f"Largest group: {largest}/32",
            f"Clusters: {nc}",
            f"Rg: {a['fixed_rg_m'][index] * 1e6:.2f} µm",
            f"Δx: {a['streamwise_displacement_m'][index] * 1e6:.2f} µm",
        ]
        for n, line in enumerate(lines):
            draw.text((left, 193 + 47 * n), line, font=style.font(24), fill="#dce8f5")
        for top, title, values, ymax, color in [
            (482, "Radius of gyration (µm)", a["fixed_rg_m"] * 1e6, 15, "#4c9cdc"),
            (750, "Largest connected group", a["largest_size"], 32, "#b1c5d7"),
        ]:
            xx0, xx1, yy0, yy1 = left + 32, 1888, top + 45, top + 205
            draw.text((left, top), title, font=style.font(22), fill="#bfcddb")
            draw.line([(xx0, yy0), (xx0, yy1), (xx1, yy1)], fill="#6d7d8d", width=1)
            coords = np.column_stack([xx0 + t / t[-1] * (xx1 - xx0), yy1 - values / ymax * (yy1 - yy0)])
            draw.line([tuple(q) for q in coords], fill="#304354", width=2)
            if index > 0:
                draw.line([tuple(q) for q in coords[: index + 1]], fill=color, width=3)
            for tick in [0, 300, 600]:
                xp = xx0 + tick / 600 * (xx1 - xx0)
                draw.text((xp, yy1 + 18), str(tick), font=style.font(18), fill="#95a8ba", anchor="mm")
            draw.text((xx0 - 8, yy0), str(ymax), font=style.font(18), fill="#95a8ba", anchor="rm")
            draw.text((xx0 - 8, yy1), "0", font=style.font(18), fill="#95a8ba", anchor="rm")
            xp = xx0 + (xx1 - xx0) / 2
            draw.line([(xp, yy0), (xp, yy1)], fill="#536d85", width=1)
            pt = coords[index]
            draw.ellipse([*(pt - 4), *(pt + 4)], fill=color)
            draw.text(
                ((xx0 + xx1) / 2, yy1 + 47), "Time (ms)", font=style.font(18), fill="#95a8ba", anchor="mm"
            )
        return im

    chosen = [0, 200, 599, 600, 900, 1200]
    stem = "Video_F5_acoustic_release_flow_60fps"
    if preview:
        for index in chosen:
            frame(index).save(out / (stem + f"_preview_{index:04d}.png"))
    else:
        start = perf_counter()
        writer = imageio_ffmpeg.write_frames(
            str(out / (stem + ".mp4")),
            (W, H),
            fps=60,
            codec="libx264",
            pix_fmt_in="rgb24",
            pix_fmt_out="yuv420p",
            macro_block_size=1,
            output_params=["-crf", "18", "-preset", "medium", "-movflags", "+faststart"],
        )
        writer.send(None)
        for index in range(len(t)):
            im = frame(index)
            writer.send(np.asarray(im))
            if index in chosen:
                im.save(out / (stem + f"_frame_{index:04d}.png"))
            if index % 200 == 0:
                print(f"F5 video {index}/{len(t)}", flush=True)
        writer.close()
        np.savez_compressed(
            out / "FRAME_MAP.npz",
            video_frame=np.arange(len(t)),
            source_frame=np.arange(len(t)),
            physical_time_s=t,
            particle_id=d["particle_id"],
            radii_m=d["radii_m"],
            ultrasound_on=np.arange(len(t)) < 600,
        )
        write(
            out / "VIDEO_IDENTITY.json",
            dict(
                frames=len(t),
                fps=60,
                size=[W, H],
                switch_frame=600,
                physical_interval_s=0.0005,
                render_wall_s=perf_counter() - start,
                local_scene=scene[3],
                true_radius_m=style.R,
                position_interpolation=False,
                display_position_or_radius_magnification=False,
                camera_fixed=True,
                background_rgb=[0, 0, 0],
                PZT_enabled_on_screen_matches_drive=True,
                actual_positions_SHA256=sha(out.parent / "science_v1/main/trajectory.h5"),
            ),
        )
    scene[0].close()


def main(label, mode):
    source = BASE / label
    out = BASE / "visualization"
    out.mkdir(exist_ok=True)
    assert js(source / "RESULT.json")["status"] == "COMPLETE"
    d = read(source / "main/trajectory.h5")
    cfg = js(source / "main/CONFIG.json")
    assert len(d["time_s"]) == 1201 and d["active"].all()
    np.testing.assert_array_equal(d["particle_id"], np.arange(32))
    analysis = {name: dict(np.load(source / (name + "_ANALYSIS.npz"))) for name in ["main", "A", "B"]}
    if mode == "figures":
        figures(out, d, cfg, analysis)
    else:
        video(out, d, cfg, analysis["main"], mode == "preview")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", default="science_v1")
    ap.add_argument("--mode", choices=["figures", "preview", "video"], required=True)
    args = ap.parse_args()
    main(args.label, args.mode)
