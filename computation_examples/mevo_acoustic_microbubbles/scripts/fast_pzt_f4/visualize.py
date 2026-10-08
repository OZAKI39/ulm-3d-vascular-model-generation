"""Original-mesh, actual-radius scientific figures and uniformly sampled true-state video."""

import argparse
import hashlib
import subprocess
from time import perf_counter
import numpy as np
import pyvista as pv
from PIL import Image
import imageio_ffmpeg
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from .common import BASE, ROOT, write, js
from .analyze import read, R

FIG = BASE / "figures"
VIS = BASE / "visualization"
COLORS = {"parent": "#2669a8", "junction": "#974b79"}
plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.labelcolor": "#222933",
        "text.color": "#222933",
        "xtick.color": "#222933",
        "ytick.color": "#222933",
        "pdf.fonttype": 42,
    }
)


def mesh():
    m = pv.read(BASE / "inputs/labelled_lumen_m.vtp")
    return m, m.extract_cells(m["ModelFaceID"] == 1).extract_surface(algorithm="dataset_surface")


def axes(p, bounds, unit, white, font=20):
    p.show_bounds(
        bounds=bounds,
        grid="back",
        location="outer",
        ticks="outside",
        n_xlabels=3,
        n_ylabels=3,
        n_zlabels=3,
        font_size=font,
        color="#3b4550" if white else "#b5c1d0",
        all_edges=False,
        xtitle="x (" + unit + ")",
        ytitle="y (" + unit + ")",
        ztitle="z (" + unit + ")",
        use_2d=False,
        bold=False,
        padding=0,
    )


def bubbles(p, x, r, white):
    sphere = pv.Sphere(radius=r, theta_resolution=26, phi_resolution=22)
    actors = []
    for pos in x:
        a = p.add_mesh(
            sphere,
            color="#226fae" if white else "#388ddd",
            smooth_shading=True,
            ambient=0.23,
            diffuse=0.72,
            specular=0.5,
            specular_power=30,
            reset_camera=False,
        )
        a.position = pos
        actors.append(a)
    return actors


def global_scene(p, white=True, initial=True):
    whole, wall = mesh()
    wall.points = wall.points.astype(float) * 1000
    p.set_background("white" if white else "#080d15")
    box = np.array(js(BASE / "GEOMETRY_FIELD_OVERLAP.json")["usable_center_bounds_m"]).ravel() * 1000
    inside = wall.clip_box(box, invert=False)
    p.add_mesh(wall, color="#ab7476", opacity=0.14, smooth_shading=True, reset_camera=False)
    p.add_mesh(inside, color="#bd555a", opacity=0.36, smooth_shading=True, reset_camera=False)
    p.add_mesh(
        pv.Box(bounds=box),
        color="#8e6aa9",
        style="wireframe",
        line_width=1.3,
        opacity=0.7,
        reset_camera=False,
    )
    for tag, label in [(2, "INLET"), (3, "OUTLET A"), (4, "OUTLET B")]:
        cap = whole.extract_cells(whole["ModelFaceID"] == tag).extract_surface(algorithm="dataset_surface")
        rim = cap.extract_feature_edges(
            boundary_edges=True, feature_edges=False, manifold_edges=False, non_manifold_edges=False
        )
        rim.points = rim.points.astype(float) * 1000
        p.add_mesh(rim, color="#923e49", line_width=2, reset_camera=False)
        center = np.array(cap.center) * 1000
        p.add_point_labels(
            [center + [0, 0, 1.6]],
            [label],
            font_size=23,
            text_color="#303741" if white else "#d4deeb",
            shape=None,
            show_points=False,
            always_visible=True,
            reset_camera=False,
        )
    if initial:
        for region in ["parent", "junction"]:
            x = np.load(BASE / ("inputs/" + region + "_initial.npz"))["positions_m"].reshape(-1, 3) * 1000
            bubbles(p, x, R * 1000, white)
            center = x.mean(axis=0)
            # ROI location marker explicitly represents a zoom region, not a bubble.
            p.add_mesh(
                pv.Cube(center=center, x_length=0.6, y_length=0.6, z_length=0.6),
                style="wireframe",
                color=COLORS[region],
                line_width=2,
                reset_camera=False,
            )
            p.add_point_labels(
                [center + [0, 0, 1.3]],
                [region.title() + " ROI"],
                font_size=22,
                text_color=COLORS[region],
                shape=None,
                show_points=False,
                always_visible=True,
                reset_camera=False,
            )
    bounds = [-14, 11, -7, 7, -1.8, 1.8]
    axes(p, bounds, "mm", white)
    center = np.array([-1.5, 0, 0])
    p.camera_position = [center + np.array([10, -22, 31]), center, [0, 0, 1]]
    p.enable_parallel_projection()
    p.camera.parallel_scale = 11.3
    return p


def local_scene(p, region, states, white=True, extent=None):
    _, wall = mesh()
    origin = np.array(js(BASE / "GEOMETRY_FIELD_OVERLAP.json")["initialization"][region]["wall_anchor_m"])
    coords = (states - origin) * 1e6
    center = (coords.min(axis=(0, 1)) + coords.max(axis=(0, 1))) / 2
    if extent is None:
        extent = max(42, float(np.ptp(coords, axis=(0, 1)).max() + 12))
    bounds = np.array([center - extent / 2, center + extent / 2]).T.ravel()
    wall.points = (wall.points.astype(float) - origin) * 1e6
    wall = wall.clip_box(bounds, invert=False)
    p.set_background("white" if white else "#080d15")
    p.add_mesh(
        wall,
        color="#c76c70",
        opacity=0.30 if white else 0.24,
        smooth_shading=False,
        specular=0.3,
        ambient=0.3,
        reset_camera=False,
    )
    actors = bubbles(p, coords[0], R * 1e6, white)
    axes(p, bounds, "um", white, font=20)
    p.camera_position = [center + extent * np.array([0.85, -1.1, 1.5]), center, [0, 0, 1]]
    p.enable_parallel_projection()
    p.camera.parallel_scale = extent * 0.93
    return (
        actors,
        origin,
        dict(
            bounds_um=bounds,
            origin_m=origin,
            camera=list(p.camera_position),
            parallel_scale_um=p.camera.parallel_scale,
        ),
    )


def panels(region, states):
    p = pv.Plotter(off_screen=True, window_size=(900, 700), lighting="three lights", border=False)
    aa, origin, meta = local_scene(p, region, states)
    p.enable_anti_aliasing("ssaa")
    images = []
    for state in states:
        for a, x in zip(aa, (state - origin) * 1e6):
            a.position = x
        p.render()
        images.append(p.screenshot(return_img=True))
    p.close()
    return images, meta


def savefig(fig, name):
    for ext in ["png", "pdf"]:
        fig.savefig(FIG / (name + "." + ext), dpi=230, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def figures():
    FIG.mkdir(exist_ok=True)
    meta = {}
    p = pv.Plotter(off_screen=True, window_size=(1500, 900), lighting="three lights", border=False)
    global_scene(p)
    p.enable_anti_aliasing("ssaa")
    p.render()
    global_img = p.screenshot(return_img=True)
    p.close()
    fig = plt.figure(figsize=(13, 9), layout="constrained")
    gs = fig.add_gridspec(2, 2, height_ratios=[1.35, 1])
    ax = fig.add_subplot(gs[0, :])
    ax.imshow(global_img)
    ax.axis("off")
    ax.set_title("a  Original Y lumen, open terminal rims and stored-field intersection", loc="left")
    for j, region in enumerate(["parent", "junction"]):
        x = np.load(BASE / ("inputs/" + region + "_initial.npz"))["positions_m"].reshape(-1, 3)
        imgs, meta[region] = panels(region, x[None])
        ax = fig.add_subplot(gs[1, j])
        ax.imshow(imgs[0])
        ax.axis("off")
        ax.set_title(chr(98 + j) + "  " + region.title() + " initialization | R = 1.78 µm", loc="left")
    fig.suptitle("F4-1 | Original Y vessel and valid acoustic region", fontsize=17)
    fig.supxlabel(
        "Red: WALL within stored field   |   Pale: unsupported vessel   |   Purple: field-domain boundary\nOpen rims are terminals. Distal daughters lie outside the field. Local coordinates: µm from recorded ROI origins.",
        fontsize=10,
    )
    savefig(fig, "Fig_F4_1_geometry_and_field")
    data = {
        r: {w: read(BASE / "aggregation_v1" / (r + "_wall_" + w)) for w in ["off", "on"]}
        for r in ["parent", "junction"]
    }
    fig = plt.figure(figsize=(13, 13), layout="constrained")
    gs = fig.add_gridspec(4, 3, height_ratios=[1.05, 1.05, 0.8, 0.8])
    plot_data = {}
    for k, region in enumerate(["parent", "junction"]):
        d = data[region]["on"]
        indexes = [0, int(np.argmin(abs(d["time_s"] - 0.02))), len(d["time_s"]) - 1]
        # Intermediate state is at the first saved frame with >=16 connected IDs.
        hit = np.flatnonzero(d["largest_size"][:, 0] >= 16)
        if len(hit):
            indexes[1] = int(hit[0])
        states = d["positions_m"][indexes, 0]
        imgs, meta[region + "_aggregation"] = panels(region, states)
        for j, (idx, img) in enumerate(zip(indexes, imgs)):
            ax = fig.add_subplot(gs[k, j])
            ax.imshow(img)
            ax.axis("off")
            ax.set_title(
                f"{region.title()} | t = {d['time_s'][idx] * 1000:.1f} ms; Nmax = {d['largest_size'][idx, 0]}",
                fontsize=11,
            )
        plot_data[region + "_snapshots_m"] = states
        plot_data[region + "_snapshot_time_s"] = d["time_s"][indexes]
    axs = [fig.add_subplot(gs[i, j]) for i in [2, 3] for j in range(3)]
    for region in ["parent", "junction"]:
        for wall in ["off", "on"]:
            d = data[region][wall]
            t = d["time_s"] * 1000
            ls = "-" if wall == "on" else "--"
            label = region.title() + " / " + wall.upper()
            arrays = [
                d["largest_size"][:, 0],
                d["global_rg_m"][:, 0] * 1e6,
                d["wall_gap_m"][:, 0].min(axis=1) * 1e6,
                (d["wall_gap_m"][:, 0] <= 0.01 * R).sum(axis=1),
                d["hard_contact_pair_count"][:, 0],
                d["wall_stats"][:, 0, :, 3].sum(axis=1) * 1e6,
            ]
            for ax, y in zip(axs, arrays):
                ax.plot(t, y, ls, color=COLORS[region], lw=1.1, label=label)
    labels = [
        "Largest cluster size",
        "Radius of gyration (µm)",
        "Minimum WALL gap (µm)",
        "Bubbles at WALL (h ≤ 0.01R)",
        "Hard-contact pairs",
        "Cumulative WALL correction (µm)",
    ]
    for ax, label in zip(axs, labels):
        ax.set(xlabel="Time (ms)", ylabel=label)
        ax.grid(alpha=0.15)
    axs[2].axhline(0, color=".5", lw=0.6)
    axs[0].legend(fontsize=8, frameon=False)
    fig.suptitle(
        "F4-2 | Actual near-wall aggregation\nOriginal WALL geometry • actual sphere radius • fixed cameras within each ROI",
        fontsize=16,
    )
    fig.supxlabel(
        "WALL OFF uses identical initial states with unconfined F3 motion. Negative wall gaps are comparator penetrations.\nContact counts use saved states; cumulative corrections sum all physical-step wall passes. GEOMETRIC_WALL_ONLY / EXPLORATORY",
        fontsize=10,
    )
    savefig(fig, "Fig_F4_2_near_wall_aggregation")
    np.savez_compressed(FIG / "Fig_F4_2_data.npz", **plot_data)
    c = {n: read(BASE / "controls_v1" / ("wall_on_" + n)) for n in ["baseline", "switched"]}
    states = np.stack(
        [
            c["baseline"]["positions_m"][0, 0],
            c["baseline"]["positions_m"][-1, 0],
            c["switched"]["positions_m"][-1, 0],
        ]
    )
    imgs, meta["control"] = panels("parent", states)
    fig = plt.figure(figsize=(13, 10), layout="constrained")
    gs = fig.add_gridspec(3, 3, height_ratios=[1.2, 0.85, 0.85])
    for j, (img, title) in enumerate(
        zip(imgs, ["Same initial state | 600 ms", "Unchanged PZT | 720 ms", "Y+ phase +30° | 720 ms"])
    ):
        ax = fig.add_subplot(gs[0, j])
        ax.imshow(img)
        ax.axis("off")
        ax.set_title(title)
    axs = [fig.add_subplot(gs[i, j]) for i in [1, 2] for j in range(3)]
    for n, color, ls in [("baseline", "#276da7", "-"), ("switched", "#b65a3b", "--")]:
        d = c[n]
        t = d["time_s"] * 1000
        centroid = (d["global_centroid_m"][:, 0] - d["global_centroid_m"][0, 0]) * 1e6
        axs[0].plot(t, np.linalg.norm(centroid, axis=1), ls, color=color, label=n.title())
        axs[1].plot(t, d["global_rg_m"][:, 0] * 1e6, ls, color=color)
        for j in range(3):
            axs[2].plot(t, d["global_eigenvalues_m2"][:, 0, j] * 1e12, ls, color=color, lw=1)
        axs[3].plot(t, d["wall_gap_m"][:, 0].min(axis=1) * 1e6, ls, color=color)
        axs[4].plot(t, d["hard_contact_pair_count"][:, 0], ls, color=color)
    contrast = c["switched"]["global_centroid_m"][:, 0] - c["baseline"]["global_centroid_m"][:, 0]
    axs[5].plot(c["baseline"]["time_s"] * 1000, np.linalg.norm(contrast, axis=1) * 1e9, color="#66518e")
    for ax, label in zip(
        axs,
        [
            "Centroid displacement (µm)",
            "Radius of gyration (µm)",
            "Shape eigenvalues (µm²)",
            "Minimum WALL gap (µm)",
            "Hard-contact pairs",
            "Control centroid contrast (nm)",
        ],
    ):
        ax.set(xlabel="Time after control (ms)", ylabel=label)
        ax.grid(alpha=0.15)
        ax.ticklabel_format(axis="y", useOffset=False)
    axs[0].legend(frameon=False, fontsize=9)
    fig.suptitle("F4-3 | PZT control from the same near-wall cluster", fontsize=17)
    fig.supxlabel(
        "Weak computed response: endpoint centroid contrast 13.9 nm; relative Rg contrast 0.0063%.\nNo further WALL contact; WALL ON/OFF trajectories coincide. Robust shape control is not established. EXPLORATORY.",
        fontsize=10,
    )
    savefig(fig, "Fig_F4_3_PZT_control")
    np.savez_compressed(FIG / "Fig_F4_3_data.npz", positions_m=states, centroid_contrast_m=contrast)
    write(
        FIG / "RENDER_GEOMETRY.json",
        dict(
            regions=meta,
            sphere_radius_m=R,
            wall="Original finite triangles, clipped for local view only; no smoothing/remeshing of physics",
            initial_marker="0.6 mm outline boxes label ROI locations; not bubbles",
            source_mesh_sha256=hashlib.sha256(
                (BASE / "inputs/labelled_lumen_m.vtp").read_bytes()
            ).hexdigest(),
        ),
    )
    print("Three F4 figures rendered", flush=True)


def video(preview=False):
    VIS.mkdir(exist_ok=True)
    a = read(BASE / "aggregation_v1/parent_wall_on")
    b = read(BASE / "controls_v1/wall_on_switched")
    np.testing.assert_array_equal(a["positions_m"][-1], b["positions_m"][0])
    x = np.concatenate([a["positions_m"][:, 0], b["positions_m"][1:, 0]])
    t = np.concatenate([a["time_s"], 0.6 + b["time_s"][1:]])
    nc = np.concatenate([a["largest_size"][:, 0], b["largest_size"][1:, 0]])
    idx = np.arange(0, len(t), 5)
    np.testing.assert_allclose(np.diff(t[idx]), 0.001, rtol=0, atol=1e-14)
    size = (1920, 1080)
    p = pv.Plotter(off_screen=True, shape=(1, 2), window_size=size, lighting="three lights", border=False)
    p.subplot(0, 0)
    global_scene(p, white=False, initial=False)
    p.camera.parallel_scale = 20.0
    ga = bubbles(p, x[0] * 1000, R * 1000, False)
    roi = x[0].mean(axis=0) * 1000
    p.add_mesh(
        pv.Cube(center=roi, x_length=0.6, y_length=0.6, z_length=0.6),
        style="wireframe",
        color="#87aac8",
        line_width=2,
        reset_camera=False,
    )
    p.add_text(
        "Original Y vessel, global view\nRed: supported WALL; purple: stored field\nOutline: location of the local view",
        position=(20, 970),
        font_size=15,
        color="#e0e9f5",
    )
    p.add_text(
        "GEOMETRIC_WALL_ONLY / EXPLORATORY\nNo flow; R = 1.78 um; N = 32\nDaughters extend beyond the stored field",
        position=(20, 20),
        font_size=14,
        color="#b5c1d0",
    )
    p.subplot(0, 1)
    aa, origin, meta = local_scene(p, "parent", x, white=False)
    p.camera.parallel_scale *= 1.12
    meta["parallel_scale_um"] = p.camera.parallel_scale
    title = p.add_text("", position=(20, 970), font_size=18, color="#edf3fd")
    p.add_text(
        "Local view, actual sphere radius\nSame physical time as global view\n1 ms physical time per video frame",
        position=(20, 20),
        font_size=14,
        color="#bac6d8",
    )
    p.enable_anti_aliasing("ssaa")
    p.render()
    if preview:
        for i in [0, 3000, len(t) - 1]:
            for actor, pos in zip(aa, (x[i] - origin) * 1e6):
                actor.position = pos
            for actor, pos in zip(ga, x[i] * 1000):
                actor.position = pos
            title.SetInput(
                f"t = {t[i] * 1000:.1f} ms; Nmax = {nc[i]}\n"
                + ("Y+ phase +30 deg" if t[i] >= 0.6 else "Fixed PZT control"),
            )
            p.render()
            p.screenshot(VIS / f"preview_{i:04d}.png")
        p.close()
        return
    output = VIS / "F4_Y_vessel_aggregation_and_PZT_control_60fps.mp4"
    writer = imageio_ffmpeg.write_frames(
        str(output),
        size,
        fps=60,
        codec="libx264",
        pix_fmt_in="rgb24",
        pix_fmt_out="yuv420p",
        output_params=["-crf", "20", "-preset", "medium", "-movflags", "+faststart"],
        macro_block_size=1,
    )
    writer.send(None)
    start = perf_counter()
    for frame, i in enumerate(idx):
        for actor, pos in zip(aa, (x[i] - origin) * 1e6):
            actor.position = pos
        for actor, pos in zip(ga, x[i] * 1000):
            actor.position = pos
        title.SetInput(
            f"t = {t[i] * 1000:.1f} ms; Nmax = {nc[i]}\n"
            + ("Y+ phase +30 deg" if t[i] >= 0.6 else "Fixed PZT control"),
        )
        p.render()
        img = p.screenshot(return_img=True)
        writer.send(np.ascontiguousarray(img[:, :, :3]))
        if frame in [0, 300, 600, 720]:
            Image.fromarray(img).save(VIS / f"frame_{frame:04d}.jpg", quality=90)
        if frame % 120 == 0:
            print(f"F4 video {frame}/{len(idx)}", flush=True)
    writer.close()
    p.close()
    sources = ["aggregation_v1/parent_wall_on", "controls_v1/wall_on_switched"]
    write(
        VIS / "VIDEO_IDENTITY.json",
        dict(
            sources={
                n: hashlib.sha256((BASE / n / "trajectory.h5").read_bytes()).hexdigest() for n in sources
            },
            frame_count=len(idx),
            fps=60,
            physical_interval_s=0.001,
            duration_s=len(idx) / 60,
            render_wall_s=perf_counter() - start,
            radius_m=R,
            local_geometry=meta,
            camera_motion=False,
            position_interpolation=False,
            particle_scaling=False,
            time_mapping="Saved state every fifth frame; duplicate continuation initial state removed; control switches at physical 600 ms",
        ),
    )
    np.savez_compressed(
        VIS / "FRAME_MAP.npz",
        video_frame=np.arange(len(idx)),
        physical_time_s=t[idx],
        combined_source_index=idx,
        source_run=np.where(idx <= 3000, 0, 1),
        source_frame_index=np.where(idx <= 3000, idx, idx - 3000),
        particle_id=np.arange(32),
    )
    print("F4 video complete", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--video", action="store_true")
    ap.add_argument("--preview", action="store_true")
    args = ap.parse_args()
    if args.video or args.preview:
        video(args.preview)
    else:
        figures()
    write(
        BASE / "OFFLINE_RENDER_IDENTITY.json",
        dict(
            commit=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
            source_sha256={
                str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                for p in [ROOT / "scripts/fast_pzt_f4/visualize.py", ROOT / "scripts/fast_pzt_f4/analyze.py"]
            },
            pyvista=pv.__version__,
            vtk=pv.vtk_version_info,
            scope="Offline rendering only; original execution identities remain in each run",
        ),
    )
