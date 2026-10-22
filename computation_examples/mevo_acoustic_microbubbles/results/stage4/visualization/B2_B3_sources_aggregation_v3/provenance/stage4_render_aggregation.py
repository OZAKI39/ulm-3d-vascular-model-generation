"""Actual PZT/field slices and recorded-state microbubble aggregation microscopy."""

import argparse
import csv
import json
from pathlib import Path
import shutil
import subprocess
import time

import h5py
import imageio.v2 as imageio
import imageio_ffmpeg
import matplotlib.colors as colors
from matplotlib import font_manager
from matplotlib.cm import ScalarMappable
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d.art3d import Line3DCollection, Poly3DCollection
import numpy as np
from PIL import Image, ImageDraw, ImageFont, ImageOps
import yaml

from scripts.stage4 import render_motion as base
from mevo_mb.bubbles.metrics import aggregation_metrics

ROOT = base.ROOT
FPS = 60
SAVED_INDICES = list(range(11)) + [15, 20, 30, 40, 50, 60, 70, 80, 100, 125, 150, 200, 250, 300]
FRAME_HOLDS = [120] * 11 + [24] * 13 + [180]
VIDEO_STATE = np.repeat(np.arange(len(SAVED_INDICES)), FRAME_HOLDS)
HALF_UM = 20.0
FIELD_BOUNDS = np.array([[-3.5, 11.5], [-7.0, 7.0], [-5.0, 5.0]])
MICRO_BOUNDS = np.array([[-HALF_UM, HALF_UM]] * 3)
BLUE, MINT = "#339dff", "#66e0c2"


def sphere_faces():
    theta = np.linspace(0, 2 * np.pi, 25)
    phi = np.linspace(0, np.pi, 13)
    p, t = np.meshgrid(phi, theta, indexing="ij")
    v = np.stack([np.sin(p) * np.cos(t), np.sin(p) * np.sin(t), np.cos(p)], axis=-1)
    faces = np.stack([v[:-1, :-1], v[:-1, 1:], v[1:, 1:], v[1:, :-1]], axis=2).reshape(-1, 4, 3)
    normal = faces.mean(axis=1)
    normal /= np.linalg.norm(normal, axis=1)[:, None]
    light = np.array([-0.4, -0.6, 1.0])
    light /= np.linalg.norm(light)
    intensity = 0.40 + 0.60 * np.maximum(0, normal @ light)
    return faces, intensity


def load():
    data, protected = base.load_data(np.array(SAVED_INDICES) * 0.02)
    gold = json.loads(
        (ROOT / "results/stage4/performance_refactor/gold_baseline/GOLD_BASELINE_MANIFEST.json").read_text()
    )
    for name, d in data.items():
        run = ROOT / d["result_id"]
        cfg = yaml.safe_load((run / "config.yaml").read_text())
        metric_path = run / "metrics.json"
        assert base.sha(metric_path) == gold["files"][f"cases/{name}/metrics.json"]["sha256"]
        protected[str(metric_path)] = base.sha(metric_path)
        recorded = json.loads(metric_path.read_text())["rows"]
        field_path = ROOT / cfg["effective_acoustic_field"]["path"]
        assert base.sha(field_path) == cfg["pzt_control"]["field_sha256"]
        protected[str(field_path)] = base.sha(field_path)
        with h5py.File(field_path) as f:
            d["axes"] = [f[f"axes/{axis}_m"][:] * 1000 for axis in "xyz"]
            d["P"] = f["pressure_Pa"][:]
            assert f.attrs["frequency_Hz"] == 2.5e6
        assert d["P"].shape == tuple(len(axis) for axis in d["axes"])
        assert np.abs(d["P"]).max() <= 39000.0001
        d["config"], d["field_path"] = cfg, field_path
        d["anchor"] = int(np.argmin(d["ids"]))  # Fixed identity; never switches between saved states.
        d["states"] = []
        for k in SAVED_INDICES:
            t, x = d["t"][k], d["x"][k]
            assert np.isclose(t, k * 0.02)
            active = t < d["exit_time"]
            active_indices = np.flatnonzero(active)
            stats, labels, _, _ = aggregation_metrics(x[active], d["radii"][active], cfg["metrics"])
            for key in [
                "cluster_count",
                "largest_cluster_size",
                "median_nearest_neighbor_m",
                "near_contact_fraction",
            ]:
                assert np.isclose(stats[key], recorded[k][key], rtol=1e-12, atol=1e-15)
            local = (x - x[d["anchor"]]) * 1e6
            roi = active & np.all(np.abs(local) + d["radii"][:, None] * 1e6 <= HALF_UM, axis=1)
            group = np.zeros(len(x), dtype=bool)
            if active[d["anchor"]]:
                anchor_label = labels[np.flatnonzero(active_indices == d["anchor"])[0]]
                group[active_indices[labels == anchor_label]] = True
            else:
                roi[:] = False
            assert np.all(roi[group]), "Tracked component must fit the microscope view"
            # The co-moving view is a rigid translation, preserving every pair distance.
            delta = local[:, None] - local[None, :]
            original_delta = (x[:, None] - x[None, :]) * 1e6
            np.testing.assert_allclose(delta, original_delta, atol=1e-10, rtol=0)
            distance = np.linalg.norm(x[active][:, None] - x[active][None, :], axis=-1)
            np.fill_diagonal(distance, np.inf)
            d["states"].append(
                dict(
                    t=float(t),
                    source_state=k,
                    active_mask=active,
                    active_particles=int(active.sum()),
                    anchor_active=bool(active[d["anchor"]]),
                    A=int(np.count_nonzero((~active) & (d["outcome"] == 1))),
                    B=int(np.count_nonzero((~active) & (d["outcome"] == 2))),
                    local=local,
                    roi=roi,
                    group=group,
                    stats=stats,
                    anchor_cluster_size=int(group.sum()),
                    visible=int(roi.sum()),
                    minimum_gap_m=float(
                        np.min(distance - d["radii"][active][:, None] - d["radii"][active][None, :])
                    ),
                )
            )
    return data, protected


def axes_style(ax, bounds, labels, zoom=1.0):
    ax.set_proj_type("ortho")
    ax.view_init(elev=24, azim=-68)
    ax.set(xlim=bounds[0], ylim=bounds[1], zlim=bounds[2])
    ax.set_box_aspect(np.ptp(bounds, axis=1), zoom=zoom)
    for axis, label in zip([ax.xaxis, ax.yaxis, ax.zaxis], labels, strict=True):
        axis.set_label_text(label, color="white", fontsize=15)
        axis.labelpad = 4
        axis.line.set_color("#42556c")
        axis.set_pane_color((0, 0, 0, 1))
        axis._axinfo["grid"]["color"] = base.GRID
        axis._axinfo["grid"]["linewidth"] = 0.55
    ax.tick_params(colors="#d5e1eb", labelsize=12, pad=0)


def source_and_field(ax, d):
    norm = colors.Normalize(0, 39)
    cmap = base.plt.get_cmap("magma")
    x, y, z = d["axes"]
    iz, ix = int(np.argmin(abs(z))), int(np.argmin(abs(x)))
    assert abs(z[iz]) < 1e-10 and abs(x[ix]) < 1e-10
    take = (x >= FIELD_BOUNDS[0, 0]) & (x <= FIELD_BOUNDS[0, 1])
    X, Y = np.meshgrid(x[take], y, indexing="ij")
    ax.plot_surface(
        X,
        Y,
        np.zeros_like(X),
        facecolors=cmap(norm(abs(d["P"][take, :, iz]) / 1000)),
        rstride=1,
        cstride=1,
        shade=False,
        linewidth=0,
        antialiased=False,
        alpha=0.82,
        zorder=1,
    )
    Y, Z = np.meshgrid(y, z, indexing="ij")
    ax.plot_surface(
        np.zeros_like(Y),
        Y,
        Z,
        facecolors=cmap(norm(abs(d["P"][ix]) / 1000)),
        rstride=1,
        cstride=1,
        shade=False,
        linewidth=0,
        antialiased=False,
        alpha=0.60,
        zorder=2,
    )
    wall = base.wall_sections(d["surface"])
    ax.add_collection3d(Line3DCollection(wall, colors="#d8e2ed", linewidths=0.7, alpha=0.68, zorder=3))
    for element in d["config"]["pzt"]:
        center = np.array(element["position_m"]) * 1000
        normal = np.array(element["normal"], dtype=float)
        normal /= np.linalg.norm(normal)
        direction = np.eye(3)[int(np.argmin(abs(normal)))]
        u = np.cross(normal, direction)
        u /= np.linalg.norm(u)
        v = np.cross(normal, u)
        a = np.linspace(0, 2 * np.pi, 97)
        rim = center + element["size_m"][0] * 500 * (np.cos(a)[:, None] * u + np.sin(a)[:, None] * v)
        strong = element["amplitude_Pa"] > 10000
        color = "#ffab45" if strong else "#dec09a"
        ax.add_collection3d(
            Poly3DCollection(
                [rim],
                facecolors=color,
                edgecolors=color,
                linewidths=1.3,
                alpha=0.14 if not strong else 0.24,
                zorder=4,
            )
        )
        ax.quiver(
            *center,
            *normal,
            length=0.75,
            normalize=True,
            color=color,
            linewidth=1.0,
            arrow_length_ratio=0.25,
            zorder=6,
        )
        label_position = center + normal * 0.3
        ax.text(
            *label_position,
            element["id"].replace("_plus", "+").replace("_minus", "-").upper(),
            color=color,
            fontsize=15,
            zorder=8,
        )
    return norm, cmap


def scene(name, d):
    fig = base.plt.figure(figsize=(19.2, 10.8), dpi=100, facecolor="black")
    fig.text(0.035, 0.946, f"{name}  |  Acoustic field & bubble aggregation", color="white", fontsize=28)
    timer = fig.text(0.965, 0.946, "", ha="right", color=MINT, fontsize=24)
    fig.text(
        0.035,
        0.902,
        "2.5 MHz  |  Flow + primary + secondary force  |  256 bubbles",
        color="#bed0e0",
        fontsize=17,
    )
    fig.text(0.035, 0.848, "PZT sources + pressure amplitude", color="white", fontsize=21)
    fig.text(0.610, 0.848, "Aggregation microscope", color="white", fontsize=21)
    micro_status = fig.text(0.610, 0.810, "", color="#bed0e0", fontsize=15)
    elements = d["config"]["pzt"]
    for row, pair in enumerate([elements[:2], elements[2:]]):
        text = "     ".join(
            f"{e['id'].replace('_plus', '+').replace('_minus', '-').upper()}: {e['amplitude_Pa'] / 1000:.2f} kPa / {np.rad2deg(e['phase_rad']):.0f}°"
            for e in pair
        )
        fig.text(0.035, 0.812 - row * 0.027, text, color="#dec09a", fontsize=13)
    global_ax = fig.add_axes(
        [0.002, 0.285, 0.566, 0.487], projection="3d", facecolor="black", computed_zorder=False
    )
    axes_style(global_ax, FIELD_BOUNDS, ["X (mm)", "Y (mm)", "Z (mm)"], zoom=1.30)
    global_ax.set_xticks([0, 5, 10])
    global_ax.set_yticks([-5, 0, 5])
    global_ax.set_zticks([-4, 0, 4])
    norm, cmap = source_and_field(global_ax, d)
    pop = global_ax.scatter(
        [], [], [], c=BLUE, edgecolors="#c1e3ff", linewidths=0.22, s=13, depthshade=False, zorder=10
    )
    locator = global_ax.scatter(
        [], [], [], facecolors="none", edgecolors=MINT, linewidths=1.2, s=100, depthshade=False, zorder=11
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
        [0.584, 0.265, 0.394, 0.512], projection="3d", facecolor="black", computed_zorder=False
    )
    axes_style(micro_ax, MICRO_BOUNDS, ["ΔX (µm)", "ΔY (µm)", "ΔZ (µm)"], zoom=1.13)
    for axis in [micro_ax.xaxis, micro_ax.yaxis, micro_ax.zaxis]:
        axis.set_ticks([-10, 0, 10])
    particles = Poly3DCollection([], edgecolors="none", antialiased=False, zsort="average", zorder=8)
    micro_ax.add_collection3d(particles)
    label = fig.text(0.610, 0.218, "", color="white", fontsize=20)
    fig.text(0.610, 0.175, "Equilibrium radius: 1.78 µm  |  True spatial scale", color="#b8cde0", fontsize=14)
    bar_ax = fig.add_axes([0.045, 0.231, 0.448, 0.015])
    colorbar = fig.colorbar(
        ScalarMappable(norm=norm, cmap=cmap), cax=bar_ax, orientation="horizontal", ticks=[0, 10, 20, 30, 39]
    )
    colorbar.ax.tick_params(colors="white", labelsize=12, pad=3)
    colorbar.outline.set_visible(False)
    fig.text(
        0.045,
        0.180,
        "|P| [kPa]  ·  z = 0 and x = 0 slices  ·  Arrows: PZT normals",
        color="#c2d2e1",
        fontsize=14,
    )
    branch_text = fig.text(0.045, 0.137, "", color="#bed0e0", fontsize=14)
    stats_text = fig.text(0.035, 0.105, "", color="white", fontsize=17)
    fig.text(
        0.035,
        0.037,
        "Saved states only  |  Detailed playback: 0–200 ms  |  Spatial clusters, no fusion",
        color="#a9bed2",
        fontsize=13,
    )
    fig.text(
        0.965,
        0.037,
        "Overview markers enlarged; microscope radii unchanged",
        ha="right",
        color="#a9bed2",
        fontsize=12,
    )
    handles = [
        Line2D([], [], color="none", marker="o", markerfacecolor=BLUE, markersize=6, label="Microbubbles"),
        Line2D(
            [], [], color="none", marker="o", markerfacecolor=MINT, markersize=6, label="Tracked bubble #0"
        ),
    ]
    fig.legend(
        handles=handles,
        loc="center",
        bbox_to_anchor=(0.775, 0.135),
        frameon=False,
        ncol=2,
        fontsize=12,
        labelcolor="white",
        columnspacing=1.4,
    )
    dynamic = [pop, locator, exits, particles, timer, label, stats_text, branch_text, micro_status]
    for artist in dynamic:
        artist.set_visible(False)
    fig.canvas.draw()
    checks = {}
    for axis_name, ax in [("field", global_ax), ("microscope", micro_ax)]:
        boxes = {}
        for axis in [ax.xaxis, ax.yaxis, ax.zaxis]:
            bbox = axis.label.get_window_extent(fig.canvas.get_renderer())
            assert bbox.x0 > 3 and bbox.x1 < 1917 and bbox.y0 > 3 and bbox.y1 < 1077, (axis_name, bbox)
            boxes[axis.label.get_text()] = list(bbox.bounds)
        checks[axis_name] = dict(axis_label_bounds=boxes)
    center = np.zeros(3)
    uv = base.screen_points(micro_ax, np.vstack([center, center + np.eye(3)]))
    singular = np.linalg.svd((uv[1:] - uv[0]).T, compute_uv=False)
    assert np.isclose(singular[0], singular[1], rtol=1e-10)
    checks["microscope"]["pixels_per_um"] = float(singular.mean())
    checks["microscope"]["equilibrium_diameter_pixels"] = float(3.56 * singular.mean())
    background = fig.canvas.copy_from_bbox(fig.bbox)
    for artist in dynamic:
        artist.set_visible(True)
    template, intensity = sphere_faces()

    def render(k):
        state = d["states"][k]
        fig.canvas.restore_region(background)
        x = d["x"][state["source_state"]] * 1000
        active = state["active_mask"]
        pop._offsets3d = tuple(x[active].T)
        exits._offsets3d = tuple(x[~active].T)
        locator._offsets3d = tuple(x[[d["anchor"]]].T) if state["anchor_active"] else ([], [], [])
        indices = np.where(state["roi"])[0]
        meshes, facecolors = [], []
        for i in indices:
            meshes.append(template * d["radii"][i] * 1e6 + state["local"][i])
            rgb = np.array(colors.to_rgb(MINT if i == d["anchor"] else BLUE))
            facecolors.append(np.c_[intensity[:, None] * rgb, np.ones(len(intensity))])
        particles.set_visible(bool(len(indices)))
        if len(indices):
            vertices = np.concatenate(meshes)
            particles.set_verts(vertices)
            particles.set_facecolor(np.concatenate(facecolors))
            screen = base.screen_points(micro_ax, vertices.reshape(-1, 3))
            assert screen[:, 0].min() > 1080 and screen[:, 0].max() < 1910
            assert screen[:, 1].min() > 280 and screen[:, 1].max() < 820
        timer.set_text(f"t = {state['t'] * 1000:.0f} ms")
        if state["anchor_active"]:
            label.set_text(
                f"Local group: {state['anchor_cluster_size']} bubbles  |  {state['visible']} in view"
            )
            micro_status.set_text("40 µm view  |  Following bubble #0")
        else:
            event = next(e for e in d["events"] if e["particle_id"] == int(d["ids"][d["anchor"]]))
            label.set_text(f"#0 crossed {event['branch'][-1]} at {event['crossing_time_s']:.3f} s")
            micro_status.set_text("Tracking complete  |  Overview continues")
        assert state["A"] + state["B"] + state["active_particles"] == 256
        branch_text.set_text(
            f"A: {state['A']}    B: {state['B']}    {'Unassigned at cutoff' if state['t'] == 6 else 'Not crossed'}: {state['active_particles']}    ·    Hollow rings: crossings"
        )
        stats = state["stats"]
        stats_text.set_text(
            f"Active: {state['active_particles']} / 256   ·   Clusters: {stats['cluster_count']}   ·   Near contact: {100 * stats['near_contact_fraction']:.1f}%   ·   Median nearest distance: {stats['median_nearest_neighbor_m'] * 1e6:.2f} µm"
        )
        for ax, artist in [(global_ax, exits), (global_ax, pop), (global_ax, locator), (micro_ax, particles)]:
            if artist.get_visible():
                artist.do_3d_projection()
                ax.draw_artist(artist)
        for artist in [timer, label, stats_text, branch_text, micro_status]:
            fig.draw_artist(artist)
        return np.asarray(fig.canvas.buffer_rgba())[:, :, :3].copy()

    return fig, render, checks


def encode(out, name):
    concat = out / f"{name}_frames.txt"
    lines = []
    for k, count in enumerate(FRAME_HOLDS):
        lines.extend([f"file 'keyframes/{name}_{k:02d}.png'", f"duration {count / FPS}"])
    lines.append(f"file 'keyframes/{name}_{len(SAVED_INDICES) - 1:02d}.png'")
    concat.write_text("\n".join(lines) + "\n")
    video = out / f"stage4_{name}_sources_field_aggregation_60fps.mp4"
    cmd = [
        imageio_ffmpeg.get_ffmpeg_exe(),
        "-v",
        "error",
        "-n",
        "-f",
        "concat",
        "-safe",
        "0",
        "-i",
        str(concat),
        "-vf",
        "fps=60",
        "-frames:v",
        str(len(VIDEO_STATE)),
        "-c:v",
        "libx264",
        "-preset",
        "medium",
        "-crf",
        "18",
        "-pix_fmt",
        "yuv420p",
        "-threads",
        "4",
        "-movflags",
        "+faststart",
        str(video),
    ]
    subprocess.run(cmd, check=True)
    reader = imageio.get_reader(video)
    assert reader.count_frames() == len(VIDEO_STATE) and reader.get_meta_data()["fps"] == FPS
    reader.close()
    return video


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--preview-only", action="store_true")
    args = parser.parse_args()
    out = args.output.resolve()
    assert not list(out.glob("*.mp4")), "Use a new output directory to preserve existing exports"
    (out / "keyframes").mkdir(parents=True, exist_ok=True)
    (out / "provenance").mkdir(exist_ok=True)
    font = Path("/mnt/c/Windows/Fonts/arial.ttf")
    font_manager.fontManager.addfont(str(font))
    base.plt.rcParams.update(
        {"font.family": "Arial", "font.weight": "normal", "font.size": 14, "axes.labelweight": "normal"}
    )
    data, protected = load()
    protected[str(font)] = base.sha(font)
    for path in [
        Path(__file__),
        ROOT / "scripts/stage4/render_motion.py",
        ROOT / "scripts/stage3/render_motion.py",
        ROOT / "src/mevo_mb/bubbles/metrics.py",
    ]:
        shutil.copy2(path, out / "provenance" / f"{path.parent.name}_{path.name}")
    checks, videos = {}, {}
    started = time.perf_counter()
    for name, d in data.items():
        fig, render, checks[name] = scene(name, d)
        for k in [0, 2, 10, 19, 24] if args.preview_only else range(len(SAVED_INDICES)):
            imageio.imwrite(out / "keyframes" / f"{name}_{k:02d}.png", render(k))
            print(
                json.dumps(dict(case=name, saved_state=k, elapsed_s=time.perf_counter() - started)),
                flush=True,
            )
        base.plt.close(fig)
        if not args.preview_only:
            video = encode(out, name)
            videos[name] = dict(
                path=video.name,
                sha256=base.sha(video),
                frames=len(VIDEO_STATE),
                fps=60,
                resolution=[1920, 1080],
                duration_s=len(VIDEO_STATE) / FPS,
            )
            print(json.dumps(dict(case=name, video="COMPLETE")), flush=True)
    (out / "SCREEN_CHECKS.json").write_text(json.dumps(checks, indent=2) + "\n")
    if args.preview_only:
        return
    comparison = out / "stage4_B2_B3_sources_field_aggregation_60fps.mp4"
    subprocess.run(
        [
            imageio_ffmpeg.get_ffmpeg_exe(),
            "-v",
            "error",
            "-n",
            "-i",
            str(out / videos["B2"]["path"]),
            "-i",
            str(out / videos["B3"]["path"]),
            "-filter_complex",
            "[0:v][1:v]hstack=inputs=2[v]",
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
            "-movflags",
            "+faststart",
            str(comparison),
        ],
        check=True,
    )
    videos["B2_B3"] = dict(
        path=comparison.name,
        sha256=base.sha(comparison),
        frames=len(VIDEO_STATE),
        fps=60,
        resolution=[3840, 1080],
        duration_s=len(VIDEO_STATE) / FPS,
    )
    print(json.dumps(dict(comparison="COMPLETE")), flush=True)
    with (out / "frame_map.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "case",
                "frame",
                "video_time_s",
                "saved_state",
                "physical_time_s",
                "local_group",
                "all_cluster_count",
                "near_contact_fraction",
                "median_nearest_neighbor_m",
                "original_saved_state",
                "anchor_active",
                "active_particles",
                "A",
                "B",
            ]
        )
        for name, d in data.items():
            for frame, k in enumerate(VIDEO_STATE):
                s = d["states"][k]
                w.writerow(
                    [
                        name,
                        frame,
                        frame / 60,
                        k,
                        s["t"],
                        s["anchor_cluster_size"],
                        s["stats"]["cluster_count"],
                        s["stats"]["near_contact_fraction"],
                        s["stats"]["median_nearest_neighbor_m"],
                        s["source_state"],
                        s["anchor_active"],
                        s["active_particles"],
                        s["A"],
                        s["B"],
                    ]
                )
    assert all(base.sha(path) == digest for path, digest in protected.items())
    manifest = dict(
        status="PASS",
        videos=videos,
        renderer_sha256=base.sha(__file__),
        protected_inputs=protected,
        source_runs={name: d["result_id"] for name, d in data.items()},
        fields={
            name: dict(
                path=str(d["field_path"].relative_to(ROOT)),
                sha256=base.sha(d["field_path"]),
                frequency_Hz=2.5e6,
                display="steady |P| in kPa, fixed 0–39 kPa; full-resolution XY z=0 and YZ x=0 slices within stored ROI",
                source_elements=d["config"]["pzt"],
            )
            for name, d in data.items()
        },
        microscope=dict(
            anchor_selection="Fixed lowest particle ID #0 in both cases; an illustrative local neighborhood, accompanied by full-population statistics",
            half_width_um=HALF_UM,
            transform="Translation only: x_i(t)-x_0(t); equal axis scale",
            radius="Actual equilibrium R0=1.78 um; no glyph enlargement; no radial oscillation simulated",
            temporal_sampling="Original saved states only: 0–200 ms in 20 ms increments held 2 video seconds each, then selected states through 6 s held .4 video seconds, final held 3 seconds. No position interpolation.",
            original_saved_indices=SAVED_INDICES,
            video_frames_per_state=FRAME_HOLDS,
            spatial_scope="All active spheres fully inside the 40 um cube shown; full component around #0 fits while it is active. Microscope ends at its recorded crossing. Overview conserves all 256 as active or hollow crossing markers.",
            definition=dict(cluster_center_distance_m=1e-5, near_contact_surface_gap_m=5e-7),
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
        near_contact_scope="Spatial clustering, not fusion. Existing model validity limitations remain. Twenty-millisecond saved states cannot resolve microsecond contacts; no new trajectory or acoustic computation.",
        screen_checks=checks,
        elapsed_s=time.perf_counter() - started,
    )
    (out / "RENDER_MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n")
    sections = []
    for name in ["B2", "B3"]:
        sections.append(
            f'<h2>{name}</h2><video controls playsinline preload="metadata" poster="keyframes/{name}_10.png" src="{videos[name]["path"]}"></video>'
        )
    html = """<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>B2/B3 声源、声场与微泡聚集</title><style>body{background:#070a0f;color:#eef4ff;font:17px/1.7 Arial,sans-serif;max-width:1600px;margin:30px auto;padding:0 24px}h1,h2{font-weight:400}p{color:#b7cee2}video{width:100%;background:black}a{color:#66e0c2}</style><h1>声源、声场与微泡聚集</h1><p>左侧：实际四阵元位置、直径和法向，正式 B2/B3 复声压场的两个三维截面，统一 0–39 kPa 色标。显示稳态声压幅值 |P|；2.5 MHz 的载波周期为 0.4 µs，不与微泡的毫秒时钟混用。</p><p>右侧：跟随初始编号最小的微泡 #0 的 40 µm 立方视野，仅作坐标平移，微泡按实际平衡半径 1.78 µm 绘制。初始孤立微泡逐步形成可分辨的泡群；近接触不代表融合。下方聚集指标统计当前仍在域内的活动微泡，分流计数始终覆盖全部 256 个微泡。#0 出流后，该局部视图结束，全景继续播放。</p><p>覆盖完整 0–6 s。前 0–200 ms 按原始 20 ms 间隔逐状态慢放，每个状态停留 2 秒；之后抽取原始保存状态加快播放，最终状态停留 3 秒。全片约 30 秒，未插值微泡位置。全程近接触、近壁和有限雷诺数适用性限制仍然保留。</p>"""
    html += "\n".join(sections)
    html += f'<h2>B2 / B3 同步对照（3840 × 1080）</h2><video controls playsinline preload="metadata" src="{comparison.name}"></video>'
    html += '<p><a href="RENDER_MANIFEST.json">来源、半径、群聚定义及逐状态核验</a> · <a href="frame_map.csv">逐帧时间映射</a></p></html>'
    (out / "OPEN_RESULTS.html").write_text(html, encoding="utf-8")
    sheet = Image.new("RGB", (1440, 900), "black")
    draw = ImageDraw.Draw(sheet)
    caption_font = ImageFont.truetype(str(font), 24)
    for row, name in enumerate(["B2", "B3"]):
        for col, k in enumerate([0, 2, 24 if name == "B2" else 17]):
            with Image.open(out / "keyframes" / f"{name}_{k:02d}.png") as im:
                # Show the microscope, not a downsampled full scene.
                crop = ImageOps.contain(im.crop((1270, 325, 1710, 695)), (480, 400), Image.Resampling.LANCZOS)
                sheet.paste(crop, (col * 480 + (480 - crop.width) // 2, row * 450))
                draw.text(
                    (col * 480 + 20, row * 450 + 409),
                    f"{name}  |  {data[name]['states'][k]['t'] * 1000:.0f} ms  |  {data[name]['states'][k]['visible']} bubbles",
                    font=caption_font,
                    fill="white",
                )
    sheet.save(out / "aggregation_contact_sheet.png")
    delivery = {
        str(p.relative_to(out)): dict(bytes=p.stat().st_size, sha256=base.sha(p))
        for p in sorted(out.rglob("*"))
        if p.is_file() and p.name not in ["DELIVERY_SHA256.json", "render.log"]
    }
    (out / "DELIVERY_SHA256.json").write_text(json.dumps(delivery, indent=2) + "\n")
    print(json.dumps(dict(status="PASS", output=str(out))), flush=True)


if __name__ == "__main__":
    main()
