"""Finite physical spheres, real WALL, fixed camera and actual last-state sampling."""

import argparse
import json
from pathlib import Path
import imageio.v2 as imageio
import numpy as np
import pyvista as pv
import vtk
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scripts.stage5a_r2.compare import actual_sample
from scripts.stage5a_r2.validation_figures import local
from scripts.stage5a_r2.seal import ROOT, BASE
from scripts.stage4l.dense_probe import sha, write


def wall_path():
    relative = "results/stage4/simvascular/mesh/20261002T141441884393Z/SV_MESH/mesh-surfaces/WALL.vtp"
    choices = [
        ROOT / relative,
        Path("/home/lzy/projects/computation_examples/mevo_acoustic_microbubbles") / relative,
        Path("/workspace/mevo_acoustic_microbubbles") / relative,
    ]
    return next(p for p in choices if p.exists())


def glyph(plotter, positions, radius, color):
    data = pv.PolyData(positions)
    source = vtk.vtkSphereSource()
    source.SetRadius(radius)
    source.SetThetaResolution(24)
    source.SetPhiResolution(20)
    mapper = vtk.vtkGlyph3DMapper()
    mapper.SetInputData(data)
    mapper.SetSourceConnection(source.GetOutputPort())
    mapper.ScalingOff()
    mapper.ScalarVisibilityOff()
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    prop = actor.GetProperty()
    prop.SetColor(*pv.Color(color).float_rgb)
    prop.SetAmbient(0.18)
    prop.SetDiffuse(0.8)
    prop.SetSpecular(0.32)
    prop.SetSpecularPower(24)
    plotter.add_actor(actor)
    return data


def decorate(plotter, wall, bounds, origin, dark, title):
    color = "#dadde4" if dark else "#282c33"
    plotter.set_background("#080b12" if dark else "white")
    w = wall.copy()
    w.points = (w.points - origin) * 1e6
    # Clip by six linear planes. A composite box scalar can miss a large WALL
    # triangle whose vertices all lie outside this microscopic field of view.
    patch = w
    for axis in range(3):
        normal = np.eye(3)[axis]
        lo = np.zeros(3)
        lo[axis] = bounds[2 * axis]
        hi = np.zeros(3)
        hi[axis] = bounds[2 * axis + 1]
        patch = patch.clip(normal=normal, origin=lo, invert=False)
        patch = patch.clip(normal=normal, origin=hi, invert=True)
    if patch.n_cells:
        plotter.add_mesh(
            patch,
            color="#a15050",
            opacity=0.24 if dark else 0.18,
            smooth_shading=True,
            ambient=0.4,
            specular=0.18,
        )
    plotter.add_text(title, position="upper_left", font_size=13, color=color, name="title")
    center = np.array([(bounds[i] + bounds[i + 1]) / 2 for i in [0, 2, 4]])
    span = max(bounds[i + 1] - bounds[i] for i in [0, 2, 4])
    axes_grid(plotter, bounds, color, dark)
    plotter.camera_position = [center + span * np.array([0.9, -1.55, 0.8]), center, [0, 0, 1]]
    plotter.enable_parallel_projection()
    plotter.camera.parallel_scale = span * 0.86
    # Standard alpha rendering is reliable for this single cropped wall sheet.
    # The local software EGL depth-peeling path suppresses translucent actors.


def axes_grid(plotter, bounds, color, dark):
    lo = np.array(bounds)[::2]
    hi = np.array(bounds)[1::2]
    span = max(hi - lo)
    segments = []
    for axis in range(3):
        for value in np.linspace(lo[axis], hi[axis], 3):
            # A floor grid plus a back grid; no surrounding box.
            if axis in [0, 1]:
                a = lo.copy()
                b = lo.copy()
                a[axis] = b[axis] = value
                other = 1 - axis
                b[other] = hi[other]
                segments.extend([a, b])
            if axis in [0, 2]:
                a = lo.copy()
                b = lo.copy()
                a[1] = b[1] = hi[1]
                a[axis] = b[axis] = value
                other = 2 - axis
                b[other] = hi[other]
                segments.extend([a, b])
    plotter.add_lines(np.array(segments), color="#465160" if dark else "#b8bec6", width=0.7, connected=False)
    points = []
    labels = []
    for axis in range(3):
        a = lo.copy()
        b = lo.copy()
        b[axis] = hi[axis]
        if axis == 1:
            a[0] = b[0] = hi[0]
        plotter.add_lines(np.array([a, b]), color=color, width=1.2)
        for fraction in [0, 0.5, 1]:
            if not dark and axis > 0 and fraction == 0:
                continue  # shared corners otherwise collide at print-size labels
            p = a + (b - a) * fraction
            if axis == 0:
                p[1] -= 0.035 * span
                p[2] -= 0.035 * span
            elif axis == 1:
                p[0] += 0.085 * span
                p[2] -= 0.025 * span
            else:
                p[0] -= 0.055 * span
            points.append(p)
            labels.append(f"{lo[axis] + fraction * (hi[axis] - lo[axis]):.1f}")
        p = (a + b) / 2
        if axis == 0:
            p[1] -= 0.16 * span
            p[2] -= 0.08 * span
        elif axis == 1:
            p[0] += 0.15 * span
            p[2] -= 0.08 * span
        else:
            p[0] -= 0.19 * span
        points.append(p)
        labels.append("XYZ"[axis] + " (um)")
    plotter.add_point_labels(
        np.array(points),
        labels,
        font_size=11 if dark else 18,
        text_color=color,
        font_family="arial",
        shape=None,
        show_points=False,
        always_visible=True,
        shadow=False,
        margin=0,
    )


def context(report):
    arrays = [
        dict(np.load(BASE / "comparisons" / f"{report['label']}_{kind}.npz"))
        for kind in ["reference", "candidate"]
    ]
    origin = arrays[0]["positions_m"][0].mean(0)
    positions = np.concatenate([d["positions_m"].reshape(-1, 3) for d in arrays])
    assert np.isfinite(positions).all()
    radius = float(arrays[0]["radius_m"])
    lower = ((positions.min(0) - origin) - 4 * radius) * 1e6
    upper = ((positions.max(0) - origin) + 4 * radius) * 1e6
    center = (lower + upper) / 2
    extent = max(upper - lower)
    bounds = np.column_stack([center - extent / 2, center + extent / 2]).ravel().tolist()
    return arrays, origin, bounds, radius


def main(video=False, requested_duration=None):
    candidates = [
        json.loads(p.read_text())
        for p in (BASE / "comparisons").glob("*_legacy.json")
        if "refined_" not in p.stem
    ]
    if requested_duration is not None:
        candidates = [r for r in candidates if abs(r["compared_duration_s"] - requested_duration) <= 1e-12]
    selected = {
        case: max(
            [r for r in candidates if r["label"].startswith(case)],
            key=lambda r: (r["compared_duration_s"], "_v4a_" in r["label"], "_v4_" in r["label"]),
        )
        for case in ["B2", "B3"]
    }
    out = BASE / "figures"
    wallfile = wall_path()
    wall = pv.read(wallfile)
    provenance = dict(
        wall=str(wallfile),
        wall_sha256=sha(wallfile),
        script_sha256=sha(Path(__file__)),
        orchestration_checkout_at_render_time=(ROOT / "SOURCE_COMMIT_R2.txt").read_text().strip(),
        cases={},
        geometry="Actual finite spheres at physical radius, uniform m-to-µm units, no group remapping or particle enlargement; real open WALL patch.",
        view="Fixed parallel camera and equal XYZ scales shared by old/new, fixed common origin per case.",
    )
    for case, report in selected.items():
        candidate_identity = json.loads((local(report["candidate"]).parent / "IDENTITY.json").read_text())
        arrays, origin, bounds, radius = context(report)
        indices = np.rint(np.linspace(0, len(arrays[0]["time_s"]) - 1, 4)).astype(int)
        plotter = pv.Plotter(shape=(2, 4), window_size=(2560, 1280), off_screen=True, border=False)
        for row, (d, name) in enumerate(zip(arrays, ["Frozen Stage4L", "R2 S1/S0 V6"])):
            for col, index in enumerate(indices):
                plotter.subplot(row, col)
                elapsed = (d["requested_time_s"][index] - d["requested_time_s"][0]) * 1e3
                title = f"{case} | {name}\nΔt = {elapsed:.3f} ms; n = {len(d['particle_ids'])}\nSample age = {d['sample_age_s'][index] * 1e6:.3g} µs"
                decorate(plotter, wall, bounds, origin, False, title)
                glyph(plotter, (d["positions_m"][index] - origin) * 1e6, radius * 1e6, "#2873ad")
        image = plotter.screenshot(return_img=True)
        plotter.close()
        fig, ax = plt.subplots(figsize=(16, 8))
        ax.imshow(image)
        ax.axis("off")
        fig.suptitle(
            f"R7 | {case}: actual motion of fixed original IDs, radius {radius * 1e6:g} µm", fontweight="bold"
        )
        fig.subplots_adjust(left=0, right=1, bottom=0.035, top=0.96)
        fig.text(
            0.01,
            0.008,
            f"Common time since t = 3 s. Comparison: {report['status']}. Physical radii; no coordinate interpolation. R2 source {candidate_identity['source_commit'][:12]}.\n"
            f"Old: {Path(report['reference']).parent.name}; R2: {Path(report['candidate']).parent.name}",
            fontsize=7,
        )
        for suffix in ["png", "pdf"]:
            fig.savefig(out / f"R7_ACTUAL_MOTION_{case}.{suffix}", dpi=200, bbox_inches="tight")
        plt.close(fig)
        provenance["cases"][case] = dict(
            report=report,
            origin_m=origin.tolist(),
            bounds_um=bounds,
            selected_indices=indices.tolist(),
            candidate_identity=candidate_identity,
        )
    if video:
        duration = min(r["compared_duration_s"] for r in selected.values())
        assert duration >= 0.001 - 1e-12, "No >=1 ms common comparison: MP4 gate not met"
        # Uniform clock at <=2 µs; show the latest actual accepted state. Sample
        # ages remain visible. Video time is a uniform slow-motion clock.
        count = int(np.ceil(duration / 2e-6)) + 1
        clock = 3.0 + np.linspace(0, duration, count)
        plotter = pv.Plotter(shape=(2, 2), window_size=(1920, 1080), off_screen=True, border=False)
        layers = []
        archive = {}
        maxage = {}
        for row, (case, report) in enumerate(selected.items()):
            arrays, origin, bounds, radius = context(report)
            members = arrays[0]["particle_ids"]
            for col, (kind, name) in enumerate(
                [("reference", "Frozen Stage4L"), ("candidate", "R2 S1/S0 V6")]
            ):
                path = local(report[kind])
                d = actual_sample(path, clock, members)
                plotter.subplot(row, col)
                decorate(plotter, wall, bounds, origin, True, f"{case} | {name}")
                data = glyph(plotter, (d["positions_m"][0] - origin) * 1e6, radius * 1e6, "#2999e2")
                layers.append((row, col, case, name, d, data, origin))
                key = f"{case}_{kind}"
                archive[key + "_positions_m"] = d["positions_m"]
                archive[key + "_actual_time_s"] = d["actual_time_s"]
                archive[key + "_particle_ids"] = members
                maxage[key] = float(d["sample_age_s"].max())
        target = out / "R7_OLD_NEW_ACTUAL_MOTION_60fps.mp4"
        with imageio.get_writer(
            target,
            fps=60,
            codec="libx264",
            quality=8,
            macro_block_size=2,
            ffmpeg_params=["-pix_fmt", "yuv420p", "-movflags", "+faststart"],
        ) as writer:
            for i, at in enumerate(clock):
                for row, col, case, name, d, data, origin in layers:
                    data.points = (d["positions_m"][i] - origin) * 1e6
                    plotter.subplot(row, col)
                    plotter.add_text(
                        f"{case} | {name} | n = {len(d['particle_ids'])}\nt = {at:.6f} s; Δt = {(at - 3) * 1e3:.3f} ms\nActual-state age = {d['sample_age_s'][i] * 1e6:.3f} µs",
                        position="upper_left",
                        font_size=13,
                        color="#dadde4",
                        name="title",
                    )
                writer.append_data(plotter.screenshot(return_img=True))
        plotter.close()
        np.savez_compressed(out / "R7_VIDEO_ACTUAL_FRAMES.npz", requested_time_s=clock, **archive)
        snapshot = out / "R7_VIDEO_RENDER_SOURCE.py"
        snapshot.write_bytes(Path(__file__).read_bytes())
        provenance["video"] = dict(
            path=str(target),
            sha256=sha(target),
            fps=60,
            frames=count,
            physical_duration_s=duration,
            max_sample_age_s=maxage,
            coordinate_interpolation=False,
            render_script_sha256=sha(snapshot),
            render_source_snapshot=str(snapshot),
            frame_clock="uniform; last actual accepted state held when no newer accepted state exists",
            actual_motion=all(
                np.max(abs(v[1:] - v[:1])) > 0 for k, v in archive.items() if k.endswith("positions_m")
            ),
        )
    if not video and (out / "R7_DATA_SOURCES.json").exists():
        previous = json.loads((out / "R7_DATA_SOURCES.json").read_text())
        if "video" in previous:
            provenance["video"] = dict(
                previous["video"],
                render_script_sha256=previous.get("video", {}).get(
                    "render_script_sha256", previous["script_sha256"]
                ),
                render_source_snapshot=str(out / "R7_VIDEO_RENDER_SOURCE.py"),
                note="Existing actual-motion video retained. This pass only increases static figure tick-label sizes; video rendering code and data are unchanged.",
            )
    write(out / "R7_DATA_SOURCES.json", provenance)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--video", action="store_true")
    p.add_argument(
        "--duration", type=float, help="Select an already completed comparison duration in seconds."
    )
    a = p.parse_args()
    main(a.video, a.duration)
