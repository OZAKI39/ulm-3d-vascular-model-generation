"""Primary transport with fixed physical centers and explicit display-only sampling."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import shutil
import time

import numpy as np
import pyvista as pv
from matplotlib import font_manager

from mevo_mb.coupling.bifurcation import BifurcationGeometry
from scripts.stage4l import render_motion as previous

ROOT, BASE, style = previous.ROOT, previous.BASE, previous.style
sha, write = previous.sha, previous.write


def display_samples(d, run, protected):
    """Anchor saved positions exactly; interpolate IDs, never connected components.

    Intermediate positions are drawing aids, not additional solver output. Small
    chord-wall intersections are corrected only for these unsaved display points.
    """
    branch_path = run / "BRANCH_METRICS.json"
    protected[str(branch_path)] = sha(branch_path)
    branches = json.loads(branch_path.read_text())
    times = np.arange(1201) * 0.005
    times[::4] = d["t"]
    exit_times = np.full(len(d["ids"]), np.inf)
    for event in branches["crossings"]:
        exit_times[d["id_index"][event["particle_id"]]] = event["crossing_time_s"]
    sampled, left, right, alpha = style.base.sample_case(d | dict(events=branches["crossings"]), times)
    sampled[::4] = d["x"]
    active = times[:, None] < exit_times[None]
    native_active = np.stack([s["active_mask"] for s in d["states"]])
    assert np.array_equal(active[::4], native_active)
    radii = np.broadcast_to(d["radii"], active.shape)
    geometry = BifurcationGeometry(d["config"]["geometry"])
    geometry.bind_mesh_surface(pv.read(d["surface"]))
    intermediate = active.copy()
    intermediate[::4] = False
    before = sampled[intermediate].copy()
    adjusted, _ = geometry.project(before, radii[intermediate])
    sampled[intermediate] = adjusted
    correction = np.zeros(active.shape)
    correction[intermediate] = np.linalg.norm(adjusted - before, axis=1)
    # Reject large invented excursions instead of concealing a bad interpolation.
    assert correction.max() < 10e-6
    assert np.array_equal(sampled[::4], d["x"])
    gap, _, valid = geometry.wall_gap_and_normal(sampled[active], radii[active])
    assert valid.all() and gap.min() >= -1.001e-12
    assert np.all(gap + radii[active] > 0)
    common = active[:-1] & active[1:]
    steps = np.linalg.norm(np.diff(sampled, axis=0), axis=2)
    checks = dict(
        saved_position_max_error_m=0.0,
        saved_states_exact=True,
        all_display_active_centers_inside=True,
        all_display_active_samples=int(active.sum()),
        minimum_surface_gap_m=float(gap.min()),
        center_coordinate_magnification=1.0,
        group_membership_affects_particle_position=False,
        interpolated_wall_adjustment_samples=int(np.count_nonzero(correction)),
        max_interpolated_wall_adjustment_m=float(correction.max()),
        max_display_frame_displacement_m=float(steps[common].max()),
        inactive_to_active_transitions=int(((~active[:-1]) & active[1:]).sum()),
        raster_cross_dissolve=False,
        scope="Saved scientific positions are unchanged. Wall-corrected linear ID interpolation is display-only, not contact dynamics or newly computed trajectories.",
    )
    assert checks["inactive_to_active_transitions"] == 0
    return dict(
        time_s=times,
        positions_m=sampled,
        active=active,
        particle_id=d["ids"],
        source_left=left,
        source_right=right,
        coordinate_fraction=alpha,
        wall_adjustment_m=correction,
    ), checks


def label_for_frame(d, source, positions, active, project, pixels_per_mm):
    # The count is from an actual saved state. Only its leader follows display motion.
    group = d["clusters"][source][0].copy()
    ix = np.array([d["id_index"][pid] for pid in group["particle_ids"]], dtype=int)
    ix = ix[active[ix]]
    if not len(ix):
        group["count"] = 0
        return group
    center = positions[ix].mean(axis=0)
    uv = project(positions[ix])
    cuv = project(center[None])[0]
    group["screen_center_top_origin"] = [float(cuv[0]), float(1080 - cuv[1])]
    group["label_radius_px"] = float(
        np.max(
            np.linalg.norm(uv - cuv, axis=1) + previous.MAGNIFICATION * d["radii"][ix] * 1000 * pixels_per_mm
        )
    )
    return group


def viewer(out, videos):
    content = """<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Stage 4L · Corrected transport</title><style>body{background:#070a0f;color:#eef4ff;font:17px/1.7 Arial,sans-serif;max-width:1500px;margin:30px auto;padding:0 24px}h1,h2{font-weight:400}p{color:#b7cee2}video{width:100%;background:black}a{color:#66e0c2}</style>
<h1>完整 6 s 输运 · 真实位置显示</h1><p>保留黑底、单一三维视角、PZT、紫色声场、蓝色球体及泡数标注。画面内不显示 B2/B3。1920 × 1080，60 fps，0.30× 匀速，每段约 20 秒。</p>
<p>取消群内坐标放大：所有原始保存时刻的泡心与原始数据完全相同，分组变化只影响数量标注，不移动微泡。球体半径仍作为可见性符号放大，超出真实血管投影轮廓的符号部分被裁剪；放大符号的重叠不表示微泡融合或实际重叠。</p>
<p>原数据每 20 ms 保存一次。帧间按同一微泡 ID 做显示插值，不再对两幅画面做渐变。少量穿过弯曲网格边界的插值点仅作显示边界修正，原始保存点完全不动。这些中间画面不构成新的数值轨迹或接触动力学证据。</p>
<p>两段均为新润滑模型的主力场结果，不含次级力。没有重新运行科学计算，也没有修改原始 HDF5。真实空间分布比旧版 60 倍群内放大更紧凑。</p>"""
    for key, v in videos.items():
        c = v["source"]["final_branch_counts"]
        content += f'<h2>{v["case"]}</h2><video controls playsinline preload="metadata" poster="keyframes/{key}_000000.png" src="{v["path"]}"></video>'
        content += f"<p>出口 A：{c['N_A']}；出口 B：{c['N_B']}；观察截止未分配：{c['N_unassigned']}。</p>"
    content += '<p><a href="DISPLAY_GEOMETRY_CHECKS.json">位置与血管边界检查</a> · <a href="DELIVERY_VERIFICATION.json">视频检查</a> · <a href="RENDER_MANIFEST.json">来源</a></p></html>'
    (out / "OPEN_RESULTS.html").write_text(content)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--case", choices=["B2", "B3", "all"], default="all")
    parser.add_argument("--preview-only", action="store_true")
    args = parser.parse_args()
    out = args.output.resolve()
    assert not out.exists(), "Preserve previous deliveries."
    (out / "keyframes").mkdir(parents=True)
    (out / "provenance").mkdir()
    runs, protected = previous.selected_runs()
    protected[str(previous.REFERENCE)] = sha(previous.REFERENCE)
    protected[str(style.FONT)] = sha(style.FONT)
    source_paths = [
        Path(__file__),
        Path(previous.__file__),
        Path(style.__file__),
        Path(style.original.__file__),
        Path(style.playback.__file__),
        Path(style.base.__file__),
        ROOT / "scripts/stage3/render_motion.py",
        ROOT / "src/mevo_mb/coupling/bifurcation.py",
        ROOT / "src/mevo_mb/bubbles/metrics.py",
        ROOT / "scripts/stage4l/check_motion.py",
    ]
    sources = {str(p): sha(p) for p in source_paths}
    snapshots = {}
    for i, p in enumerate(source_paths):
        target = out / "provenance" / f"{i:02d}_{p.name}"
        shutil.copy2(p, target)
        snapshots[str(target)] = sha(target)
    font_manager.fontManager.addfont(str(style.FONT))
    style.base.plt.rcParams.update({"font.family": "Arial", "font.weight": "normal", "font.size": 14})
    started = time.perf_counter()
    videos, screens, geometry_checks, groups, rows = {}, {}, {}, {}, []
    cases = ["B2", "B3"] if args.case == "all" else [args.case]
    for case in cases:
        key = case + "_primary"
        d, _, info = previous.load(case, "primary", runs["primary"][case], protected)
        samples, checks = display_samples(d, runs["primary"][case], protected)
        d["display_geometry"] = "physical_centers"
        d["scope_label"] = "Normal lubrication  |  Primary only  |  0.30\u00d7 playback"
        fig, render, screens[key] = style.scene(
            case, d, sphere_mode=True, magnification=previous.MAGNIFICATION
        )
        mask = np.asarray(render.lumen_mask)
        outside = mask == 0
        background = np.asarray(render.static_image)[outside]
        render.lumen_mask.save(out / f"{key}_LUMEN_MASK.png")
        render.static_image.save(out / f"{key}_BACKGROUND.png")
        np.savez_compressed(out / f"{key}_DISPLAY_SAMPLES.npz", **samples)
        checks.update(
            frames_checked_for_glyph_clipping=0,
            outside_lumen_glyph_pixels=0,
            maximum_display_wall_adjustment_pixels=checks["max_interpolated_wall_adjustment_m"]
            * 1000
            * screens[key]["pixels_per_mm"],
        )

        def frame_image(frame):
            k = frame // 4
            x, active = samples["positions_m"][frame], samples["active"][frame]
            im = render(k, positions_m=x, active_mask=active)
            assert np.array_equal(np.asarray(im)[outside], background)
            checks["frames_checked_for_glyph_clipping"] += 1
            label = label_for_frame(d, k, x, active, render.screen_points, screens[key]["pixels_per_mm"])
            return style.overlay(im, d, k, samples["time_s"][frame], label_group=label)

        selected = [0, 1, 2, 40, 384, 385, 388, 600, 608, 609, 612, 1200]
        if args.preview_only:
            for frame in selected:
                frame_image(frame).save(out / "keyframes" / f"{key}_{frame:06d}.png")
            geometry_checks[key] = checks
            style.base.plt.close(fig)
            print(
                json.dumps(dict(preview=case, geometry=checks, elapsed_s=time.perf_counter() - started)),
                flush=True,
            )
            continue
        path = out / f"stage4L_{case}_primary_physical_centers_60fps.mp4"
        with style.playback.writer(path) as encoder:
            for frame, t in enumerate(samples["time_s"]):
                im = frame_image(frame)
                encoder.append_data(np.asarray(im))
                if frame in selected:
                    im.save(out / "keyframes" / f"{key}_{frame:06d}.png")
                left, right = frame // 4, min(frame // 4 + 1, 300)
                rows.append(
                    [
                        key,
                        frame,
                        frame / 60,
                        float(t),
                        left,
                        right,
                        float(d["t"][left]),
                        float(d["t"][right]),
                        0.0,
                        (frame % 4) / 4,
                    ]
                )
                if frame % 240 == 0:
                    print(
                        json.dumps(
                            dict(video=key, frame=frame, total=1201, elapsed_s=time.perf_counter() - started)
                        ),
                        flush=True,
                    )
        style.base.plt.close(fig)
        geometry_checks[key] = checks
        info["sampling"] = (
            "Saved positions exact; linear same-ID display sampling, with wall correction only at unsaved display times. No raster blending or group-coordinate scaling."
        )
        videos[key] = dict(
            path=path.name,
            fps=60,
            frames=1201,
            duration_s=1201 / 60,
            resolution=[1920, 1080],
            sha256=sha(path),
            kind="primary",
            case=case,
            source=info,
            playback_speed=0.3,
            display_samples=f"{key}_DISPLAY_SAMPLES.npz",
            lumen_mask=f"{key}_LUMEN_MASK.png",
        )
        groups[key] = [
            dict(
                time_s=s["t"],
                source_state=s["source_state"],
                stats=s["stats"],
                visible_particle_ids=d["ids"][s["active_mask"]].tolist(),
                groups=[
                    {k: v for k, v in g.items() if k not in ["display_positions_m", "display_radii_m"]}
                    for g in d["clusters"][i]
                ],
            )
            for i, s in enumerate(d["states"])
        ]
    assert all(sha(p) == digest for p, digest in (protected | sources).items())
    write(out / "DISPLAY_GEOMETRY_CHECKS.json", geometry_checks)
    write(out / "SCREEN_CHECKS.json", screens)
    if args.preview_only:
        return
    write(out / "STATE_GROUPS.json", groups)
    with (out / "frame_map.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "video",
                "frame",
                "video_time_s",
                "display_time_s",
                "source_left",
                "source_right",
                "source_left_time_s",
                "source_right_time_s",
                "image_blend_fraction",
                "coordinate_fraction",
            ]
        )
        writer.writerows(rows)
    write(
        out / "RENDER_MANIFEST.json",
        dict(
            status="RENDERED_PENDING_DELIVERY_CHECK",
            videos=videos,
            protected_inputs=protected,
            renderer_sources=snapshots,
            original_source_sha256=sources,
            display_geometry="physical_centers",
            screen_checks=screens,
            case_identifiers_in_video=False,
            coordinate_magnification=1.0,
            sphere_symbol_magnification=previous.MAGNIFICATION,
            coordinate_interpolation=True,
            raster_cross_dissolve=False,
            scientific_data_modified=False,
            simulations_launched=False,
            geometry_checks=geometry_checks,
            temporal_scope="Interpolation and small mesh-boundary corrections are drawing aids only. No interpolated point is represented as a solver state. All 301 stored positions are preserved bitwise.",
            spatial_scope="One shaded visibility sphere per active ID at its physical/display-interpolated center. The inflated sphere layer is clipped to the projection of the exact vessel mesh. Clusters affect labels only.",
            elapsed_s=time.perf_counter() - started,
        ),
    )
    viewer(out, videos)
    (out / "REPRODUCE.txt").write_text(
        f"cd {ROOT}\nPYTHONPATH=src:. /home/lzy/projects/computation_examples/_shared/python_environment/bin/python -m scripts.stage4l.render_transport --output <new-directory>\nPYTHONPATH=src:. /home/lzy/projects/computation_examples/_shared/python_environment/bin/python -m scripts.stage4l.check_motion <new-directory>\n"
    )
    print(json.dumps(dict(status="RENDERED", output=str(out))), flush=True)


if __name__ == "__main__":
    main()
