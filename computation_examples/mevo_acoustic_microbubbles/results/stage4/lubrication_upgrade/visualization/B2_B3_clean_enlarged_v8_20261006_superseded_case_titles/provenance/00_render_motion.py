"""Render completed Stage 4L data with the established v7 single-scene style.

This module reads recorded trajectories only. Primary and dense-checkpoint
experiments remain separate videos with different, constant playback speeds.
"""

from __future__ import annotations

import argparse
from collections import OrderedDict
import csv
import html
import json
from pathlib import Path
import shutil
import time

import h5py
from matplotlib import font_manager
import numpy as np
from PIL import Image

from scripts.stage4 import render_single_overview as style
from mevo_mb.bubbles.metrics import aggregation_metrics

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "results/stage4/lubrication_upgrade"
REFERENCE = Path(
    "/home/lzy/projects/computation_examples/mevo_acoustic_microbubbles/results/stage4/"
    "visualization/B2_B3_clean_enlarged_v7/stage4_B2_spheres_uniform_60fps.mp4"
)
FPS = 60
MAGNIFICATION = 60.0
sha = style.base.sha


def write(path, value):
    def convert(x):
        return x.tolist() if isinstance(x, np.ndarray) else x.item()

    path.write_text(json.dumps(value, indent=2, allow_nan=False, default=convert) + "\n")


def selected_runs():
    """Use the reviewed final evidence, never the latest arbitrary run directory."""
    branch = BASE / "comparison_legacy/20261005T232125840331Z_branches/SUMMARY.json"
    reports = {str(branch): sha(branch)}
    primary = json.loads(branch.read_text())["primary_comparison"]
    result = {"primary": {}, "dense_probe": {}}

    def project_path(value):
        prefix = "/workspace/mevo_acoustic_microbubbles_stage4l/"
        relative = Path(value.removeprefix(prefix))
        assert not relative.is_absolute() and ".." not in relative.parts, value
        return ROOT / relative

    for case in ["B2", "B3"]:
        result["primary"][case] = project_path(primary[case]["candidate"])
        p = BASE / "comparison_legacy/paired_5ms_actual_states_v5" / f"{case}.json"
        reports[str(p)] = sha(p)
        result["dense_probe"][case] = project_path(json.loads(p.read_text())["lubrication"])
    return result, reports


def resolve_input(path, expected_sha):
    relative = str(path).removeprefix("/workspace/mevo_acoustic_microbubbles/")
    alternatives = [
        relative,
        relative.removeprefix("results/stage4/performance_refactor/gold_baseline/inputs/"),
    ]
    for root in [
        ROOT,
        Path("/home/lzy/projects/computation_examples/mevo_acoustic_microbubbles"),
    ]:
        for relative in alternatives:
            candidate = root / relative
            if candidate.is_file() and sha(candidate) == expected_sha:
                return candidate
    raise FileNotFoundError(f"No checksum-matching local input: {path}")


def load(case, kind, run, protected):
    identity_path = run / ("IDENTITY.json" if kind == "dense_probe" else "RUN_IDENTITY.json")
    identity = json.loads(identity_path.read_text())
    cfg = identity["config"]
    assert cfg["particles"]["translation_backend"] == "lubrication_coupled"
    assert cfg["particles"]["secondary_enabled"] == (kind == "dense_probe")
    summary_path = run / "SUMMARY.json"
    summary = json.loads(summary_path.read_text())
    assert summary["status"] == "COMPLETE", run
    trajectory = next(run.rglob("trajectories.h5"))
    archive = BASE / "comparison_legacy/final_stage4l_20261005/ARTIFACT_INDEX.json"
    records = json.loads(archive.read_text())["raw_trajectories"]
    expected = next(r["sha256"] for r in records if r["path"] == str(trajectory.relative_to(ROOT)))
    assert sha(trajectory) == expected
    for p in [identity_path, summary_path, trajectory, archive]:
        protected[str(p)] = sha(p)
    with h5py.File(trajectory, "r") as f:
        assert f.attrs["complete"] and f.attrs["translation_backend"] == "lubrication_coupled"
        times = f["time_s"][:]
        ids, radii = f["particle_id"][:], f["radii_m"][:]
        assert np.all(np.diff(times) > 0)
        if kind == "dense_probe":
            assert f.attrs["all_accepted_substeps_recorded"]
            assert np.isclose(times[0], 3.0) and np.isclose(times[-1], 3.005)
            # Every video frame uses an actual accepted state at/before its uniform clock.
            display_times = np.linspace(times[0], times[-1], 1001)
            selected = np.searchsorted(times, display_times, side="right") - 1
            source_indices, inverse = np.unique(selected, return_inverse=True)
            timeline = [(int(i), int(i), 0.0, float(t)) for i, t in zip(inverse, display_times, strict=True)]
        else:
            np.testing.assert_allclose(times, np.arange(301) * 0.02, atol=1e-12, rtol=0)
            source_indices = np.arange(len(times))
            timeline = list(style.uniform_timeline())
        x = f["positions_m"][source_indices]
        active = f["active"][source_indices]
    assert len(ids) == len(set(ids.tolist())) == 256
    assert x.shape == (len(source_indices), len(ids), 3) and np.isfinite(x).all()
    assert np.allclose(radii, 1.78e-6, rtol=0, atol=1e-20)
    inputs = identity.get("input_sha256", identity.get("inputs_sha256"))
    field_source = cfg["effective_acoustic_field"]["path"]
    mesh_source = cfg["geometry"]["mesh_surface"]
    field_path = resolve_input(field_source, inputs[field_source])
    surface = resolve_input(mesh_source, inputs[mesh_source])
    protected[str(field_path)], protected[str(surface)] = sha(field_path), sha(surface)
    assert protected[str(field_path)] == cfg["pzt_control"]["field_sha256"]
    with h5py.File(field_path, "r") as f:
        axes = [f[f"axes/{axis}_m"][:] * 1000 for axis in "xyz"]
        pressure = f["pressure_Pa"][:]
        assert f.attrs["frequency_Hz"] == 2.5e6
    assert pressure.shape == tuple(len(a) for a in axes)
    assert abs(pressure).max() <= 39000.0001
    states = []
    recorded = None
    if kind == "primary":
        metric_path = run / "metrics.json"
        recorded = json.loads(metric_path.read_text())["rows"]
        protected[str(metric_path)] = sha(metric_path)
    for k, index in enumerate(source_indices):
        stats, _, _, _ = aggregation_metrics(x[k, active[k]], radii[active[k]], cfg["metrics"])
        if recorded is not None:
            for key in [
                "cluster_count",
                "largest_cluster_size",
                "median_nearest_neighbor_m",
                "near_contact_fraction",
            ]:
                assert np.isclose(stats[key], recorded[index][key], rtol=1e-12, atol=1e-15), (
                    case,
                    index,
                    key,
                )
        states.append(
            dict(
                t=float(times[index]),
                source_state=int(index),
                active_mask=active[k],
                active_particles=int(active[k].sum()),
                stats=stats,
            )
        )
    d = dict(
        t=times[source_indices],
        x=x,
        ids=ids,
        radii=radii,
        states=states,
        config=cfg,
        id_index={int(pid): i for i, pid in enumerate(ids)},
        axes=axes,
        P=pressure,
        surface=surface,
        field_path=field_path,
        result_id=str(run.relative_to(ROOT)),
    )
    if kind == "dense_probe":
        d.update(
            scene_title=f"{case}  |  Acoustic field & bubble aggregation",
            scene_subtitle=f"2.5 MHz  |  Flow + primary + secondary force  |  {int(active[0].sum())} active bubbles",
            scope_label="Normal lubrication  |  5 ms probe from legacy t = 3 s",
            clock_precision=3,
        )
    else:
        d.update(
            scene_title=f"{case}  |  Acoustic field & bubble transport",
            scene_subtitle="2.5 MHz  |  Flow + primary force  |  256 bubbles",
            scope_label="Normal lubrication  |  Primary only  |  0.30\u00d7 playback",
        )
    style.prepare_clusters(d)
    ages = np.array([t - d["t"][left] for left, _, _, t in timeline])
    assert ages.min() >= -1e-12
    info = dict(
        run=str(run),
        trajectory=str(trajectory),
        trajectory_sha256=expected,
        kind=kind,
        case=case,
        secondary_enabled=cfg["particles"]["secondary_enabled"],
        time_start_s=float(times[0]),
        time_end_s=float(times[-1]),
        source_saved_states=len(times),
        rendered_source_states=len(source_indices),
        initial_active=int(active[0].sum()),
        final_active=int(active[-1].sum()),
        max_last_accepted_age_s=float(ages.max()) if kind == "dense_probe" else None,
        final_branch_counts={
            k: summary.get("branch_metrics", {}).get(k) for k in ["N_A", "N_B", "N_unassigned"]
        },
        sampling="actual last-accepted states; no interpolation"
        if kind == "dense_probe"
        else "all 301 saved states; raster cross-dissolves only; no coordinate interpolation",
    )
    return d, timeline, info


def viewer(out, videos):
    page = """<!doctype html><html lang="zh-CN"><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>Stage 4L · B2/B3</title>
<style>body{background:#070a0f;color:#eef4ff;font:17px/1.7 Arial,sans-serif;max-width:1500px;margin:30px auto;padding:0 24px}h1,h2{font-weight:400}p{color:#b7cee2}video{width:100%;background:black}a{color:#66e0c2}</style>
<h1>Stage 4L · 新润滑模型 B2/B3 动画</h1>
<p>沿用 v7 单一三维场景、黑底、蓝色球体、PZT 声源及紫色稳态声场。1920 × 1080，60 fps；每段全程匀速，不改变开头速度，不循环填充片段。</p>
<p>每个活动微泡对应一个圆球；半径和群内相对坐标均放大 60 倍，群体质心保持实际位置，保留真实逐泡排布和距离/半径比。此为可见性放大，不代表毫米级团簇尺寸或微泡融合。聚集阈值仍为中心距离 10 µm。</p>
<p>两组实验独立展示。聚集探针从旧模型 t=3 s 的密集团簇快照出发，仅运行新模型 5 ms；它不表示新模型从初始条件形成了同样的团簇。完整 6 s 视频仅含主力场；本轮尚无完整 6 s 新次级力结果。近接触物理局限和硬接触抖动仍存在，画面平稳不代表物理精度已验证。</p>"""
    for key, v in videos.items():
        title = "5 ms 聚集探针（含次级力）" if v["kind"] == "dense_probe" else "完整 6 s 输运（仅主力场）"
        explanation = (
            "t=3000–3005 ms，约 16.7 秒视频；每帧取该时刻之前最近的实际接受状态，未插值粒子位置。"
            if v["kind"] == "dense_probe"
            else "t=0–6000 ms，约 20 秒视频，0.30× 匀速；301 个保存状态之间使用原版画面渐变，未重建未保存的接触轨迹。"
        )
        page += f"<h2>{v['case']} · {title}</h2><p>{explanation}</p>"
        page += f'<video controls playsinline preload="metadata" poster="keyframes/{key}_000000.png" src="{html.escape(v["path"])}"></video>'
        if v["kind"] == "primary":
            counts = v["source"]["final_branch_counts"]
            page += f"<p>出口 A：{counts['N_A']}；出口 B：{counts['N_B']}；观察截止未分配：{counts['N_unassigned']}。</p>"
    page += '<p><a href="RENDER_MANIFEST.json">来源与显示方式</a> · <a href="frame_map.csv">逐帧时间映射</a> · <a href="DELIVERY_VERIFICATION.json">视频核验</a></p>'
    page += "<script>const vs=[...document.querySelectorAll('video')];vs.forEach(v=>v.addEventListener('play',()=>vs.forEach(o=>{if(o!==v)o.pause()})));</script></html>"
    (out / "OPEN_RESULTS.html").write_text(page)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--kind", choices=["all", "dense_probe", "primary"], default="all")
    parser.add_argument("--case", choices=["all", "B2", "B3"], default="all")
    parser.add_argument("--preview-only", action="store_true")
    args = parser.parse_args()
    out = args.output.resolve()
    assert not out.exists(), "Use a new directory; existing deliveries are preserved."
    (out / "keyframes").mkdir(parents=True)
    (out / "provenance").mkdir()
    runs, protected = selected_runs()
    protected[str(REFERENCE)] = sha(REFERENCE)
    protected[str(style.FONT)] = sha(style.FONT)
    source_files = [
        Path(__file__),
        Path(style.__file__),
        Path(style.original.__file__),
        Path(style.playback.__file__),
        Path(style.base.__file__),
        ROOT / "scripts/stage3/render_motion.py",
        ROOT / "src/mevo_mb/bubbles/metrics.py",
        ROOT / "scripts/stage4l/check_motion.py",
    ]
    sources = {str(p): sha(p) for p in source_files}
    for i, p in enumerate(source_files):
        shutil.copy2(p, out / "provenance" / f"{i:02d}_{p.name}")
    font_manager.fontManager.addfont(str(style.FONT))
    style.base.plt.rcParams.update({"font.family": "Arial", "font.weight": "normal", "font.size": 14})
    videos, checks, rows, states = {}, {}, [], {}
    started = time.perf_counter()
    kinds = ["dense_probe", "primary"] if args.kind == "all" else [args.kind]
    cases = ["B2", "B3"] if args.case == "all" else [args.case]
    for kind in kinds:
        for case in cases:
            key = f"{case}_{kind}"
            d, timeline, info = load(case, kind, runs[kind][case], protected)
            fig, render, checks[key] = style.scene(case, d, sphere_mode=True, magnification=MAGNIFICATION)
            checkpoints = sorted(set([0, 1, 2, 10, 40, len(timeline) // 2, len(timeline) - 1]))
            cache = OrderedDict()

            def source_image(k):
                if k not in cache:
                    cache[k] = render(k)
                image = cache[k]
                cache.move_to_end(k)
                while len(cache) > 2:
                    cache.popitem(last=False)
                return image

            def frame_image(frame):
                left, right, alpha, t = timeline[frame]
                a = source_image(left)
                im = Image.blend(a, source_image(right), alpha) if alpha else a.copy()
                return style.overlay(im, d, left, t)

            if args.preview_only:
                for frame in checkpoints:
                    frame_image(frame).save(out / "keyframes" / f"{key}_{frame:06d}.png")
                style.base.plt.close(fig)
                print(
                    json.dumps(
                        dict(
                            preview=key,
                            source=info,
                            elapsed_s=time.perf_counter() - started,
                        )
                    ),
                    flush=True,
                )
                continue
            path = out / f"stage4L_{case}_{kind}_spheres_uniform_60fps.mp4"
            with style.playback.writer(path) as encoder:
                for frame, (left, right, alpha, t) in enumerate(timeline):
                    im = frame_image(frame)
                    encoder.append_data(np.asarray(im))
                    if frame in checkpoints:
                        im.save(out / "keyframes" / f"{key}_{frame:06d}.png")
                    rows.append(
                        [
                            key,
                            frame,
                            frame / FPS,
                            t,
                            d["states"][left]["source_state"],
                            d["states"][right]["source_state"],
                            d["t"][left],
                            d["t"][right],
                            alpha,
                        ]
                    )
                    if frame % 240 == 0:
                        print(
                            json.dumps(
                                dict(
                                    video=key,
                                    frame=frame,
                                    total=len(timeline),
                                    elapsed_s=time.perf_counter() - started,
                                )
                            ),
                            flush=True,
                        )
            style.base.plt.close(fig)
            videos[key] = dict(
                path=path.name,
                fps=FPS,
                frames=len(timeline),
                duration_s=len(timeline) / FPS,
                resolution=[1920, 1080],
                sha256=sha(path),
                kind=kind,
                case=case,
                source=info,
                playback_speed=(info["time_end_s"] - info["time_start_s"]) / ((len(timeline) - 1) / FPS),
            )
            states[key] = [
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
    write(out / "SCREEN_CHECKS.json", checks)
    if args.preview_only:
        return
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
            ]
        )
        writer.writerows(rows)
    write(out / "STATE_GROUPS.json", states)
    manifest = dict(
        status="RENDERED_PENDING_DELIVERY_CHECK",
        videos=videos,
        protected_inputs=protected,
        renderer_sources=sources,
        screen_checks=checks,
        style_reference=str(REFERENCE),
        structure_magnification=MAGNIFICATION,
        coordinate_interpolation=False,
        scientific_data_modified=False,
        simulations_launched=False,
        scene="Same v7 camera, dimensions, grid, PZT, purple field, sphere lighting and cluster counts; no zoom panel or footer.",
        geometry="X_display = centroid + 60*(X_real-centroid); R_display = 60*R_real. Exact current connected components and IDs; no synthetic packing.",
        scope="Separate complete primary-only 6 s runs and 5 ms secondary-on probes from legacy t=3 s snapshots. No new full-secondary continuation is implied.",
        temporal_scope="Dense: actual last accepted states at a uniform 5 us display clock. Primary: original v7 4-frame raster cross-dissolves between saved 20 ms states. No particle-coordinate filtering or interpolation.",
        elapsed_s=time.perf_counter() - started,
    )
    write(out / "RENDER_MANIFEST.json", manifest)
    viewer(out, videos)
    (out / "REPRODUCE.txt").write_text(
        f"cd {ROOT}\nPYTHONPATH=src:. /home/lzy/projects/computation_examples/_shared/python_environment/bin/python "
        f"-m scripts.stage4l.render_motion --output <new-directory> --kind {args.kind} --case {args.case}\n"
        "PYTHONPATH=src:. /home/lzy/projects/computation_examples/_shared/python_environment/bin/python "
        "-m scripts.stage4l.check_motion <new-directory>\n"
    )
    print(json.dumps(dict(status="RENDERED", output=str(out))), flush=True)


if __name__ == "__main__":
    main()
