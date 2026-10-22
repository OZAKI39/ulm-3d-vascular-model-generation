"""Full-vessel movies with every recorded bubble in similarity-scaled group glyphs."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from pathlib import Path
import shutil
import time

import imageio.v2 as imageio
import numpy as np
from PIL import ImageDraw
import pyvista as pv

from scripts.stage4l import render_motion as previous
from scripts.stage4l import render_morphology as spheres
from scripts.stage4l import magnified_scene as context
from scripts.stage4l import magnified_geometry as layout

ROOT, write, sha = previous.ROOT, previous.write, previous.sha
SCALE, FPS, FRAMES = 24.0, 60, 1201


def write_gzip(path, value):
    with gzip.open(path, "wt") as stream:
        json.dump(value, stream, separators=(",", ":"), allow_nan=False)


def annotate(image, d, k, samples, groups, project):
    draw = ImageDraw.Draw(image)
    draw.text(
        (134, 390),
        f"t = {d['t'][k] * 1000:.1f} ms",
        font=context.font(30),
        fill=previous.style.MINT,
        anchor="la",
    )
    group = groups[k][0]
    if group["count"] > 1:
        ix = np.array([d["id_index"][pid] for pid in group["particle_ids"]])
        uv, scale, *_ = project(samples["positions_m"][k, ix])
        center = uv.mean(axis=0)
        radius = np.max(np.linalg.norm(uv - center, axis=1) + samples["radii_m"][ix] * scale)
        target = center + np.array([radius + 65, -radius - 65])
        target = np.clip(target, [200, 650], [context.WIDTH - 750, context.HEIGHT - 290])
        direction = target + [-16, 20] - center
        direction /= np.linalg.norm(direction)
        start = center + direction * (radius + 10)
        label = f"Largest: {group['count']} bubbles"
        draw.line([tuple(start), tuple(target + [-16, 20])], fill=previous.style.BLUE, width=2)
        box = draw.textbbox(tuple(target), label, font=context.font(23), anchor="lm")
        draw.rounded_rectangle((box[0] - 10, box[1] - 8, box[2] + 10, box[3] + 8), radius=6, fill="#080e15")
        draw.text(tuple(target), label, font=context.font(23), fill="#c1e3ff", anchor="lm")
    return image


def viewer(out, videos):
    body = """<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Microbubble morphology ×24</title><style>body{background:#070a0f;color:#eef4ff;font:17px/1.7 Arial,sans-serif;max-width:1600px;margin:30px auto;padding:0 24px}h1,h2{font-weight:400}p{color:#bbcede}video{width:100%;background:black}a{color:#66e0c2}</style><h1>完整血管全景 · 逐泡聚集形态 ×24</h1><p>两组均采用已完成 Stage 4L 的 0–6 s 主力场数据。黑底、紫色声场、蓝色球体、英文坐标和泡数标注；画面无 B2/B3 字样。四个 PZT 沿用论文图 v6 的灰蓝薄圆盘、透明边缘和朝外反光，声源前表面与实际配置一致，0.35 mm 背板为显示示意。</p><p>按用户选定的显示方案，每个真实泡群的球径和全部群内位移统一放大 24 倍。保留所有活动 ID、每个球体、相对排列、方向及半径比例；不合并成大球，不重排、不随机补点、不画聚集环。真实半径为 1.78 µm，显示球径为 85.44 µm；画面已标注 Morphology ×24 (display)。</p><p>完整血管使用原始开放端口壁面，无减面或裁短。为使放大泡群仍完整位于血管内部，并避免不同泡群的显示球体互相穿插，必要时将整群一致平移；因此显示球心位置、群间距离和距壁距离不代表物理坐标。真实泡数和聚集判定仍直接来自原数据。蓝色球体作为可见显示层叠加于透明上下文，球体相互间采用真实球面深度遮挡。</p><p>3840×2160 / 60 fps / 0.30×，约 20 s，所有 301 个保存状态均使用。原数据间隔 20 ms，每个状态固定显示 4 帧，末态 1 帧；没有插值或混合聚集形态。原始保存状态之间的形变、分组改变及显示群体位置调整仍可见，60 fps 编码不等于 60 fps 的独立模拟状态。</p><p>原始数据和科学结论未改动，没有启动新模拟。显示放大不证明局部物理精度；约 10 pm 原求解器球对容差保留。完整显示映射、逐泡投影和核验记录可下载检查。</p>"""
    for case, video in videos.items():
        body += f'<h2>{case} · 完整 6 s 主力场</h2><video controls playsinline preload="metadata" poster="keyframes/{case}_300.png" src="{video["path"]}"></video>'
    body += '<p><a href="DISPLAY_CHECKS.json">显示几何检查</a> · <a href="RENDER_MANIFEST.json">来源与方法</a> · <a href="DELIVERY_VERIFICATION.json">交付核验</a></p></html>'
    (out / "OPEN_RESULTS.html").write_text(body)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--preview-only", action="store_true")
    parser.add_argument("--case", choices=["B2", "B3", "all"], default="all")
    args = parser.parse_args()
    out = args.output.resolve()
    assert not out.exists(), "Preserve earlier deliveries"
    for folder in ["keyframes", "provenance"]:
        (out / folder).mkdir(parents=True)
    runs, protected = previous.selected_runs()
    sources = [
        Path(__file__),
        Path(context.__file__),
        Path(layout.__file__),
        Path(spheres.__file__),
        Path(previous.__file__),
        ROOT / "scripts/stage4l/check_magnified.py",
        ROOT / "tests/test_stage4l_magnified_geometry.py",
    ]
    source_hash = {str(p): sha(p) for p in sources}
    for i, p in enumerate(sources):
        shutil.copy2(p, out / "provenance" / f"{i:02d}_{p.name}")
    reference = Path(
        "/home/lzy/projects/computation_examples/mevo_acoustic_microbubbles/results/stage4/publication_figures/PZT_control_reference_v6"
    )
    for name in ["paper_pzt_reference.py", "PROVENANCE.json", "Fig1a_PZT_control.png"]:
        p = reference / name
        protected[str(p)] = sha(p)
        shutil.copy2(p, out / "provenance" / f"PZT_v6_{name}")
    protected[str(previous.style.FONT)] = sha(previous.style.FONT)
    videos, checks, rows = {}, {}, []
    start = time.perf_counter()
    for case in ["B2", "B3"] if args.case == "all" else [args.case]:
        d, _, info = previous.load(case, "primary", runs["primary"][case], protected)
        wall_path = d["surface"].parent / "mesh-surfaces/WALL.vtp"
        wall = pv.read(wall_path).triangulate()
        protected[str(wall_path)] = sha(wall_path)
        exterior = pv.read(d["surface"])
        nodes = {int(pid): i for i, pid in enumerate(exterior["GlobalNodeID"])}
        assert np.array_equal(wall.points, exterior.points[[nodes[int(pid)] for pid in wall["GlobalNodeID"]]])
        native = spheres.audit_native(d, info)
        samples, groups, report = layout.prepare(d, SCALE)
        np.savez_compressed(out / f"{case}_DISPLAY_SPHERES.npz", **samples)
        np.savez_compressed(
            out / f"{case}_RECORDED_SPHERES.npz",
            positions_m=d["x"],
            radii_m=d["radii"],
            active=samples["active"],
            time_s=d["t"],
            particle_id=d["ids"],
        )
        write_gzip(out / f"{case}_GROUP_TRANSFORMS.json.gz", groups)
        background, project, scene = context.scene(d, wall, SCALE)
        background.save(out / f"{case}_BACKGROUND.png")
        checks[case] = dict(native=native, display=report, scene=scene)
        selected = [0, 10, 50, 96, 97, 100, 150, 152, 153, 200, 250, 300]
        states = selected if args.preview_only else range(301)
        path = out / f"stage4L_{case}_full_vessel_morphology_x24_4k60.mp4"
        encoder = (
            None
            if args.preview_only
            else imageio.get_writer(
                path,
                fps=FPS,
                codec="libx264",
                pixelformat="yuv420p",
                macro_block_size=None,
                ffmpeg_log_level="error",
                output_params=[
                    "-preset",
                    "medium",
                    "-crf",
                    "12",
                    "-threads",
                    "4",
                    "-profile:v",
                    "high",
                    "-level:v",
                    "5.2",
                    "-g",
                    "120",
                    "-movflags",
                    "+faststart",
                ],
            )
        )
        projections = {}
        try:
            for k in states:
                active = samples["active"][k]
                im, evidence = spheres.analytic_spheres(
                    background,
                    None,
                    samples["positions_m"][k, active],
                    samples["radii_m"][active],
                    d["ids"][active],
                    supersampling=4,
                    projector=project,
                )
                assert evidence["submitted_ids"] == d["ids"][active].tolist()
                projections[str(k)] = dict(
                    source=k, largest_cluster_ids=groups[k][0]["particle_ids"], **evidence
                )
                im = annotate(im, d, k, samples, groups, project)
                if k in selected:
                    im.save(out / "keyframes" / f"{case}_{k:03d}.png")
                for frame in range(k * 4, min(k * 4 + 4, FRAMES)):
                    rows.append([case, frame, frame / FPS, frame * 0.005, k, float(d["t"][k])])
                    if encoder is not None:
                        encoder.append_data(np.asarray(im))
                if k % 25 == 0:
                    print(
                        json.dumps(
                            dict(case=case, state=k, total=301, elapsed_s=time.perf_counter() - start)
                        ),
                        flush=True,
                    )
        finally:
            if encoder is not None:
                encoder.close()
        write_gzip(out / f"{case}_SPHERE_PROJECTIONS.json.gz", projections)
        videos[case] = dict(
            path=path.name,
            source=info,
            fps=FPS,
            frames=FRAMES,
            resolution=[context.WIDTH, context.HEIGHT],
            playback_speed=0.3,
            sha256=sha(path) if path.exists() else None,
        )
        write(out / "DISPLAY_CHECKS.json", checks)
    assert all(sha(p) == h for p, h in (protected | source_hash).items())
    with (out / "frame_map.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            ["case", "frame", "video_time_s", "requested_time_s", "source_index", "source_time_s"]
        )
        writer.writerows(rows)
    write(
        out / "RENDER_MANIFEST.json",
        dict(
            status="PREVIEW" if args.preview_only else "RENDERED_PENDING_CHECK",
            videos=videos,
            protected_inputs=protected,
            source_sha256=source_hash,
            scientific_data_modified=False,
            simulations_launched=False,
            display_similarity_scale=SCALE,
            particle_interpolation=False,
            particle_merging=False,
            particle_filtering=False,
            group_rotation=False,
            case_identifiers_in_video=False,
            rigid_display_boundary_translation=True,
            physical_state_cadence_s=0.02,
            rigid_display_intergroup_separation=True,
            intergroup_clearance_m=5e-6,
            sphere_supersampling=4,
            elapsed_s=time.perf_counter() - start,
        ),
    )
    viewer(out, videos)


if __name__ == "__main__":
    main()
