"""Read-only audit separating recorded rearrangement from rendering and rigid motion."""

import argparse
import json
from pathlib import Path
import shutil
import zipfile

import h5py
import numpy as np
from scipy.spatial.distance import pdist, squareform

from scripts.stage4.render_motion import ROOT, sha


def distribution(values):
    a = np.asarray(values, dtype=float)
    return dict(
        minimum=float(a.min()),
        median=float(np.median(a)),
        p95=float(np.percentile(a, 95)),
        maximum=float(a.max()),
    )


def audit(output):
    output.mkdir(parents=True, exist_ok=True)
    source = ROOT / "results/stage4/visualization/B2_B3_spheres_uniform_v6/RENDER_MANIFEST.json"
    manifest = json.loads(source.read_text())
    protected = {str(source): sha(source)}
    result = dict(
        status="READ_ONLY_AUDIT_COMPLETE",
        interval_s=[3, 6],
        original_state_interval_s=0.02,
        cases={},
        protected_inputs=protected,
    )
    for name in ["B2", "B3"]:
        run = ROOT / manifest["source_runs"][name]
        paths = [
            run / "particles/trajectories.h5",
            run / "particles/SUMMARY.json",
            run / "PHYSICS_BANNER.json",
            run / "logs/source_snapshot.zip",
        ]
        for path in paths:
            protected[str(path)] = sha(path)
        trajectory = paths[0]
        assert sha(trajectory) == manifest["protected_inputs"][str(trajectory)]
        with h5py.File(trajectory) as h:
            t, x, ids, radii = (h[k][:] for k in ["time_s", "positions_m", "particle_id", "radii_m"])
            saved_dt = h["accepted_dt_s"][:]
            forces = {
                k: h[k][:]
                for k in [
                    "primary_force_N",
                    "secondary_force_N",
                    "flow_velocity_m_s",
                    "free_total_velocity_m_s",
                    "total_velocity_m_s",
                ]
            }
            assert h.attrs["complete"] and not h.attrs["all_accepted_substeps_recorded"]
        assert np.all(radii == radii[0]) and np.isclose(radii[0], 1.78e-6)
        states = manifest["states"][name]
        id_index = {int(pid): i for i, pid in enumerate(ids)}
        by_members = [{tuple(sorted(g["particle_ids"])): g for g in s["overview_groups"]} for s in states]
        summary = json.loads(paths[1].read_text())
        banner = json.loads(paths[2].read_text())
        # Read the frozen run source, not a newly edited numerical implementation.
        source_evidence = {}
        with zipfile.ZipFile(paths[3]) as archive:
            for member in [
                "src/mevo_mb/bubbles/tracker.py",
                "src/mevo_mb/bubbles/contact.py",
                "src/mevo_mb/coupling/bifurcation_experiment.py",
            ]:
                lines = archive.read(member).decode().splitlines()
                source_evidence[member] = [
                    dict(line=i + 1, text=line.strip())
                    for i, line in enumerate(lines)
                    if any(
                        token in line
                        for token in [
                            "free = flow +",
                            "constrain_pair_velocities(",
                            "corrected, wall, contacts",
                            "geometric_dt =",
                            "Geometric nonpenetration only",
                            "amount = max(0.",
                            "x[i] -=",
                            "x[j] +=",
                        ]
                    )
                ]
        case = dict(
            result_id=manifest["source_runs"][name],
            radius_um=float(radii[0] * 1e6),
            individual_shape="Fixed spherical glyphs; no nonspherical deformation or fusion in this model",
            model_banner=banner,
            frozen_source_evidence=source_evidence,
            full_run_numerical_summary={
                k: summary[k]
                for k in ["accepted_steps", "minimum_dt_s", "contact_pair_events", "maximum_projection_m"]
            },
            full_run_summary_scope="Counters and maximum correction cover 0–6 s, not specifically 3–6 s",
            saved_last_substep_s=distribution(saved_dt[150:]),
            timestep_scope="Last accepted substep attached to each stored frame, not all accepted substeps",
            cohorts=[],
        )
        for group in states[150]["overview_groups"]:
            if group["count"] < 2:
                continue
            members = tuple(sorted(group["particle_ids"]))
            indices = [id_index[pid] for pid in members]
            present = [k for k in range(301) if members in by_members[k]]
            matched = [
                k for k in range(150, 300) if members in by_members[k] and members in by_members[k + 1]
            ]
            if not matched:
                continue
            records = []
            for k in matched:
                a, b = x[k, indices], x[k + 1, indices]
                ca, cb = a.mean(0), b.mean(0)
                ac, bc = a - ca, b - cb
                u, _, vt = np.linalg.svd(ac.T @ bc)
                rotation = u @ np.diag([1.0, 1.0, np.linalg.det(u @ vt)]) @ vt
                pa, pb = pdist(ac), pdist(bc)
                touching_a = pa <= 2 * radii[0] + 0.5e-6
                touching_b = pb <= 2 * radii[0] + 0.5e-6
                union = np.count_nonzero(touching_a | touching_b)
                records.append(
                    dict(
                        from_s=float(t[k]),
                        to_s=float(t[k + 1]),
                        centroid_displacement_um=float(np.linalg.norm(cb - ca) * 1e6),
                        relative_particle_rms_um=float(
                            np.sqrt(np.mean(np.sum((bc - ac) ** 2, axis=1))) * 1e6
                        ),
                        after_best_rigid_fit_rms_um=float(
                            np.sqrt(np.mean(np.sum((ac @ rotation - bc) ** 2, axis=1))) * 1e6
                        ),
                        same_ID_pair_distance_change_rms_um=float(np.sqrt(np.mean((pa - pb) ** 2)) * 1e6),
                        unordered_pair_distribution_change_rms_um=float(
                            np.sqrt(np.mean((np.sort(pa) - np.sort(pb)) ** 2)) * 1e6
                        ),
                        near_neighbor_turnover_fraction=float(
                            1 - np.count_nonzero(touching_a & touching_b) / union
                        )
                        if union
                        else 0.0,
                        source_radius_of_gyration_um=float(np.sqrt(np.mean(np.sum(ac**2, axis=1)))) * 1e6,
                        unlabeled_cloud_nearest_rms_um=float(
                            np.sqrt(
                                np.mean(
                                    np.min(np.linalg.norm(ac[:, None] - bc[None, :], axis=-1), axis=1) ** 2
                                )
                            )
                        )
                        * 1e6,
                    )
                )
            valid_frames = [k for k in range(150, 301) if members in by_members[k]]
            force_summary = {
                key: distribution(np.linalg.norm(a[valid_frames][:, indices], axis=-1))
                for key, a in forces.items()
            }
            nearest = []
            for k in valid_frames:
                matrix = squareform(pdist(x[k, indices]))
                np.fill_diagonal(matrix, np.inf)
                nearest.extend(matrix.min(axis=1) * 1e6)
            case["cohorts"].append(
                dict(
                    particle_ids=list(members),
                    count=len(members),
                    first_exact_component_s=float(t[present[0]]),
                    states_present_in_3to6s=len(valid_frames),
                    adjacent_intervals=len(matched),
                    members_unchanged_entire_3to6s=len(valid_frames) == 151,
                    centroid_at_3s_mm=(np.array(group["center_m"]) * 1000).tolist(),
                    metrics={
                        key: distribution([r[key] for r in records])
                        for key in records[0]
                        if key not in ["from_s", "to_s"]
                    },
                    nearest_neighbor_um=distribution(nearest),
                    stored_instantaneous_fields=force_summary,
                    force_scope="Magnitudes at original stored states; velocities are instantaneous constrained/free values, not finite-difference average motion",
                    adjacent_state_records=records,
                )
            )
        result["cases"][name] = case
    result["findings"] = dict(
        recorded="Largest B2/B3 cohorts keep identical membership throughout 3–6 s, while relative coordinates change in the original files. Rigid translation/rotation alone cannot explain those changes.",
        shape_vs_identity="Same-ID particle displacement is not whole-outline displacement. The largest clouds remain comparably compact: median unordered pair-distance distribution changes are only about 0.12/0.10 um, while centered nearest-cloud distances are about 1.93/1.92 um for B2/B3. Internal identity rearrangement is much larger than the change in overall extent.",
        display="60x within-group magnification makes micrometre rearrangements conspicuous. 20 ms recorded states appear at 15 states/s under constant 0.30x playback, cross-dissolved into 60 fps. No special rendering branch starts at 3 s.",
        model="Cycle-averaged forces continue driving independent spheres subject to geometric nonpenetration. No adhesion/fusion model makes a permanently rigid aggregate.",
        caution="Near-contact dynamics are explicitly unvalidated. Strong near-contact forcing and repeated geometric projections can contribute numerical rearrangement/jitter; this audit does not establish that they are the unique cause or that the movement is physically accurate.",
        not_resolved="No all-substep motion/projection record or newly completed timestep convergence is available for this late-time motion in this audit. Real physical dynamics cannot be separated quantitatively from numerical artifacts here.",
        action="No trajectory, force, timestep, radius or group shape was altered or smoothed. Only layout changed in the animation.",
    )
    assert all(sha(path) == digest for path, digest in protected.items())
    (output / "SHAPE_AUDIT.json").write_text(json.dumps(result, indent=2) + "\n")
    rows = []
    for name, case in result["cases"].items():
        for cohort in case["cohorts"][:2]:
            mm = cohort["metrics"]
            rows.append(
                f"<tr><td>{name}</td><td>{cohort['count']}</td><td>{cohort['states_present_in_3to6s']}/151</td>"
                f"<td>{mm['centroid_displacement_um']['median']:.3f}</td>"
                f"<td>{mm['relative_particle_rms_um']['median']:.2f}</td>"
                f"<td>{mm['after_best_rigid_fit_rms_um']['median']:.2f}</td>"
                f"<td>{mm['unlabeled_cloud_nearest_rms_um']['median']:.2f}</td></tr>"
            )
    text = """<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>3–6 s 泡群形状变化核查</title><style>body{font:17px/1.7 Arial,sans-serif;max-width:1150px;margin:32px auto;padding:0 24px;background:#080c12;color:#e5edf5}h1,h2{font-weight:400}table{border-collapse:collapse;width:100%}td,th{padding:10px;border-bottom:1px solid #344455;text-align:left}a{color:#66e0c2}</style><h1>为什么 3000 ms 后泡群频繁换形？</h1><p>变化已经存在于原始保存坐标中，并非每帧重新随机摆球。B2 最大泡群的 75 个成员、B3 最大泡群的 91 个成员在 3–6 s 完全不变；球半径也始终是 1.78 µm。发生变化的是这些球之间的相对位置。</p><p>下表统计相邻 20 ms 保存状态，数值均为这一时间段内的中位数。先去掉质心平移，再做逐 ID 的最佳刚体旋转配准；仍很大的残差说明变化不能仅用整体转动解释。</p><table><tr><th>案例</th><th>泡数</th><th>同一群体出现状态数</th><th>质心位移 µm</th><th>群内相对位移 RMS µm</th><th>去除刚体运动后 RMS µm</th><th>忽略身份的点云近邻 RMS µm</th></tr>"""
    text += (
        "".join(rows)
        + """</table><p>这里约 10 µm 的数值是同一微泡的内部位移，不是整个泡群轮廓摆动 10 µm。忽略身份后，最大泡群的相邻点云近邻差异中位数约为 B2 1.93 µm、B3 1.92 µm；无序两两距离分布变化仅约 0.12/0.10 µm。也就是说，整体紧密程度相近，内部成员和表面排布反复变化。点云近邻指标只去掉平移，不去掉转动，也不等同于轮廓误差。</p><p>B2 的最终 75 泡成员组合首次出现于 1.98 s，B3 的 91 泡组合首次出现于 2.82 s。约 3 s 后，最大的密集群体已形成，内部运动更容易在放大的球体画面中被注意到。</p><h2>为何看起来特别明显</h2><p>动画对群内结构和球半径共同做了 60 倍显示放大，微米级变化因此变得醒目。原始文件每 20 ms 才保存一次；在恒定 0.30× 播放时，每秒只有 15 个原始状态，60 fps 的短渐变会在相邻形状之间呈现短暂重影。这与视频编码缺帧是两回事，3000 ms 处没有切换速度、随机排布或更换形状的渲染分支。</p><h2>目前能作出的物理解读</h2><p>模型把微泡作为各自受流场、周期平均主/次级声辐射力驱动的球体，接触只施加几何不穿透，没有黏附、融合或把泡群锁成刚体的机制。原始结果确实持续重排，但不能据此认定真实微泡群会以相同频率变形。</p><p>原运行明确标记近接触动力学未验证；强近接触作用与反复几何投影可能带来数值抖动。全程最大单次位置修正约为 B2 10.75 µm、B3 11.85 µm（这是 0–6 s 全程统计，不能直接归因到某个 3 s 后的画面）。本次仅审计已有数据，没有全部接受子步记录，也没有为这段运动新增时间步收敛实验，因此不能定量分离真实运动和数值误差。</p><p>新版仅调整布局，没有为了消除换形而锁定、平滑或修改科学轨迹。</p><p><a href="SHAPE_AUDIT.json">完整数值、来源哈希及逐状态记录</a></p></html>"""
    )
    (output / "SHAPE_AUDIT.html").write_text(text, encoding="utf-8")
    shutil.copy2(__file__, output / "audit_cluster_shape.py")
    print(
        json.dumps(
            {
                name: dict(count=d["cohorts"][0]["count"], metrics=d["cohorts"][0]["metrics"])
                for name, d in result["cases"].items()
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    audit(parser.parse_args().output.resolve())
