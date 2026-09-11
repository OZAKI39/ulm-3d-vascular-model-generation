"""CPU-only interpretation of closed native failure snapshots; never launches Mirheo."""
import csv
import json
import time
from pathlib import Path

import numpy as np

from common import read, write, sha
from cloud_geometry import read_off, order_vertices
from rbc_checks import frame_check

REPORT_JOB_ID = 'cloud-rbc-full-20260911T000225Z'

def _state(path, shift, gap):
    rows = np.atleast_1d(np.genfromtxt(path, names=True, delimiter=',', dtype=None, encoding='utf-8'))
    xyz = np.column_stack([rows[k] for k in ('x', 'y', 'z')]) + shift
    old = np.column_stack([rows[k] for k in ('old_x', 'old_y', 'old_z')]) + shift
    velocity = np.column_stack([rows[k] for k in ('vx', 'vy', 'vz')])
    finite = all(np.isfinite(rows[k]).all() for k in rows.dtype.names)
    speeds = np.linalg.norm(velocity, axis=1)
    outside = (xyz[:, 2] < -1e-5) | (xyz[:, 2] > gap + 1e-5)
    result = dict(file=path.name, sha256=sha(path), particles=len(rows),
                  unique_ids=len(np.unique(rows['id'])), all_numeric_fields_finite=bool(finite),
                  bounds_effective= [xyz.min(0).tolist(), xyz.max(0).tolist()],
                  old_positions_bounds_effective=[old.min(0).tolist(), old.max(0).tolist()],
                  max_speed=float(speeds.max()), max_speed_id=int(rows['id'][speeds.argmax()]),
                  max_displacement_from_old_positions=float(np.linalg.norm(xyz-old, axis=1).max()),
                  z_wall_violations=int(outside.sum()), z_wall_violation_ids=rows['id'][outside].tolist())
    return result, rows['id'], xyz, velocity


def analyze_failures(archive, destination, checks):
    """Preserve failed-step phase and refuse membrane membership on invalid geometry."""
    start = time.monotonic()
    root, destination = Path(archive), Path(destination)
    spec = read(root/'full_spec.json')
    ref, faces = read_off(root/'reference.off')
    cfg = spec['config']; length = cfg['geometry']['periodic_length']; gap = cfg['geometry']['gap']
    # Native trace stores rank-local coordinates. This frozen task has one
    # spatial rank and domain (L,L,gap+2*pad); remove pad after local->global.
    shift = np.array([length/2, length/2, gap/2])
    edges = np.unique(np.sort(np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]]), axis=1), axis=0)
    reference_lengths = np.linalg.norm(ref[edges[:, 1]]-ref[edges[:, 0]], axis=1)
    failures = []
    for counter in sorted((root/'simulation').glob('bounce_counts_*.csv')):
        phase = counter.stem.removeprefix('bounce_counts_')
        with counter.open() as stream:
            for row in csv.DictReader(stream):
                if int(row['coarse_count']) <= int(row['coarse_capacity']) and int(row['fine_count']) <= int(row['fine_capacity']):
                    continue
                prefix = counter.parent/f"bounce_{phase}_{row['pv']}_{row['locality']}_{row['step']}"
                item = dict(phase=phase, native_step_zero_based=int(row['step']), native_time_label=float(row['time']),
                            observation='after integration, before this bouncer resolves collisions; not a successful run return or a scheduled beforeForces frame',
                            counter=row, fine_candidates='NOT_EVALUATED' if int(row['fine_count']) < 0 else int(row['fine_count']),
                            coordinate_conversion=dict(native='rank-local', effective_shift=shift.tolist(), periodic_axes=['x','y'], wall_axis='z'),
                            confirmed_self_intersections=None, area_relative_drift=None, volume_relative_drift=None)
                membrane_path = Path(str(prefix)+'_membrane.csv')
                if not membrane_path.exists():
                    item['geometry_status']='SNAPSHOT_MISSING'; failures.append(item); continue
                state, ids, xyz, velocity = _state(membrane_path, shift, gap)
                item['membrane'] = state
                vertices = order_vertices(ids, xyz, len(ref))
                delta = vertices[edges[:,1]]-vertices[edges[:,0]]
                delta[:,:2] -= length*np.rint(delta[:,:2]/length)
                lower_bound = np.linalg.norm(delta, axis=1)/reference_lengths*cfg['mirheo_membrane']['x0']
                item['edge_length_evidence'] = dict(max_WLC_extension_lower_bound=float(lower_bound.max()),
                    edges_with_lower_bound_at_least_one=int((lower_bound >= 1).sum()),
                    definition='Each edge minimized independently over X/Y images. Lower bound only when a globally consistent periodic surface cannot be constructed; not a reconstructed valid membrane.')
                try:
                    checked = frame_check(ids, xyz, velocity, faces, ref, spec)
                    item.update(geometry_status='CHECKED', geometry=checked['geometry'],
                                confirmed_self_intersections=checked['intersections']['confirmed_count'],
                                area_relative_drift=checked['area_relative_drift'], volume_relative_drift=checked['volume_relative_drift'],
                                hard_failures=checked['hard_failures'])
                except ValueError as error:
                    item.update(geometry_status='INVALID_PERIODIC_GEOMETRY', geometry_error=str(error),
                                hard_failures=['INVALID_PERIODIC_GEOMETRY'],
                                unavailable_geometry='Self-intersection and A/V at this failure snapshot are indeterminate; no valid periodic unwrapping.')
                if state['z_wall_violations']: item['hard_failures'].append('RAW_Z_MEMBRANE_WALL_CROSSING')
                if (lower_bound >= 1).any(): item['hard_failures'].append('WLC_EXTENSION_LOWER_BOUND_OUTSIDE_DOMAIN')
                fluid_path = Path(str(prefix)+'_fluid.csv')
                if fluid_path.exists(): item['fluid'] = _state(fluid_path, shift, gap)[0]
                item['failure_membership'] = dict(status='NOT_CHECKED', reason='Invalid membrane prevents reliable inside/outside classification; only this outer-PV failure snapshot exists at this call.',
                    strict_impermeability='NOT_VERIFIED', reclassification_performed=False)
                failures.append(item)
    stats = []
    for group in checks['raw_native_statistics']:
        rows = group['rows']
        if not rows: continue
        stats.append(dict(file=group['file'], samples=len(rows), first=rows[0], last=rows[-1],
                          particle_count_values=sorted({int(r['num_particles']) for r in rows}),
                          kBT_min=min(float(r['kBT']) for r in rows), kBT_max=max(float(r['kBT']) for r in rows),
                          max_speed=max(float(r['maxv']) for r in rows), definition=group['definition']))
    resource_rows=[json.loads(line) for line in (root/'resources.jsonl').read_text().splitlines()]
    sampled_pids=sorted({p['pid'] for sample in resource_rows for p in sample['processes']})
    rank_pids=sorted(r['pid'] for r in checks['actual_mpi_ranks'])
    missing_pids=sorted(set(rank_pids)-set(sampled_pids))
    resource_coverage=dict(status='INCOMPLETE' if missing_pids else 'REGISTERED_RANKS_SAMPLED',
                           sampled_pids=sampled_pids,actual_rank_pids=rank_pids,unsampled_rank_pids=missing_pids,
                           rank_cpu_time_and_rss='UNMEASURED' if missing_pids else 'SEE_RAW_SAMPLES',
                           rss_scope='Only observed runner/mpirun processes; not total task memory when rank PIDs are missing. GPU readings are whole-device; cgroup memory covers the entire container.')
    result = dict(failure_snapshots=failures, native_statistics_summary=stats,resource_sampling_coverage=resource_coverage,
                  sampling_limit='Scheduled beforeForces geometry and native afterIntegration Stats have distinct phases. Failure snapshots are separate observations; never inserted into the trajectory or treated as completed endpoints.',
                  root_cause='Native coarse-candidate overflow is confirmed. A severe velocity/geometry blow-up is already present before collision resolution at the failing call. These records do not isolate the original instability trigger, nor prove a hardware, integrator, membrane-force or shared-bouncer root cause.',
                  original_archive_unmodified=True, new_solver_runs=0, cpu_analysis_wall_s=time.monotonic()-start)
    write(destination/'postmortem.json', result)
    return result


def report_text(root, checks, verified, postmortem):
    """Case-specific interpretation of this one authorized failed attempt."""
    root=Path(root); spec=read(root/'full_spec.json'); execution=read(root/'execution.json')
    if root.name != REPORT_JOB_ID:
        raise ValueError('POSTMORTEM_CASE_REPORT_JOB_MISMATCH')
    frames=read(root/'frames.json')['frames']; geometry=checks['geometry']; members=checks['membership']
    timings=checks['timings']; failed=postmortem['failure_snapshots']; config=spec['config']
    lines=[
        '# Mirheo 云端完整单红细胞验证：实际结果',
        '没有完成 Γ=4。准备阶段发生原生碰撞候选溢出并退出，未进入剪切；本轮没有重试。',
        f"CLOUD_RBC_RUN_COMPLETE = {checks['CLOUD_RBC_RUN_COMPLETE']}\nCLOUD_RBC_NUMERICAL_SCREEN = {checks['CLOUD_RBC_NUMERICAL_SCREEN']}\nRESULTS_RETURN_VERIFIED = PASS",
        f"JOB_ID：{root.name}\n本地原始归档：{root}\n云端任务：/workspace/bloodflow/cloud_runs/{root.name}/",
        f"构建 mirheo-sm120-20731713865ae510；实际计算库 {spec['native_library']}；SHA-256 {spec['native_library_sha256']}。两个 MPI rank 的实际加载记录一致，单 GPU、sm_120、原单精度。",
        f"修复 local_before_halo_v2；运行方式每阶段一次连续 u.run；本次顶层 spec 的 bouncer_policy={spec['bouncer_policy']}，同一个 bounce_back bouncer 绑定两侧液体。嵌套 config.repair 的 independent_per_pv_control 是保留的历史配置字段，worker 实际读取顶层值。原生候选容量固定，未扩容。",
        f"从零初始化：液体 110592（outer 109950、inner 642），粒子质量 1；膜 642 顶点 / 1280 三角面。有效域 24³、X/Y 周期、Z 壁面。dt={spec['dt']}，实际 float32 dt 见碰撞原始记录；准备上限 {spec['prep_steps']} 步、正式剪切计划 {spec['steps']} 步（γ̇=0.02，Γ=4）。WLC/Kantor、ks={spec['ks']}、kb={spec['kb']}，DPD kBT=1；完整参数在 actual_parameters.json 与 full_spec.json。",
        f"准备未通过：u.run(60000) 未正常返回，成功返回的准备步数未知（null）。定期完整帧 {len(frames)}，最后准备步标记 {checks['last_saved_progress_lower_bound']['relaxation']}、名义 t*=22.5，只是进度下界。正式剪切 0 步、Γ=0；4000 步壁粒子准备不计入 RBC 正式应变。无成功 completion，也没有准备末态或剪切末态。",
        f"定期保存膜帧：面积最大相对漂移 {100*geometry['maximum_area_relative_drift']:.6f}%，体积 {100*geometry['maximum_volume_relative_drift']:.6f}%；确证非相邻自交帧 {geometry['confirmed_intersection_frames']}，近接触帧 {geometry['uncertain_contact_frames']}，未见定期帧穿墙、退化或 NaN/Inf。该结论仅覆盖这些帧，不包括报错瞬间快照。",
    ]
    for item in failed:
        row=item['counter']; membrane=item.get('membrane',{}); fluid=item.get('fluid',{}); edge=item.get('edge_length_evidence',{})
        lines += [
            f"原生失败：准备阶段零基步标签 {item['native_step_zero_based']}、原生时间标签 {item['native_time_label']:.12g}；{row['pv']}/{row['locality']} 粗候选 {row['coarse_count']} > 容量 {row['coarse_capacity']}。已保存 {row['stored_coarse']} 个候选对，其内部重复 {row['stored_coarse_duplicates']}；不能据此判断被截断部分。fine_count=-1 表示该失败调用尚未执行细筛，不能写成零碰撞。",
            f"故障快照已出现严重膜损坏：{membrane.get('z_wall_violations')} 个顶点的实际 Z 越过壁面，最小有效 Z={membrane.get('bounds_effective',[[None]*3])[0][2]}。膜最大速度 {membrane.get('max_speed'):.6g}，outer 液体最大速度 {fluid.get('max_speed'):.6g}。周期展开失败（{item.get('geometry_error')}）；至少 {edge.get('edges_with_lower_bound_at_least_one')} 条边的最短周期像 WLC 延伸值已≥1，最大下界 {edge.get('max_WLC_extension_lower_bound'):.6g}。",
            '故障快照位于该步积分之后、此 bouncer 碰撞处理之前，未当作正常完成帧。Z 为非周期轴，越墙证据不依赖 X/Y 展开。由于无法构造一致闭合周期曲面，该瞬间自交数、物理 A/V 和成员分类均标为无法可靠判定，不用定期帧的零自交替代。',
        ]
    lines += [
        f"成员检查仅覆盖真正初态：内外各 96 个按排序 ID 等间距抽取的探针，共 {members['tested_point_frames']} 点次；确证不符 {members['confirmed_mismatches']}，近膜不确定 {members['near_surface_uncertain']}。没有故障前完整探针轨迹，不能证明全程或所有粒子不可渗透。strict_impermeability=NOT_VERIFIED。未做事后成员修正，静壁到动壁交接 NOT_REACHED。",
        '速度证据：最后定期帧（beforeForces，45000）膜最大速度约 43.31，之前定期帧最大约 5.99；原生 Stats 在其 afterIntegration 样本 t*=22.5005 记录 inner 最大速度约 94.26、kBT≈5.658，膜最大速度≈59.15。原生 kBT=m·mean(|v|²)/3 包含整体流动动能；不是扣除流速的温度，也不是开尔文。',
        '原生 Stats 各保存 46 行，记录时液体总数/质量不变。局部流场保存 45 个可读、有限的 8³ 原生分箱场；尚无正式剪切响应。因调用未返回，vertices.csv、moments.csv、profiles.csv、local_flow.csv、timings.csv 只有表头；它们不代表零物理量。实际膜轨迹来自原生 HDF5/XMF，派生 geometry.csv 有 46 行，原生统计与局部场保留在 simulation 下。',
        '原因范围：已确认候选容量溢出，且失败调用的碰撞处理之前已存在速度与几何爆发。现有数据不足以确定最初触发来自积分、膜力、共享 bouncer 状态或其他机制；不能仅把增大容量视为修复，也不把该结果归因为 GPU 硬件。',
        f"求解进程退出码 {execution['exit_code']}；墙钟 {execution['solver_process_wall_s']:.6f} 秒（1800 秒授权的一次尝试，已消耗 1/1 次）。包括初始化与必要准备，不等于 GPU 活跃时间或 Vast 租赁费用；剩余额度不自动授权第二次运行。setup 已完成 {timings['setup_s']:.6f} 秒，其中壁粒子准备 {timings['wall_preparation_s']:.6f} 秒，不能重复相加。phase_timings 的 relaxation_s=0 是未返回调用的占位值，实际准备阶段耗时未单独闭合计量，不能解释为零耗时。云端事后 CPU 分析 {checks['analysis_cpu_wall_s']:.6f} 秒另计；本地事后分析时间记录于 postmortem.json。",
        f"资源覆盖不完整：{checks['resources']['samples']} 次进程采样只包含 runner 7455 与 mpirun 7458，漏掉实际 MPI rank 7495/7496。因此所见 RSS 峰值 {checks['resources']['peak_registered_process_rss_bytes']} 字节仅为这两个管理进程，不能当作整个求解任务峰值；rank 的 CPU 时间与 RSS 未测。实际 rank 身份由各自 rank_0.json/rank_1.json 证实。GPU 为全设备采样，观测显存峰值 571 MiB；cgroup 为全容器。事后只读复核四个登记 PID 均已不存在。该采样遗漏仅记录为限制，本轮没有为补测重跑。",
        f"回传清单通过：{verified['files']} 文件、{verified['bytes']} 字节（不含清单/校验收据自身）；SHA-256 {verified['manifest_sha256']}。RESULTS_READY、清单、JSON/CSV 与原始 HDF5 读取核查通过。远端 report_zh.md / numerical_checks.json 中回传 PENDING 是封包前状态；本地 LOCAL_ARCHIVE_VERIFIED.json 和本派生报告给出实际 PASS。",
        '原始云端封包及本地归档保持不变；本地补充分析、HTML 和浏览器证据位于同级 -review 目录。浏览器实际检查见 browser_final/browser_check.json（仅该文件及对应 HTML 哈希一致且状态 PASS 时才算通过），人工验收 PENDING。',
        'PHYSICAL_MODEL_VALIDATION：既有材料 NOT_MATCHED；本轮材料标定、空间/时间收敛及 HemoCell 比较 NOT_TESTED。qualified_speedup=null。',
        '没有重新编译、安装依赖、修改物理参数/原生修复、运行 HemoCell、推送 Git、删除云端数据或停止 Vast 实例。',
        f"查看/下载不会启动求解：\ncd /home/lzy/projects/cloud_compute\npython3 -B cloud_run.py status {root.name}\npython3 -B cloud_run.py fetch {root.name}\npython3 -B cloud_run.py verify {root.name}\npython3 -B cloud_run.py rbc-review {root.name}\n原 rbc-full --preflight-only 和 rbc-full --new --execute 已实现；本轮授权与启动收据已消费，不能用执行命令重复本轮任务。",
    ]
    return '\n\n'.join(lines)+'\n'
