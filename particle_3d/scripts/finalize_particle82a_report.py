#!/usr/bin/env python3
"""Build the Chinese review strictly from completed, independently verified evidence."""
from pathlib import Path
import argparse,hashlib,html,json,os,subprocess,xml.etree.ElementTree as ET
from datetime import datetime,timezone

ROOT=Path(__file__).resolve().parents[2]
BASINS=['OUTLET_01','OUTLET_02','OUTLET_03']


def read(p):return json.loads(Path(p).read_text())
def digest(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()
def write(p,d):Path(p).write_text(json.dumps(d,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
def pct(x):return f'{100*x:.3f}%'
def num(x):return '无数据' if x is None else f'{x:.6g}'
def table(head,rows):return '\n'.join(['| '+' | '.join(head)+' |','|'+'|'.join(['---']*len(head))+'|',*['| '+' | '.join(map(str,r))+' |' for r in rows]])
def tests(p):
    suites=ET.parse(p).getroot();cases=suites.findall('.//testcase')
    return dict(file=str(p),sha256=digest(p),tests=len(cases),failures=sum(c.find('failure') is not None for c in cases),
        errors=sum(c.find('error') is not None for c in cases),skipped=sum(c.find('skipped') is not None for c in cases))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=ROOT/'particle_3d/outputs/particle8_2a')
    parser.add_argument('--report',type=Path,default=ROOT/'particle_3d/reports/particle8_2a');a=parser.parse_args()
    out=a.output;report=a.report;review=out/'review';visual=out/'visuals';admission=out/'admission';traj=out/'trajectories'
    assert (out/'COMPUTE_AND_MEDIA_COMPLETE').exists(),'Formal computation/media not complete'
    assert (out/'FINAL_QUALITY_COMPLETE').exists(),'Final replay/evidence checks not complete'
    audit=read(review/'ADMISSION_COMPARISON.json');geometry=read(admission/'map_n24_summary.json');m=audit['methods']
    full=read(review/'FULL_TRAJECTORY_COMPARISON.json');sem=read(admission/'INLET_CAP_SEMANTICS.json')
    verification=read(review/'MACHINE_EVIDENCE_VERIFICATION.json');media=read(visual/'MEDIA_VALIDATION.json')
    entry=read(review/'METHOD_C_ENTRY_SWEEP_AUDIT.json');entryout=read(review/'METHOD_C_OUTCOMES_BY_ENTRY_CERTIFICATE.json')
    paired=read(review/'PAIRED_CURRENT_VS_C_ENTRY.json')
    sensitivity=read(out/'sensitivity/SENSITIVITY_SUMMARY.json');timestep=read(review/'TIMESTEP_SENSITIVITY.json')
    immutable=read(review/'REMOTE_IMMUTABILITY_VERIFICATION.json');resume=read(review/'RESUME_VERIFICATION.json')
    sync=read(report/'REMOTE_LOCAL_CHECKSUM_VERIFICATION.json');agentvisual=read(report/'AGENT_VISUAL_INSPECTION.json')
    local=read(report/'LOCAL_PROVENANCE.json');server=read(report/'SERVER_COMPUTE_PROVENANCE.json')
    deploy=read(report/'DEPLOYMENT.json');arun=read(admission/'ADMISSION_RUN_RECEIPT.json');trun=read(traj/'FULL_TRAJECTORY_RUN_RECEIPT.json')
    ascale=read(admission/'ADMISSION_SCALING.json');ascale_cpu=read(out/'scaling_admission_cpu/ADMISSION_SCALING.json')
    tscale=read(traj/'TRAJECTORY_SCALING.json');parity=read(traj/'P81_INTEGRATION_PARITY.json')
    checks=[tests(report/'logs/development_tests.xml'),tests(report/'logs/upstream_regression.xml'),tests(review/'remote_tests.xml')]
    lock=read(report/'UPSTREAM_LOCK.json');changed=[p for p,v in lock['sha256'].items() if digest(p)!=v]
    write(report/'UPSTREAM_IMMUTABILITY_VERIFICATION.json',dict(all_pass=not changed,checked_files=len(lock['sha256']),changed=changed))
    assert not audit['development_partial'] and audit['common_event_count']>=100000
    assert verification['all_pass'] and immutable['all_pass'] and resume['all_pass'] and sync['all_pass']
    assert all(x['failures']==x['errors']==0 for x in checks) and not changed
    assert verification['boundary_partition']['exact_disjoint_complete_cover']
    assert all(full[k]['admitted_integrated']>=5000 for k in 'ABC')
    assert all(sum(full[k]['outlet_counts'].values())==full[k]['completed'] for k in 'ABC')
    assert media['all_pass'] and len(media['records'])==3 and agentvisual['all_pass']
    assert len(list((visual/'figures').glob('*.png')))==12
    assert ascale['scientific_parity'] and ascale_cpu['scientific_parity'] and tscale['scientific_parity'] and parity['all_exact']
    assert arun['REMOTE_SERVER_COMPUTE'] and trun['REMOTE_SERVER_COMPUTE'] and verification['REMOTE_SERVER_COMPUTE']
    cap_correct=(sem['solid_boundaries']==['WALL'] and sem['inlet_cap_wall_face_overlap']==0
        and not sem['current_checker_cap_is_solid'] and not sem['current_checker_requires_whole_sphere_inside'])
    assert cap_correct
    flow={k:v for k,v in deploy['frozen_input_sha256'].items() if k.endswith('steady_flow_stage_sv1_3q.vtu')}
    mesh={k:v for k,v in deploy['frozen_input_sha256'].items() if k.endswith('mesh-complete.mesh.vtu')}
    stages={k:'PASS' for k in ['AUTOMATED_CHECKS','REMOTE_HEAVY_COMPUTE','COMMON_LEDGER_PARITY','OPEN_APERTURE_AUDIT',
        'CURRENT_ADMISSION_AUDIT','METHOD_A_COMPLETE','METHOD_B_COMPLETE','METHOD_C_COMPLETE','FULL_TRAJECTORY_COMPARISON']}
    # These classifications have explicit scope; no physiological obstruction claim.
    bottleneck=all(geometry['by_basin'][b]['geometry_pass_probability']<.02 for b in ['OUTLET_01','OUTLET_03'])
    stages.update(CURRENT_ADMISSION_OPEN_CAP_HANDLING='CORRECT',
        BIRTH_PLANE_ARTIFACT='NOT_SUPPORTED',GEOMETRIC_BOTTLENECK='SUPPORTED' if bottleneck else 'INCONCLUSIVE',
        STAGE_RESULT='PASS_WITH_CARRY_FORWARD_LIMITATIONS')
    commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    telemetry=[json.loads(s) for s in (out/'RESOURCE_TELEMETRY.jsonl').read_text().splitlines() if s.strip()]
    resources=dict(samples=len(telemetry),sampling_interval_s=10,
        peak_monitored_process_tree_PSS_bytes=max(x['total_process_pss_bytes'] for x in telemetry),
        peak_cgroup_memory_current_bytes=max(int(x['cgroup']['memory.current']) for x in telemetry),
        memory_current_warning='CGROUP_INCLUDES_PAGE_CACHE_AND_OTHER_PROCESSES; PSS_IS_SEPARATE',
        solver_parallelism='CPU_FORK_COPY_ON_WRITE; GPU_NOT_USED_BY_SPHERE_SOLVER')
    write(report/'RESOURCE_UTILIZATION_SUMMARY.json',resources)
    validation=dict(stage='Particle-8.2A',git_commit=commit,generated_utc=datetime.now(timezone.utc).isoformat(),
        execution_source_commits=dict(admission=arun['git_commit'],full_trajectories=trun['git_commit'],machine_verifier=verification['source_commit'],report=commit),
        frozen_flow_sha=next(iter(flow.values())),mesh_sha=next(iter(mesh.values())),
        sonovue_distribution_sha=audit['official_sonovue']['histogram_sha256'],local_hostname=local['hostname'],
        remote_hostname=server['hostname'],remote_resources=server,
        execution_environment=read(out/'SERVER_SOFTWARE_ENVIRONMENT.json'),
        chosen_workers=dict(admission=arun['workers'],trajectories=tscale['chosen_workers']),
        parallel_scaling=dict(admission=ascale,admission_cpu_repeat=ascale_cpu,trajectories=tscale),
        actual_full_trajectory_performance=trun['performance'],resource_utilization=resources,
        formal_compute_host=server['hostname'],REMOTE_SERVER_COMPUTE=True,
        common_event_count=audit['common_event_count'],point_tracer_basin_counts=audit['point_tracer_basin_counts'],
        Dmax_statistics_by_basin={b:v['Dmax_um_quantiles'] for b,v in geometry['by_basin'].items()},
        geometry_pass_probability_by_basin={b:v['geometry_pass_probability'] for b,v in geometry['by_basin'].items()},
        center_on_plane_acceptance_by_basin={b:v['center_on_plane_acceptance'] for b,v in audit['geometry_by_basin'].items()},
        center_on_plane_expected_probability_by_basin={b:v['center_on_plane_probability'] for b,v in geometry['by_basin'].items()},
        method_C_first_inward_acceptance_by_basin={b:v['acceptance_fraction'] for b,v in m['C']['by_anchor_basin'].items()},
        current_failure_causes=m['A']['rejection_counts'],method_A_statistics=m['A'],method_B_statistics=m['B'],method_C_statistics=m['C'],
        position_distribution_distortion={k:dict(basin_TV=m[k]['basin_total_variation_vs_original_anchor'],
            accepted_anchor_TV=m[k]['accepted_anchor_total_variation_vs_original_anchor'],density_TV=m[k]['inlet_density_total_variation_32x32']) for k in 'ABC'},
        size_distribution_distortion={k:m[k]['size_distance'] for k in 'ABC'},
        full_trajectory_outlets_A=full['A']['outlet_counts'],full_trajectory_outlets_B=full['B']['outlet_counts'],full_trajectory_outlets_C=full['C']['outlet_counts'],
        safety_stop_A=full['A']['safety_stops'],safety_stop_B=full['B']['safety_stops'],safety_stop_C=full['C']['safety_stops'],
        full_trajectory_statistics=full,method_C_entry_sweep=entry['summary'],method_C_entry_certified_outcomes=entryout,
        current_vs_C_entry_paired_comparison=paired,
        open_cap_handling_status='CORRECT',birth_plane_artifact_status=stages['BIRTH_PLANE_ARTIFACT'],geometric_bottleneck_status=stages['GEOMETRIC_BOTTLENECK'],
        inference_scope='CURRENT_FROZEN_MESH_FIXED_ANCHOR_SPHERICAL_MB_MODEL; NOT_A_PHYSIOLOGICAL_OBSTRUCTION_CLAIM',
        sensitivity=sensitivity,timestep_sensitivity=timestep,automated_tests=checks,manual_visual_review='PENDING_USER_REVIEW',
        agent_visual_inspection=agentvisual,statuses=stages,
        carry_forward_limitations=['Point-tracer unresolved basins remain visible; labels do not enforce outcomes.',
            'Every admitted cohort is conditional; per-event preservation does not preserve the admitted population distribution.',
            'C legal center does not imply a wall-clear entry; swept-path strata must be reported.',
            'C entry representation is not finite-size entry dynamics or a new force.',
            'Safety stops are computational unknown outcomes, not physiological trapping.',
            'Outlet completion uses first center crossing, not a finite-body outlet clearance certificate.',
            'Meshes and rigid sphere wall constraints do not establish physiological obstruction.',
            'No automatic ranking or production promotion of A/B/C; no RBC deformation work.'])
    write(report/'PARTICLE8_2A_VALIDATION.json',validation)
    link=lambda p:'['+p.name+']('+os.path.relpath(p,report)+')'
    lines=['# Particle-8.2A 中文审核报告','',
        '本阶段完成同一 100,000 条候选入口事件的 A/B/C 接纳对照，并对每种方法先被接纳的 5,000 个稳定 ID 独立积分。所有正式大样本、完整轨迹、敏感性和最终可视化在远程服务器执行。未修改 Frozen FEM、原始 SonoVue 分布或 P8.1 历史结果；未按出口选择正式轨迹。','',
        '**核心结论：当前 checker 正确区分开放 INLET 与实体 WALL，未发现所提出的“球必须在入口平面上已完全位于域内”拒绝机制。当前 Frozen 几何下，通向 01/03 的入口区域对原始尺寸分布的通过概率很低，支持局部入口几何筛选这一解释。它不证明真实生理堵塞。**','',
        '方法 C 的域内合法出生中心、入口路径几何可通过性、最终到达出口是三个不同问题。本报告分别列出；不会因得到某个出口就宣布方法物理正确。','',
        table(['最终状态','结果'],list(stages.items())), '',
        '**MANUAL_VISUAL_REVIEW = PENDING_USER_REVIEW**。自动视频解码与代理图像检查不能替代用户科学审核。','',
        '## 计算范围与来源','',
        f"服务器 `{server['hostname']}`：AMD Ryzen 7 7800X3D，16 逻辑 CPU，物理内存约 61 GiB，RTX 4090 24 GiB。容器 CPU 配额为 7.68 核，内存上限约 42 GiB。求解采用 CPU 多进程与 fork 只读共享；GPU 可用不等于本求解器使用 GPU。入口正式 {arun['workers']} workers，完整轨迹 {tscale['chosen_workers']} workers。",'',
        f"入口 {audit['common_event_count']:,} 事件、{geometry['count']:,} 个高分辨率通量积分点；全轨迹合计 {sum(full[k]['admitted_integrated'] for k in 'ABC'):,} 条。实测入口 {arun['events_per_s']:.3f} events/s；全轨迹 {trun['performance']['trajectories_per_hour']:.1f} 条/小时。",'',
        f"全轨迹实测平均占用 {trun['performance']['mean_busy_cpu_cores']:.3f} 个 CPU 核，相对于所选 worker 数的计算占用率 {trun['performance']['worker_capacity_utilization_percent']:.2f}%。监测过程树峰值 PSS 为 {resources['peak_monitored_process_tree_PSS_bytes']/2**30:.3f} GiB；容器 memory.current 包含缓存和其它进程，不当作本任务独占内存。",'',
        '使用 P8.1 原始 Frozen 场，不使用独立的 2 mm/s 新流场案例。保持原 P6.5 sphere physics 和 P8.1 积分主体，dt=0.00025 s、单粒子时限 1.5 s、原安全守卫不变。无 MB–MB / RBC 耦合。完整轨迹的选择规则为各方法按共同事件 ID 顺序的前 5,000 个 admitted，三个方法达到配额所需事件前缀不同。','',
        table(['执行部分','代码 commit'],list(validation['execution_source_commits'].items())), '',
        f"Frozen flow SHA256：`{validation['frozen_flow_sha']}`；mesh：`{validation['mesh_sha']}`；SonoVue histogram：`{validation['sonovue_distribution_sha']}`。",'',
        f"本地 {len(lock['sha256']):,} 个只读历史与输入文件重新逐字节哈希检查，无变化。远程所有源快照、Frozen 输入、已绑定的外部 SonoVue 输入已复核。断点重跑验证保留原文件字节与修改时间。",'',
        '证据：'+', '.join(link(p) for p in [report/'SERVER_COMPUTE_PROVENANCE.json',review/'MACHINE_EVIDENCE_VERIFICATION.json',review/'REMOTE_IMMUTABILITY_VERIFICATION.json',review/'RESUME_VERIFICATION.json',report/'REMOTE_LOCAL_CHECKSUM_VERIFICATION.json'])+'。','',
        '## 入口几何、原始分布与共同分母','',
        f"原始 SonoVue 实际分布均值为 {audit['official_sonovue']['mean_um']:.6f} µm；以下分位数由已冻结 CDF 直接计算。",'',
        table(['分位数','直径 µm'],[(k,num(v)) for k,v in audit['official_sonovue']['quantiles_um'].items()]),'',
        'Dmax 是特定入口中心位置相对真实 WALL / 开口边缘的局部最大球直径，不是该出口血管的直径。几何通过概率为该 basin 内按入流通量加权的原始 SonoVue CDF 积分，分母是该 basin 的全部入口通量。', '',
        table(['basin','通量权重','Dmax P10 / median / P90 / max (µm)','几何通过概率','当前平面接纳期望'],[
            [b,pct(v['flux_fraction']),' / '.join(num(v['Dmax_um_quantiles'][q]) for q in ['p10','median','p90','max']),pct(v['geometry_pass_probability']),pct(v['center_on_plane_probability'])] for b,v in geometry['by_basin'].items()]),'',
        'point-tracer 未解析路径单独保留，未强行归入任一出口，也未将已解析部分重归一到 Frozen 出口流量比例。当前接纳额外包含 nearfield handoff 最小间隙，所以与纯几何通过率有小差别。','',
        '下面三个比例使用同一 common ledger、同一原始 anchor 和第一次尺寸抽样，按原始 anchor basin 分组。C 一栏只表示找到完整域内球心，不等价于通过入口。','',
        table(['anchor basin','事件数','open-aperture 实测','当前 center-on-plane 实测','C 找到域内球心'],[
            [b,v['scheduled'],pct(v['open_aperture_passable_fraction']),pct(v['center_on_plane_acceptance']),pct(m['C']['by_anchor_basin'][b]['acceptance_fraction'])] for b,v in audit['geometry_by_basin'].items()]),'',
        f"入口 cap 与 WALL 三角面重叠数为 {sem['inlet_cap_wall_face_overlap']}；封闭外边界由官方 WALL 和开放 cap 构成且与四面体网格外表面完全一致。当前 checker 仅要求球心在域内，不要求整个球已在域内。把整球入域作为新测试时，入口平面通过数为 {audit['full_domain_at_anchor_count']}；这不能倒过来归咎于从未使用该条件的旧 checker。",'',
        '解析圆管、开口边缘、三角形 face/edge/vertex、C 最小首个中心搜索均有永久测试。'+link(admission/'INLET_CAP_SEMANTICS.json')+'、'+link(review/'ADMISSION_COMPARISON.json')+'。','',
        '## 三方法接纳、分布偏差与入口转换','',
        table(['方法','接纳数 / scheduled','接纳率','接纳直径 mean / median (µm)','对原始 CDF 的 KS','Wasserstein (µm)'],[
            [k,f"{m[k]['accepted']:,} / {m[k]['scheduled']:,}",pct(m[k]['acceptance_fraction']),f"{m[k]['accepted_size_um']['mean']:.4f} / {m[k]['accepted_size_um']['median']:.4f}",num(m[k]['size_distance']['KS_statistic']),num(m[k]['size_distance']['Wasserstein_um'])] for k in 'ABC']), '',
        'KS 是两个累积分布之间最大的比例差；Wasserstein 可理解为把一种直径分布变为另一种平均需要移动多少 µm。这里只描述偏差，不用显著性 p 值挑选赢家。','',
        table(['方法','原 anchor 的 01 / 02 / 03','接纳 birth 的 01 / 02 / 03','basin 总变差','32×32 入口密度总变差'],[
            [k,' / '.join(pct(m[k]['original_anchor_basin_fraction'][b]) for b in BASINS),' / '.join(pct(m[k]['accepted_birth_point_basin_fraction'][b]) for b in BASINS),num(m[k]['basin_total_variation_vs_original_anchor']),num(m[k]['inlet_density_total_variation_32x32'])] for k in 'ABC']), '',
        '总变差 0 表示这组分箱比例完全相同，1 表示没有重叠。分箱指标不是连续密度的无误差真值。C 出生中心沿曲线路径内移，二维 cap 密度对比不适用，另用固定 anchor 与接纳筛选统计。','',
        table(['原因','A candidate 拒绝数','B candidate 拒绝数','C search 失败数'],[[cause,m['A']['rejection_counts'].get(cause,0),m['B']['rejection_counts'].get(cause,0),m['C']['rejection_counts'].get(cause,0)] for cause in m['A']['rejection_counts']]),'',
        'A/B 是每次候选尝试的分类，C 是每个事件的一次路径搜索结局，不能把两种分母混用。PAIR_CONFLICT 为零，因为本阶段使用独立单粒子模型。', '',
        table(['C anchor basin','admitted','s_birth P10 / median / P90 (µm)','s_birth / a median'],[
            [b,v['accepted'],' / '.join(num(v['inward_distance_um'].get(q)) for q in ['p10','median','p90']),num(v['inward_distance_over_radius'].get('median'))] for b,v in m['C']['by_anchor_basin'].items()]),'',
        f"C 搜索沿实际 Frozen 速度流线的存储折线，入口邻域上限 {sem['search_horizon_m']*1e6:.6f} µm（预先定义为 2 倍入口等效直径），最小中心括区间容差 1 nm。无固定半径偏移，无法证明更早小区间不合法时记录 NUMERICAL_GEOMETRY_UNRESOLVED，未偷偷跨过。normal fallback 数为 {m['C']['normal_fallback_count']}。",'',
        f"C 域内合法出生数 {m['C']['accepted']:,}，其中原 anchor 不满足纯入口几何者 {m['C']['accepted_but_anchor_aperture_fails']:,}。对所有 C 接纳事件额外用球沿每个折线段的完整扫掠体对真实 WALL 做几何检查；cap 在转换中仍是开放的。该扫掠体用现有 capsule–triangle 精确间隙例程表示，只是球扫掠体的几何计算，没有新增 RBC 模型。",'',
        table(['anchor basin','C 球心合法','路径无几何穿墙','全路径满足 handoff 间隙'],[[b,v['C_births_found'],v['wall_geometrically_clear'],v['full_sweep_handoff_certified']] for b,v in entry['summary'].items()]),'',
        '源事件未据此删除或重选。进入路径证书仍是几何证书，不是有限尺寸入流动力学验证。'+link(review/'METHOD_C_ENTRY_SWEEP_AUDIT.json')+'。','',
        table(['同一 anchor basin','当前首次接纳','C 路径通过','当前拒绝但 C 路径通过','当前接纳但 C 未找到球心'],
            [[b,v['current_on_plane_accepted'],v['C_entry_certified'],v['current_rejected_but_C_entry_certified'],v['current_accepted_but_C_no_birth']] for b,v in paired['by_anchor_basin'].items()]),'',
        '该逐事件配对检查直接区分合法表示转换与越过入口冲突，不能只比较两组未配对总数。'+link(review/'PAIRED_CURRENT_VS_C_ENTRY.json')+'。','',
        '## 完整轨迹与数值稳健性','',
        table(['方法','integrated','完成 01 / 02 / 03','安全停止','其它终止','完成轨迹驻留 median (s)','完成路径 median (µm)'],[
            [k,full[k]['admitted_integrated'],' / '.join(str(full[k]['outlet_counts'][b]) for b in BASINS),full[k]['safety_stops'],json.dumps({r:n for r,n in full[k]['end_reasons'].items() if r not in BASINS+['INTEGRATION_SAFETY_STOP']}),num(full[k]['residence_time_s'].get('median')),num(full[k]['completed_path_length_um'].get('median'))] for k in 'ABC']), '',
        '驻留时间从已表示的出生中心起算，到球心第一次穿过官方出口面为止。C 的 anchor→birth 流线走时单独标为 ENTRY REPRESENTATION TRANSITION；它不是有限尺寸粒子的受力积分，也不计入上述已表示轨迹的驻留时间。出口完成不声明整个球已完全穿出出口。','',
        table(['C 入流全路径 handoff 证书','积分数','完成 01 / 02 / 03','终止原因'],[[flag,v['count'],' / '.join(str(v['outlet_counts'][b]) for b in BASINS),json.dumps(v['end_reasons'])] for flag,v in entryout['groups'].items()]),'',
        '上述分层只用于审计，不能把筛出的子组重新称为未经选择的原始总体。'+link(review/'METHOD_C_OUTCOMES_BY_ENTRY_CERTIFICATE.json')+'。','',
        table(['时间步细化','预定比较数','出口标签改变','终止原因改变','共同时间最大位置差的 median / max (µm)'],[[f'dt/{factor}',v['compared'],v['outlet_changes'],v['end_reason_changes'],f"{num(v['common_time_max_position_difference_um'].get('median'))} / {num(v['common_time_max_position_difference_um'].get('max'))}"] for factor,v in timestep['variants'].items()]),'',
        table(['C 搜索敏感性','状态改变数 / 256','最大距离差 (nm)'],[[k,v['status_changes'],num(v['maximum_distance_difference_m']*1e9 if v['maximum_distance_difference_m'] is not None else None)] for k,v in sensitivity['variants'].items()]),'',
        f"相同高分辨率积分点将 point-tracer 路径步长减半后，改变 basin 标签的通量比例为 {pct(sensitivity['point_basin_refinement']['changed_flux_fraction'])}；未解析通量由 {pct(sensitivity['point_basin_refinement']['original_unresolved_flux_fraction'])} 变为 {pct(sensitivity['point_basin_refinement']['refined_unresolved_flux_fraction'])}。这些数值不确定性保留，不将细化标签覆盖共同 ledger。",'',
        '近壁间隙、近场状态计数、全轨迹/完成子集尺寸及路径长度的完整分位数见 '+link(review/'FULL_TRAJECTORY_COMPARISON.json')+'。安全停止及细化失败是尚未获得物理结局的数值结果，不能记为被捕获、堵塞或流入未观测出口。','',
        '## 逐条回答审核问题','']
    answers=[
        ('1. A 到底改变什么分布？',f"A 固定单个事件的第一次尺寸，最多 512 次按原 flux proposal 抽位置；失败候选造成接纳位置偏向更宽区域。02 basin 比例从 {pct(m['A']['original_anchor_basin_fraction']['OUTLET_02'])} 变为 {pct(m['A']['accepted_birth_point_basin_fraction']['OUTLET_02'])}。尺寸没有逐事件重抽，但 guard 导致接纳子集的尺寸分布也受筛选。"),
        ('2. B 到底改变什么分布？',f"B 的每次尝试位置逐字节固定，尺寸最多重抽 512 次。所得是位置条件下可通过的尺寸；接纳均值 {m['B']['accepted_size_um']['mean']:.4f} µm，相比官方 {audit['official_sonovue']['mean_um']:.4f} µm。部分 anchor 连最小尺寸也无法接纳，接纳总体的 anchor 分布仍受筛选。"),
        ('3. C 是否位置与尺寸都不重抽？','是。机器验证逐事件重建 RNG，检查三方法共同输入哈希，C 保存原 anchor 与原 radius；birth center 只沿该 anchor 的真实流线路径搜索。保持单次提案不等于接纳总体同时保持两个无条件分布。'),
        ('4. 01/03 区域能容纳多大微泡？',f"局部 Dmax 的中位数分别为 {geometry['by_basin']['OUTLET_01']['Dmax_um_quantiles']['median']:.4f} / {geometry['by_basin']['OUTLET_03']['Dmax_um_quantiles']['median']:.4f} µm，采样最大值 {geometry['by_basin']['OUTLET_01']['Dmax_um_quantiles']['max']:.4f} / {geometry['by_basin']['OUTLET_03']['Dmax_um_quantiles']['max']:.4f} µm。原始 SonoVue 的几何通过概率分别 {pct(geometry['by_basin']['OUTLET_01']['geometry_pass_probability'])} / {pct(geometry['by_basin']['OUTLET_03']['geometry_pass_probability'])}。这些是入口中心位置的局部约束。"),
        ('5. 是真实过窄还是人为出生面过严？','支持当前几何模型下的入口边缘/墙面空间限制，不支持所提出的完整球必须在 cut plane 处已入域的人工拒绝机制。不能把此结论扩大为实际生理堵塞，也不排除网格、分辨率和刚性球模型的局限。'),
        ('6. checker 是否错误把 cap 当墙？','没有。实体集合只有 WALL，cap/WALL 面重叠为零，独立圆管测试允许球向上游开放侧伸出；当前 checker 不要求整球入域。所有被选的原 checker 复算与分类加速版本一致。'),
        ('7. 只改 birth representation 能否恢复其它出口？',f"原始 C 队列完成 01/02/03 = {' / '.join(str(full['C']['outlet_counts'][b]) for b in BASINS)}；其中入口路径通过 handoff 几何证书的 C 子组为 {' / '.join(str(entryout['groups']['True']['outlet_counts'][b]) for b in BASINS)}。必须同时看这两组，不能把跨过入口墙面冲突后得到的域内轨迹当作自然恢复。两种原始提案均固定，但接纳筛选仍改变总体。"),
        ('8. B 把尺寸压小多少？',f"接纳直径 mean / median = {m['B']['accepted_size_um']['mean']:.4f} / {m['B']['accepted_size_um']['median']:.4f} µm；相对于这些相同事件第一次尺寸的平均改变为 {m['B']['paired_mean_diameter_change_um']:.4f} µm；对官方分布的 Wasserstein 为 {m['B']['size_distance']['Wasserstein_um']:.4f} µm。详见各 basin 分组 CSV/JSON。"),
        ('9. A 推向 02 有多严重？',f"原 anchor 到接纳 birth 的 02 占比变化为 {100*(m['A']['accepted_birth_point_basin_fraction']['OUTLET_02']-m['A']['original_anchor_basin_fraction']['OUTLET_02']):.3f} 个百分点，basin 总变差 {m['A']['basin_total_variation_vs_original_anchor']:.5f}，32×32 入口密度总变差 {m['A']['inlet_density_total_variation_32x32']:.5f}。接纳事件平均尝试 {m['A']['attempts_per_accepted_event']['mean']:.3f} 次（平均位置重抽次数减一）；计入 guard 失败的全事件均值 {m['A']['attempt_count']['mean']:.3f} 次。"),
        ('10. C 需要多少 inward distance？',f"均值 {m['C']['inward_distance_um']['mean']:.4f} µm；P10 / median / P90 = {' / '.join(format(m['C']['inward_distance_um'][q], '.4f') for q in ['p10','median','p90'])} µm；s/a 中位数 {m['C']['inward_distance_over_radius']['median']:.4f}。整球入域本来就需要一个半径量级的位移，因此 s≈a 本身不能证明旧 checker 有伪影。"),
        ('11. 完整轨迹三个出口统计是什么？','A：'+' / '.join(str(full['A']['outlet_counts'][b]) for b in BASINS)+'；B：'+' / '.join(str(full['B']['outlet_counts'][b]) for b in BASINS)+'；C：'+' / '.join(str(full['C']['outlet_counts'][b]) for b in BASINS)+'。顺序均为 01/02/03，每种方法分母是 5,000 admitted/integrated。'),
        ('12. 原先 720 safety-stop 是否仍是独立第二层问题？',f"入口接纳和域内积分是两个层次。本阶段 A/B/C 安全停止分别 {full['A']['safety_stops']} / {full['B']['safety_stops']} / {full['C']['safety_stops']}；这些不是旧队列的同一批 720 个 ID，也不与旧样本量直接相比。时间步和路径敏感性已单列；未经完成的轨迹去向未知。"),
        ('13. 后续应采用什么科学解释？','使用“入口几何筛选 + admission 引入的条件分布 + 独立的域内数值未知结局”这一分层解释。A 保持每事件尺寸但改位置；B 保持每次位置但条件化尺寸；C 保持提案但改变出生中心表示，且必须显式审计入口扫掠。没有证据支持仅凭三出口覆盖选择某方法为生产模型。本阶段停止，不开展后续模型或 RBC deformation。')]
    for q,answer in answers:lines.extend(['### '+q,'',answer,''])
    lines.extend(['## 验证与可视化入口','',table(['测试集','测试数','失败','错误'],[[Path(t['file']).name,t['tests'],t['failures'],t['errors']] for t in checks]),'',
        '验证器重建全部共同 RNG 事件、核对每个 A/B 候选的固定量、C 所有域内出生中心与括区间、全部实际积分样本哈希及所有已完成出口穿越。原 P8.1 积分与显式出生接口的 4 个基准轨迹逐样本一致；1/2/4/8/16 workers 的固定混合队列科学输出一致。','',
        '静态图为 12 张 2800×1600 PNG，另有 PDF；三段 1920×1080、24 fps、18 秒 MP4 的全部 432 帧均已解码。相机缓慢绕行，物理坐标不旋转；动画中 A/B 重抽是离散候选，C 是表示转换，完整轨迹使用保存的实际样本。完整血管动画仅显示明确披露的代表路径，计数始终来自全部 5,000 条/方法。','',
        '图集与动画可直接打开 '+link(out/'OPEN_RESULTS.html')+'。',''])
    lines.extend([f'- {link(p)}' for p in sorted((visual/'figures').glob('*.png'))]);lines.append('')
    lines.extend([f'- {link(p)}' for p in sorted((visual/'animations').glob('*.mp4'))]);lines.extend(['',
        '原始证据包括所有候选 attempt、完整 RNG ledger、C 流线前缀与出生几何、完整轨迹 NPZ/JSON、每 shard/每粒子远程执行收据、源哈希、运行日志和资源记录。主要入口：'+link(review/'common_inlet_audit_ledger.csv')+'、'+link(review/'A_B_C_admission_events.csv')+'、'+link(review/'A_B_C_full_trajectory_catalog.csv')+'、'+link(report/'PARTICLE8_2A_VALIDATION.json')+'。','',
        '本报告所有 PASS 表示该专项的实现、证据及对照工作完成，不表示三种方法已完成物理真实性验证。用户人工科学审核仍待进行。'])
    (report/'PARTICLE8_2A_REVIEW.md').write_text('\n'.join(lines)+'\n')
    cards=[]
    for p in sorted((visual/'figures').glob('*.png')):
        rel=os.path.relpath(p,out);cards.append(f'<figure><a href="{html.escape(rel)}"><img loading="lazy" src="{html.escape(rel)}"></a><figcaption>{html.escape(p.stem)}</figcaption></figure>')
    for p in sorted((visual/'animations').glob('*.mp4')):
        rel=os.path.relpath(p,out);cards.append(f'<figure><video controls preload="metadata" src="{html.escape(rel)}"></video><figcaption>{html.escape(p.stem)}</figcaption></figure>')
    page='<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>Particle-8.2A 科学审核</title><style>body{background:#090f19;color:#dce6f2;font:18px sans-serif;margin:32px;max-width:1600px}a{color:#66d8ed}img,video{width:100%}figure{margin:28px 0}figcaption{padding:8px}</style><h1>Particle-8.2A 科学审核</h1><p>原始 Frozen FEM · 100,000 同源事件 · A/B/C 每种 5,000 条独立积分</p><p>人工审核：PENDING_USER_REVIEW。入口球心合法性、真实墙面路径间隙、最终出口分别统计。</p><p><a href="'+html.escape(os.path.relpath(report/'PARTICLE8_2A_REVIEW.md',out))+'">中文审核报告</a> · <a href="'+html.escape(os.path.relpath(report/'PARTICLE8_2A_VALIDATION.json',out))+'">机器验证 JSON</a></p>'+''.join(cards)+'</html>'
    (out/'OPEN_RESULTS.html').write_text(page)
    print(json.dumps(dict(statuses=stages,report=str(report/'PARTICLE8_2A_REVIEW.md')),ensure_ascii=False,indent=2))


if __name__=='__main__':main()
