"""Assemble the Chinese scientific review from computed evidence, never invented values."""
from pathlib import Path
import json,hashlib,html,subprocess
import numpy as np
ROOT=Path(__file__).resolve().parents[1];D=ROOT/'data/remote_results'
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def fmt(v):return '无定义/无样本' if v is None else f'{v:.5g}'
def main():
    s=read(D/'statistics.json');e=read(D/'eulerian_summary.json');host=read(D/'compute_provenance.json');done=read(D/'compute_complete.json')
    post=read(D/'post_verification.json');inv=read(ROOT/'data/trajectory_inventory.json');ind=read(ROOT/'data/independent_verification.json')
    whole=s['unweighted']['whole'];near=s['unweighted']['near_wall'];far=s['unweighted']['far_field'];branch=s['unweighted']['branch'];valid=s['validity']
    tests=read(ROOT/'data/test_results.json');visual=read(ROOT/'data/visualization_validation.json')
    extra=read(D/'supplementary_statistics.json');delivery=read(ROOT/'data/delivery_verification.json')
    inspection=read(ROOT/'data/FIGURE_INSPECTION.json')
    status='PASS' if done['status']==post['status']==ind['status']==delivery['status']==inspection['status']=='PASS' and tests['passed'] and visual.get('video',{}).get('frames')==432 else 'FAIL'
    valid_fraction=valid['certified_classic_valid_fraction']
    result=dict(stage='Shear-Lift Magnitude Audit for Finite-Size Microbubble Transport',status=status,git_commit=inv['git_commit'],
        flow_sha={k:v['flow_sha'] for k,v in e['flows'].items()},trajectory_dataset_sha=inv['dataset_sha'],compute_host=host,
        sample_count=dict(trajectory_states=done['sample_count'],accepted_intervals=post['accepted_interval_rows'],initial_diagnostics=post['initial_diagnostic_rows'],
                          eulerian_by_flow={k:v['count'] for k,v in e['flows'].items()}),
        viscosity=dict(value=.00345312,unit='Pa s'),density=dict(value=1056,unit='kg/m^3'),radius_quantiles=e['radius_quantiles'],
        shear_rate_statistics=whole['shear_s_inv'],slip_statistics=whole['slip_speed_m_s'],endpoint_misaligned_slip_statistics=whole['endpoint_slip_speed_m_s'],
        Re_p_statistics=whole['Re_p'],Re_G_statistics=whole['Re_G'],candidate_lift_statistics=whole['CANDIDATE_SAFFMAN_MAGNITUDE'],
        drag_statistics=whole['drag_N'],lift_drag_ratio_statistics=whole['lift_drag_ratio'],lubrication_statistics=near['lubrication_N'],
        lift_lubrication_ratio_statistics=near['lift_lubrication_ratio'],saffman_valid_fraction=valid_fraction,
        near_wall_invalid_fraction=valid['near_wall_invalid_fraction'],validity_detail=valid,
        formula_source=['10.1017/S0022112065000824','10.1017/S0022112068999990','10.1017/S0022112001007145','10.1103/PhysRevFluids.6.104309'],
        formula_convention=dict(coefficient_radius=6.46,coefficient_diameter=1.615,size='RADIUS',slip='FLUID_MINUS_PARTICLE_AT_FORCE_EVALUATION_POSITION',
            shear='sqrt(2 E:E); PROJECT_3D_MAGNITUDE_PROXY_NOT_A_GENERAL_SAFFMAN_THEOREM',Re_p='a*slip/nu',Re_G='a*a*shear/nu',
            direction='slip CROSS curl; NORMALIZED_ILLUSTRATION_ONLY',wall_length='(a+h)/sqrt(nu/shear)',
            zero_guard='velocity_tol=512*6*eps*(1+zeta/self)*max(norm(u),norm(v)); denominator tol = self*velocity_tol for drag, max(self,zeta)*velocity_tol for lubrication; numerator tol = Saffman_prefactor*velocity_tol; NUMERICAL_NOT_PHYSICAL_THRESHOLD'),
        rigid_particle_formula_role='ORDER_OF_MAGNITUDE_ONLY',sonovue_applicability_status='NOT_A_VALIDATED_LIPID_SHELL_LIFT_MODEL',
        conclusion='MAGNITUDE_CANNOT_BE_TRUSTED_BECAUSE_MODEL_VALIDITY_FAILS',
        current_lift_model='NOT_INCLUDED',formal_solver_modified=False,trajectory_recomputed=False,
        all_region_statistics=s['unweighted'],time_weighted_statistics=s['time_weighted'],counts=s['counts'],zero_guard=s['zero_guard'],
        supplementary_statistics=extra,delivery_verification=delivery,figure_inspection=inspection,
        extrema=post['extrema'],eulerian_audit=e['flows'],nearwall_geometry_parity=ind['geometry_parity_samples'],
        protected_file_count=ind['protected_file_count'],tests=tests,video=visual['video'],
        field2_slip_dependent_lift='NOT EVALUABLE WITHOUT PARTICLE SLIP',
        next_stage_started=False)
    (ROOT/'SHEAR_LIFT_MAGNITUDE_VALIDATION.json').write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
    def summary(metric,group=whole):
        r=group[metric];return f"中位数 {fmt(r['median'])}，P99 {fmt(r['P99'])}，最大 {fmt(r['max'])}"
    md=['本阶段没有把剪切升力加入轨迹方程，只评估其可能量级和适用性。','',
        '# Shear-Lift Magnitude Audit：中文审核报告','',
        f"本轮审计状态 **{status}**。这表示只读计算、数据检查和可视化通过，不表示 Saffman 已被验证为本项目微泡的真实受力模型。最终科学结论为 **C：MAGNITUDE_CANNOT_BE_TRUSTED_BECAUSE_MODEL_VALIDITY_FAILS**。",'',
        '## 先用普通语言解释五个量','',
        '**Shear（剪切）**：相邻位置的血流速度变化有多快。本轮从微泡位置的 FEM 速度梯度计算，不拿壁面剪切应力 WSS 替代。WSS 是壁面上的切向应力；微泡局部剪切率是流体内部的位置量，两者的单位与取样位置不同。','',
        '**Slip（滑移）**：同一时间、同一位置下，血流速度与微泡速度之差。微泡随背景流走时，滑移接近零；受到近壁阻力或不穿透约束时，两者才可能明显不同。','',
        '**Lift（升力）**：这里指可能让颗粒偏离原流动方向的横向水动力。经典 Saffman 理论涉及小但非零的惯性；本轮只算候选量级，没有把它加入运动。','',
        '**Drag（阻力）**：现有模型中，血流与球体的速度差通过 6πμa 转换成的 Stokes 阻力贡献。','',
        '**Lubrication（润滑阻力）**：微泡接近固壁时，窄间隙内流体对法向运动的额外阻碍。本项目只含 leading 法向阻力与既有正则化。本报告比较的是等价水动力贡献；不是额外启动惯性力积分，也不把不穿透约束反力冒充润滑力。','',
        '## 1. 当前轨迹里已有 shear lift 吗？','',
        '**没有。** 逐函数核对了 `SavedTrajectoryStepper → Particle65Stepper → v1_trial → assemble_v1 → solve_resistance`。右端只有 Stokes self resistance 乘背景自由速度；矩阵增加球壁法向近场阻力，随后处理不穿透约束。没有 Saffman、壁面惯性升力、变形升力、Magnus、Brownian、重力/浮力或声辐射力项。旋转跟随局部涡量不等于已有 Magnus 力。详见 [源码核对](SOURCE_MODEL_AUDIT.md)。','',
        '## 2. 真实轨迹中的剪切率有多大？','',
        f"旧流场全体保存状态的局部 γ_E：{summary('shear_s_inv')} s⁻¹。远场模型状态的中位数为 {fmt(far['shear_s_inv']['median'])} s⁻¹，近壁为 {fmt(near['shear_s_inv']['median'])} s⁻¹，分叉描述性 ROI 为 {fmt(branch['shear_s_inv']['median'])} s⁻¹。定义 γ_E=√(2E:E)，E=(∇u+∇uᵀ)/2。三维场中的这一标量不是对任意局部流型都严格成立的 Saffman 剪切参数。",'',
        '## 3. 真实滑移有多大？','',
        f"与正式力求值位置对齐后：{summary('slip_speed_m_s')} m/s。远场中位数 {fmt(far['slip_speed_m_s']['median'])} m/s，近壁中位数 {fmt(near['slip_speed_m_s']['median'])} m/s。",'',
        f"如果直接把保存的区间速度与该行终点流速相减，会得到中位数 {fmt(whole['endpoint_slip_speed_m_s']['median'])} m/s。这包含显式步进的空间/时间错配，不能直接解释成额外物理滑移。每行同时保存这个原始诊断；主要力比较采用该速度实际求值的前一接受位置，并核对时间和位移恒等式。{post['initial_diagnostic_rows']:,} 个初始行是自由流诊断，时间权重为零。",'',
        '## 4. 经典 Saffman 候选量级有多大？','',
        f"半径形式 F_c=6.46 μa²|u−v|√(γ_E/ν)，全体状态 {summary('CANDIDATE_SAFFMAN_MAGNITUDE')} N。近壁 {summary('CANDIDATE_SAFFMAN_MAGNITUDE',near)} N。所有字段都叫 `CANDIDATE_SAFFMAN_MAGNITUDE`，不是正式微泡升力。",'',
        '系数已经用原文半径定义、权威公式和独立归一化转换交叉核对；1965 摘要中的 81.2 不能直接照抄。直径形式系数是 1.615。文献访问限制也如实列在 [文献审计](SHEAR_LIFT_LITERATURE_AUDIT.md)，没有把无法读到的全文写成已核验。','',
        '## 5. 相对 Stokes drag 是多少？','',
        f"在分母有数值分辨率的 {whole['lift_drag_ratio']['count']:,} 个状态中，F_c/F_drag：{summary('lift_drag_ratio')}。Stokes drag 本身 {summary('drag_N')} N。",'',
        f"全体状态中，有 {s['zero_guard']['whole']['1']:,} 个被标为 `BOTH_FORCE_MAGNITUDES_NEAR_ZERO`，其比值不定义为 0，也不加分母下限制造大比值。绝对力仍保留。远场很多差值只在浮点舍入级别。比值在非零滑移极限可化为 (6.46/6π)√Re_G，但这不授权把真实 0/0 填成该值。代表轨迹图用灰色虚线单独表示这个极限。",'',
        '## 6. 相对 wall lubrication 是多少？','',
        f"近壁可定义的 {near['lift_lubrication_ratio']['count']:,} 个状态中，F_c/F_lub：{summary('lift_lubrication_ratio',near)}。近壁润滑贡献本身 {summary('lubrication_N',near)} N。",'',
        '采用正式模型已有的 ζ_n=w(h/a)6πμa²/max(h,h_lower)，然后 F_lub=|ζ_n(v·n)|。间隙、法向和速度都在同一个求值位置。进入 handoff 后，法向速度可能近零，此时有限阻力系数并不保证非零润滑贡献；不穿透约束可以承担剩余平衡。数值零值预算分别通过候选力、Stokes 和润滑算子传播，避免大 ζ_n 把法向速度舍入误差放大后误认为真实润滑力。无定义比值保留为 NaN；原绝对力没有截断或修改。','',
        f"近壁共有 {s['counts']['near_wall']:,} 状态，其中 {extra['ratio_bins']['near_wall']['lift_lubrication_ratio']['undefined']:,} 个润滑分母无法可靠分辨。可定义比值中仅 {extra['ratio_bins']['near_wall']['lift_lubrication_ratio']['gt_1']} 个大于 1。最大比值对应候选力 {post['extrema']['lift_lubrication_ratio']['candidate_N']:.5g} N、润滑贡献 {post['extrema']['lift_lubrication_ratio']['lubrication_N']:.5g} N、Stokes 阻力 {post['extrema']['lift_lubrication_ratio']['drag_N']:.5g} N；较大的润滑比不等于升力主导全部平衡。该点位于近壁 handoff 分叉区，同样超出经典公式的无限域适用范围。",'',
        '数值预算倍数 0.1、1、10 的单独敏感性复核已在服务器完成，结果见下表。这里的条件数只对应 self-scaled 单球/单法向壁面块，不构成对多接触约束和几何误差的严格总误差界。','',
        '| 数值零值预算倍数 | 可定义润滑比数 | 中位数 | P99 | 最大 |','|---|---:|---:|---:|---:|']
    for factor,q in extra['numerical_guard_multiplier_sensitivity'].items():
        md.append(f"| {factor} | {q['count']:,} | {fmt(q['median'])} | {fmt(q['P99'])} | {fmt(q['max'])} |")
    md += ['',
        '## 7. 哪些位置最大？','']
    for metric,name in [('CANDIDATE_SAFFMAN_MAGNITUDE','候选升力最大'),('lift_drag_ratio','候选/drag 最大'),('lift_lubrication_ratio','候选/润滑最大')]:
        q=post['extrema'][metric];xyz=np.array(q['x_m'])*1e6
        md.append(f"- {name}：{q['dataset']} #{q['id']}，t={q['time_s']:.6g} s，位置 ({xyz[0]:.4f}, {xyz[1]:.4f}, {xyz[2]:.4f}) µm；h/a={q['h_over_a']:.5g}，壁距/剪切长度={q['wall_margin']:.5g}；候选/drag/润滑绝对值依次 {q['candidate_N']:.5g}/{q['drag_N']:.5g}/{q['lubrication_N']:.5g} N。")
    md+=['','这些是本次已有保存数据中的最大值，不是整个连续空间的严格极值。完整的入口、主干中心、分叉前、分叉、各端口邻域以及数据集分组表见 [statistics.csv](data/remote_results/statistics.csv)。ROI 是明确记录的描述性空间窗口，不是已验证的解剖分区，组间允许重叠。','',
        '| 区域 | 状态数 | 候选力 P99 / 最大 (N) | Stokes P99 / 最大 (N) | 润滑 P99 / 最大 (N) |',
        '|---|---:|---:|---:|---:|']
    for region in ['inlet_vicinity','trunk_center','prebranch','branch','outlet_01_vicinity','outlet_02_vicinity','outlet_03_vicinity','far_field','near_wall','handoff']:
        q=s['unweighted'][region]
        md.append(f"| {region} | {s['counts'][region]:,} | "+' | '.join(f"{fmt(q[k]['P99'])} / {fmt(q[k]['max'])}" for k in ['CANDIDATE_SAFFMAN_MAGNITUDE','drag_N','lubrication_N'])+' |')
    md += ['',
        '## 8. Saffman 有效条件在哪些区域成立？','',
        f"采用半径 Re_p=a|u−v|/ν、Re_G=a²γ_E/ν。Re_p {summary('Re_p')}；Re_G {summary('Re_G')}。然而小 Reynolds 数只满足部分要求：经典无限域推导还需要壁面远于剪切扰动长度 L_G=√(ν/γ_E)。本次 (a+h)/L_G 的最大值只有 {fmt(whole['wall_margin']['max'])}，全部小于 1；因此连“壁面位于扰动外部”的宽松必要次序也不满足。可确认的经典无限域有效覆盖率为 **{100*valid_fraction:.2f}%**。",'',
        '逐状态保留 Re_p/√Re_G、√Re_G、(a+h)/L_G 三个连续裕度。0.03、0.1、0.3 的渐近分离参数仅作描述性灵敏度筛查，不是文献给出的通用有效/无效阈值；三种筛查都无法消除这里的壁面问题。三维非简单剪切和壳层边界条件还带来额外限制。背景是冻结的连续介质 FEM 场；这些独立 MB 轨迹没有解析 RBC–MB 相互作用，因此本审计也不能排除悬浮细胞造成的其他侧向输运。','',
        '## 9. near-wall 已超出经典公式适用范围吗？','',
        f"是，本次近壁状态中 {100*valid['near_wall_invalid_fraction']:.2f}% 的壁面距离小于剪切扰动长度。并且代码里的 FAR_FIELD 只是表示 h/a≥0.05、既有近场修正关闭，并不代表相对于 Saffman 扰动长度处于无限域。不能因为“没有启用 lubrication”就宣布 Saffman 有效。文献中的单一平壁模型也不能直接用于曲壁、多分叉和纳米间隙。",'',
        '## 10. 能直接作为 SonoVue 正式模型吗？','',
        '**不能。** 经典对象是无滑移刚性球；SonoVue 是磷脂壳层包覆气泡。洁净剪切自由气泡、污染气泡及具有膜流变的壳层气泡具有不同界面条件。洁净气泡低 Re 的 4/9 对比仅说明界面会改变结果，不是 SonoVue 可直接套用的校正系数。此次角色保持 `ORDER_OF_MAGNITUDE_ONLY`；没有建立或验证壳层侧向迁移模型。','',
        '## 11. 下一步是否值得开发？','',
        '值得把“受限血管中的 bubble-specific lateral migration”列为后续研究候选，首先明确壳层界面条件、曲壁/多壁作用和近壁验证基准。当前证据不支持直接把无限域 Saffman 项加入正式求解器，也不能仅因候选/drag 比值较小就宣称所有横向迁移机制都可忽略。本阶段到此停止，未启动新模型开发。','',
        '## 数据覆盖、计算和两版流场','',
        f"原始清点 {len(inv['records']):,} 个数组记录，按“数组 SHA256 + 半径”识别 {len(inv['records'])-inv['unique_files']} 个重复内容，合并统计 {inv['unique_files']:,} 个唯一记录；其中有保存状态的轨迹 {post['initial_diagnostic_rows']:,} 条，共 {done['sample_count']:,} 状态。空的入口未接纳记录保留在清单，不计为空间样本。未完成路径全部纳入；没有只筛成功出口。",'',
        '| 数据集 | 原记录数 | 已有终止情况 |', '|---|---:|---|']
    for ds,c in inv['counts'].items():md.append(f"| {ds} | {sum(c.values())} | "+'；'.join(f'{k}: {v}' for k,v in c.items())+' |')
    md+=['','“1249 条”是 P81 到达 Outlet 02 的成功条数，不是其全部记录。PPT 的 500 条由 B 入场策略筛选尺寸，不能代表未经入口选择的 SonoVue 人群。P82A 仅使用现存 formal/A 数组；历史 guard/timestep/cache/worker benchmarks、点示踪线和理想管 RBC 动画不混作新的正式 MB 样本。Outlet 01 没有完成 MB 轨迹；Outlet 03 只有 PPT 中已有的 3 条，未补造路线。','',
        f"正式后处理主机 `{host['hostname']}`，{host['workers']} 个 CPU 工作进程，未使用 GPU 计算；运行耗时 {done['elapsed_s']:.1f} s。主机、CPU、内存、GPU 查询原文在 compute_provenance.json。首次启动因清单外依赖漏打包被哈希校验阻止，补齐原文件后重新运行；未绕过冻结检查。",'',
        f"受保护的 {ind['protected_file_count']:,} 个文件前后 SHA256 一致。所有已审计轨迹在远程也逐文件复核，改动数为 0；重建间隙与保存值最大差 {post['maximum_gap_reconstruction_error_m']:.5g} m，另有 {ind['geometry_parity_samples']} 个小样本与原始精确几何和润滑政策独立比对。",'',
        '| 流场 | Eulerian 单元数 | 剪切率中位数 (s⁻¹) | P99 (s⁻¹) | 最大值 (s⁻¹) | 滑移相关力 |','|---|---:|---:|---:|---:|---|']
    for k,v in e['flows'].items():
        q=v['shear_statistics'];md.append(f"| {k} | {v['count']:,} | {fmt(q['median'])} | {fmt(q['P99'])} | {fmt(q['max'])} | "+('对应旧正式轨迹另行统计' if '0P353' in k else 'NOT EVALUABLE WITHOUT PARTICLE SLIP')+' |')
    md+=['','两版场分别从实际冻结 VTU 读取，没有把旧速度乘 5.67。Eulerian 每个四面体取重心；P1 梯度在该单元内为常数，另给体积加权统计。新场只输出剪切、Re_G 和候选力/滑移前因子，未伪造新场微泡滑移。表面图的颜色来自相邻单元局部剪切率，不是 WSS；三维显示复用固定 Z 轴、放大视野、旁置端口箭头、透明文字和渐变标签。','',
        '全状态统计是对本次保存样本的描述；自适应步长会增加近壁状态数量。因此 JSON 同时提供 dt 加权表，仍不把这些独立轨迹的加权结果当作生理浓度、停留概率或注入后的药代分布。','',
        '## 全数据统计（状态等权；SI）','',
        '| 量 | 有限值数 | 中位数 | P90 | P95 | P99 | 最大 |','|---|---:|---:|---:|---:|---:|---:|']
    for metric in ['shear_s_inv','slip_speed_m_s','Re_p','Re_G','CANDIDATE_SAFFMAN_MAGNITUDE','drag_N','lift_drag_ratio','lubrication_N','lift_lubrication_ratio']:
        group=near if metric in ['lubrication_N','lift_lubrication_ratio'] else whole;q=group[metric]
        md.append('| '+metric+(' (near wall)' if group is near else '')+f" | {q['count']:,} | "+' | '.join(fmt(q[k]) for k in ['median','P90','P95','P99','max'])+' |')
    md+=['','初始行尚未完成阻力平衡，不宜把其重建润滑贡献作为已实现动力学受力。单独排除这些行后的完整统计见 [accepted_interval_statistics.csv](data/remote_results/accepted_interval_statistics.csv)。接受区间候选/润滑比中位数为 '+fmt(extra['initial_and_accepted_statistics']['accepted_intervals']['lift_lubrication_ratio']['median'])+'，最大值不变。','',
        '### 描述性比值分箱','',
        '**DESCRIPTIVE ONLY**：下列 0.01、0.1、1 仅用于展示分布，不是经验证的是否加入模型的阈值；无定义值单独计数。','',
        '| 比值口径 | <0.01 | [0.01,0.1) | [0.1,1] | >1 | 无定义 |','|---|---:|---:|---:|---:|---:|']
    for group,key in [('whole','lift_drag_ratio'),('near_wall','lift_lubrication_ratio')]:
        q=extra['ratio_bins'][group][key]
        md.append('| '+group+' / '+key+' | '+' | '.join(f"{q[k]:,}" for k in ['lt_0p01','from_0p01_to_0p1','from_0p1_to_1','gt_1','undefined'])+' |')
    md+=['','## 正式 SonoVue 粒径敏感性','', '| 分位 | 直径 (µm) | 半径 (µm) |','|---|---:|---:|']
    for k,q in e['radius_quantiles'].items():md.append(f"| {k} | {q['diameter_um']:.8f} | {q['radius_m']*1e6:.8f} |")
    md+=['','这些分位从冻结经验直径 CDF 直接求逆，然后除以 2 得半径。固定滑移/剪切时，候选量级随 a²、|slip| 和 √γ 增长；经典公式本身没有壁距变量，因此壁距改变的是适用性和既有润滑贡献。敏感性 CSV 明示 `HYPOTHETICAL_ONLY`；它没有生成新的正式轨迹，也没有把该尺度图当成新流场实测力。','',
        '## 输出与可视化审核','',
        '[全部图片与旋转动画](OPEN_RESULTS.html)；[逐状态、CSV 和计算记录](data/remote_results/)；[Validation JSON](SHEAR_LIFT_MAGNITUDE_VALIDATION.json)。', '',
        f"永久测试：{tests['summary']}。图像人工检查记录在 data/FIGURE_INSPECTION.json。旋转视频 {visual['video']['frames']} 帧、{visual['video']['fps']} fps、{visual['video']['duration_s']} s；这是相机播放时间，不是新增物理模拟时长。少量微泡点取自既有轨迹；画面标注 Candidate Saffman magnitude only / Not included in trajectory dynamics。所有图中候选方向为归一化示意，不宣称可预测 SonoVue 的真实偏转。",'',
        '**停止状态：** 正式求解器未修改，正式轨迹未重新计算，下一阶段未启动。','']
    (ROOT/'SHEAR_LIFT_MAGNITUDE_REVIEW.md').write_text('\n'.join(md))
    gallery=['<!doctype html><meta charset="utf-8"><title>Shear-lift audit</title><style>body{font:17px system-ui;background:#eef2f6;color:#21344a;max-width:1400px;margin:35px auto}img,video{width:100%;background:white}section{margin:28px 0}a{color:#20648a}h2{font-size:20px}</style><h1>Shear-Lift Magnitude Audit</h1><p>Candidate magnitude only · Not included in trajectory dynamics</p><p><a href="SHEAR_LIFT_MAGNITUDE_REVIEW.md">中文审核报告</a> · <a href="SHEAR_LIFT_MAGNITUDE_VALIDATION.json">Validation JSON</a> · <a href="SHEAR_LIFT_LITERATURE_AUDIT.md">文献审计</a></p><video controls loop preload="metadata" src="animations/shear_lift_audit_rotating.mp4"></video>']
    for f in sorted((ROOT/'figures').glob('*.png')):gallery.append(f'<section><h2>{html.escape(f.stem.replace("_"," "))}</h2><a href="figures/{f.name}"><img loading="lazy" src="figures/{f.name}"></a></section>')
    (ROOT/'OPEN_RESULTS.html').write_text('\n'.join(gallery))
    print(status)

if __name__=='__main__':main()
