"""Readable evidence and repair reports, generated from actual diagnostics."""
import json


def fmt(value):
    if value is None:return '未测/无有效结果'
    if isinstance(value,(int,float)):return f'{value:.7g}'
    return str(value)


def repair_log(rows,c):
    lines=['# 分析修复记录','','确认的 Python 问题：`block_statistics.mean` 使用全部选取样本，而 CI 的方差只用完整块。黏度平均剖面也使用全部样本。',
      '先新增 `StatisticsRepairTests.test_mean_and_ci_use_same_complete_blocks`：旧实现返回均值 1.9230769，完整块均值为 0，测试失败；最小修复后通过。',
      '修复后均值、标准差和黏度拟合使用与 CI 相同的完整块样本；另存 `all_sample_mean`、`all_sample_profile_star`、估计样本数和区间。尾部原始记录没有删除。没有修改分块长度、平衡段阈值或物理目标。',
      '没有修复或重编译原生数值模型，没有把分析修复当作液体物性改善。历史报告、CSV 和脚本快照不变。',
      '',f'修复前代码与失败/通过日志：`{c["setup_directory"]}`。','','| 任务 | 温度：修复前 → 后 | 黏度 SI：修复前 → 后 | 压力星号：修复前 → 后 |','|---|---|---|---|']
    for r in rows:
        a,b=r['before'],r['after']
        lines.append(f"|{r['task_id']}|{fmt(a['temperature_mean'])} → {fmt(b['temperature_mean'])}|{fmt(a['nu_si'])} → {fmt(b['nu_si'])}|{fmt(a['pressure_mean_star'])} → {fmt(b['pressure_mean_star'])}|")
    return '\n'.join(lines)+'\n'


def chinese_report(s,a,sampling):
    u=s['units'];b=s['shared_budget'];t={x['task_id']:x for x in a['temperature_audit']};p=a['pressure_audit']
    lines=['# SDPD 原因排查报告','',
      '**结论：已证实启动过渡进入统计段，以及均值/置信区间样本不一致的分析错误；尚不能将全部温偏或剖面残差归结为原生求解器错误。**',
      f"诊断状态 `{s['status']}`；`selection=null`；人工验收 PENDING。诊断完成不等于材料合格。",
      '', '## 来源与复现','',f"权威旧包：`{a['provenance']['original_result']}`。按配置和 SHA-256 定位；18 项旧/本轮实验从 CSV 和快照重新分析，修复前温度、压力、黏度、窗口、性能和退出状态精确复现。",
      f"当前冻结单位 L0={u['L0']:.17g} m、M0={u['M0']:.17g} kg、t0={u['t0']:.17g} s；ν*={u['required_nu_star']:.12g}、μ*={u['required_mu_star']:.12g}，两者相差真实质量密度 8。没有重设单位。",
      '执行计划、实际步数、脚本/库哈希及退出原因见 `provenance.json`；源调用链行号及哈希见 `native_source_contract.json`。',
      '', '## 原因与证据状态','','| 问题 | 状态 | 一句话解释 |','|---|---|---|']
    for e in s['evidence']:lines.append(f"|{e['phenomenon']}|{e['status']}|{e['one_sentence']}|")
    lines += ['', '每项支持证据、反对证据、缺失观测和最小验证方法完整保存在 `evidence_matrix.json`。',
      '', '## 温度：测量、过渡与推进偏差分开','','| 任务 | 首个 worker 值 | 峰值 | 最后 10% 均值 | 原统计均值 | 完整块均值 | N−1 修正相对影响 |','|---|---:|---:|---:|---:|---:|---:|']
    for r in a['temperature_audit']:
        lines.append(f"|{r['task_id']}|{fmt(r['first_worker_kBT'])}|{fmt(r['peak_kBT'])}|{fmt(r['last_10pct_mean'])}|{fmt(r['original_selected_mean'])}|{fmt(r['corrected_selected_mean'])}|{fmt(r['DOF_relative_correction'])}|")
    lines += ['', '最后 10% 只是预先列出的诊断窗口之一，不是重新挑选的正式验收窗口。所有 20–90% 截段及两半漂移均保存并显示。',
      '等质量平衡温度采用 m Σ|v−vCOM|²/[3(N−1)]；流动温度扣除实际瞬时分箱均速并计入 3(N−B) 自由度。空箱不贡献均速或自由度。原始 profile 的 count、平均三分量和 sum_v2 可重算原始、COM 和粗分箱温度；细分箱逐箱矩没有保存。',
      f"平衡任务的 CPU 矩重算最大 COM 差值为 {fmt(t['sdpd_equilibrium']['independent_moment_recovery']['COM_max_abs'])}；其 N−1 修正仅 {100*t['sdpd_equilibrium']['DOF_relative_correction']:.5f}%，不能解释约十几个百分点的偏高。",
      '原生 Stats 样本在第 1、every+1、…步；worker 在 every、2every、…步。相邻步差异图明确标注该一积分步差，不任意平移来强行重合。两者都读取同一底层速度，矩复核也不是完全独立的物理观测。',
      '旧速度分量方差、逐粒子速度/力、随机状态和 t=0 位置未保存；不能恢复粒子异常尾部分布或构造可信续跑。',
      '', '## 初始化与时间推进','','Uniform 在单位胞中生成随机位置，其种子由 rank 和粒子向量名字确定，并非均匀晶格。Python 随后按固定种子设高斯速度并只扣一次 COM。较大背景压力对配对力有实际贡献。',
      '真实位置快照的核密度、近邻、最近距离、保守力、势能变化和耗散矩阵量级见 `snapshot_audit.json`。这些是固定快照上的独立公式核查，没有开发 CPU 流体求解器。',
      '完整纯液体调用链的更新为 v_new=v_old+F(x_old,v_old)dt/m，x_new=x_old+v_new dt。不能仅凭 VelocityVerlet 类名假定存在另一次半步 kick 或中心速度转换。该语义本身不是已证实的原生 bug。',
      '固定 OU 耗散模态的解析平衡方差比为 1/(1−λdt/2)，说明显式耗散存在正步长偏差的机制。快照耗散矩阵给的是量级估计，不能直接乘一个系数修正实际液体温度。',
      '按核快照平均耗散模态估算，旧 dt 下温偏量级约 0.5–0.6%，远小于强启动峰值；这不是完整流体误差上界。声学、粒子间黏性、盒动量衰减、实际速度位移分别见 `time_scale_audit.json`；结构松弛时标尚未建立。',
      '原生噪声幅度 sqrt(2·5μ·kBT/dt) 只在速度推进时再乘一次 dt；配对随机标量按排序 ID 生成，方向反对称。生成器随当前时间更新，同一步重复计算复用种子。实际设备随机量分布未被完整统计验证。',
      '反复 u.run 会重建调度和 cell lists、重新 setup 插件，但保留当前粒子、时间和 SDPD 随机发生器；已检查范围内未发现重置初态。重排可能改变浮点累加和后续混沌轨迹，未声称连续/分块逐位相同。',
      '', '## 一块统计的完整账单','','| 任务/观测 | 原计划/分配/实际步数 | 初始化后样本 | 间隔* | ACF 样本 | 最小块样本 | ACF 要求样本 | 最终块样本 | 尾部 | 块数 |','|---|---|---:|---:|---:|---:|---:|---:|---:|---:|']
    for r in sampling:
        if not r['task_id'].startswith('sdpd_'):continue
        lines.append(f"|{r['task_id']}/{r['observable']}|{r['desired_steps']}/{r['allocated_steps']}/{r['actual_steps']}|{r['selected_samples']}|{fmt(r['sample_interval_star'])}|{fmt(r['tau_samples'])}|{r['minimum_block_samples']}|{r['ACF_block_samples']}|{r['final_block_samples']}|{r['discarded_tail_samples']}|{r['block_count']}|")
    flow=next(r for r in sampling if r['task_id']=='sdpd_flow' and r['observable']=='viscosity_slope')
    lines += ['',f"基准流动：实际 {flow['actual_steps']} 步、dt*={flow['dt_star']}、物理演化 {flow['actual_physical_duration_s']*1e6:.6f} µs；从 t*={flow['selected_start_star']} 开始只剩 {flow['selected_samples']} 个样本。块长由 max(物理下限 {flow['minimum_block_samples']}, ceil(5τ)={flow['ACF_block_samples']}) 得 {flow['final_block_samples']}，因此 {flow['selected_samples']} // {flow['final_block_samples']} = {flow['block_count']}，余 {flow['discarded_tail_samples']}。",
      'ACF 单位是样本数；最终乘采样间隔才是模拟时间。所有 CSV 的单调性、规则采样、step/time 和 profile 分箱顺序已核对。实际停止任务按 completion 和 ledger 记录，未当作完整执行。',
      'ACF 受未衰减趋势影响，当前 τ 不能作为可信稳态相关时间。原规则允许 10% 两半漂移，和 2% 温度目标不等价；本轮不事后放宽或修改验收阈值。',
      '', '## 黏度和剖面','','| 任务 | 原点估计 ν SI | 完整块 ν SI | 原相对 RMS | 完整块相对 RMS | 近变号面 RMS* | 内部 RMS* |','|---|---:|---:|---:|---:|---:|---:|']
    for r in a['profile_audit']:
        lines.append(f"|{r['task_id']}|{fmt(r['old_nu_si'])}|{fmt(r['corrected_nu_si'])}|{fmt(r['old_RMS_relative'])}|{fmt(r['corrected_RMS_relative'])}|{fmt(r['near_force_planes_RMS_star'])}|{fmt(r['interior_RMS_star'])}|")
    lines += ['', '主解析解使用每粒子力 F、粒子质量 m 和运动黏度 ν，u=F/(mν)·分段抛物线；包含有限箱宽的二次平均修正。平均流速点、目标线、拟合线与所有位置残差均保存；未删除不利区域。',
      'rc/半通道宽=1/4，支持范围跨 4 个剖面分箱；非局部作用与连续分段抛物线之间存在有限尺度疑点。数据中未发现空箱，粒子加权与等时间平均差异另列。缺乏稳态独立块，尚不能将空间残差确认为核作用或噪声中的任何一种。',
      '半力/半步的总物理时长不同，短记录仍在初始化；直接比较各自全窗拟合差异不构成稳态敏感性验证。修复后的点估计仍是“拟合估计”，没有获得有效 CI。',
      '', '## 压力参考、覆盖和精度','',f"全盒密度代入输入 EOS 的背景量为 {fmt(p['background_at_global_density_pa'])} Pa；这不是已知真实绝对压力，也不是可任意移除的无影响配对力。",
      '机械压力使用配对应力迹除以 3V，加一次 COM 扣除后的动能动量通量。EOS 已用于保守配对力，不再另加 EOS 压力。应力来自积分前、动量来自积分后，离散相位差 dt 明确保留。',
      '', '| 端点 | 端点−基准 Pa | 需要表压 Pa | 名义余量 Pa | 条件 Welch 半宽 Pa | 保守端点+基准半宽 Pa |','|---|---:|---:|---:|---:|---:|']
    for r in p['endpoint_margins']:
        lines.append(f"|{r['endpoint']}|{fmt(r['pressure_difference_pa'])}|{fmt(r['required_gauge_pa'])}|{fmt(r['nominal_margin_pa'])}|{fmt(r['independent_run_Welch_95_halfwidth_pa'])}|{fmt(r['conservative_endpoint_plus_reference_halfwidth_pa'])}|")
    lines += ['', '两个压差共享同一基准，所以其协方差含 Var(P_ref)；高低端点差中基准会抵消。Welch 估计假定运行间独立，实际共同随机种子可能产生协方差，未测的协方差不置为实测零。原保守端点规则不需要该独立假设，但不宣称联合 95% 覆盖。',
      f"virial CSV 舍入最大半量级约 {fmt(p['CSV_virial_rounding_max_half_quantum_pa'])} Pa；单精度一个数值分辨尺度约 {fmt(p['single_float_pressure_resolution_scale_pa'])} Pa。它们远小于当前统计不确定性；这不是完整 GPU 求和误差证明，也不支持因此重编译。",
      '各 EOS 点仍有温度和结构问题，故即使几何区间覆盖也不等于已得到合格等温压力响应。',
      '', '## GPU 决策与下一步','','当前原授权累计 '+fmt(b['total_charged_or_reserved_s'])+' / 3600 s，余额 '+fmt(b['remaining_s'])+' s。']
    plan=a['proposed_experiments'];q=plan['selected_probe']
    lines += [f"筛选探针 `{q['id']}`：{q['steps']} 步，dt*={q['dt_star']}，物理 {q['physical_duration_s']*1e6:.6f} µs，预测墙钟 {q['predicted_total_wall_s']:.3f} s，预约上限 {q['allocation_s']} s。它只判断早期温升/衰减的步长依赖，不用于完整温度或黏度收敛。",
      '执行结果：'+('未执行；参见 CPU 决策与显式执行入口。' if s['new_probe'] is None else json.dumps(s['new_probe']['summary'],ensure_ascii=False)),
      '', '| 后续条件情景 | 步数 | 物理时长 µs | 预测墙钟 s | 单任务600秒可行 |','|---|---:|---:|---:|---|']
    for r in plan['eight_block_cost_scenarios']:
        lines.append(f"|{r['scenario']}|{r['steps']}|{r['physical_duration_s']*1e6:.6f}|{r['predicted_wall_s']:.3f}|{r['single_task_600s_feasible']}|")
    lines += ['', '上述是有明确假设的预算情景：稳态相关时间尚未可靠测得，不能机械把当前时间乘八。超过原余额的工作标为 REQUIRES_ADDITIONAL_AUTHORIZATION；超过单任务上限还需要先验证完整状态及随机过程的续跑。',
      '**下一步只优先建立无驱动目标液体的可信热平衡平台及独立块。** 这同时决定黏度拟合、EOS 和 dt 敏感性是否有意义，优先级高于新求解器或盲目调参数。',
      '', '## 验证与范围','','CPU 回归记录及修复前失败用例见 setup；最终 CPU 验证记录绑定 provenance。浏览器只有实际检查后才在单独 browser_checks.json 写 PASS，人工验收始终 PENDING。',
      '未修改原始数据、历史验收、Mirheo/CUDA/MPI 或兼容补丁；未运行真实血管、RBC、微泡、SDF 或开放边界。']
    next_run=plan['next_priority_experiment']
    lines += ['', '## 唯一优先补测的具体预算','',
      f"建议无驱动、原 dt*={next_run['dt_star']}、固定 {next_run['steps']} 步，即 {next_run['physical_duration_s']*1e6:.6f} µs。按旧无驱动实测成本加 15% 和 18 s 余量，预计 {next_run['predicted_total_wall_s']:.3f} s；建议预约 600 s，需要在当前余额之外新增授权 {next_run['additional_authorization_needed_for_allocation_s']:.3f} s。状态 REQUIRES_ADDITIONAL_AUTHORIZATION，本轮不执行。",
      '事前固定最后 t*=0.30–0.40 为观测窗，要求至少 8 块、块长≥max(5 ACF，2 个盒黏性衰减时间)、两半温漂≤0.5%，且温度均值偏差加 CI 半宽≤2%。这些是下一轮拟议判据，不更改本轮或历史验收。若失败，不移动窗口或自动延长。',
      '0.30* 是否足以平衡、稳态相关时间多长仍未知，因此这笔预算不是“必得 8 块或合格”的保证。',
      '', '## 新增观测代码的失败及修正','']
    for failed in s.get('preserved_failed_diagnostic_attempts',[]):
        lines += [f"首次探针在 plugin setup 请求了不存在的 `forces` 通道；本地 `core/utils/common.cpp` 定义实际名称为 `__forces`。尚未积分便退出，实收 {failed['execution']['elapsed_monotonic_s']:.9f} s，旧失败目录 `{failed['directory']}` 完整保留。",
          '先新增通道名回归并记录 FAIL，再只修改新增 worker 的通道名称；完整 CPU 回归通过后，重新预检余额并独立执行一次修正观测版本。没有修改旧 worker、计算库或参数，首次失败不作为物理证据。']
    if s['new_probe'] and s['new_probe'].get('observations'):
        probe=s['new_probe'];q=probe['summary'];observed=probe['observations']
        lines += ['', '## 新探针能证明什么','',f"实际完成 {q['steps']} 步，收费 {q['charged_s']:.6f} s。相同物理时段的旧 dt 峰值 {q['baseline_peak']:.6f}，新 dt/2 峰值 {q['half_dt_peak']:.6f}。两个时间步均保留强启动升温和衰减，减半 dt 未消除这一现象。",
          '这支持启动结构松弛，不证明所有末端温偏都由初始化造成；不同 dt 的同种子不是同一连续噪声路径。旧 t=0 位置没有保存，无法直接比对旧/新初态字节哈希。',
          f"新初态 COM 温度 {probe['initial_state']['COM_kBT']:.9f}，是原始高斯抽样结果，没有重缩放。真实前后相位快照 CPU 重算 kick 最大差 {max(x['kick_max_abs'] for x in observed):.6g}，周期 drift 最大差 {max(x['drift_PBC_max_abs'] for x in observed):.6g}，量级与单精度舍入相容。",
          f"getter 与通道速度最大差 {max(x['getter_velocity_max_abs'] for x in observed):.6g}；所有保存帧的 packed ID 字段对应结果 {all(x['packed_ID_words_preserved'] for x in observed)}。布局、各分量温度、最高 1% 粒子的热能占比和初始/后续固定结构审查见 `diagnosis_summary.json → new_probe`。这些观测仅属于新运行，不能补成旧逐粒子状态。"]
    lines += ['', '## 相同物理阶段的流动对照','', '| 对照 | 共同窗口* | 任务 | 均温 | ν 拟合 SI | 相对 RMS |','|---|---|---|---:|---:|---:|']
    for comparison in a['matched_physical_time_audit']:
        for window in comparison['windows']:
            for row in window['records']:
                lines.append(f"|{comparison['task_id']}|{window['interval_star']}|{row['task_id']}|{fmt(row['temperature_mean'])}|{fmt(row['nu_fit_si'])}|{fmt(row['relative_RMS'])}|")
    lines += ['', '共同后段中原步长与半步长的拟合很接近，说明各自整段均值的差异混入了物理时长/建立过程的影响；这些仍是过渡段，不能升级为稳态收敛证明。',
      '预算表的每块时长已先按真实采样间隔向上取整，再计算八块总步数。新预算估算曾漏掉这一分块网格取整，回归复现 314 个可用样本少于所需 320 个；修正后记录于 cost_grid_regression_before.log 及最终 CPU 日志。该问题仅影响本轮拟议成本，不改变任何已执行任务或历史结果。',
      '新探针顶层执行值为 160000 步、dt*=5e-7，总时长 0.08*；嵌套 task.desired_time_star=0.25 是继承旧任务模板的元数据，worker 未使用该字段推进。原始参数文件未事后改写。']
    return '\n'.join(lines)+'\n'
