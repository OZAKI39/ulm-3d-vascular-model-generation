# SDPD 原因排查报告

**结论：已证实启动过渡进入统计段，以及均值/置信区间样本不一致的分析错误；尚不能将全部温偏或剖面残差归结为原生求解器错误。**
诊断状态 `DIAGNOSED_WITH_LIMITATIONS`；`selection=null`；人工验收 PENDING。诊断完成不等于材料合格。

## 来源与复现

权威旧包：`/home/lzy/projects/mirheo_starter/data/fluid_model_comparison/result_7cc1db5bcaec2c16`。按配置和 SHA-256 定位；18 项旧/本轮实验从 CSV 和快照重新分析，修复前温度、压力、黏度、窗口、性能和退出状态精确复现。
当前冻结单位 L0=4.9999999999999998e-07 m、M0=1.6499999999999997e-17 kg、t0=3.1655773384195756e-05 s；ν*=414.057515865、μ*=3312.46012692，两者相差真实质量密度 8。没有重设单位。
执行计划、实际步数、脚本/库哈希及退出原因见 `provenance.json`；源调用链行号及哈希见 `native_source_contract.json`。

## 原因与证据状态

| 问题 | 状态 | 一句话解释 |
|---|---|---|
|温度测量|RULED_OUT_WITHIN_TESTED_SCOPE|保存的分箱矩可独立重算 worker 的原始、COM 和粗分箱温度，未发现能解释大幅升温的质量或自由度错误。|
|温度偏高|CONFIRMED|选入统计的轨迹仍包含初始化后的降温过程，所报均值不能当作已证明稳定的平台温度。|
|初始化结构与背景压力|SUPPORTED_BUT_NOT_CONFIRMED|真实初态温度约为目标值，dt 与 dt/2 均出现近似的强温升与衰减；结合核密度和保守势能下降，更支持初始结构松弛，完整因果能量收支尚未测得。|
|时间推进|SUPPORTED_BUT_NOT_CONFIRMED|本地纯液体调用链实际采用一次 kick-drift；固定耗散模态会有正温度步长偏差，但不能把简化估计直接套到液体。|
|流动温度扣除|RULED_OUT_WITHIN_TESTED_SCOPE|平衡任务也升温，粗细分箱差异很小；漏扣平均 Poiseuille 流速不足以解释当前整体温升。|
|只有一个统计块|CONFIRMED|预算缩短轨迹、初始化截段及漂移信号的长 ACF 联合留下一个完整块；块数整除本身没有算错。|
|均值与置信区间|CONFIRMED|旧均值包含不完整尾部，而区间只用完整块；已通过先失败后通过的回归测试修复，并保留全窗口均值。|
|平衡段判据|CONFIRMED|现有截段只搜索 20–50% 且允许 10% 两半漂移，这不足以证明满足 2% 温度要求的稳态。|
|黏度与剖面不符|SUPPORTED_BUT_NOT_CONFIRMED|拟合系数只描述一个投影；过渡和不足统计使剖面剩余形状无法可靠归因，有限作用范围也尚未排除。|
|压力覆盖|CONFIRMED|端点减基准的名义余量小于其统计不确定性；参考压力共享产生相关误差，移动零点不会增加覆盖能力。|
|重复 u.run|RULED_OUT_WITHIN_TESTED_SCOPE|本地调用链保留物理状态、累计时间和 SDPD 随机发生器；未发现每次采样重新初始化液体。|
|可观测性限制|INSUFFICIENT_EVIDENCE|旧位置快照不足以恢复逐粒子速度、力、随机状态或可靠续跑。|

每项支持证据、反对证据、缺失观测和最小验证方法完整保存在 `evidence_matrix.json`。

## 温度：测量、过渡与推进偏差分开

| 任务 | 首个 worker 值 | 峰值 | 最后 10% 均值 | 原统计均值 | 完整块均值 | N−1 修正相对影响 |
|---|---:|---:|---:|---:|---:|---:|
|equilibrium|1.026191|1.119522|1.023738|1.026023|1.026086|0.0002442002|
|high_gamma_equilibrium|1.064029|1.151048|1.052168|1.051766|1.051781|0.0002442002|
|sdpd_probe|2.917892|9.389353|9.365473|未测/无有效结果|未测/无有效结果|0.0002442002|
|sdpd_equilibrium|2.917892|9.389353|1.089941|1.146279|1.146279|0.0002442002|
|sdpd_flow|2.901968|9.217657|1.036301|1.061956|1.069992|0.0002442002|
|sdpd_half_dt|2.856025|9.237247|1.099109|1.181124|1.181124|0.0002442002|
|sdpd_half_force|2.901957|9.219006|1.05076|1.089712|1.089712|0.0002442002|
|sdpd_eos_low|3.029424|10.39643|1.104079|1.187201|1.187201|0.0002464268|
|sdpd_eos_high|3.059031|10.08444|1.100458|1.178119|1.178119|0.0002339729|

最后 10% 只是预先列出的诊断窗口之一，不是重新挑选的正式验收窗口。所有 20–90% 截段及两半漂移均保存并显示。
等质量平衡温度采用 m Σ|v−vCOM|²/[3(N−1)]；流动温度扣除实际瞬时分箱均速并计入 3(N−B) 自由度。空箱不贡献均速或自由度。原始 profile 的 count、平均三分量和 sum_v2 可重算原始、COM 和粗分箱温度；细分箱逐箱矩没有保存。
平衡任务的 CPU 矩重算最大 COM 差值为 3.552714e-15；其 N−1 修正仅 0.02442%，不能解释约十几个百分点的偏高。
原生 Stats 样本在第 1、every+1、…步；worker 在 every、2every、…步。相邻步差异图明确标注该一积分步差，不任意平移来强行重合。两者都读取同一底层速度，矩复核也不是完全独立的物理观测。
旧速度分量方差、逐粒子速度/力、随机状态和 t=0 位置未保存；不能恢复粒子异常尾部分布或构造可信续跑。

## 初始化与时间推进

Uniform 在单位胞中生成随机位置，其种子由 rank 和粒子向量名字确定，并非均匀晶格。Python 随后按固定种子设高斯速度并只扣一次 COM。较大背景压力对配对力有实际贡献。
真实位置快照的核密度、近邻、最近距离、保守力、势能变化和耗散矩阵量级见 `snapshot_audit.json`。这些是固定快照上的独立公式核查，没有开发 CPU 流体求解器。
完整纯液体调用链的更新为 v_new=v_old+F(x_old,v_old)dt/m，x_new=x_old+v_new dt。不能仅凭 VelocityVerlet 类名假定存在另一次半步 kick 或中心速度转换。该语义本身不是已证实的原生 bug。
固定 OU 耗散模态的解析平衡方差比为 1/(1−λdt/2)，说明显式耗散存在正步长偏差的机制。快照耗散矩阵给的是量级估计，不能直接乘一个系数修正实际液体温度。
按核快照平均耗散模态估算，旧 dt 下温偏量级约 0.5–0.6%，远小于强启动峰值；这不是完整流体误差上界。声学、粒子间黏性、盒动量衰减、实际速度位移分别见 `time_scale_audit.json`；结构松弛时标尚未建立。
原生噪声幅度 sqrt(2·5μ·kBT/dt) 只在速度推进时再乘一次 dt；配对随机标量按排序 ID 生成，方向反对称。生成器随当前时间更新，同一步重复计算复用种子。实际设备随机量分布未被完整统计验证。
反复 u.run 会重建调度和 cell lists、重新 setup 插件，但保留当前粒子、时间和 SDPD 随机发生器；已检查范围内未发现重置初态。重排可能改变浮点累加和后续混沌轨迹，未声称连续/分块逐位相同。

## 一块统计的完整账单

| 任务/观测 | 原计划/分配/实际步数 | 初始化后样本 | 间隔* | ACF 样本 | 最小块样本 | ACF 要求样本 | 最终块样本 | 尾部 | 块数 |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
|sdpd_equilibrium/temperature|250000/172000/172000|431|0.0002|115.6362|40|579|579|431|0|
|sdpd_equilibrium/pressure|250000/172000/172000|430|0.0002|7.221263|40|37|40|30|10|
|sdpd_flow/temperature|400000/270400/270400|677|0.0002|101.9736|40|510|510|167|1|
|sdpd_flow/pressure|400000/270400/270400|676|0.0002|21.55116|40|108|108|28|6|
|sdpd_flow/viscosity_slope|400000/270400/270400|677|0.0002|82.44347|40|413|413|264|1|
|sdpd_half_dt/temperature|700000/302800/302800|379|0.0002|109.4531|40|548|548|379|0|
|sdpd_half_dt/pressure|700000/302800/302800|378|0.0002|12.98206|40|65|65|53|5|
|sdpd_half_dt/viscosity_slope|700000/302800/302800|379|0.0002|112.4502|40|563|563|379|0|
|sdpd_half_force/temperature|350000/221200/219800|550|0.0002|116.4499|40|583|583|550|0|
|sdpd_half_force/pressure|350000/221200/219800|549|0.0002|2.039702|40|11|40|29|13|
|sdpd_half_force/viscosity_slope|350000/221200/219800|550|0.0002|98.18969|40|491|491|59|1|
|sdpd_eos_low/temperature|200000/155600/155600|390|0.0002|119.0403|40|596|596|390|0|
|sdpd_eos_low/pressure|200000/155600/155600|389|0.0002|7.524926|40|38|40|29|9|
|sdpd_eos_high/temperature|200000/155600/155600|390|0.0002|112.3576|40|562|562|390|0|
|sdpd_eos_high/pressure|200000/155600/155600|389|0.0002|4.96159|40|25|40|29|9|

基准流动：实际 270400 步、dt*=1e-06、物理演化 8.559721 µs；从 t*=0.1352 开始只剩 677 个样本。块长由 max(物理下限 40, ceil(5τ)=413) 得 413，因此 677 // 413 = 1，余 264。
ACF 单位是样本数；最终乘采样间隔才是模拟时间。所有 CSV 的单调性、规则采样、step/time 和 profile 分箱顺序已核对。实际停止任务按 completion 和 ledger 记录，未当作完整执行。
ACF 受未衰减趋势影响，当前 τ 不能作为可信稳态相关时间。原规则允许 10% 两半漂移，和 2% 温度目标不等价；本轮不事后放宽或修改验收阈值。

## 黏度和剖面

| 任务 | 原点估计 ν SI | 完整块 ν SI | 原相对 RMS | 完整块相对 RMS | 近变号面 RMS* | 内部 RMS* |
|---|---:|---:|---:|---:|---:|---:|
|sdpd_flow|3.099559e-06|3.021927e-06|0.1222802|0.1371241|0.05502405|0.03317508|
|sdpd_half_dt|2.658362e-06|2.658362e-06|0.1862352|0.1862352|0.09837808|0.05779115|
|sdpd_half_force|3.239939e-06|3.209188e-06|0.3051304|0.3133272|0.06564061|0.03965872|

主解析解使用每粒子力 F、粒子质量 m 和运动黏度 ν，u=F/(mν)·分段抛物线；包含有限箱宽的二次平均修正。平均流速点、目标线、拟合线与所有位置残差均保存；未删除不利区域。
rc/半通道宽=1/4，支持范围跨 4 个剖面分箱；非局部作用与连续分段抛物线之间存在有限尺度疑点。数据中未发现空箱，粒子加权与等时间平均差异另列。缺乏稳态独立块，尚不能将空间残差确认为核作用或噪声中的任何一种。
半力/半步的总物理时长不同，短记录仍在初始化；直接比较各自全窗拟合差异不构成稳态敏感性验证。修复后的点估计仍是“拟合估计”，没有获得有效 CI。

## 压力参考、覆盖和精度

全盒密度代入输入 EOS 的背景量为 3793.679 Pa；这不是已知真实绝对压力，也不是可任意移除的无影响配对力。
机械压力使用配对应力迹除以 3V，加一次 COM 扣除后的动能动量通量。EOS 已用于保守配对力，不再另加 EOS 压力。应力来自积分前、动量来自积分后，离散相位差 dt 明确保留。

| 端点 | 端点−基准 Pa | 需要表压 Pa | 名义余量 Pa | 条件 Welch 半宽 Pa | 保守端点+基准半宽 Pa |
|---|---:|---:|---:|---:|---:|
|sdpd_eos_low|-87.22629|-13.70063|73.52566|59.31711|88.9672|
|sdpd_eos_high|137.1762|132.2045|4.971612|48.35275|73.25181|

两个压差共享同一基准，所以其协方差含 Var(P_ref)；高低端点差中基准会抵消。Welch 估计假定运行间独立，实际共同随机种子可能产生协方差，未测的协方差不置为实测零。原保守端点规则不需要该独立假设，但不宣称联合 95% 覆盖。
virial CSV 舍入最大半量级约 0.0003215941 Pa；单精度一个数值分辨尺度约 0.0003916184 Pa。它们远小于当前统计不确定性；这不是完整 GPU 求和误差证明，也不支持因此重编译。
各 EOS 点仍有温度和结构问题，故即使几何区间覆盖也不等于已得到合格等温压力响应。

## GPU 决策与下一步

当前原授权累计 3492.919 / 3600 s，余额 107.0807 s。
筛选探针 `equilibrium_half_dt_trend_corrected_channel`：160000 步，dt*=5e-07，物理 2.532462 µs，预测墙钟 233.955 s，预约上限 260 s。它只判断早期温升/衰减的步长依赖，不用于完整温度或黏度收敛。
执行结果：{"status": "COMPLETED_FIXED_PHYSICAL_WINDOW", "charged_s": 195.86299424999743, "steps": 160000, "physical_duration_s": 2.5324618707356607e-06, "baseline_peak": 9.389353015164602, "half_dt_peak": 9.376462297026295, "peak_relative_difference": -0.0013729080286455364, "window_means": [{"interval_star": [0.0, 0.02], "baseline_dt": {"samples": 100, "mean": 6.894225821874479, "min": 2.9178915803546857, "max": 9.389353015164602, "complete_interval": true, "CI95": null, "CI_reason": "Single drifting trajectory; intervals are not independent repeats."}, "probe_half_dt": {"samples": 100, "mean": 6.874757968008514, "min": 2.8700215836531027, "max": 9.376462297026295, "complete_interval": true, "CI95": null, "CI_reason": "Single drifting trajectory; intervals are not independent repeats."}}, {"interval_star": [0.02, 0.04], "baseline_dt": {"samples": 100, "mean": 3.0063554741775445, "min": 2.1495815048963234, "max": 4.3317413908253535, "complete_interval": true, "CI95": null, "CI_reason": "Single drifting trajectory; intervals are not independent repeats."}, "probe_half_dt": {"samples": 100, "mean": 3.017870966706486, "min": 2.1684938124456745, "max": 4.2925475635248445, "complete_interval": true, "CI95": null, "CI_reason": "Single drifting trajectory; intervals are not independent repeats."}}, {"interval_star": [0.04, 0.06], "baseline_dt": {"samples": 100, "mean": 1.8349907845575104, "min": 1.590816112551783, "max": 2.203597310427226, "complete_interval": true, "CI95": null, "CI_reason": "Single drifting trajectory; intervals are not independent repeats."}, "probe_half_dt": {"samples": 100, "mean": 1.824156833125243, "min": 1.550728168750045, "max": 2.1541247657778024, "complete_interval": true, "CI95": null, "CI_reason": "Single drifting trajectory; intervals are not independent repeats."}}, {"interval_star": [0.06, 0.08], "baseline_dt": {"samples": 100, "mean": 1.4406487576764402, "min": 1.2983949454281096, "max": 1.640153006524971, "complete_interval": true, "CI95": null, "CI_reason": "Single drifting trajectory; intervals are not independent repeats."}, "probe_half_dt": {"samples": 100, "mean": 1.4331029306042984, "min": 1.2857183278443913, "max": 1.6043981173464865, "complete_interval": true, "CI95": null, "CI_reason": "Single drifting trajectory; intervals are not independent repeats."}}], "moment_COM_max_abs_error": 3.552713678800501e-15, "interpretation": "Both matched unforced runs retain a large early heating peak and decay. Halving dt does not remove that transient; initialization/structure relaxation is supported. Stable late temperature bias is not tested.", "qualification_capable": false, "same_noise_path": false, "direct_old_initial_state_hash_available": false, "remaining_own_processes": []}

| 后续条件情景 | 步数 | 物理时长 µs | 预测墙钟 s | 单任务600秒可行 |
|---|---:|---:|---:|---|
|optimistic_relaxation_floor|364000|11.522702|523.391|True|
|observed_drifting_ACF_upper_scenario|960800|30.414867|1352.011|False|

上述是有明确假设的预算情景：稳态相关时间尚未可靠测得，不能机械把当前时间乘八。超过原余额的工作标为 REQUIRES_ADDITIONAL_AUTHORIZATION；超过单任务上限还需要先验证完整状态及随机过程的续跑。
**下一步只优先建立无驱动目标液体的可信热平衡平台及独立块。** 这同时决定黏度拟合、EOS 和 dt 敏感性是否有意义，优先级高于新求解器或盲目调参数。

## 验证与范围

CPU 回归记录及修复前失败用例见 setup；最终 CPU 验证记录绑定 provenance。浏览器只有实际检查后才在单独 browser_checks.json 写 PASS，人工验收始终 PENDING。
未修改原始数据、历史验收、Mirheo/CUDA/MPI 或兼容补丁；未运行真实血管、RBC、微泡、SDF 或开放边界。

## 唯一优先补测的具体预算

建议无驱动、原 dt*=1e-06、固定 400000 步，即 12.662309 µs。按旧无驱动实测成本加 15% 和 18 s 余量，预计 582.124 s；建议预约 600 s，需要在当前余额之外新增授权 492.919 s。状态 REQUIRES_ADDITIONAL_AUTHORIZATION，本轮不执行。
事前固定最后 t*=0.30–0.40 为观测窗，要求至少 8 块、块长≥max(5 ACF，2 个盒黏性衰减时间)、两半温漂≤0.5%，且温度均值偏差加 CI 半宽≤2%。这些是下一轮拟议判据，不更改本轮或历史验收。若失败，不移动窗口或自动延长。
0.30* 是否足以平衡、稳态相关时间多长仍未知，因此这笔预算不是“必得 8 块或合格”的保证。

## 新增观测代码的失败及修正

首次探针在 plugin setup 请求了不存在的 `forces` 通道；本地 `core/utils/common.cpp` 定义实际名称为 `__forces`。尚未积分便退出，实收 1.869825645 s，旧失败目录 `/home/lzy/projects/mirheo_starter/runs/sdpd_diagnostics/thermal_cause_20260909/equilibrium_half_dt_trend` 完整保留。
先新增通道名回归并记录 FAIL，再只修改新增 worker 的通道名称；完整 CPU 回归通过后，重新预检余额并独立执行一次修正观测版本。没有修改旧 worker、计算库或参数，首次失败不作为物理证据。

## 新探针能证明什么

实际完成 160000 步，收费 195.862994 s。相同物理时段的旧 dt 峰值 9.389353，新 dt/2 峰值 9.376462。两个时间步均保留强启动升温和衰减，减半 dt 未消除这一现象。
这支持启动结构松弛，不证明所有末端温偏都由初始化造成；不同 dt 的同种子不是同一连续噪声路径。旧 t=0 位置没有保存，无法直接比对旧/新初态字节哈希。
新初态 COM 温度 0.997714377，是原始高斯抽样结果，没有重缩放。真实前后相位快照 CPU 重算 kick 最大差 4.36401e-07，周期 drift 最大差 1.19208e-07，量级与单精度舍入相容。
getter 与通道速度最大差 0；所有保存帧的 packed ID 字段对应结果 True。布局、各分量温度、最高 1% 粒子的热能占比和初始/后续固定结构审查见 `diagnosis_summary.json → new_probe`。这些观测仅属于新运行，不能补成旧逐粒子状态。

## 相同物理阶段的流动对照

| 对照 | 共同窗口* | 任务 | 均温 | ν 拟合 SI | 相对 RMS |
|---|---|---|---:|---:|---:|
|sdpd_half_dt|[0, 0.0757]|sdpd_flow|3.356441|2.13853e-06|0.2927575|
|sdpd_half_dt|[0, 0.0757]|sdpd_half_dt|3.349713|2.127309e-06|0.2933159|
|sdpd_half_dt|[0.0757, 0.1514]|sdpd_flow|1.182664|2.657524e-06|0.180161|
|sdpd_half_dt|[0.0757, 0.1514]|sdpd_half_dt|1.181124|2.658362e-06|0.1862352|
|sdpd_half_force|[0, 0.1099]|sdpd_flow|2.700583|2.256251e-06|0.2512869|
|sdpd_half_force|[0, 0.1099]|sdpd_half_force|2.700016|2.452084e-06|0.5385464|
|sdpd_half_force|[0.1099, 0.2198]|sdpd_flow|1.089986|2.937197e-06|0.1398032|
|sdpd_half_force|[0.1099, 0.2198]|sdpd_half_force|1.089712|3.239939e-06|0.3051304|

共同后段中原步长与半步长的拟合很接近，说明各自整段均值的差异混入了物理时长/建立过程的影响；这些仍是过渡段，不能升级为稳态收敛证明。
预算表的每块时长已先按真实采样间隔向上取整，再计算八块总步数。新预算估算曾漏掉这一分块网格取整，回归复现 314 个可用样本少于所需 320 个；修正后记录于 cost_grid_regression_before.log 及最终 CPU 日志。该问题仅影响本轮拟议成本，不改变任何已执行任务或历史结果。
新探针顶层执行值为 160000 步、dt*=5e-7，总时长 0.08*；嵌套 task.desired_time_star=0.25 是继承旧任务模板的元数据，worker 未使用该字段推进。原始参数文件未事后改写。
