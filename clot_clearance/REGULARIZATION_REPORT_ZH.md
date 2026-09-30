# 能量正则化开发报告

本轮沿用现有 `/home/lzy/projects/clot_clearance`，已完成新模块、18 组候选/敏感性计算、验证、报告和无箭头动画。**没有获得达到首次断键/脱落的主算例，也没有证明网格收敛。** 未为了得到碎裂画面提高载荷或降低断键阈值。

## 状态与验收边界

| 状态 | 结论 |
|---|---|
| SOFTWARE VERIFIED | 48 项当前测试通过；18 组主研究的数值/图/身份审计通过；旧强制碎裂重放数组逐位一致。限于已测试软件行为，不代表物理验证。|
| ENERGY REGULARIZED | 仅指本轮定义的循环张开 coupon 实测能量标定和不可逆损伤层；不是任意加载路径下已证明的客观断裂能。|
| NOT CONVERGED | 中/细档平均损伤比粗档下降 56.37% / 67.44%，损伤耗散下降 65.19% / 76.08%；没有失效/脱落事件供比较。|
| NOT_ASSESSED_NO_DETACHMENT | 新模型没有脱落物，不能判定碎片分辨率充分，也不能将 0/0 粉碎指标写成零。|
| FRAGMENT_RESOLUTION_INADEQUATE | 只对重新诊断的旧强制碎裂算例成立：P_singleton=78.899%，P_lowrank=82.569%。其载荷不同，不能作为同条件模型优劣比较。|

提示词软件验收第 5 项“粗、中档至少到首次脱落”未达到；含首次断键的 dt 窗口、断裂局域性、有限尺寸碎片占比与输运收敛亦未获得证据。按提示词要求报告 NOT CONVERGED，不为通过这些指标继续调弱材料。

## 保留与身份

工程没有 Git 仓库，使用源码快照与 SHA-256 确认身份；未发现适用 AGENTS.md。受保护的 627 个旧源码/配置/结果/动画文件全部未变，兼容链接保留。旧 mild、forced-failure 数据未覆盖，BraVa 旧微泡批次及其他停止任务未启动。每个新算例保存配置、标定、源码快照及其哈希。

证据：[保护核对](/home/lzy/projects/clot_clearance/verification/regularization/FINAL_PRESERVATION.json)、[旧结果逐数组重放核对](/home/lzy/projects/clot_clearance/verification/regularization/LEGACY_REPRODUCIBILITY.json)、[48 项测试](/home/lzy/projects/clot_clearance/verification/regularization/FINAL_TESTS_SCOPED.log)、[新文件与全部输出绝对路径索引](/home/lzy/projects/clot_clearance/verification/regularization/SOURCE_OUTPUT_INDEX.md)。第一次全目录 pytest 误收集了受保护的测试副本，报同名模块冲突，保留在 FINAL_TESTS.log；指定 `tests` 后 48 项通过，未删除快照。40 条提示来自 VTK/NumPy 数组 shape API 弃用，不是数值测试失败。

原冻结入口保留原样；新兼容入口 `pd_clot.regularization.legacy` 对 D_break<1 发出警告并保存 `forced_failure_verification=true`。新科学候选强制 D_break=1。旧代码中的 Q_ref 等字段仅因复制配置而保留，新损伤模式不使用它们。

## 几何、场与位置

直管 R=1.25 mm、长度 6 mm、流量 18 mL/min、黏度 0.00345312 Pa·s、流体/固体密度 1056 kg/m³。矩形壁面血栓尺寸 1.25×0.75×0.75 mm，总体积 0.703125 mm³；G=1000 Pa、K=50000 Pa，沿用原 NOSB 材料与稳定化实现。论文原件见 [PDF](/home/lzy/projects/clot_clearance/references/1-s2.0-S2666496826000567-main.pdf)，原核验说明见 [PAPER_READING.md](/home/lzy/projects/clot_clearance/references/PAPER_READING.md)；本轮不是原论文公式的完全复现。

中心现在为 **(-0.825, 0, -0.725) mm**，位于 clot 左面 x=-0.625 mm 上游 0.2 mm，与 clot 的 y/z 中点齐平。这是新计算的场参数，不是移动旧动画图标。白色空心菱形仅标记制造场的局部化中心；它不是顶部向下的点力、实际微泡表面或牵引方向。

`ManufacturedStreamingProvider` 使用 A_y=-UL exp(-r²/(2L²)) 的旋度，局部 u=U exp(-r²/(2L²))(-dz/L,0,dx/L)，叠加原 Poiseuille u/p/grad；U=0.01 m/s，L=0.35 mm，时间系数 1+0.9 sin(2π·500t)，局部压力幅度为 0。局部场解析散度为零，具有涡旋结构和空间衰减。牵引统一由 t=[-pI+μ(grad u+grad uᵀ)]n 计算，脱落后使用同一场的速度。这是人为构造的验证场，没有求解微泡边界或 Navier–Stokes 方程；1 MHz/2 µm 兼容配置字段不构成已求解的超声或泡半径效应。

## A. 断裂能与损伤律

用户授权使用演示值，选取 **Gc_demo=0.01 J/m²**，不是实验标定值。每档使用相同 0.5 mm 立方 coupon：两半刚性块位移控制循环张开、对称弱平面，仅跨平面键按同一损伤律演化，达到 D=1 才失效；没有手动切键。投影裂面面积 2.5×10⁻⁷ m²（不是两个面的面积之和）。标定通过标量求根调整 ψcrit，使实际 NOSB 储能下降累计/裂面面积达到 Gc_demo。以独立六点 Gauss 力–位移积分核验外功。低 C_E 扫参允许延续相同幅度增量最多三倍步数，直到跨平面键全部失效；粗、中、细基准都保持原提前完成的路径。

| Case | h (mm) | Particles | Initial bonds | ψcrit (Pa) | Gc measured (J/m²) | Gc error (%) | First failure / detach N |
|---|---|---|---|---|---|---|---|
| coarse | 0.125 | 360 | 12592 | 0.012520992 | 0.009999937912 | 0.000620878 | Not observed / not observed |
| medium | 0.083333333 | 1215 | 51911 | 0.0345811901 | 0.009999840214 | 0.001597862 | Not observed / not observed |
| fine | 0.0625 | 2880 | 135180 | 0.0479243093 | 0.009999737837 | 0.002621631 | Not observed / not observed |


coupon 外功相对残差范围约 2.15×10⁻¹⁵–2.85×10⁻¹⁴；目标耗散约 2.5×10⁻⁹ J。配置容差为 5%，实测误差远小于该容差。W_cycle,ij=0.5 E_young a_ij² V_ij，a 是循环键伸长的半极差；V_ij 从初始加权邻域分配并满足求和等于材料体积。Wcrit,ij=ψcrit V_ij。ΔD=ΔN C_E[max(W_cycle/Wcrit−threshold,0)]^m_E；基准 C_E=0.001、m_E=1、threshold=0，无非零疲劳门槛。损伤不可逆，D_break=1。spacing/horizon/material/Gc/损伤参数不一致时拒绝复用标定。

这一能量驱动是方向应变能代理，不是 NOSB 自由能的精确逐键分解。指定 coupon 的能量一致不能推出主流场加载下的网格客观性。每个宏步实际储能降为损伤耗散；按 ΔD×驱动分配到各键的累计耗散是诊断估计。每次真正失败会写 failure_ledger.jsonl 和 fracture_energy_ledger.jsonl（失效周期、之前损伤、驱动、累计耗散、h、horizon、支撑秩）；本轮主算例没有失败，所以两个失败账本为空，这是正确结果。

## B. 网格敏感性与体积

三档 horizon/h=3.01、同一物理几何/材料/流场/Gc、25,000 代表周期；每档单独标定。所有颗粒保留，最终均连接固定基底，全部支撑秩为 3。未观测到首次失效或脱落的 N 在 JSON 为 null，不是 N=0。

| Case | Attached (mm³) | Detached (mm³) | Resolved (mm³) | Debris (mm³) | Singleton (mm³) | Rank 0/1 detached (mm³) | Damage dissipation (J) |
|---|---|---|---|---|---|---|---|
| coarse | 0.703125 | 0 | 0 | 0 | 0 | 0 | 2.112127e-17 |
| medium | 0.703125 | 0 | 0 | 0 | 0 | 0 | 7.35266331e-18 |
| fine | 0.703125 | 0 | 0 | 0 | 0 | 0 | 5.05228765e-18 |


**固定层离散是明确的混杂因素。** 统一物理规则为 z<origin_z+0.125 mm，中心恰位于上边界者排除（1e-10 h 舍入容差）；粗/中/细固定体积分数为 1/6、1/9、1/6。因此尽管连续几何相同，边界约束积分并非相同，不能将所有差异归因于损伤律。10% 仅为描述性比较容差，不是数学收敛证明。脱落体积都为零不能证明碎片或断裂收敛。

## C. 周期跳跃与时间步长

| Case | ΔN | dt (µs) | Mean D | Δ mean D (%) | Δ dissipation (%) | Max displacement (µm) | Final signed energy residual (%) |
|---|---|---|---|---|---|---|---|
| coarse | 1000 | 1 | 0.000104395759 | 0 | 0 | 0.26126723 | -0.01266413 |
| cycle_500 | 500 | 1 | 0.000104110926 | -0.272839 | -0.143144 | 0.26127049 | -0.01174284 |
| cycle_250 | 250 | 1 | 0.000103968606 | -0.409167 | -0.214566 | 0.26127212 | -0.01025135 |
| dt_half | 1000 | 0.5 | 0.000104397221 | 0.00140048 | 0.00125719 | 0.26126722 | -0.003166068 |
| medium | 1000 | 1 | 4.55492604e-05 | -56.3687 | -65.1883 | 0.28245005 | -0.01935095 |
| fine | 1000 | 1 | 3.39874662e-05 | -67.4436 | -76.0796 | 0.26347362 | -0.02623106 |


ΔN=1000/500/250 不改变材料或重新标定；在相同代表 N 比较，并检查损伤单调。由于没有事件，首次断键/脱落的周期偏移无法计算，不能写成“零偏移”。dt/2 覆盖完整 25,000 周期，平均损伤差约 0.0014%，但没有包含断裂的时间窗，因此只能支持当前未断裂响应的时间步敏感性判断。

周期跳跃的实际模拟代表力学时长分别为 0.052/0.102/0.202 s，包含初始暖启动；改变 ΔN 也改变实际积分时长。25,000/500=50 s 不能代替当前已模拟的流固运动时间，外功不能凭空乘 ΔN。后续碎片输运的物理时间映射仍需要独立解决。

## D. 局域性

参考坐标中以中心半径 R=0.525 mm 定义影响区域，键用参考中点。没有冻结/删除远端键来限制损伤。

| Case | Damage inside R (%) | Material initially inside R (%) | Damage-weighted distance (mm) | Fixed volume (%) |
|---|---|---|---|---|
| coarse | 50.0597 | 17.7778 | 0.5917913 | 16.6667 |
| medium | 50.7 | 16.7901 | 0.5751262 | 11.1111 |
| fine | 52.7891 | 17.5 | 0.570932 | 16.6667 |


粗档约 17.78% 材料所在区域内集中了 50.06% 的损伤，可描述为损伤对该区域偏聚；仍有约一半损伤在外部，不能称为全部局限于近源区域。断键比例和脱落来源比例因没有对应事件而为 null。

## E. 能量账本

| Final coarse quantity | J |
|---|---|
| elastic_energy_J | 9.32368459013e-14 |
| stabilization_energy_J | 6.46298672971e-15 |
| kinetic_energy_J | 3.25552193412e-17 |
| external_work_J | 2.16972305732e-13 |
| damping_dissipation_J | 1.17191318963e-13 |
| damage_dissipation_estimate_J | 2.1121269959e-17 |
| transport_work_J | 0 |
| numerical_energy_residual_J | -2.74776553105e-17 |


残差定义 R=U_elastic+U_stabilization+K−E_initial+D_damping+D_damage−W_surface−W_transport。外功沿真实代表积分路径梯形积分，阻尼/输运按实际动能改变量记账。粗档末态相对残差 -0.012664%，dt/2 为 -0.003166%；小残差是代码账本证据，不能证明局部力学或 CFD 物理精度。全时程见 energy_history.png/PDF 和 history.csv。

## F–G. 碎片、低秩、输运

诊断阈值为至少 8 颗粒、体积≥1.5625×10⁻¹¹ m³、至少 12 条内键、至少 80% 体积支撑秩≥2。它们是数值分辨要求，不是物理血栓尺寸阈值。单颗粒单独分类，不合并、不删除；其他小连通块为 UNDER_RESOLVED_DEBRIS。P_singleton>20% 或 P_lowrank>30% 会标记 FRAGMENT_RESOLUTION_INADEQUATE 并发出警告；低秩 intrinsic fallback 仅是带诊断的数值延续。新算例 P 值均未定义，所有脱落、resolved/debris/singleton/低秩体积与各类清除量均为零，没有可报告的新碎片尺寸分布。

旧结果重新诊断得到 P_singleton=78.899%、P_lowrank=82.569%，展示颗粒数、体积、等效直径三个分布。新旧加载强度、位置和损伤律均不同，不能把“新算例没有碎片”宣称为粉碎问题已改善。

可选 `transport.mode=fragment_drag` 对 resolved 连通块计算质量、体积、质心、等效球直径、迎风面积、平均同场流速和 Re，使用低 Re Stokes 阻力及 Schiller–Naumann 修正，Re≥1000 使用 Cd=0.44 的球形代理。每半步冻结系数的解析松弛，按质量分配总冲量，保持内部相对速度。参考公式：[OpenFOAM Foundation 源码](https://cpp.openfoam.org/v12/SchillerNaumann_8C_source.html)。未实现力矩、不适用于定量预测不规则可变形血栓。低分辨碎屑仍用旧粒子松弛且分开标记。主算例默认保留 relaxation，因无脱落未触发输运。

独立预分类连通块测试对照 ODE，dt=10/5 µs 的质心速度误差为 1.13×10⁻⁸/5.65×10⁻⁹ m/s，守恒检查通过；这是模块验证，没有在主模型中人为制造碎片。resolved/under-resolved 清除量按首次越线时类别记录且相加等于总量，越过清除平面不能称为临床溶栓。

## 参数敏感性

| Run | Gc demo (J/m²) | m_E | C_E | U (m/s) | Final mean D | Max D | First failure / detach N |
|---|---|---|---|---|---|---|---|
| surface_Gc_0.005_U_0.005 | 0.005 | 1.0 | 0.001 | 0.005 | 5.36094e-05 | 0.0002399158 | Not observed |
| surface_Gc_0.005_U_0.02 | 0.005 | 1.0 | 0.001 | 0.02 | 0.0008426697 | 0.003843757 | Not observed |
| surface_Gc_0.02_U_0.005 | 0.02 | 1.0 | 0.001 | 0.005 | 1.301727e-05 | 5.824638e-05 | Not observed |
| surface_Gc_0.02_U_0.02 | 0.02 | 1.0 | 0.001 | 0.02 | 0.0002044748 | 0.0009310573 | Not observed |
| sweep_C_E_0.0005 | 0.01 | 1.0 | 0.0005 | 0.01 | 8.28509e-05 | 0.0003757651 | Not observed |
| sweep_C_E_0.002 | 0.01 | 1.0 | 0.002 | 0.01 | 0.0001315396 | 0.0005966815 | Not observed |
| sweep_Gc_demo_J_m2_0.005 | 0.005 | 1.0 | 0.001 | 0.01 | 0.0002113513 | 0.0009589593 | Not observed |
| sweep_Gc_demo_J_m2_0.02 | 0.02 | 1.0 | 0.001 | 0.01 | 5.131418e-05 | 0.000232709 | Not observed |
| sweep_m_E_0.8 | 0.01 | 0.8 | 0.001 | 0.01 | 0.001426735 | 0.005333286 | Not observed |
| sweep_m_E_1.2 | 0.01 | 1.2 | 0.001 | 0.01 | 7.726234e-06 | 4.258811e-05 | Not observed |
| sweep_streaming_velocity_scale_m_s_0.005 | 0.01 | 1.0 | 0.001 | 0.005 | 2.648195e-05 | 0.000118501 | Not observed |
| sweep_streaming_velocity_scale_m_s_0.02 | 0.01 | 1.0 | 0.001 | 0.02 | 0.0004160712 | 0.001895644 | Not observed |


Gc、m_E、C_E 改变时按新值重新执行 coupon；速度改变时复用同一能量标定。因此 m_E/C_E 图是**保持目标 coupon Gc 的联合标定比较**，不是固定 ψcrit 的单参数偏导。额外四角计算与基准/OAT 构成 3×3 Gc–U 网格，response_surface 显示全部九个采样点；没有连续拟合、择优保留或为碎裂调参。所有 18 组均无失效或脱落，最大粒子损伤约 0.00533（m_E=0.8）。每组 hash、事件、体积、能量残差保存在 [CSV](/home/lzy/projects/clot_clearance/verification/regularization/comparisons/study_summary.csv)。

## H. 导入器与最小下一步

`ResolvedStreamingFieldProvider` 支持 VTU/legacy VTK 的线性四面体原单元 P1 插值，以及显式声明 convex_hull 的 VTP/CSV 采样插值。u/p 当前要求节点值。体数据缺 grad 时按 P1 一致计算；平面数据必须给完整三维 grad，不能推断法向导数。拒绝缺单位、NaN、越界、不可覆盖的表面、过大散度及不一致梯度，保存观察范围与可选范围阈值诊断；只声明 outward_solid 并不构成法向几何方向的证明。当前为静态快照，不支持跨时间插值，非凸表面/孔洞不适用 convex_hull 近似。

四种格式的独立仿射测试通过；SI 最大误差 u=8.67×10⁻¹⁹ m/s、p=1.14×10⁻¹³ Pa、grad=7.22×10⁻¹⁶ s⁻¹。另一个弱均匀合成 VTU 已通过实际 PD runner 一宏步接入审计，见 [接口报告](/home/lzy/projects/clot_clearance/verification/regularization/resolved_import_with_runner/IMPORT_VERIFICATION.json) 和 [runner 审计](/home/lzy/projects/clot_clearance/verification/regularization/resolved_import_with_runner/runner_case/AUDIT.json)。它们不是实际微泡 CFD。

最小下一步：取得一份覆盖血栓全部受载表面及预期运动范围的**实际微泡解析体网格 VTU 快照**（节点 u、p，注明压力基准和坐标配准，最好同时提供 grad），给出以下显式元数据，在复制的新配置内设 `streaming.provider="resolved_field"`、`field_path`、`metadata_path`，先做导入诊断和一宏步对照；无需修改 PD 求解器。静态平均场的循环驱动含义仍需重新定义，不能把静态导入自动解释为超声疲劳。

```json
{
  "length_unit": "m", "velocity_unit": "m/s", "pressure_unit": "Pa", "time_unit": "s",
  "normal_convention": "outward_solid",
  "arrays": {"velocity": "velocity", "pressure": "pressure", "velocity_gradient": "velocity_gradient"},
  "validation": {"maximum_relative_divergence": 0.05, "maximum_relative_gradient_inconsistency": 0.25}
}
```


材料、疲劳律均未实验/声致溶栓标定；制造场不是包膜微泡 CFD；尚无双向流体–PD 耦合；表面面积/法向、固定层离散、稳定化与低秩延续仍是精度限制。完成软件开发不等于完成所要求的科学验收。

## 动画、图片与复现

[本地可视化入口](/home/lzy/projects/clot_clearance/visualization/streaming_regularized_demo/index.html) · [MP4](/home/lzy/projects/clot_clearance/visualization/streaming_regularized_demo/regularized_damage_transport.mp4) · [GIF](/home/lzy/projects/clot_clearance/visualization/streaming_regularized_demo/regularized_damage_transport.gif) · [QC](/home/lzy/projects/clot_clearance/visualization/streaming_regularized_demo/QC.json)。视频 1920×1080、60 fps、751 帧、12.5167 s，黑底、Arial 常规字重、单损伤主图、无箭头/大蓝泡/Force 子图/N 顶标、固定相机与真实位移。源标记静止仅表示空间位置，不用其亮度伪造载荷相位。色标固定 0–0.0005，黄色不是 D=1；粒子显示损伤是邻域键完整度推导的值，断键阈值作用在键 D 上。

26 个原始状态每隔 30 显示帧严格对齐，位置只做线性插值、损伤保持前一宏步原值；固定点逐位不变，全部粒子保留。运动最大仅约 0.261 µm，某些相邻解码帧相同并不代表编码掉帧；60 fps 时间戳恒定，损伤宏步跳变也没有被虚构中间科学状态掩盖。四阶段图用 N=0/8000/17000/25000，C/D 明确写未到断裂/输运，局部流线均为代表峰值相位的说明叠图。

- [cycle_jump_convergence PNG](/home/lzy/projects/clot_clearance/visualization/streaming_regularized_demo/figures/cycle_jump_convergence.png) · [PDF](/home/lzy/projects/clot_clearance/visualization/streaming_regularized_demo/figures/cycle_jump_convergence.pdf)
- [energy_history PNG](/home/lzy/projects/clot_clearance/visualization/streaming_regularized_demo/figures/energy_history.png) · [PDF](/home/lzy/projects/clot_clearance/visualization/streaming_regularized_demo/figures/energy_history.pdf)
- [four_stage_summary PNG](/home/lzy/projects/clot_clearance/visualization/streaming_regularized_demo/figures/four_stage_summary.png) · [PDF](/home/lzy/projects/clot_clearance/visualization/streaming_regularized_demo/figures/four_stage_summary.pdf)
- [fracture_energy_calibration PNG](/home/lzy/projects/clot_clearance/visualization/streaming_regularized_demo/figures/fracture_energy_calibration.png) · [PDF](/home/lzy/projects/clot_clearance/visualization/streaming_regularized_demo/figures/fracture_energy_calibration.pdf)
- [fragment_size_distribution PNG](/home/lzy/projects/clot_clearance/visualization/streaming_regularized_demo/figures/fragment_size_distribution.png) · [PDF](/home/lzy/projects/clot_clearance/visualization/streaming_regularized_demo/figures/fragment_size_distribution.pdf)
- [fragment_size_distribution_detail PNG](/home/lzy/projects/clot_clearance/visualization/streaming_regularized_demo/figures/fragment_size_distribution_detail.png) · [PDF](/home/lzy/projects/clot_clearance/visualization/streaming_regularized_demo/figures/fragment_size_distribution_detail.pdf)
- [localization_metrics PNG](/home/lzy/projects/clot_clearance/visualization/streaming_regularized_demo/figures/localization_metrics.png) · [PDF](/home/lzy/projects/clot_clearance/visualization/streaming_regularized_demo/figures/localization_metrics.pdf)
- [manufactured_field PNG](/home/lzy/projects/clot_clearance/visualization/streaming_regularized_demo/figures/manufactured_field.png) · [PDF](/home/lzy/projects/clot_clearance/visualization/streaming_regularized_demo/figures/manufactured_field.pdf)
- [mesh_convergence PNG](/home/lzy/projects/clot_clearance/visualization/streaming_regularized_demo/figures/mesh_convergence.png) · [PDF](/home/lzy/projects/clot_clearance/visualization/streaming_regularized_demo/figures/mesh_convergence.pdf)
- [parameter_sensitivity PNG](/home/lzy/projects/clot_clearance/visualization/streaming_regularized_demo/figures/parameter_sensitivity.png) · [PDF](/home/lzy/projects/clot_clearance/visualization/streaming_regularized_demo/figures/parameter_sensitivity.pdf)
- [response_surface PNG](/home/lzy/projects/clot_clearance/visualization/streaming_regularized_demo/figures/response_surface.png) · [PDF](/home/lzy/projects/clot_clearance/visualization/streaming_regularized_demo/figures/response_surface.pdf)
- [time_step_sensitivity PNG](/home/lzy/projects/clot_clearance/visualization/streaming_regularized_demo/figures/time_step_sensitivity.png) · [PDF](/home/lzy/projects/clot_clearance/visualization/streaming_regularized_demo/figures/time_step_sensitivity.pdf)
- [topology PNG](/home/lzy/projects/clot_clearance/visualization/streaming_regularized_demo/figures/topology.png) · [PDF](/home/lzy/projects/clot_clearance/visualization/streaming_regularized_demo/figures/topology.pdf)

复现命令（新目录存在则换名；不会覆盖未完成结果）：

```bash
cd /home/lzy/projects/clot_clearance
export PYTHONDONTWRITEBYTECODE=1
export NUMBA_CACHE_DIR="$PWD/build_regularization/numba_cache"
export MPLCONFIGDIR="$PWD/build_fragmentation/matplotlib"
export PYTHONPATH="$PWD"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
CL_PY=/home/lzy/projects/computation_examples/_shared/python_environment/bin/python

# Coupon: use a new output directory; units are m and J/m².
"$CL_PY" -B -m pd_clot.regularization.coupon --config configs/streaming_regularized_demo_coarse.json --spacing 0.000125 --Gc-demo 0.01 --output verification/regularization/replay_coupon
# Demo using that new calibration.
"$CL_PY" -B -m pd_clot.regularization.runner --config configs/streaming_regularized_demo_coarse.json --calibration verification/regularization/replay_coupon/CALIBRATION.json --output results/streaming_regularized_demo_replay/coarse
# Medium/fine coupons are independently recalibrated by the mesh command.
"$CL_PY" -B -m pd_clot.regularization.campaign --study mesh --output-root results/streaming_regularized_demo_replay
# Cycle tests retain the baseline material and coarse calibration.
"$CL_PY" -B -m pd_clot.regularization.campaign --study cycle --output-root results/streaming_regularized_demo_replay
"$CL_PY" -B -m pd_clot.regularization.campaign --study timestep --output-root results/streaming_regularized_demo_replay
"$CL_PY" -B -m pd_clot.regularization.campaign --study sweep --output-root results/streaming_regularized_demo_replay
"$CL_PY" -B -m pd_clot.regularization.campaign --study sweep_surface --output-root results/streaming_regularized_demo_replay
# Current tests only; preserved snapshot directories are not a test suite.
"$CL_PY" -B -m pytest tests -q -p no:cacheprovider
"$CL_PY" -B scripts/verify_regularized_interfaces.py --mode drag --output verification/regularization/drag_replay
"$CL_PY" -B scripts/verify_regularized_interfaces.py --mode import --output verification/regularization/import_replay
# Presentation replay requires a fresh output directory.
"$CL_PY" -B scripts/render_regularized_damage.py --run results/streaming_regularized_demo/coarse --output visualization/streaming_regularized_demo_replay
"$CL_PY" -B scripts/check_regularized_damage_v2.py --output visualization/streaming_regularized_demo_replay
"$CL_PY" -B -m pd_clot.regularization.figures --output visualization/streaming_regularized_demo_replay/figures
```

