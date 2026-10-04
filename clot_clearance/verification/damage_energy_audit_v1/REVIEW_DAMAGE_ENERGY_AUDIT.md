# 固定坐标损伤更新的储能机制审计

审计日期：2026-10-01 UTC。工作目录：`/home/lzy/projects/temp_storage/clot_clearance_github_sync_20260930/clot_clearance`。
基线分支：`clot-clearance/energy-regularization-20260930`；Git HEAD：`f0e891686b9b107510ab3292aa27ce83638b3339`。基线还包含上一轮已存在但未提交的机械平衡/侵蚀文件，其身份通过审计前 SHA-256 清单保护，不能仅用 HEAD 代替。

**审计已完成，35 个已保存跳变的增能来源已定位：重构项使本构储能增加，幅度超过 remaining 降低及稳定化能量下降之和。生产模型未修改，耗散性仍未建立。**

## 1. 为什么做这一轮

上一轮已经在当前开发算例上完成机械平衡验证、ADR 与小算例参考解对照，以及侵蚀慢时间步一致性检查。但 0.02、0.01、0.005 s 三档慢步长分别有 5、10、20 次损伤更新，在粒子位置不变时都使储能增加。这不是可用启动瞬态解释的现象。本轮只调查原因，不求新的平衡解、不推进侵蚀、不生成替代物理结果。

已审阅指定的七个力学/侵蚀文件及机械平衡报告；旧报告中的待论文阶段与 `adr_completion` 完成阶段作了区分。未发现适用的 AGENTS.md。实现前的代数审查见 [PRE_AUDIT.md](PRE_AUDIT.md)，预声明容差与抽样计划见 [PREDECLARED_PLAN.json](provenance/PREDECLARED_PLAN.json)。

历史保护结果：审计前记录的 **3706 个文件全部未变**，包括要求保护的四个 results/verification 目录。见 [PRESERVATION_CHECK.json](provenance/PRESERVATION_CHECK.json)。没有修改生产 solver、配置、材料、损伤率、ADR 容差或历史 JSON/CSV/NPZ/图。

新增代码仅有 8 个文件：`pd_clot/regularization/damage_energy_audit.py`；`scripts/run_damage_energy_audit.py`、`analyze_damage_energy_audit.py`、`check_damage_reconstruction_identity.py`、`visualize_damage_energy_audit.py`；三个指定的新测试文件。其余新增成果均在本审计目录。

测试：原有 **112 passed**，新增 **18 passed**，全量 **130 passed / 0 failed / 0 skipped**。全量耗时 36.37 s，60 个 warning 中 57 个是既有 VTK/NumPy 弃用提示，3 个来自原有测试有意保留的未平衡/能量未建立警告。见 [全量日志](tests/FULL_TESTS.log) 与 [新增测试日志](tests/NEW_TESTS_FIRST.log)。新增检查覆盖分项和、精确分解、未变形零能、均匀缩放、历史回放、差分稳定性、输入只读及生产结果等价。

硬件为 i7-13700HX，12 物理核/24 逻辑核，RAM 7.62 GiB；按可用 RAM 自动限定进程池，预留至少 20% RAM，每 worker 的 BLAS/OMP/MKL 线程数均为 1。35 个回放使用 2 workers，约 5.95 s；106 个独立机制任务使用 3 workers，约 29.71 s。软件版本、可用内存、磁盘与分配规则均见 [HARDWARE.json](provenance/HARDWARE.json) 及两个 `*_BATCH.json`。没有为 GPU 改写内核。

## 2. 用一句话解释当前问题

**血栓在损伤更新前已达到当前判据下的静力平衡，但程序降低部分键的完整性时，也改变了对局部变形的测量权重；位置虽然没动，重新估计的局部变形却让弹性能略微上升。**

全文统一使用 `ΔU = U_before − U_after`。正值表示储能下降；负值表示储能上升。没有取绝对值或裁掉负值。图的负半轴保留原始符号。

## 3. 当前 PD 能量是怎么组成的

粒子 i 的储能是两部分之和：

```text
U = Σ_i [ V_i r_i W(F_i) + V_i H_i ]
W(F) = μ/2 (F:F − 3) − μ log J + λ/2 (log J)²
λ = bulk − 2μ/3，J = det F
r_i = Σ_j g_ij V_j / Σ_j V_j，g_ij = 1 − D_ij
H_i = 1/2 [α μ r_i / tr(K_i)] Σ_j w_ij g_ij V_j |η_ij − F_i ξ_ij|²
```

第一部分是剩余材料的本构能量，第二部分是抵抗非仿射变形的稳定化能量。r 的平均不含 influence w，K 和 F 的加权积分含 w；二者不能混同。本构应力进入力态前也乘 r。这里不把外力势能、动能或假定的损伤耗散加入 stored energy。

当前 360 粒子开发算例的材料为 μ=1000 Pa、bulk=50000 Pa、α=3，12592 条键，60 个固定粒子。所有本轮状态保持三维满秩。上述公式按满秩状态写出；审计 API 仍记录真实 support rank，不把降维输出伪装成三维。

API 直接复用 `shape_data`、`prepare_supported_shape`、`evaluate_supported` 和 `potential`。逐粒子输出 r、rank、K/特征值/条件数、F、J、本构/稳定化能量密度及焦耳值，并提供全局三项能量。没有复制新的 mechanics solver。

## 4. fixed-x damage 时程序内部发生了什么

```text
ξ = X_j − X_i，η = x_j − x_i
K_i = Σ_j w g V_j ξ⊗ξ
A_i = Σ_j w g V_j η⊗ξ
F_i = A_i K_i^−1
```

损伤更新后，g、r、K、K 逆、trace、F、J、应力与稳定化残差都重新计算。x 固定只保证 η 固定，并不保证加权拟合 F 固定。

从当前公式可直接推得以下有限增量恒等式：

```text
F1 − F0 = { Σ_j w (g1−g0) V_j [η−F0 ξ]⊗ξ } K1^−1
```

这说明权重变化作用在旧状态的非仿射残差上；残差的方向、损伤在哪些键上以及新的 K 逆共同决定 F 的变化。它不是只看 damage 均值就能确定的量。额外使用全部 35 个保存状态做了该恒等式的代数核对，最大分量误差为 2.88584e−15，原始误差见 [RECONSTRUCTION_IDENTITY.json](comparisons/RECONSTRUCTION_IDENTITY.json)。这是审计后的数学交叉检查，没有新增或放宽通过容差。

严格仿射时 η=F_macro ξ，残差为零，任意保持满秩的权重变化都应保持 F。全键同比例缩放时，分子与 K 同步缩放，也保持 F。这两种情形与真实非仿射平衡态上的局部损伤不同。

## 5. controlled fixtures 结果

36 粒子小算例共 60 个固定坐标试验，全部不运行侵蚀或重新平衡。包括：未变形；类单轴伸长 `diag(1.01,0.9975,0.9975)`；各向同性伸长 `1.005 I`；压缩 `diag(0.99,1.0025,1.0025)`；剪切 `F_xy=0.01`；在类单轴状态上对中心附近粒子施加 `0.002h × (1,0.3,−0.2)` 的非仿射位置扰动。

每种状态检查 c=1、0.99、0.95、0.9、0.8 的全键 integrity，以及幅值 0.01 的单键、粒子星形、对称局部、空间连贯、随机稀疏损伤。

| 检查 | 结果 | 含义 |
|---|---:|---|
| 未变形态最大绝对储能 | 2.38235e−24 J | 远小于预声明 1e−20 J；无 HIGH_SEVERITY_FINDING |
| 仿射 F 相对宏观 F 的最大分量误差 | 1.44329e−15 | 正确再现仿射变形 |
| 仿射稳定化最大能量 | 3.43273e−38 J | 接近舍入零 |
| 仿射组超过数值噪声的增能情况 | 0 | 本组软化方向正确 |
| 均匀缩放最大相对 F 变化 | 1.58364e−15 | 与 F 不变的代数预测一致 |
| 正能状态最大相对 U(c)−cU(1) 误差 | 6.72429e−13 | 小于预声明 1e−7 |

必须披露一个汇总问题：第一份 `MECHANISM_SUMMARY.json` 的 `uniform_monotonic_with_roundoff` 为 false，因为汇总脚本把正能态的相对误差判据也用于零能态。三条负值均来自未变形组，分别约 −4.77e−25、−1.95e−25、−1.07e−24 J。原始数据和该 false 标记全部保留。依据开始前已声明的零态绝对容差，零态仍为零；其余 25 个正能态均匀缩放案例全部单调。解释与分类见 [UNIFORM_SCALING_INTERPRETATION.json](comparisons/UNIFORM_SCALING_INTERPRETATION.json)，没有修改任何容差、数据或测试。

非仿射扰动确实激活了重构/稳定化，但这组扰动下所有测试的总能量差仍非负。例如粒子星形损伤的软化、重构、稳定化、总差分别为 `+2.87213e−14、−2.40444e−15、+1.38379e−16、+2.64552e−14 J`。**非仿射本身不是产生总增能的充分条件**。

## 6. uniform vs heterogeneous damage

从代码先推导：全键 g 乘 c 时 K、A、r、trace 都乘 c；F 保持不变；r/trace 不变；因此本构和稳定化储能都乘 c。当前满秩正能状态下，均匀降完整性给出储能下降。

另在同一个真实基线首跳固定 x 上，比较七种损伤空间分布，各用 1e−5、1e−4、1e−3、1e−2 四个最大幅值，随机种子 20261001。实际分布用已保存的 ΔD 向量除以其最大值，只有幅值重标，不调用侵蚀更新。

| 损伤空间分布 | ΔU，幅值 1e−5（J） | ΔU，幅值 1e−2（J） | 四档符号 |
|---|---:|---:|---|
| uniform | +9.96349e−19 | +9.96345e−16 | 均下降 |
| one bond | +6.89188e−22 | +6.88031e−19 | 均下降 |
| particle star | +1.00530e−20 | +1.00474e−17 | 均下降 |
| symmetric patch | +5.56079e−20 | +5.55384e−17 | 均下降 |
| coherent local patch | +7.62737e−20 | +7.62201e−17 | 均下降 |
| random sparse | +4.97952e−20 | +4.97428e−17 | 均下降 |
| actual erosion pattern | **−6.42152e−20** | **−6.44831e−17** | **均增加** |

同最大幅值的各分布，其损伤总量不同。为避免把数量差误当空间效应，每个案例还保存了“相同 ΣΔD 的均匀损伤”对照；这些对照的 ΔU 全部为正。它们见各 `fixtures/pattern_study/*/SUMMARY.json`。

结论是：当前异常与**实际损伤权重和非仿射局部变形的配合**有关；不能写成“只要不均匀损伤就一定增能”。[图：PNG](figures/uniform_vs_local_damage.png) / [PDF](figures/uniform_vs_local_damage.pdf)。

## 7. actual 35 damage jumps replay

仅读取 `results/streaming_erosion_equilibrium_v1/{case}/damage_jump_NNN/state.npz` 的 x、D_before、D_after，与各旧 AUDIT.json 对照。逐跳源文件 SHA-256、差值、容差单独保存。预声明容差为 `1e−28 J + 1e−12 × max(|U_before|,|U_after|)`，超过即停止后续审计。

**35/35 精确复现；新旧总能量差最大差异为 0 J，before/after 同样一致。** 本次没有用重算轨迹替换历史。三档仍分别为 5/5、10/10、20/20 增能。

主表：[damage_jump_audit.csv](actual_run_replay/damage_jump_audit.csv) / [JSON](actual_run_replay/damage_jump_audit.json)。每行包含 case、exposure time、macro step、before/after、本构/稳定化分项、精确分解、残差、最大 K/F 变化和抽样负 Y 数。

每跳另外检查 ΔD 最大的 10 条键，共 350 次“键-状态”样本，在所用差分判据下全部得到稳定负 Y。它们可能跨跳重复，不是 350 条不同的键，也不是全模型普查。

## 8. energy decomposition

对每个粒子逐项计算：

```text
softening      = V (r0−r1) W0
reconstruction = V r1 (W0−W1)
stabilization  = Ustab0 − Ustab1
total          = softening + reconstruction + stabilization
```

35 次真实跳变的最大全局闭合残差为 **5.43582e−29 J**。若连同受控/模式/α 试验的全部 135 个已保存前后状态一起计，最大全局残差为 **5.11414e−26 J**，来自各向同性伸长的 c=0.99 案例（U_before≈3.93557e−10 J），仍在原预声明容差内。逐粒子和全局残差都已保存；这是代数闭合证据，不是物理精度证明。

以下每一行对应独立的 0.1 s 暴露慢步长试验，单位均为 J：

| 慢步长 | 软化 | 重构 | 稳定化 | 总差 |
|---|---:|---:|---:|---:|
| 0.02 s，5 跳 | +2.885594e−17 | −5.215009e−17 | +2.276637e−18 | −2.101751e−17 |
| **0.01 s，10 跳，基线** | **+2.886065e−17** | **−5.214787e−17** | **+2.276748e−18** | **−2.101047e−17** |
| 0.005 s，20 跳 | +2.886300e−17 | −5.214670e−17 | +2.276803e−18 | −2.100690e−17 |

把三档 35 行作审计账目求和：软化 `+8.657960e−17`，重构 `−1.564447e−16`，稳定化 `+6.830187e−18`，总差 `−6.303488e−17 J`。**此合计不是一次处理过程的总耗散或一条更长轨迹**。

基线首跳的分解最直观：`+2.88184e−18 −5.21436e−18 +2.27543e−19 = −2.10498e−18 J`。材料剩余量减少本来降能，但 W(F) 上升占上风。[图：PNG](figures/energy_decomposition_actual_jumps.png) / [PDF](figures/energy_decomposition_actual_jumps.pdf)。

五种反事实已在首跳与末跳单独保存，全部标明 **DIAGNOSTIC ONLY / NOT A CANDIDATE MODEL / NOT USED IN PRODUCTION**。首跳为：

| 只读反事实 | ΔU（J） | 边界 |
|---|---:|---|
| 实际 formulation | −2.104981e−18 | 真实历史审计 |
| 固定旧 F/W，仅减本构 r；稳定化保持旧值 | +2.881840e−18 | 隔离纯软化；不是自洽的新力学状态 |
| 旧 K 逆与 trace，新 A 与 r | −3.261289e−13 | 失去仿射再现，不能当修复方案 |
| 只看本构 | −2.332524e−18 | 移除稳定化记账后仍增能 |
| 只看稳定化 | +2.275431e−19 | 稳定化本身下降 |

冻结 K 并不等于冻结 F。此反事实的巨大增能恰好表明：保留旧算子而单独削减 numerator 会引入完全不同的运动学误差，不能直接推广为 production model。

## 9. K/F reconstruction diagnostics

全部 35 跳中，最大相对 F 变化为 `8.81199e−9`，最大相对 K 变化为 `1.20144e−4`；最低 support rank 始终为 3，最大 K 条件数约 **2.68876**。没有观察到 rank loss 或病态 K，但小幅 F 变化依然足以影响当前很小的弹性储能。

以基线首跳的 360 个粒子计算 Spearman 相关系数：

| 特征 | 与 signed total ΔU | 与 reconstruction ΔU |
|---|---:|---:|
| 相对 F 变化 | −0.7664 | −0.8464 |
| 相对 K 变化 | −0.1429 | −0.5106 |
| K 条件数变化 | −0.2439 | −0.3863 |
| K 各向异性变化 | −0.3475 | −0.4995 |

F 变化与增能的关联较强，条件数变化关联更弱。因此“损伤改变重构”有直接代数/分解依据，而“K 病态导致错误”没有被数据支持。相关性本身不证明边界、source 或条件数的因果作用。

首跳 total 增能最大的粒子是 **93、87、123**。93 号位于参考坐标 `(-0.3125,0.0625,-0.6625) mm`，属于 SOURCE_NEAR；其 total 差 `−4.03965e−20 J`，重构差 `−5.71539e−20 J`，K 条件数仅由约 1 增至 1.0000585。前 20 名及坐标/分解见 [TOP_ENERGY_INCREASE_PARTICLES.json](comparisons/TOP_ENERGY_INCREASE_PARTICLES.json)。

[F 图 PNG](figures/F_reconstruction_vs_energy.png) / [PDF](figures/F_reconstruction_vs_energy.pdf)；[K 图 PNG](figures/K_condition_vs_energy.png) / [PDF](figures/K_condition_vs_energy.pdf)。所有旧/新 K 特征值、F、J、r 及逐粒子分量在各 `particle_audit.npz` 中，标量表在 `particles.csv`。

## 10. thermodynamic force Y diagnostics

定义诊断量 `Y_k = −∂U/∂D_k |_x`。它的单位是 J / unit D，不是一个空间力矢量，也没有被用作新的损伤律。负 Y 表示在这个 fixed-x 状态增加该键损伤会增加储能。

每条键使用 ε=1e−6、3e−6、1e−5、3e−5、1e−4；D 两侧可用时用中央差分，边界用前向/后向。先逐粒子做能量差再求和，避免无关粒子的能量抵消。预声明规则是从小到大寻找首个稳定的连续三档，极差不超过 `1e−22 + 0.03×|median|`；符号仅在中位值绝对值大于“极差+绝对底值”时认定。该规则是数值稳定诊断，不是严格误差上界。

| 状态与抽样范围 | 样本数 | 符号已分辨 | 未分辨 | 稳定负 Y | 占全部样本 |
|---|---:|---:|---:|---:|---:|
| 小仿射，全 550 键 | 550 | 550 | 0 | 0 | 0% |
| 小平衡态，全 550 键 | 550 | 513 | 37 | 123 | 22.36% |
| 基线首跳，分层并集 | 951 | 906 | 45 | 300 | **31.55%** |
| 基线末跳，分层并集 | 957 | 918 | 39 | 302 | **31.56%** |

基线并集含 ΔD 最大 top100、全部接触前 5 个增能热点粒子的键、源附近 100、固定层附近 100、确定性随机 100；去重后得到以上样本量。它偏重热点，因此**这些比例不能推断全部 12592 键的总体比例**。若只在已分辨样本中计数，首跳是 300/906=33.11%，不是表中的 31.55%。未分辨样本没有被归为零或正值。

已分辨 Y 的 min/median/max（J / unit D）：小仿射 `7.98128e−14 / 9.28111e−14 / 1.07441e−13`；小平衡态 `−6.33094e−16 / 8.41304e−18 / 1.68083e−15`；基线首跳 `−1.14300e−16 / 7.46106e−18 / 1.75261e−16`；末跳 `−1.14440e−16 / 7.32725e−18 / 1.75127e−16`。未收敛者保留 `THERMODYNAMIC_FORCE_NUMERICALLY_UNRESOLVED`。

首跳 Y 与实际 ΔD 的 Spearman 相关为 −0.5011，与 F 对损伤的敏感性为 −0.2647。最负的键为 4960（连接 58、93 号粒子），中点 `(-0.375,0.0625,-0.600) mm`，Y≈−1.14300e−16，而该次实际 ΔD=0。这说明“潜在增能敏感键”不等于“该跳实际损伤最多的键”。[前 100 个稳定负 Y 键](bond_sensitivity/baseline_first/TOP_NEGATIVE_BONDS.json) 保存位置、不确定性、实际 ΔD 与 F 敏感性。

[Y 分布 PNG](figures/thermodynamic_force_distribution.png) / [PDF](figures/thermodynamic_force_distribution.pdf)；额外提供 [ε 稳定性 PNG](figures/finite_difference_stability.png) / [PDF](figures/finite_difference_stability.pdf)，展示稳定正、稳定负及未分辨的代表曲线，不作平滑。

## 11. stabilization contribution

35/35 跳的稳定化能量差为正，表示该部分能量下降。基线累计下降 `2.27675e−18 J`，抵消约 **9.78%** 的本构增能量；它是次要的抵消项，不能称为增能主因。只看本构仍累计增能 `2.32872e−17 J`。

首跳固定 x、旧/新 D 不变，只在审计副本改变 α：

| α | 与当前值比例 | total ΔU（J） |
|---:|---:|---:|
| 0 | 0× | −2.332524e−18 |
| 1.5 | 0.5× | −2.218752e−18 |
| 3 | 1× | −2.104981e−18 |
| 6 | 2× | −1.877438e−18 |

四档都增能；较大的 α 在此状态抵消更多增能。不能据此选一个“修复值”。这不是不同 α 下重新求平衡的比较。

仿射组稳定化接近零，对 α 不敏感；非仿射扰动组稳定化被激活，但其 total ΔU 在四档均为正。因此本轮证据不支持“非仿射 fixture 中总增能问题必然更强”。[图 PNG](figures/stabilization_sensitivity.png) / [PDF](figures/stabilization_sensitivity.pdf)，所有输出醒目标注 DIAGNOSTIC ONLY。

## 12. surface/fixed-region association

按已有几何定义：固定层附近为 z 不高于最上固定粒子中心加 h；source 附近为距配置中的源不超过已有 influence radius；free surface 复用已有 weighted-support 判据且排除固定粒子；interior 为上述区域之外。源仍为 `(-0.825,0,-0.725) mm`，位于上游左侧，不移动到 clot 顶部。

粒子互斥分类优先级为 `FIXED_BASE_NEAR > SOURCE_NEAR > FREE_SURFACE > INTERIOR > OTHER`；原始重叠布尔掩码也保存。首跳统计：

| 互斥粒子类别 | 数量 | 平均 reconstruction ΔU（J） | 平均 ||ΔF|| | 平均条件数变化 | total ΔU 求和（J） |
|---|---:|---:|---:|---:|---:|
| INTERIOR | 80 | −2.61554e−20 | 3.93240e−9 | 2.96286e−5 | −1.36581e−18 |
| FREE_SURFACE | 114 | −1.14073e−20 | 2.5835e−9 | 5.9497e−5 | −3.76401e−19 |
| FIXED_BASE_NEAR | 120 | −6.26823e−21 | 1.96325e−9 | 8.125e−6 | **+1.46891e−19** |
| SOURCE_NEAR | 46 | −2.32459e−20 | 4.18604e−9 | 5.9946e−5 | −5.09661e−19 |
| OTHER | 0 | 不适用 | 不适用 | 不适用 | 0 |

键的 Y 使用端点区域掩码：interior 要求两个端点都在内区，其他类别任一端点满足即可，因此类别可重叠。首跳各类稳定负 Y / 该类全部抽样键：interior 32/107=29.91%；surface 179/436=41.06%；fixed-near 74/332=22.29%；source-near 172/560=30.71%。这些是抽样相关，不可把重叠数相加，也不可与互斥粒子均值当作同一总体。

内部与源附近的平均重构增能较大；固定层附近的 total 求和反而为正。表面相关性值得后续检查，但**不能把异常归因为固定底层，也不能只归因为表面支撑截断**。见 [粒子表](comparisons/PARTICLE_CATEGORIES.csv)、[键表](comparisons/BOND_CATEGORIES.csv)、[相关系数](comparisons/CORRELATIONS.json)。

[空间图 PNG](figures/spatial_energy_anomaly_map.png) / [PDF](figures/spatial_energy_anomaly_map.pdf) 用真实保存 x、1× 位移比例，显示上游源星标和固定粒子空心圈，无箭头和微泡球体。所有 8 张图均为黑底、Arial 常规字重、英文标注，每张同时导出 PNG/PDF；字体来源和数据路径见 [FIGURE_MANIFEST.json](figures/FIGURE_MANIFEST.json)。

## 13. open-source / literature comparison

本地机制审计完成后才调查原文与源码。完整证据、公式位置、访问限制及固定行号见 [LITERATURE_REVIEW.md](literature/LITERATURE_REVIEW.md)。

Peridigm 固定 commit `f25c4928e850b2aa0012bc5d3b8b1ec5c78c8903`、BSD-3-Clause。经典 correspondence 的 K/F 函数存在零值 bondDamage 占位，不能把它当成活动损伤重构示范；bond-associated material 则在损伤增长时更新含 damage 的重构权重。[经典函数](https://github.com/peridigm/peridigm/blob/f25c4928e850b2aa0012bc5d3b8b1ec5c78c8903/src/materials/correspondence.cxx#L810-L920)、[bond-associated 调用](https://github.com/peridigm/peridigm/blob/f25c4928e850b2aa0012bc5d3b8b1ec5c78c8903/src/materials/Peridigm_BondAssociatedCorrespondenceMaterial.cpp#L1066-L1110)。没有在所查实现中找到对本项目损伤/能量形式的固定 x 单调性保证；也没有运行 Peridigm 对照算例，不能宣称两者数值等价。

[Silling et al. 2007](https://doi.org/10.1007/s10659-007-9125-1) 的原始 correspondence 公式解释 K/F/force 的构造。已读取作者上传全文的相关章节；出版社/UNL PDF 下载受限的记录保留。稳定化对照读取 [Silling 的完整 SAND2016-6710 报告](https://www.osti.gov/servlets/purl/1457611)，并核对 [2017 期刊文章身份](https://doi.org/10.1016/j.cma.2017.03.043)。该报告分析对形变扰动的稳定性，与本轮对 damage 方向的能量单调性是不同问题；不据此给本项目背书。

[Rausch 2021 官方全文](https://par.nsf.gov/servlets/purl/10317555) 的式 (11)–(14) 明确区分运动学、不可逆局部变量、非局部场及总自由能。本轮没有移植其材料参数或模型；仅用它说明非负耗散必须依据完整的自由能和演化关系检查。其非局部惩罚能也提醒我们，任意手工更新内部变量并不自动满足总能量单调性。详细的偏导解释明确标为本审计推导，不伪装成论文原文演化方程。

当前证据不足以把“出现重构增能”定性为代码 bug，也不足以说换用成熟开源代码就一定消除它。

## 14. 哪些解释被排除

| 解释 | 本轮判断及范围 |
|---|---|
| 只是 startup transient / 尚未达到机械平衡 | 不能解释这些保存的已平衡状态上的固定坐标跳变 |
| 只是旧能量日志记录错了 | 35/35 的 before/after 和 signed ΔU 精确回放，不支持 |
| 分解漏了一项导致假增能 | 三项闭合到 5.44e−29 J，不支持 |
| stabilization 是主要增能源 | 全 35 次该项下降；去掉记账后本构仍增能，不支持 |
| uniform degradation 同样产生有意义的增能 | 正能态符合解析缩放；未变形零态的负舍入已单独披露 |
| K 已奇异/降秩 | 全部支持 rank=3，最大条件数约 2.69，不支持 |
| 任意 local damage 或任意非仿射状态都会增能 | 多个直接反例；不支持 |
| fixed base 是已证实原因 | 固定附近并非最严重；没有隔离边界效应的因果实验，不支持 |

未排除的范围包括：变形场与局部损伤耦合的更一般规律、离散化/支撑布局依赖、未分辨的 Y 数值精度、网格客观性及真实血栓物理。全局闭合、测试通过和图的平滑都不证明局部物理准确。

## 15. 当前最可能机制

**SUPPORTED：** 在本轮全部 35 个真实 fixed-x jumps，损伤改变加权重构 F，使 W(F) 上升；其负 reconstruction contribution 超过正 softening 和正 stabilization contribution，精确解释总增能。实际 damage 空间分布在四档幅值下都复现该方向，且有稳定负 Y。`ENERGY INCREASE SOURCE: IDENTIFIED FOR AUDITED STATES`。

**PARTIALLY SUPPORTED：** 非仿射残差与损伤空间方向的配合是重要条件；source/interior 热点、表面截断和 K 几何变化与异常存在不同程度关联。代数恒等式支持权重/残差耦合，但本轮没有通过独立边界、负载或网格实验证明这些位置因素各自的因果份额。

**NOT SUPPORTED：** 单纯 r 软化导致增能、稳定化主导增能、K 奇异、固定层主导、任意非均匀损伤都增能、全局代码实现错误已被证明、能量问题已被修复。

## 16. 当前是否能把 D 叫 thermodynamic damage

当前 D 可以继续作为本 prototype 的损伤/侵蚀状态变量，但**尚不能称为已建立热力学耗散保证的 damage variable**。诊断 Y 有稳定负值，当前演化也没有从完整自由能和耗散约束推导；本轮不直接使用“违反第二定律”表述。

```text
DAMAGE ENERGY AUDIT: COMPLETED
ENERGY INCREASE SOURCE: IDENTIFIED (audited fixed-x states)
CURRENT DAMAGE VARIABLE: THERMODYNAMIC DISSIPATION NOT ESTABLISHED
MECHANICAL EQUILIBRIUM: REMAINS VERIFIED FOR CURRENT DEVELOPMENT CASE
EROSION SLOW-STEP CONSISTENCY: REMAINS PASSED FOR CURRENT DEVELOPMENT CASE
PHYSICAL CALIBRATION: NOT PERFORMED
MESH OBJECTIVITY: STILL OPEN
WALL ATTACHMENT: STILL OPEN
FAILURE / DETACHMENT: NOT OBSERVED IN THE HISTORICAL CASES
REAL MICROBUBBLE CFD: NOT USED
PRODUCTION DAMAGE LAW: UNCHANGED
```

审计完成只指任务证据齐备，不代表所有导数已分辨、所有 formulation 已查遍或物理已校准。保留先前有限范围的软件/平衡结论，也保留其原有科学限制。

## 17. 下一轮候选方案

仅提供以下 **两类待评审候选**，本轮均未实现：

1. **从完整离散自由能与不可逆演化共同定义损伤。** 借鉴 Rausch 等的能量、驱动力、非局部项分开定义的框架，先明确本项目运动学与 damage 的依赖，再检查共轭力及耗散约束。候选门槛是能量/力一致性、仿射再现、不可逆性和当前 35 状态的 signed energy/Y 检查；不能只改一个退化函数或移植血栓参数。
2. **对照已有 bond-associated correspondence 的重构方案。** Peridigm 提供固定版本的可审查实现，可作为独立小算例研究对象；其抗零能模态或权重更新方式不等于非负损伤耗散保证，仍需相同审计。冻结旧 K 的本轮反事实不是此候选，也不能转入生产。

本轮在交付这些证据、图和报告后停止。不自动进入模型修复，不提交/推送新分支，不恢复旧批次。
