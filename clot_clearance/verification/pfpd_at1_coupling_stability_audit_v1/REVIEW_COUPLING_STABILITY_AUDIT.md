# PFPD_AT1_COUPLING_STABILITY_AUDIT_V1 人工审核报告

本轮完成的是原因诊断。当前证据支持 **COUPLED_BRANCH_LOSS_OR_LOCALIZATION_INSTABILITY**：失败载荷附近找到的一阶平衡状态具有可行负曲率；最后一个已接受状态也已有负曲率。未找到失败载荷下通过曲率检查的局部极小值。实际迭代确实有非光滑分支切换，约 1.93 只能解释为与观测方向吻合的**条件冻结分支**特征值，不能称为实际映射的普通谱半径。这些结论限于当前小立方体候选能量，不证明真实血栓裂纹失稳。

工作区：`/home/lzy/projects/temp_storage/clot_clearance_github_sync_20260930/clot_clearance`。Git 根目录是其父目录；分支 `clot-clearance/energy-regularization-20260930`，HEAD `f0e891686b9b107510ab3292aa27ce83638b3339`。HEAD 之外原有的大量未跟踪文件也纳入保护。适用祖先与工作目录中未发现 AGENTS.md。原文件 SHA256、链接和最终核对见 [PRE_RUN_MANIFEST](provenance/PRE_RUN_MANIFEST.json)、[POST_RUN_MANIFEST](provenance/POST_RUN_MANIFEST.json)、[PRESERVATION_CHECK](provenance/PRESERVATION_CHECK.json)：保护 15,710 个普通文件、241 个符号链接；原有文件修改 0、缺失 0、链接改动 0、production 修改 0、历史结果修改 0。

## 1. 上一轮到底哪里失败

固定位置时的 AT1 标量更新、启动阈值、不可逆性和能量最小化已验证；失败发生在力学与相场交替求解之后。原标准加载—卸载—重载路径停在第 010 步，工程应变 0.014；上一接受状态是第 009 步、应变 0.013。加载后尚未走到卸载和重载。增加交替迭代次数的另一条原路径更早停在第 008 步、应变 0.012。不能用前面步骤通过代替整条路径通过。

原单调拉伸路径接受到应变 0.04，最大相场约 0.939325；下一应变 0.05 的提议约 0.961546 触发原有 `sc=0.95` 高损伤门限。这与本轮两条路径的纯拉伸域退出是不同停止原因。

所有分析先读取原 `accepted_state.npz`、逐次机械/相场状态和失败日志，随后才执行诊断调用；没有先重新生成失败数据。完整身份、载荷、位移、上一接受相场来源、迭代号、残差、能量和 Green 应变见 [FAILURE_STATE_INDEX](state_index/FAILURE_STATE_INDEX.json)。

## 2. 本轮区分哪两类问题

可以把能量想象成地形。第一类情况是谷底仍在，求解算法却找不到；第二类是原来走的路线变成了山脊，附近一个方向向下，另一个方向向上，停在原地也不稳定。仅看“力接近零”无法区分谷底与山脊。

还有第三种困难：不可逆更新像一个带折角的开关，不同方向的小扰动会走不同公式。本轮同时检查这个折角；不因为一条方向曲线很直，就宣称所有方向都有同一个普通导数。

## 3. Fixed-point map 是什么

新增只读 API 位于 [coupling_stability_audit.py](../../pd_clot/pfpd_candidate/coupling_stability_audit.py)。`Phi(s_input)` 依次复用原 ADR 力学平衡、原 BAQP 驱动能 `Y`、原 AT1 更新，始终 `h=1`。每次返回位形、力残差、驱动能、更新分支、能量、最小 Green 应变和 ADR 次数。位移边界保持给定值。

原交替算法的内循环使用 `s_output=max(s_input,s_trial)`。诊断 Phi 必须保留这个**随输入移动的内循环下界**，才能复现原行为。联合能量问题则使用**上一接受加载步的相场**作为物理不可逆下界，两者不可混用。

复现样本：标准失败步的迭代 0、1、6、12、23，以及 extra 路径的迭代 0、6、20、26。9/9 通过；最大相场差和最大位形差均为 **0**。预声明容差分别为 `5e-10` 与 `1e-12 m`，未事后调整。证据：[REPRODUCTION](fixed_point/REPRODUCTION.json)。

## 4. 实际 amplification 是否稳定

对保存的相场更新向量计算 L2、最大范数的相邻比值，以及相邻更新的带符号余弦。第 0 次是加载后的初始校正，不属于之后的小扰动放大阶段；相邻比值没有人为取绝对方向。

| 原路径 | 中段 L2 放大率 | 中段最大范数放大率 | 相邻更新余弦 |
|---|---:|---:|---:|
| 标准失败步 | 1.929487–1.929768 | 1.929389–1.929767 | 0.9999999504–约 1 |
| extra 失败步 | 1.933160–1.933356 | 1.933091–1.933361 | 0.9999853570–约 1 |

这里“中段”采用相场最大增量 `1e-9 < max Δs < 1e-4` 的数据范围；完整逐次结果保留在 [SAVED_PATHS](energy_path/SAVED_PATHS.json)。正余弦趋近 1，表示逐次增长沿近似同一方向，并非正负振荡。稳定对照 A 的两个相场更新均为零，其 `0/0` 比值记为未定义。

## 5. Spectral radius 结果

27 粒子、185 根键、27 个自由位移坐标。A/B/C 分别取应变 0.009 的接受对照、应变 0.013 的最后接受状态、应变 0.014 的第 6 次失败迭代输入。完整 Jacobian 差分测试共 `3 × 3 × 185 = 1665` 列。相对 epsilon 为 `1e-6 / 3e-6 / 1e-5`，实际每列步长乘 `max(0.1,|s_b|)`。边界附近使用二阶单侧差分，其余使用中心差分。

| 状态 | 普通 Phi 谱半径 | 冻结所选分支主特征值（中心 epsilon） | 切换列数 | 分支临界键数 |
|---|---|---:|---:|---:|
| A | 1.000000；边界可行单侧映射，边际 | 1 + 0i | 0/185 | 0 |
| B | **UNRESOLVED** | 1.9315586149 + 0i，条件结果 | 185/185 | 185 |
| C | **UNRESOLVED** | 1.9297683120 + 0i，条件结果 | 185/185 | 35 |

三档 epsilon 的条件谱半径相对跨度：A 为 0，B 为 `3.73e-9`，C 为 `3.51e-8`，低于预声明 2%。数值重复性良好并不能消除分支切换。切换产生的原始割线矩阵另存为 `raw_secant_NOT_ordinary_J`，没有拿它计算并宣称普通 Jacobian 谱。

A 的所有相场变量都是物理下界强活跃变量；小的正相场扰动不能愈合，Phi 在这些方向保留输入，因此给出单位特征值。它不是 `rho<1` 的收缩例子，但在固定物理下界下的耦合能量仍可局部稳定。B/C 的 185 根键相对上一加载步下界均在内部，均不靠近 `Y=Yc` 或 `sc`。这里检测到的主要折角属于 `max(s_input,s_trial)` 的**内循环选择器**，不能混称为损伤启动阈值切换。

原始矩阵、特征值实部/虚部/模和特征向量均在 [spectrum](spectrum/)，综合判断见 [EPSILON_STABILITY](spectrum/EPSILON_STABILITY.json)。普通谱的总体状态保持 `UNRESOLVED`。

## 6. Dominant mode 是否解释实际发散

标准路径第 4–20 次迭代与 C 条件主模态的最小余弦为 `0.9999997123`；L2 放大率与条件主特征值的最大相对差为 `1.4582e-4`，即约 **0.0146%**，均满足预声明关联阈值。可表述为 `OBSERVED_DIVERGENCE_IS_CONSISTENT_WITH_DOMINANT_FIXED_POINT_MODE`，必须同时注明“条件冻结分支”。

额外沿正、反模态方向各做三档扰动：沿实际增长方向，线性预测的相对 JVP 误差约 `6.37e-7–6.36e-6`；反方向误差约 `1.09–1.11`。两边均保留了发生切换的记录。这直接说明同一个普通线性算子不能描述两侧实际映射。数据见 [MODE_ALIGNMENT](spectrum/MODE_ALIGNMENT.json)。

只计算了理论松弛表，没有运行松弛算法。对 C 的正实条件特征值，omega 为 `0.1 / 0.25 / 0.5 / 0.75 / 1` 时，条件 `lambda_relaxed` 分别约 `1.092977 / 1.232442 / 1.464884 / 1.697326 / 1.929768`，仍均大于 1。该条件结果不支持把简单减小松弛因子作为修复建议。见 [RELAXATION_THEORY_ONLY](spectrum/RELAXATION_THEORY_ONLY.json)。

## 7. Monolithic minimization 是否找到稳定解

候选目标完全复用原能量：`Pi = U_BAQP + Σ_b weight_b × 2Yc × s_b`。固定边界不参与优化；自由位移与相场联合优化。失败载荷为应变 0.014，所有起点使用加载步 009 的接受相场作为下界，上界仍为 `sc=0.95`，最终必须严格小于 sc。

位移变量以 `1.25e-6 m` 缩放，能量以 `5.3263037109375e-12 J` 缩放。四个初值分别是：上一接受状态按新边界作仿射抬升；原失败步第 0 次相场更新后；第 23 次相场更新后；上一接受状态抬升后加固定随机种子的微小正相场扰动。仿射抬升仅是初值选择，不约束自由节点继续作仿射运动。

梯度采用 `grad_x Pi=-force_free` 和 `grad_s Pi=weight×[-2(1-s)Y+2Yc]`，没有另造权重或双计数约定。优化前在 A/B/失败加载初态各检查 12 个确定性随机方向、3 档 epsilon，共 108 项。最大归一化误差分别为 `1.17e-6 / 4.55e-5 / 3.97e-7`，均低于预声明 `5e-5`；B 接近这个门槛，未隐去。归一化分母含预声明梯度尺度下限 `1e-4`，避免近驻点的相对误差被零除放大。差分方向允许探测能量的代数延拓，越域探针从未作为有效物理状态。见 [GRADIENT_GATE](monolithic/GRADIENT_GATE.json)。

| 同一失败载荷的起点 | L-BFGS-B | SLSQP |
|---|---|---|
| 上一接受状态仿射抬升 | 试探点越纯拉伸域，停止 | 试探点越纯拉伸域，停止 |
| 失败步第 0 次相场更新后 | success；一阶通过；负曲率 | success；一阶通过；负曲率 |
| 失败步第 23 次相场更新后 | 试探点越纯拉伸域，停止 | 试探点越纯拉伸域，停止 |
| 抬升后微小确定性相场扰动 | 试探点越纯拉伸域，停止 | 试探点越纯拉伸域，停止 |

因此 8 次尝试中，6 次退出当前模型域，2 次停在鞍形驻点，**0 个通过局部极小值审核**。SLSQP 从本已近乎驻定的第 0 次状态出发，仅一次函数求值就返回；L-BFGS-B 是 3 次迭代、33 次求值。不能把二者停在相近位置宣传成独立找到两个稳定解。两个控制计算中，A 同时通过一阶与曲率检查，B 一阶通过但曲率为负。

越域尝试按预声明规则在首个越域试探点终止，未给 SciPy 回溯到域内继续尝试。这是本次参考的严格边界，也限制了搜索能力：**没有找到有效极小值，不等于证明附近或远处一切极小值都不存在**。保存的最后域内点是最后一个函数求值，未必是优化器接受迭代点，更不是解。详见 [ACCEPTANCE_AND_COMPARISON](monolithic/ACCEPTANCE_AND_COMPARISON.json)。优化器接口依照 [SciPy L-BFGS-B](https://docs.scipy.org/doc/scipy/reference/optimize.minimize-lbfgsb.html) 与 [SciPy SLSQP](https://docs.scipy.org/doc/scipy/reference/optimize.minimize-slsqp.html)；所有选项已预声明。

## 8. KKT 是否满足

KKT 是“在允许方向上是否还存在一阶下降驱动力”的检查；它是局部最小值的必要条件，不能替代曲率检查。本轮独立计算 `r_s=s-clip(s-q_scaled,lower,sc)`，其中 `q_scaled=[-2(1-s)Y+2Yc]/max(Y,Yc)`，报告其最大绝对值。未依赖 SLSQP 未返回的边界乘子。原力残差仍按原力尺度归一化，门槛 `1e-9`；相场投影 KKT 门槛 `1e-8`。

| 状态 | 归一化力残差 | 投影相场 KKT | 最大 s | 最小 Green 应变 |
|---|---:|---:|---:|---:|
| 最后接受 B | 9.28e-12 | 6.38e-11 | 0.4100463935 | 约 -3.3e-18 |
| 早期迭代 C 输入 | 6.75e-16 | 2.79e-8 | 0.4918205429 | 在原域内 |
| L-BFGS-B 第 0 起点驻点 | 3.45e-11 | 2.08e-11 | 0.4918205358 | -3.70e-18 |
| SLSQP 第 0 起点驻点 | 3.61e-11 | 1.67e-16 | 0.4918205358 | -3.60e-18 |

C 输入不是联合驻点，因为相场 KKT 未通过；没有把它当作 monolithic 解。两个优化器驻点相对物理下界的最小余量约 `0.08177414`，距 sc 约 `0.45817946`，不可逆和高损伤条件均满足。它们仍被负曲率否决。

## 9. Hessian / Schur complement 说明什么

Hessian 表示局部地形弯曲。单独固定损伤求位移、固定位置求损伤，都可能像谷底；一起改变位移和损伤时仍可能出现下坡方向。本轮形成完整 212×212 梯度差分矩阵，先量测对称性误差，再对称化用于特征值计算。三档 epsilon 为 `1e-5 / 3e-6 / 1e-6`，保留全部结果。

相场对角块与 `Hss=diag(2Y×weight/E0)` 一致。A/B/C 的中心档原始对称性相对误差约 `2.74e-9 / 4.12e-9 / 4.04e-9`；Hss 解析误差约 `6.50e-11 / 8.74e-12 / 9.22e-12`，低于预声明容差。A 排除 185 个物理强下界变量，仅检查 27 个自由位移方向；B/C 检查全部 212 个方向。没有未处理的弱下界自由度。

| 状态 | 最小 reduced 特征值 | 最小 Hxx 特征值 | Hxx 条件数 | 最小 Schur 特征值 |
|---|---:|---:|---:|---:|
| A，0.009 | +0.0628263463 | +0.0628263463 | 3.4481 | 无内部相场自由度 |
| B，0.013 | -0.0434294677 | +0.0220185209 | 3.4282 | -0.1049929890 |
| C，0.014 | -0.0407522531 | +0.0163658247 | 3.4232 | -0.1217320366 |
| L-BFGS-B 驻点，0.014 | -0.0407522534 | +0.0163658249 | 3.4232 | -0.1217320367 |
| SLSQP 驻点，0.014 | -0.0407522534 | +0.0163658249 | 3.4232 | -0.1217320367 |

表中数值对应预声明缩放后的无量纲坐标和能量；数值大小不应解释为 SI 刚度。可逆坐标缩放不改变负曲率是否存在。三档 epsilon 的符号与数量级一致，远离 `1e-7×max(||H||₂,1)` 曲率分类容差。

Hxx 条件数远低于 `1e10` 限制且为正定，故用线性 `solve` 形成 `S=Hss-Hsx solve(Hxx,Hxs)`，没有显式硬求逆。Schur 为负而 Hss、Hxx 单独为正，证据指向耦合方向失稳。[SciPy eigh](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.eigh.html) 用于对称特征分解。另以 Hessian 隐函数关系重建平滑 trial 分支导数，与 Phi 差分 trial 导数的相对差 B/C 为 `1.64e-7 / 1.10e-7`，形成独立交叉核对；没有将该平滑 trial 导数冒充非光滑 Phi 导数。

负模态用 `1e-8 / 1e-7 / 1e-6` 幅度检查正反方向的 bounds 和原纯拉伸域，均找到双侧可行扰动。其能量差可小到约 `1e-26 J`，接近双精度抵消尺度，不作为单独决定性证据；主要依据是梯度一致性、矩阵块一致性、epsilon 稳定性和可行性共同支持的负曲率。稳定 A 最大模态试探幅度可能越域，更小幅度仍可行，未放宽门限。

空间图采用 B 的同一个负 Hessian 模态，左为按原键权重聚合到粒子的相场分量，右为自由节点位移幅值，分别归一化。相场呈上下相反的分量，自由中间层承担位移，提示当前夹持几何相关的对称破缺方向；不据此定位真实裂纹。原始向量与块矩阵见 [hessian](hessian/)。

## 10. Load-step subdivision 说明什么

从完全相同的已接受应变 0.013 状态出发，目标都为 0.014，仅将这一段分成 2、4、8 步。Gc、sc、能量、ADR、纯拉伸域及交替收敛容差均不变，未覆盖原 load schedule。

| 步数 | 首次尝试应变 | 已接受子步数 | 到达 0.014 | 停止原因 |
|---|---:|---:|---|---|
| 原 1 步 | 0.014 | 0 | 否 | 纯拉伸域越界 |
| 2 步 | 0.0135 | 0 | 否 | 纯拉伸域越界 |
| 4 步 | 0.01325 | 0 | 否 | 纯拉伸域越界 |
| 8 步 | 0.013125 | 0 | 否 | 纯拉伸域越界 |

各路径的终点位形、相场、最大 s、能量以及与联合参考的终点距离均不存在，因此 JSON 写 null，图中写 N/A。不同载荷的最后内循环状态保留作诊断，但没有冒充相同终点结果。联合参考本身为鞍形驻点，亦不能作为稳定真值。此次细分没有改变“首段即失败”的诊断行为，不能归因于只要减小加载步就能解决。见 [ENDPOINT_COMPARISON](continuation/ENDPOINT_COMPARISON.json)，状态始终 `STEP_SUBDIVISION_DIAGNOSTIC_ONLY`。

## 11. Pure-tension gate 是原因还是伴随现象

标准路径清楚的放大从迭代 4 已出现，首次越域在迭代 24；extra 路径分别为迭代 6 和 27。清楚放大判据是 L2 比值大于 `1.02` 且最大相场增量大于 `1e-9`。首次退出时最小 Green 应变分别为 `-1.115206881e-9` 与 `-1.773933946e-9`，原门槛仍为 `-1e-9`。

因此时间顺序支持“域越界是已经开始的放大之后出现的伴随现象”，不是已观察到的最初事件。这不是严格因果证明。每次相场后残差与机械平衡后残差分别记录：机械刚平衡时力可很小，随后改变相场仍会重新产生不平衡，不能混用两种残差。见 [TENSION_GATE_ORDERING](energy_path/TENSION_GATE_ORDERING.json)。

## 12. Energy path 说明什么

上一接受状态在其原载荷的总能量为 `7.510335337926564e-12 J`。在同一失败目标载荷下，上一状态仿射抬升的能量为 `8.015976794611827e-12 J`；失败步第一次相场更新后为 `7.945889255577194e-12 J`，第 23 次后为 `7.945888117077936e-12 J`。图 4 同时提供相对旧接受态的全局标尺和相同目标载荷下的放大细节，明确区分跨载荷基线与同载荷比较。

两个联合驻点总能量约 `7.94588925557718e-12 J`，弹性能约 `2.7067181663e-12 J`、裂纹能约 `5.2391710893e-12 J`，与第一次相场更新后几乎相同；比第 23 次的失败相场后状态**高约 `1.138499e-18 J`**。没有证据表明找到“交替算法漏掉的更低能稳定状态”。一些越域前函数求值更低能，但未满足力/KKT 与稳定性验收。

内循环后段中，相场增量增长而能量继续下降；这与沿负曲率方向离开驻点一致。机械子步与相场子步的能量变化独立记录。标准路径最大正相场能量差仅 `1.62e-27 J`，机械最大正差 `4.20e-26 J`，属于当前能量量级下的浮点抵消，未作绝对值美化；extra 对应约 `8.08e-28 / 3.88e-26 J`。固定位置相场不增能检查保留原验收标准。ADR 求得驻点不自动保证局部极小值，候选内部能下降也不能证明全局外功—能量平衡已经成立。

每个起点、终止点的 total/elastic/crack、力/KKT、最大/加权平均 s、位形与相场距离均保存在 [ACCEPTANCE_AND_COMPARISON](monolithic/ACCEPTANCE_AND_COMPARISON.json)、[REFERENCE_STATES](energy_path/REFERENCE_STATES.json) 及逐迭代历史中。

## 13. 当前根因分类

依据本轮预声明规则 B：失败载荷的联合一阶驻点具有明确、epsilon 稳定且满足当前约束的负曲率，采用 `COUPLED_BRANCH_LOSS_OR_LOCALIZATION_INSTABILITY`。这里直接证实的是**采样候选分支的负曲率**，并未定位稳定性何时丧失，更未证明某个真实物理裂纹分支。B 状态虽曾“accepted”，也不是局部能量稳定性的保证。

未采用规则 A，因为失败载荷未找到通过曲率检查的局部极小值，而且普通 Phi 谱仍未分辨。未采用规则 C 的主根因标签，因为该规则要求存在联合局部极小值；尽管非光滑切换本身已得到明确证据，这个额外条件没有满足。参考搜索受到当前纯拉伸域的严格限制，也不排除未搜索到的其它稳定状态。

```text
COUPLING_STABILITY_AUDIT: COMPLETE
FIXED_POINT_MAP: VERIFIED_REPRODUCTION
FIXED_POINT_SPECTRAL_RADIUS: UNRESOLVED
MONOLITHIC_REFERENCE: NO_VALID_LOCAL_MINIMUM_FOUND
REDUCED_HESSIAN: NEGATIVE_CURVATURE_FOUND
ROOT_CAUSE: COUPLED_BRANCH_LOSS_OR_LOCALIZATION_INSTABILITY
PURE_TENSION_GATE: FOLLOWS_DIVERGENCE
STEP_SUBDIVISION: DIAGNOSTIC_ONLY
PFPD_AT1_EVOLUTION_V1: REMAINS_NOT_VERIFIED
PRODUCTION MODEL: UNCHANGED
```

审计完成与模型验证通过是两个不同状态。全局能量平衡、完整加载—卸载—重载和当前耦合演化仍未验证。

## 14. 测试与文件保护

基线 **208 passed / 3 failed / 0 skipped**；新增诊断独立运行 **23 passed**；最终完整运行 **231 passed / 相同 3 failed / 0 skipped**，62 条原有 warning。原 208 项通过测试无回归；失败 ID、failure message 和 traceback 完整比对未变，未修改原测试、未增加 xfail/skip。

三个保留失败均来自 `tests/test_pfpd_at1_staggered.py`：

1. `test_saved_coupled_mechanics_phase_and_energy_all_converged`：接受 10 步，原路径期望 19 步。
2. `test_coupled_no_bond_heals_on_unload_and_reload`：未到达卸载/重载，接受 10/19 步。
3. `test_extra_stagger_iterations_change_no_accepted_state_materially`：可比较 8 步，原路径期望 19 步。

新增 5 个测试文件覆盖精确复现、输入和旧状态只读、生产/配置哈希、谱 epsilon 稳定性、特征对复算、切换标记、能量块独立交叉核对、目标有限性、方向梯度、力符号、物理下界与固定边界、投影 KKT、全部优化器起点、鞍点与退出拒绝、缺失终点、原始 Hessian 对称性、Hss 解析结构、约化变量、Schur 及负方向可行性。见 [COMPARISON](tests/COMPARISON.json)、[独立测试日志](tests/DIAGNOSTIC.log)、[完整测试日志](tests/FINAL_FULL.log)。完整 pytest 的退出码 1 由原三个失败产生，不能报告为全绿。

运行硬件为 i7-13700HX，12 个物理核/24 个逻辑核，RAM 约 8.18 GB；Python 3.13.11、NumPy 2.5.3、SciPy 1.18.1、Numba 0.67.0。按实际可用 RAM 扣除至少总内存 20% 与父进程预算后，ProcessPoolExecutor 选用 3 个 worker，每 worker 预算 1 GiB，BLAS/OMP/MKL 均单线程。详见 [HARDWARE](provenance/HARDWARE.json) 与 worker 预算 JSON。

诊断、联合参考、综合分析、图件、最终检查入口分别为 `scripts/run_coupling_stability_audit.py`、`scripts/run_monolithic_stability_audit.py`、`scripts/analyze_coupling_stability_audit.py`、`scripts/visualize_coupling_stability_audit.py`、`scripts/review_coupling_stability_audit.py`。全部为本轮新文件，使用原共享 Python 环境；`PYTHONDONTWRITEBYTECODE=1`、新的 NUMBA_CACHE_DIR、pytest `-p no:cacheprovider` 和新 basetemp 避免写旧缓存。计算配置先于计算保存，SHA256 为 `887184d55d4f4caedae05bd3bda058e19a303d7463878ea0433910c820337d29`。不应直接重跑覆盖本次诊断；若人工授权复算，应新建输出副本。最终自动检查见 [FINAL_REVIEW](provenance/FINAL_REVIEW.json)。

## 15. 图怎么看

全部图件为黑底、Arial 常规字重、英文、无箭头的静态 PNG 和可分享 PDF。未生成本轮未要求的新动画，也未改变科学数组来改善画面。

| 图 | 审核问题与注意点 | 文件 |
|---|---|---|
| 1 iteration_amplification | 更新是衰减还是放大；对照的零更新不伪造比值 | [PNG](figures/iteration_amplification.png) / [PDF](figures/iteration_amplification.pdf) |
| 2 fixed_point_spectrum | 单位圆之外的是条件分支特征值；B/C 普通谱未分辨 | [PNG](figures/fixed_point_spectrum.png) / [PDF](figures/fixed_point_spectrum.pdf) |
| 3 observed_vs_eigenmode | 实际增长与条件模态接近；下图显示 `1-cosine` 以看清微小差别 | [PNG](figures/observed_vs_eigenmode.png) / [PDF](figures/observed_vs_eigenmode.pdf) |
| 4 monolithic_vs_staggered_energy | 能量下降但未收敛；驻点圆点与退出前求值叉号严格区分 | [PNG](figures/monolithic_vs_staggered_energy.png) / [PDF](figures/monolithic_vs_staggered_energy.pdf) |
| 5 smallest_hessian_eigenvalue | B 已低于零；只显示采样点，未插值断言临界载荷 | [PNG](figures/smallest_hessian_eigenvalue.png) / [PDF](figures/smallest_hessian_eigenvalue.pdf) |
| 6 unstable_mode_map | 同一耦合模态的相场与位移幅值；不是物理裂纹预测 | [PNG](figures/unstable_mode_map.png) / [PDF](figures/unstable_mode_map.pdf) |
| 7 substep_diagnostic | 相同终点无法比较，所以显示 N/A；联合参考是鞍点 | [PNG](figures/substep_diagnostic.png) / [PDF](figures/substep_diagnostic.pdf) |
| 8 tension_gate_ordering | 先放大，后越域；时间先后不等于因果证明 | [PNG](figures/tension_gate_ordering.png) / [PDF](figures/tension_gate_ordering.pdf) |

## 16. 唯一下一步建议

人工审核后，开展一次**原候选能量的加载分支稳定性追踪**：在 A 与 B 之间定位最小可行耦合曲率首次变号的载荷，并同时记录不可逆选择器的方向响应，厘清已接受分支何时不再是局部谷底。保持当前物理参数、门限和模型域，用诊断结果决定后续求解策略。本轮没有实施该追踪，也没有自动加入求解器修复或扩展物理模型。
