# PFPD_BAQP_FINITE_DOMAIN_BOUNDARY_CONSISTENCY_V1

**审计完成，最终 `ROOT_CAUSE_PRIMARY=INCONCLUSIVE`。** 已定位非零仿射力的计算来源：有限域中邻点归一化与 Z/形状张量系数经 BAQP 的离散一阶变分传递；它不要求当前点自身的 family 被截断。L=0.75 mm 的 27 个 full-horizon 点确实未通过用户指定的 bulk 检验，但其中没有一个点具有完整的多跳力模板。独立能量差分和固定作者 package 都复现该残差，故不能把它直接归为本地装配 bug，也尚不能宣布为正确的物理非局部边界牵引。

分支 `clot-clearance/pfpd-progress-20261005`；起始及完成时 HEAD 为 `695e75166aa864809460ea9a1845fe72abce229c`。本轮没有提交或推送；新增代码仅用于诊断，既有 production 力学逻辑与旧科学结果未改动。

## 1. 上一轮已经解决了什么

上一轮已确认，C5/C7 平衡后第一内层位移增量未按 h 比例缩小，导致键伸长、Ψ⁺峰值和激活阈值漂移；affine Fbar/Ftilde 正确，侧向 clamp 释放没有改善热点。这些事实本轮直接使用，没有重跑此前平衡或激活 campaign。

## 2. 本轮真正要解释什么

只问：**为什么 uniform affine deformation 需要非零节点力；这些力属于正确的边界相互作用，还是有限域离散算子的问题？**

相场完全禁用。所有新诊断显式读取冻结材料 `SVK, μ=1000 Pa, bulk=666.6666666666666 Pa, λ=0`，没有调用默认 bulk=50000 的 Material。ε=0.003、δ=251.25 μm 固定。

| 试件 | L [mm] | n | h [μm] | 粒子 | 无向键 | d_min≥δ 的粒子 |
|---|---:|---:|---:|---:|---:|---:|
| C5 | 0.25 | 5 | 62.5 | 125 | 5,392 | 0 |
| C7 | 0.25 | 7 | 41.666667 | 343 | 46,059 | 0 |
| L050 | 0.50 | 9 | 62.5 | 729 | 54,360 | 0 |
| L075 | 0.75 | 13 | 62.5 | 2,197 | 196,064 | 27 |

全部使用端面半权重、棱四分之一、角八分之一的节点体积；总体积和六个面面积均守恒。新代码为 [diagnostic helper](../../pd_clot/pfpd_candidate/boundary_force_diagnostics.py)、[审计入口](../../scripts/finite_domain_boundary_audit.py)、[系数分解](../../scripts/finite_domain_normalization.py)；原 `baqp.evaluate()` 未改动。

## 3. BAQP affine energy 是否全局正确

连续参考仅由当前 SVK 给出：

\[
F=\operatorname{diag}(1,1,1+\epsilon),\quad
E_{zz}=\epsilon+\epsilon^2/2=0.0030045,
\]
\[
W=\mu E_{zz}^2=0.00902702025\ \mathrm{Pa},\qquad
P_{zz}=2\mu(1+\epsilon)E_{zz}=6.027027\ \mathrm{Pa}.
\]

\(R_{cont}=P_{zz}L^2\)、\(U_{cont}=WL^3\)、\(dU/d\epsilon=R_{cont}L\)。C5/C7 的 R_cont=3.766891875×10⁻⁷ N，L075 为 3.3902026875×10⁻⁶ N。

唯一无向键的端点权重满足恒等式：

\[
\sum_{i<j}2\bar\omega_{ij}V_iV_j
=\sum_i\frac{V_i}{\omega_{0i}}\sum_{j\in H_i}\omega_{ij}V_j
=\sum_i V_i.
\]

实际 weight-sum/L³ 分别为 1、0.9999999999999998、1、1.0000000000000002；使用重构 F 计算的 U/U_cont 为 0.999999999999938、0.999999999999973、0.999999999999991、0.999999999999931。`GLOBAL_AFFINE_ENERGY_NORMALIZATION=EXACT` 指上述代数恒等式，数值能量保留浮点重构误差。

## 4. 为什么 affine F 正确但 force 不平衡

运动学完整性控制的是 F 的值；逐点平衡控制的是 U 对各个节点独立虚位移的导数。前者和一参数仿射能量恒等式并不强制后者为零。

固定 affine P 时，某节点的力仍涉及 **邻居的 ω0、K⁻¹、Z**。即使当前点的 family 完整，邻居系数仍可能不均匀。以 L075 core 残差最大的一个点为例，精确的 self/neighbor 分解给出：

| 有向节点行贡献 | 占该点有符号 fz |
|---|---:|
| direct 自身归一化 | 约 2.9×10⁻¹³ |
| direct 邻居归一化 | 4.14365% |
| 自身 Z 项 | 约 −3.4×10⁻¹⁶ |
| 邻居 Z 项 | 95.85635% |

自身完整 family 的奇对称求和已经消去自身项，残差来自邻居系数。这是代数定位，不是根据 hotspot 图的相关性猜测。详见 [NORMALIZATION_CHANNELS.json](NORMALIZATION_CHANNELS.json)。

另做 8 个独立能量差分方向，含 core 单点、core 随机、collar 随机虚位移，每方向三个预设步长；均与 −f·δx 一致。core 最强单点方向误差低至 3.63×10⁻⁹，全部方向至少一个步长的误差 ≤1.09×10⁻⁶。见 [FIRST_VARIATION.json](FIRST_VARIATION.json)。因此非零 core 力确实是当前离散能量的梯度，并非分解 helper 制造出来的值。

## 5. 非零 force 位于哪些 z layers

以下残差统一为 **自由节点 max|f_i|/R_cont**，不同于上一轮按求解器 force floor 定义的残差：C5=0.07452309433，C7=0.02834450528，两者均 `FAIL`。

下半域的 layer Σfz/R_cont 为：

| 层 | C5 | C7 |
|---|---:|---:|
| 最外加载面 | 0.62413533 | 0.44500075 |
| 第一内层 | 0.75172934 | 0.66764675 |
| 第二内层 | 中面，约 0 | 0.32970425 |

上半域取相反符号，中面靠对称性抵消。完整精确值、L1/L2、max/mean 和层标记保存在 [C5_layers.csv](../../results/pfpd_baqp_finite_domain_boundary_consistency_v1/tables/C5_layers.csv) 与 [C7_layers.csv](../../results/pfpd_baqp_finite_domain_boundary_consistency_v1/tables/C7_layers.csv)。

[图 1：仿射节点力的层分布与 L2](../../results/pfpd_baqp_finite_domain_boundary_consistency_v1/figures/affine_internal_force_layers.png)。所有粒子的坐标、六面距离、support、ω0、K condition、各项力保存在 `tables/*_particles.csv`。

λ=0 的侧面检查中，x/y 节点力仅为数值噪声，侧面净侧向力/R_cont 不超过约 1.4×10⁻¹⁵。侧面节点可同时处于端部体积层，故它们的轴向力不能直接解释成侧面连续 traction。

## 6. 第一内层粒子的 force 来自哪些 bonds

C5 中心轴第一内层 ID 为 bottom=61、top=63；C7 为 169、173。全部 incident bonds 已按邻点所在层、|nz| 的四档和 r/δ 的四档逐项保存，没有只筛选热点键。

| 分组占第一内层有符号净 fz | C5 | C7 |
|---|---:|---:|
| 邻点在加载面 | −20.71796% | −9.72252% |
| 其余邻点合计 | 120.71796% | 109.72252% |

百分比有负值或大于 100%，是实际抵消关系，未取绝对值伪造正贡献。加载面邻居抵抗内层向中心的合力，更深内部邻居的贡献更强；底部第一内层的合力向上、顶部向下，与上一轮内层向中心松弛、端部位移增量变大的方向一致。

`*_first_layer_bonds.csv` 保存逐键 direct/Z、长度、方向和端点权重；`*_bond_groups.csv` 保存全部分组，包括同一第一内层与更深部分。不同 grouping 是同一批键的不同分箱，不能彼此相加。

## 7. direct 与 Z correction 谁贡献最大

| 第一内层 | direct | Z correction |
|---|---:|---:|
| C5 | 18.42164% | 81.57836% |
| C7 | 21.53659% | 78.46341% |

Z 是两级中较大的项。配置在计算前将“强主导”设为绝对贡献份额 ≥80%；C7 未达到，所以统一状态标签为 `MIXED`，并明确报告 Z 约占四分之三以上，避免标签掩盖数值。

direct+Z 重建原 force 的相对误差为 3.82×10⁻¹⁶（C5）、6.37×10⁻¹⁶（C7），远小于 10⁻¹²。该 Z 项本身并未证明“错误”；它是能量变分的一部分，不能因其占比大就删除或调弱。

[图 3：direct/Z 与邻点分组](../../results/pfpd_baqp_finite_domain_boundary_consistency_v1/figures/first_layer_force_decomposition.png)。

## 8. wbar endpoint normalization 的作用

按唯一键的 i/j 下标各分一半并**重新计算各自 Z**后，两个 force 场相加重建原场，C5/C7 相对误差 3.43×10⁻¹⁶ / 6.80×10⁻¹⁶。i/j 贡献约各半，但下标先后没有“边界端点”的物理含义。

进一步按几何来源分解：加载表面端点的归一化贡献，C5 为 −2.03410%，C7 为 +1.56613%；其余端点分别为 102.03410%、98.43387%。这削弱了“最外加载面端点单独主导”的解释。**C5/C7 所有 family 都被截断**，不能把“非加载面”误称为完整 bulk。

再用同 h、δ 的完整规则晶格计算唯一 ω0∞，只作线性恒等分解：ω̄=ω/ω0∞+(ω̄−ω/ω0∞)。后一个归一化超额项占 C5/C7 第一内层净力的 54.63633% / 54.28754%；固定晶格权重部分仍占约 45%。L075 core 最大残差点则为 38.59581% / 61.40419%。因此归一化不可忽略，但不是唯一来源。

此处两个部分始终相加为原模型；没有把 ω/ω0∞ 当成替代权重或跑任何修正。单独采用它会改变全局能量归一化。`WBAR_ENDPOINT_NORMALIZATION_EFFECT=INCONCLUSIVE` 指尚不能判定需要修复的“bias”，并非作用为零。

[图 4：端点及归一化分解](../../results/pfpd_baqp_finite_domain_boundary_consistency_v1/figures/endpoint_normalization_decomposition.png)。

## 9. larger specimen 是否出现 full-support bulk

按任务定义 d_min≥δ，L075 确实存在 27 个点，占 1.22895%。几何设计没有失败。L050 的中心距离 250 μm 仍小于 δ=251.25 μm，故它没有完整 family。

同时检查完整力模板：L075 没有任何点同时满足“自身 family 及全部邻居 family 完整”。再向外检查一层仍为零。这个结果解释了为什么 **几何 full-horizon 点不等于不受边界系数影响的力算子内点**。

BAQP force 读取邻点 Z，Z 又读取邻点的邻点 ω0。按当前公式可以给出最多三跳、距离不超过 3δ 的几何系数依赖上界；这是本轮代码推导，既不宣布实际最小边界层为 3δ，也不建议直接施加 3δ BC。见 [理论说明 §5](BOUNDARY_THEORY_NOTES.md)。

## 10. full-support bulk 是否通过 affine force patch test

**没有。阈值保持原样，没有用“力模板”解释把 FAIL 改为 PASS。**

| L075 统计 | full-support core | δ collar |
|---|---:|---:|
| 粒子数 | 27 | 2,170 |
| max‖f‖ [N] | 3.61981993×10⁻¹⁰ | 1.81082924×10⁻⁸ |
| RMS‖f‖ [N] | 2.94856112×10⁻¹⁰ | 7.13211831×10⁻⁹ |
| max‖f‖/R_cont | 1.06772965×10⁻⁴ | 5.34135982×10⁻³ |
| RMS‖f‖/R_cont | 8.69730040×10⁻⁵ | 2.10374393×10⁻³ |

bulk/collar 的 max 比为 0.01998984688，仅低约 50 倍，不满足至少低 10⁴ 倍，也不接近数值零。

为区分本地程序失配与模型本身，新增一个 **L075、同坐标/体积/材料/kernel** 的固定作者 package 仿射对照。调用未修改的 `RKCMaterial(monomial=:C1,epsilon=0)`；不调用时间推进。全场力相对误差 5.86×10⁻¹⁵，core 误差 2.32×10⁻¹⁴，能量误差 1.63×10⁻¹³。作者实现也存在同一 residual。[对照结果](author_L075/COMPARISON.json)

[图 5：固定 h 的大域 affine patch test](../../results/pfpd_baqp_finite_domain_boundary_consistency_v1/figures/large_domain_affine_patch_test.png)。当前能确定的是**用户定义的几何 core 检验失败**；是否属于真正不受边界影响的 bulk 缺陷，尚不能由这个只有 1.493δ 半宽的试件确认。

## 11. manufactured external force 是否局限在 boundary collar

按任务定义 f_ext,req=−f_int,affine。C5 与 L075 的所有自由节点直接平衡残差为零；这是 manufactured balance 的构造检验，不是新平衡求解，也不能单独作为物理正确性的证据。

L075 中 99.93121% 的绝对轴向节点力在 d_min<δ 区域；core 占 0.06879%。尽管大部分力在 collar，core 尾部明显高于噪声，因此标签为 **`BOUNDARY_FORCE_LOCALIZATION=MIXED`**，不称“只有 δ collar”。全局合力与转矩的归一化误差约 10⁻¹⁶ / 10⁻¹⁷。

[图 6：manufactured force 的三维位置、截面与层分布](../../results/pfpd_baqp_finite_domain_boundary_consistency_v1/figures/manufactured_external_force_map.png)。图中明确标为 manufactured force，不称物理裂纹载荷。top/bottom/side collar 在角部可能重叠，区域统计不能直接求和。

## 12. 其积分是否恢复 continuum traction

按任务要求直接计算上半域总外力：

| 试件 | R_manufactured_top/R_cont | 相对差 |
|---|---:|---:|
| C5 | 1.3758646719 | +37.58647% |
| C7 | 1.4423517512 | +44.23518% |
| L050 | 1.0982921860 | +9.82922% |
| L075 | 1.0022695587 | +0.22696% |

底部相反、大小对称。L075 单独满足预声明 1% 比较标准，但 C5/C7/L050 不满足，故本轮整体标签 **`CONTINUUM_TRACTION_RECOVERY=INCONSISTENT`**，不能用大试件掩盖原试件偏差。

然而另一项严格检查在全部几何上成立：

\[
-\sum_i f^{int}_{iz}(Z_i+L/2)=R_{cont}L.
\]

这是分布外力与仿射位移的虚功。定义上半域有符号力的重心 z̄_F，则 2R₊z̄_F=R_cont L。因此把同样虚功的力分布到物体内部时，半域总力和表面集中力不必相等。C5/C7 的 z̄_F=90.85196/86.66402 μm，小于外表面 125 μm；L075 为 374.15084 μm，接近外表面 375 μm。该代数关系对四例的误差 <4.3×10⁻¹⁴。[VIRTUAL_WORK_TRACTION.json](VIRTUAL_WORK_TRACTION.json)

| 到达观测半域合力百分比的首个层深度 [μm] | d50 | d90 | d95 | d99 |
|---|---:|---:|---:|---:|
| C5 | 62.5 | 62.5 | 62.5 | 62.5 |
| C7 | 41.6667 | 83.3333 | 83.3333 | 83.3333 |
| L075 | 0 | 62.5 | 62.5 | 62.5 |

分位数分母是**实际半域合力**，曲线纵轴仍为 /R_cont；两者没有混用。L075 有负层力和 overshoot，故这些首个 crossing 不能解释成载荷已在该深度完全收敛；其累计值此后始终处于最终值 ±1% 内所需深度为 312.5 μm=1.24378δ。[图 2：累计力](../../results/pfpd_baqp_finite_domain_boundary_consistency_v1/figures/cumulative_nonlocal_boundary_force.png)

## 13. Peridynamics.jl / Peridigm / 文献如何处理边界

定向来源和限定见 [BOUNDARY_THEORY_NOTES.md](BOUNDARY_THEORY_NOTES.md)：

- PFPD 原文声明有限厚度点层；其示例也有单层体力密度加载，未给本 cube 的普适最优层厚。
- 固定 Peridynamics.jl 拉伸 benchmark 延伸左端三个点层，右端一层按 F/ΣV 分配力密度，horizon/spacing=3.015。该 benchmark 用 CMaterial/BBMaterial，不能冒充当前 BAQP 的验证。
- 固定 Peridigm MeshConverter 按用户输入的 `Boundary Layer Thickness` 生成六面体积 node sets。
- D’Elia–Yu 的局部到非局部映射使用 interaction domain；文中 Poisson 厚 δ，而双积分 LPS 厚 2δ。这支持“应检查完整算子支持”，其 O(δ²) 定理不直接适用于本 BAQP。

因此 `LITERATURE_BOUNDARY_TREATMENT=MIXED_BY_METHOD`，不是一条“δ 厚 Dirichlet 必然正确”的统一规则。

## 14. surface-only BC 是否与当前 nonlocal operator 相容

在当前有限域能量、仅最外顶底节点固定且其余节点无外力时，affine displacement **不是驻定解**，这一事实已确定。第一内层净力的方向解释了为何放开自由节点后朝中心松弛、端部 gap 增大。

但“该位移边界问题是定义明确但不同的非局部边值问题”，与“目标物理牵引本应配额外 interaction domain”仍未区分。不能因为 affine 不平衡就认定位移约束无效，或直接把制造力当成物理 BC。C5/C7 的 δ 大于 L，内置 δ 厚上下层会重叠、可能消除全部自由内区；这种控制即使得到 affine，也不能构成有效因果证据。

## 15. boundary control 的结果

**`BOUNDARY_CONTROL_RUN=NO`，`BOUNDARY_CONTROL_EFFECT=NOT_RUN`。** 按任务 §45 CASE 2，显著 core residual 出现后先核对算子，不能直接做 BC fix。我们完成了能量差分、作者实现和邻居系数来源检查；它们尚未隔离“扩大作用范围的边界尾部”与“真正 bulk 缺陷”，所以没有进入 §35–41。

没有生成假 comparison 图，没有把 manufactured forcing 当成修正，没有在原 C5/C7 上增加新控制。`RECOMMENDED_BOUNDARY_MODEL=DEFERRED_PENDING_FORCE_STENCIL_AUDIT`；没有提出可直接投入生产的模型。

## 16. 真正 root cause

**已定位机制，尚未满足真正根因判定的排他证据。**

| 假设 | 本轮判断 |
|---|---|
| H1：正确牵引仅在 δ collar 分布 | 大部分力在 collar，但严格支持范围、原 C5/C7 半域总力均不满足；不能确认 |
| H2：surface Dirichlet 与所需非局部 BC 不相容 | 仍有可能；因 core gate 失败未开展有鉴别力的边界控制 |
| H3：BAQP bulk force patch failure | 用户定义的 d_min≥δ 检验确实 FAIL；未隔离完整力模板，不能升级为无限域 bulk 缺陷 |
| H4：wbar boundary bias | 归一化超额明显贡献，但保留其余项仍有 residual；没有获理论与控制支持的替代权重 |
| H5：Z correction assembly discrepancy | Z 较大，但作者代码逐点一致，能量变分一致；未发现本地装配差异 |
| H6：组合边界算子作用 | 当前最有力线索，尚缺完全闭合力模板与物理边界比较来完成因果区分 |

最低标准 A–E 已完成；F“两个独立证据支持同一个可排他物理根因”未满足。因此不宣布 `ROOT_CAUSE_IDENTIFIED`，最终 **`ROOT_CAUSE_PRIMARY=INCONCLUSIVE`、`ROOT_CAUSE_SECONDARY=NONE`**。不能为了给出非 INCONCLUSIVE 标签而把当前几何 core 改称真正无边界影响的 bulk。

`BAQP_OPERATOR_DEFECT_LOCALIZED_TO=NOT_ESTABLISHED`。定位到的**残差通道**主要为 neighbor ZⱼKⱼ⁻¹ 及 ω0ⱼ；“残差通道”不是“已证明需要修复的代码行”。

7 个直接相关科学测试通过：[日志](TARGETED_TESTS.log)。全部运行只是完整态计算、差分或对照，没有新的机械平衡、相场、稳定性或 CFD。主 Python 审计用 8 workers、每进程 BLAS=1，仿射阶段 4 个独立 job、差分阶段 8 个。预估并发峰值 3.583 GB < 总 RAM 的 60%；观察主阶段进程树峰值 2.422 GB，最低可用 4.201 GB，任务进程 VmSwap=0、该阶段 swap-out 增量为零。已有 swap 没有被清空。

主审计 5.765 s；作者包单线程对照 72.537 s（含加载/编译），该子树峰值 1.894 GB。分阶段资源见 [RESOURCES.json](RESOURCES.json)、[作者资源](author_L075/RESOURCES.json)、[RESOURCE_SUMMARY.json](RESOURCE_SUMMARY.json)。未采样所有报告/绘图进程的全流程瞬时合计，因此阶段峰值不冒充全任务严格最大值。GPU_NOT_USED_BY_DESIGN。

## 17. 对当前 C5/C7 activation 漂移意味着什么

漂移的直接机械链条进一步清楚了：有限域一阶变分产生朝内的第一内层合力；表面约束、自由内层和零外力条件允许内层松弛；原两级近似有限的端部 gap 被更小 h 放大为局部应变。它与相场求解器无关。

本轮没有证明哪种边界表达对应目标物理载荷，也没有选择一个修正。当前 activation brackets、接受态与科学结论保持原样，不能声称已解决 mesh objectivity。以后若改变 interaction domain 或 BC，必须在人工批准后重新审视机械参考与激活比较，不能用 C5/C7 activation 对齐来挑边界模型。

## 18. 下一步只给一个最高价值开发任务

**完全闭合 BAQP 力模板的仿射参考测试。** 固定现有 h、δ、材料与 kernel，围绕一个中心观测区补足算子所需的完整邻域链，并与同一固定 RKCMaterial 比较；分别检查一跳、两跳、三跳闭合后的系数与 affine 力。只需少量中心观测点/仿射评价，无需 fracture、平衡大计算或参数调节。这样才能区分“边界作用范围被低估”与真正的 bulk force inconsistency，再决定是否值得开发物理 volume-constraint/collar 边界。

本轮到此停止，等待人工审核。没有自动开始上述任务。

```text
GLOBAL_AFFINE_ENERGY_NORMALIZATION: EXACT
AFFINE_KINEMATICS: VERIFIED
AFFINE_FORCE_EQUILIBRIUM_C5: FAIL
AFFINE_FORCE_EQUILIBRIUM_C7: FAIL
FULL_SUPPORT_BULK_AVAILABLE: YES
FULL_SUPPORT_BULK_PATCH_TEST: FAIL
BOUNDARY_FORCE_LOCALIZATION: MIXED
CONTINUUM_TRACTION_RECOVERY: INCONSISTENT
FIRST_LAYER_FORCE_DOMINANT_TERM: MIXED
WBAR_ENDPOINT_NORMALIZATION_EFFECT: INCONCLUSIVE
LITERATURE_BOUNDARY_TREATMENT: MIXED_BY_METHOD
BOUNDARY_CONTROL_RUN: NO
BOUNDARY_CONTROL_EFFECT: NOT_RUN
ROOT_CAUSE_PRIMARY: INCONCLUSIVE
ROOT_CAUSE_SECONDARY: NONE
CURRENT_C5_C7_ACTIVATION_RESULTS: UNCHANGED
C7_STABILITY: NOT_RUN
FRACTURE: NOT_RUN
```
