# PFPD_AT1_ACTIVATION_ROOT_CAUSE_AND_FIX_V1

**结论：Phase A 完成；停止状态 `DIAGNOSIS_INCONCLUSIVE`。** 已定位激活漂移的直接数值机制：真实平衡位移在端面与第一内层之间形成近似不随 h 缩小的位移增量，细化后局部应变和 Ψ⁺ 峰值上升。仿射重构没有热点；解除端面侧向约束没有改善热点。现有证据尚不足以在有限域非局部支持、BAQP 力分配及轴向边界施加方式之间作因果区分，因此本轮没有强行开发修正。

模型为 **`PFPD_AT1_SPECTRAL_CANDIDATE_V1`**：BAQP、AT1 不可逆相场、SVK Green 应变谱分裂、g(s)=(1−s)²、当前 h(s)=1、PETSc RSLS。**不是 2604.09215v1 的逐式复现。** 论文继续作为理论/归一化/结构参考，不再是开发硬门。

## 1. 现在到底是什么问题

同一 L=250 μm、δ=251.25 μm 试件中，历史 C5/C7 的激活应变分别为 0.004497390736212 和 0.003137607455254；首个激活位置均为上下端面中心附近、d/h=0.5 的轴向键。C5 已有早期稳定性证据，C7 已有早期驻定分支，但这些均不证明激活阈值具有空间客观性。

在 s=0 时，原候选给出

\[
q_b=2(Y_c-\Psi_b^+),\qquad Y_c=0.1010025\ \mathrm{Pa}.
\]

因此先解释平衡机械场的 Ψ⁺，再讨论激活。活动键比例不是本轮主要物理量。公式来源是原 [q_fields](../../scripts/c5_activation_audit.py) 与 [AT1 候选](../../pd_clot/pfpd_candidate/at1_evolution.py)，未修改。

## 2. 为什么不再关注旧 cleanup tests

按本轮用户指令，历史缓存/manifest 的 byte-hash failures 不再阻塞科学开发。本轮没有恢复旧缓存，没有修改旧断言或清单，没有运行全仓库哈希扫描或全量 pytest；这些旧失败仍是独立历史事项，未宣称已修复。

读取了最新 activation mesh audit、C5 stability 和 C7 early branch 报告。复用 C7 ε=0.0025/0.003 的两个已接受 intact 状态；新增 C5 baseline 两点以及 C5/C7 axial-only 各两点，共 **6 次新平衡求解 + 2 个复用状态**。没有重跑历史 campaign。

## 3. C5/C7 activation 为什么值得怀疑

在相同 ε 下，材料、Gc、Yc、δ 和相场均相同；两个网格却距离 q=0 显著不同：

| 网格 | ε | baseline max Ψ⁺/Yc | axial-only max Ψ⁺/Yc |
|---|---:|---:|---:|
| C5 | 0.0025 | 0.30792575 | 0.31221585 |
| C7 | 0.0025 | 0.63384863 | 0.64776905 |
| C5 | 0.0030 | 0.44380049 | 0.44998922 |
| C7 | 0.0030 | 0.91389212 | 0.93399135 |

全部状态保持 s=0，并通过自由力 ≤1e−9、positive J、安全、固定分量精确满足和反力平衡门槛。最大自由力残差 2.755e−13，最小 J=0.99920178，最大反力平衡相对误差 4.706e−15。负 Green 主应变允许存在，未当成失败。

## 4. affine state 结果

采用与轴向位移一致的诊断场

\[
u_z^{aff}=\varepsilon(Z+L/2),\quad u_x^{aff}=u_y^{aff}=0,
\quad F^{aff}=\operatorname{diag}(1,1,1+\varepsilon).
\]

当前 λ=0、ν=0，故均匀单轴状态的侧向应变为零；不能只凭“端面侧向固定”就预设明显的 Poisson 约束效应。此 affine 场仅作为 oracle，并没有被当作已平衡状态。

ε=0.003 时，Ezz=0.0030045，Ψ⁺=0.00902702025 Pa，q≈0.1839509595 Pa，两网格相同。Fbar/Ftilde 对 exact F 的最大误差 5.11e−15；两应变所有案例 Ψ⁺ 的最大 CV=5.46e−13。各轴向距离层和方向分组均未发现 affine hotspot。

但 affine 自由力残差 C5/C7 分别为 0.0235082/0.0156735，远大于平衡门槛。**仿射重构准确，并不意味着该有限域端面加载下仿射场满足离散机械平衡。**

见 [Figure 1](../../results/pfpd_at1_activation_root_cause_and_fix_v1/figures/baseline_affine_vs_equilibrium.png)。各方向及距离层的 mean/median/P95/max/CV 保存在 `analysis/` 和结果 `tables/`。

## 5. equilibrium hotspot 结果

定义 Aψ=Ψ⁺eq/Ψ⁺aff、ΔΨ=Ψ⁺eq−Ψ⁺aff，以及最大正 Green 主应变的放大率 AE。仿射分母须 >1e−12 Yc；本次没有需要排除的键。

ε=0.003 的 baseline 峰值：

| 量 | C5 | C7 |
|---|---:|---:|
| max Ψ⁺ [Pa] | 0.04482496 | 0.09230539 |
| max Aψ | 4.965643 | 10.225455 |
| 热点 AE | 2.228328 | 3.197681 |
| 热点实际键伸长应变 | 0.00667275 | 0.00956172 |
| 第一段中心轴位移增量 [μm] | 0.417047 | 0.398405 |
| h [μm] | 62.5 | 41.666667 |
| 第一段占总端位移比例 | 55.606% | 53.121% |
| 全体粒子 non-affine RMS [μm] | 0.108686 | 0.135932 |

两级第一段位移增量仅下降约 4.5%，而 h 下降 1/3，所以局部伸长应变增长约 43.3%。上下对称的第一段占据大部分端位移，中间某些中心轴区间出现压缩。已逐 z 层保存 non-affine RMS、uz 偏差、侧向位移均值/最大值。

BAQP 中

\[
\widetilde F_b=\overline F_{ij}
+\frac{(\eta_b-\overline F_{ij}\xi_b)\otimes\xi_b}{|\xi_b|^2},
\qquad \widetilde F_b\xi_b=\eta_b.
\]

这一恒等式的数值误差 ≤8.14e−20 m。热点 Ftilde 的轴向伸长因此与求解位移的实际端点差一致；不是均匀位移场被错误算出一个 Ψ⁺ 热点。它仍不能排除 BAQP 离散能量/力分配导致非仿射解的可能性。公式直接来自未修改的 [baqp.evaluate](../../pd_clot/pfpd_candidate/baqp.py)。

当前 λ=0，Ψ⁺=μΣ⟨Ea⟩₊²，热点 Aψ 与 AE² 接近，符合“位移差集中 → 应变集中 → Ψ⁺ 增大 → q 更早过零”的机制。

见 [中心轴位移增量图](../../results/pfpd_at1_activation_root_cause_and_fix_v1/figures/centerline_displacement_gaps.png) 和 [统一空间图](../../results/pfpd_at1_activation_root_cause_and_fix_v1/figures/driver_spatial_map.png)。

## 6. hotspot 随 d/h 还是 d/δ

对 |nz|≥0.95 的轴向键，按键中点到最近端面的距离分组，以层平均 Aψ 为主指标；使用预声明的 1001 点公共线性插值网格，不拟合、不外推。

| ε | RMSE_h | RMSE_delta | RMSE_delta/RMSE_h |
|---|---:|---:|---:|
| 0.0025 | 1.674167 | 0.370351 | 0.2212 |
| 0.0030 | 1.675051 | 0.370465 | 0.2212 |

两应变均满足 RMSE_delta < 0.5 RMSE_h，故按指定规则标为 **PHYSICAL_HORIZON_SUPPORTED**。P95/max Aψ、Ψ⁺/Yc 和 q/Yc 的补充 profile 也给出更小的 δ 坐标 RMSE。

**适用范围有限：** h 坐标公共区间为 [0.5,2]，δ 坐标公共区间为 [0.124378,0.497512]；C7 最近第一层的 d/δ=0.082919 不在后一公共区间。该结果描述重叠区间内的 profile 相似度，不能证明最靠近端面的峰值收敛、有限物理宽度或 activation objectivity。两个网格 δ 相同，所以这里也不能独立区分 δ 尺度与一般固定物理长度尺度。

见 [Figure 2](../../results/pfpd_at1_activation_root_cause_and_fix_v1/figures/baseline_profile_scaling.png)。原始 profile 与每个插值比较网格均保留。

## 7. support 的作用

提取原 Family 的 neighbor count、ω0=ΣωV、K 的特征值/condition、trace(K)、1/ω0；没有引入新的校正模型。δ≈L，最内部可用粒子仍不具有完整球形 horizon，报告不称其为 full-horizon bulk。

ε=0.003、全部轴向键的 Spearman ρ：

| 指标与 Aψ 的关联 | C5 | C7 |
|---|---:|---:|
| 两端平均 weighted support | −0.31275 | −0.38223 |
| 两端最大 shape condition | +0.42721 | +0.45716 |
| d/h | −0.77224 | −0.57652 |

按事前 |ρ|≥0.7 才为两级 STRONG 的规则，**SUPPORT_ASSOCIATION=MIXED**。键共享端点，表中的统计不构成独立样本因果检验，p 值不用于判定机制。

更关键的对照：历史激活键端点在各自 z 层的 weighted support 都处于最高位；shape condition 则处于最低约 4.0%（C5）/2.04%（C7）。端面中心 K condition≈1.645/1.653，第一内层中心≈1.113/1.221。**热点不在同层支持最差或最病态的粒子上。** 这削弱了“局部 K 病态/支持量最低直接造成热点”的解释；不能因此排除整体有限域非局部支持和边界力分配。

见 [Figure 3](../../results/pfpd_at1_activation_root_cause_and_fix_v1/figures/support_and_driver.png)、`*_PARTICLE_SUPPORT.csv` 及逐端点同层比较。

## 8. full clamp 的作用

axial-only 控制仅在上下端面固定 uz。使用同一几何规则：底面中心固定 ux/uy；底面 (+L/2,0,−L/2) 点仅再固定 uy。没有硬编码粒子编号；约束对六个刚体模式的秩为 6，额外侧向 gauge 恰为三个。C5/C7 自由分量为 322/928。

除这些位移约束外，材料、δ、相场律、Gc、Yc、网格和 PETSc RSLS 设置不变。控制是 **DIAGNOSTIC CONTROL，NOT PRODUCTION BC**。

ε=0.003 对照：

| 量 | C5 full → axial-only | C7 full → axial-only |
|---|---:|---:|
| 第一层 max Ψ⁺/Yc | 0.443800 → 0.449989 | 0.913892 → 0.933991 |
| 第一层峰值 reduction | −1.3945% | −2.1993% |
| max Aψ | 4.96564 → 5.03489 | 10.22545 → 10.45034 |
| 全键 Ψ⁺ P95 [Pa] | 0.0130034 → 0.0131150 | 0.0121050 → 0.0119935 |
| 反力名义应力 [Pa] | 3.75085 → 3.74358 | 2.81739 → 2.80386 |
| non-affine RMS [μm] | 0.108686 → 0.109297 | 0.135932 → 0.136899 |

C5/C7 最大驱动的对称相对差异没有减小，反而增大约 1.003%。两级峰值变化及差异变化均小于 5%，故 **CLAMP_ASSOCIATION=WEAK**，远未达到峰值降低 ≥20% 且跨网格差异降低 ≥30% 的 STRONG 条件。

min q 仍出现在原两根上下对称轴向键，d/h=0.5，未被随机移到别处。此控制只否定“侧向夹持是主要原因”的当前证据，不否定两种 BC 共有的轴向端面加载可能参与形成边界层。

见 [Figure 4](../../results/pfpd_at1_activation_root_cause_and_fix_v1/figures/full_clamp_vs_axial_only.png)。

## 9. 最可信 root cause

**直接机制已定位，原始因果尚未唯一确定。** 当前排序如下：

| 假设 | 当前证据判断 |
|---|---|
| H2：有限域非局部支持/端面加载相容性 | 优先待判定；有限域没有 full-horizon bulk，affine 自由力不为零、平衡后边界位移增量集中；但没有支持校正的控制证据 |
| H3：BAQP 边界力分配/非仿射重构作用 | 与 H2 尚不能区分；affine completeness 通过，排除的是 affine 重构热点，未排除非仿射解中的算子作用 |
| H4：侧向 full clamp 与支持的组合 | 侧向释放效应很小，不支持其为主导的组合修正 |
| H1：端面侧向约束主导 | 本轮控制为 WEAK，不支持据此选用 BC fix |
| H5：相场缺乏有限宽度 | 本轮 s=0，不能作为当前 pre-activation 漂移的已验证解释；后激活支持未测 |
| H6 | 最终标签 **ROOT_CAUSE_PRIMARY=INCONCLUSIVE** |

不把中等 support correlation 当成因果，也不把 d/δ 的较小 RMSE 直接当成 surface correction 的依据。

## 10. 本轮开发了什么 candidate fix

实现了独立的 [component-wise BC 控制适配器](../../pd_clot/pfpd_candidate/activation_boundary_control.py)，基线选择 `FULL_CLAMP_BASELINE`；代码内 `AXIAL_ONLY_CANDIDATE` 只是可选控制模式名称。本轮没有将其晋升为 Phase B 科学修正，故 **CANDIDATE_FIX=NONE，CANDIDATE_IMPLEMENTED=NO** 指“没有选定并开发获科学证据支持的修正”，不表示控制实验未实现。

原 BAQP、spectral split、AT1、历史 solver 和历史结果没有被替换。因为 Phase A 仍 inconclusive，没有自行添加 surface correction、修改 horizon/material/phase law，亦没有运行 candidate activation search。

## 11. candidate 是否减少 C5/C7 discrepancy

**N/A，未进入 Phase B。** 两点局部线性外推仅用于检查 pre-activation mechanics 的趋势：

| 网格/BC | 两点线性预测 εact | 已知真实 baseline εact | 预测误差 |
|---|---:|---:|---:|
| C5 baseline | 0.005046736 | 0.004497391 | +12.215% |
| C7 baseline | 0.003153740 | 0.003137607 | +0.514% |
| C5 axial-only | 0.004996071 | N/A | N/A |
| C7 axial-only | 0.003115310 | N/A | N/A |

这些预测支持“共同应变机械峰值已能解释 C7 更早激活的方向”，但 C5 的远距离线性外推误差较大，不能替代真实平衡 bracket。没有把控制预测列成已测 activation。

用户指定的对称差异为 |ε7−ε5|/mean(ε7,ε5)，历史 baseline 值 **35.61974%**。此前约 −30.235% 使用 C5 作分母，两者不是矛盾。候选差异和改善等级均为 N/A / NOT_RUN。

历史激活反力为 C5 5.63844 Pa、C7 2.94749 Pa；没有候选激活反力可以比较。见 [Figure 6](../../results/pfpd_at1_activation_root_cause_and_fix_v1/figures/activation_before_after.png) 和 [Figure 7](../../results/pfpd_at1_activation_root_cause_and_fix_v1/figures/activation_discrepancy.png)，未运行部分明确标记，不画假零值。

## 12. hotspot 是否变得更合理

axial-only 控制没有改善端面第一层锁定：峰值仍是 C5 IDs 4488/5059、C7 IDs 40174/44110，均上下对称、轴向、d/h=0.5。d 分别为 31.25/20.8333 μm，d/δ 为 0.124378/0.082919。控制的热点略增强，没有提供更合理的物理支持范围证据。

这里只描述 ε=0.003 的 pre-activation 热点，不冒充新的 candidate activation 位置。

## 13. phase support 是否仍局域到 2 bonds

**本轮 NOT_ASSESSED。** 旧分支的两个活动键是历史事实；本轮没有进入候选后激活分支，没有产生新的 PERSISTS/BROADENS 结论，也未生成 `early_phase_support` 图。两个键不是裂纹，activation 也不是 fracture。

## 14. 当前还没有解决什么

尚未解释：为什么当前有限域 BAQP 能量在仅端面施加轴向位移时选择近似有限的边界位移增量；支持截断、非局部自然边界与力分配各自贡献多大；怎样用有明确依据的方法修正；激活峰值及有限相场宽度是否随细化客观。

本轮必需数值门槛全部通过：[NUMERICAL_GATES.json](NUMERICAL_GATES.json)。C5 baseline/control 接受态上 dU/dε 与 Rtop·L 的误差分别 8.08e−10/8.58e−10；C5/C7 两 BC 的零载荷力为零。5 项直接相关测试通过，无全量回归。[TARGETED_TESTS.json](TARGETED_TESTS.json)

12 physical / 24 logical cores、7.62 GiB RAM；6 个 worker、1 MPI rank、每进程单线程。有效 Phase A 约 305 s，观察进程树峰值约 3.13 GiB，最低可用 RAM 34.9%，超过 15% 保留要求。GPU_NOT_USED_BY_DESIGN。新增计算没有稳定性、C9 或参数 sweep。

## 15. 对最终 clot simulation 意味着什么

目前可以明确：阈值漂移在尚未激活相场时已存在于平衡力学场，不应通过 PETSc 活动集调参、相场续算或简单解除侧向夹持来掩盖。候选仍不能用于声称真实血栓起裂阈值已收敛，亦不能将 axial-only 作为最终血管壁边界条件。

保持当前生产实现和已有 C5/C7 科学状态，保留此次对照证据；没有本轮有效 fix，不报告 refinement sensitivity reduced。

## 16. 下一步只给一个最高价值任务

**有限域 BAQP 轴向边界载荷一致性审计。** 在同一 C5/C7 几何上，以恒应变仿射场为起点，对照离散内部力的层分布、端面反力与连续表面牵引的虚功，核对当前边界施加方式是否与现有 BAQP 及其原始边界理论一致；仅在找到明确文献/成熟实现支持后选择一个最小边界处理候选。

本轮在此停止，等待人工审核。该下一步尚未启动；不自动进入 C7 stability、C9、larger specimen、δ sweep、phase-gradient、材料扩展、侵蚀或 CFD。

```text
MODEL: PFPD_AT1_SPECTRAL_CANDIDATE_V1
BASELINE_ACTIVATION_OBJECTIVITY: NOT_ESTABLISHED
AFFINE_HOTSPOT: ABSENT
EQUILIBRIUM_HOTSPOT: PRESENT
PROFILE_SCALING: PHYSICAL_HORIZON_SUPPORTED
SUPPORT_ASSOCIATION: MIXED
CLAMP_ASSOCIATION: WEAK
ROOT_CAUSE_PRIMARY: INCONCLUSIVE
CANDIDATE_FIX: NONE
CANDIDATE_IMPLEMENTED: NO
C5_C7_ACTIVATION_DISCREPANCY_BASELINE: 0.3561974074803402
C5_C7_ACTIVATION_DISCREPANCY_CANDIDATE: N/A
CANDIDATE_EFFECT: NOT_RUN
PHASE_SUPPORT_LOCALIZATION: NOT_ASSESSED
C7_STABILITY: NOT_RUN
C9: NOT_RUN
FINAL_STATUS: DIAGNOSIS_INCONCLUSIVE
```
