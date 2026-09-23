# 残差图科学表达修订与原始日志核验

本次只读取现有求解日志并更新残差图、说明和 PPT。原速度、压力、WSS、网格及三段固定轴旋转视频保持不变。

## 主图与辅助图

- [非线性双面板主图 PNG](figures/residual_nonlinear.png) / [PDF](figures/residual_nonlinear.pdf) / [SVG](figures/residual_nonlinear.svg)：A 仅展示每步初始残差随全程演化；B 展示每步内部残差下降的数量级 `log10(start/last)`，柱高 10 即该步下降 10 个数量级。
- [逐步首末配对图 PNG](figures/residual_nonlinear_pairs.png) / [PDF](figures/residual_nonlinear_pairs.pdf) / [SVG](figures/residual_nonlinear_pairs.svg)：每条竖线只连接同一步的第一条与最后一条记录，不连接不同步的末值。
- [四面板诊断总图](figures/residual_convergence.png)：A/B 为上述两个不同过程，C/D 为原始线性求解结束残差与代表性 GMRES 历史。

图中统一使用白底和英文。不存在两条首末残差包络同时连续连接的旧画法。

## R 的含义与固定参考值

此处 R 是本次求解器用于非线性收敛判断的**缩放线性系统初始残差范数**，对应每次非线性修正关联的 KSP 历史第一项。`petsc_impl.cpp` 将它传入 `FSILS.RI.iNorm`；`Integrator.cpp` 使用它判断非线性停止，`output.cpp` 将相对量打印成 NS 行的 `Ri/R0` 和 `Ri/R1`。它不是分别报告的动量残差或质量残差，也不是以 Pa 表示的压力误差。

当前是 fresh run（`Continue_previous_simulation=false`）。全程参考量只在第一次确定：

`R_ref = R^(1,1) = 5.714650454390e-03`，其中非线性记录编号 k 从 1 开始。

所有时间步均使用同一个 `R_ref`。因此 Panel A 的每步第一条记录无需等于 1；如果除以该步自己的第一条记录，才应等于 1。Panel B 使用同一步的首末比值，公共参考量会抵消。

新图采用 KSP 监测输出的高精度初始范数重构比值。全部 167 条结果均与原 NS 文本的 `Ri/R0`、`Ri/R1` 在四位有效数字的打印误差内一致；最大相对差约 3.466e-04。没有改变原日志，也没有以平滑或插值替代迭代数据。

**最后记录的时点：**记录在本次线性修正执行后打印，但其中用于非线性判断的 R 来自本次修正开始时。它是最后一条已有非线性记录，不是最后一次更新后额外重新装配并评估的残差。图中使用 `last`，避免把它表述为新增评估得到的精确终态误差。

## 实际停止条件与容差线

配置为最少 2 次、最多 12 次非线性迭代，非线性容差为 **1e-10**。源码停止规则为：达到最大迭代次数，或者在满足最少次数后，以下任一条件成立：

1. `R/R_ref ≤ 1e-10`：相对全程固定参考的条件。
2. `R/R_step_start ≤ 1e-10`：相对本步初始参考的条件。

图中的虚线是第一项的真实非线性容差，**不是浮点下限**，也不是 PETSc 线性绝对容差。全 71 步都在第一次满足上述收敛条件的记录处结束；没有任何一步触及 12 次上限。所有末记录均满足全程参考条件，其中 19 步同时满足步内相对条件。

后期每步开始时已经低于全程容差，但仍必须执行最少 2 次迭代。因此 Panel B 后期下降量很小不表示该步未收敛，也不能要求每一步都再下降 10 个数量级。

## 第 4、23 步末值跳高的核验

| 时间步 | 非线性记录数 | 起始 R/R_ref | 最后 R/R_ref | 步内下降数量级 |
|---:|---:|---:|---:|---:|
| 3 | 4 | 1.2612e+00 | 5.8092e-14 | 13.337 |
| 4 | 3 | 1.0931e+00 | 7.7565e-11 | 10.149 |
| 10 | 3 | 9.2544e-02 | 2.7708e-12 | 10.524 |
| 20 | 3 | 3.0050e-04 | 3.6630e-14 | 9.914 |
| 22 | 3 | 8.7386e-05 | 3.5991e-14 | 9.385 |
| 23 | 2 | 4.6794e-05 | 8.5858e-11 | 5.736 |
| 50 | 2 | 7.6468e-13 | 3.1735e-14 | 1.382 |
| 71 | 2 | 4.7327e-14 | 3.1046e-14 | 0.183 |

**第 3 → 4 步：**第 3 步在第 3 条记录时，全程相对残差约为 1.412e-10，步内相对残差约为 1.120e-10，均未达到阈值，必须继续第 4 次。第 4 步在第 3 条记录时，这两项已分别约为 7.756e-11 和 7.096e-11，因此在 3 次后停止。末记录高于前一步的末记录，是少进行一次迭代后的正常停止结果。

**第 22 → 23 步：**第 22 步的第 2 条记录仍为 1.594e-10（全程参考）和 1.824e-6（步内参考），均未达到阈值，因此继续第 3 次。第 23 步的第 2 条记录为 8.586e-11（全程参考）和 1.835e-6（步内参考），全程条件已经满足，且达到最少 2 次，因此正常停止。不能将此误写成“第 23 步满足了 1e-10 的步内相对下降”。

这些步骤的 `SV13Q_BEGIN` 均记录 `reuse=1`、`rebuild_reason=REUSE`，没有发生预条件器重建；第 4、23 步的线性求解全部为 `CONVERGED_RTOL`，没有失败或恢复重试。残差记录定义在这些点也没有改变。

证据：[带原始行号的日志摘录](data/nonlinear_transition_log_excerpt.txt)、[逐非线性记录 CSV](data/nonlinear_iteration_audit.csv)、[逐步 CSV](data/nonlinear_step_audit.csv)、[跳变核验 CSV](data/nonlinear_transition_audit.csv)、[完整机器核验 JSON](NONLINEAR_RESIDUAL_AUDIT.json)。

## 低残差平台的正确解释

约 1e-14 的平台只标注为**观察到的平台**。当前没有通过高精度计算、缩放对比或其他独立实验确定其具体来源，因此不将它命名为已证实的 round-off floor、归一化下限或双精度极限。源码的近零替换逻辑也没有被本日志的初始 KSP 范数触发。

图注使用：**A lower nonlinear residual indicates better satisfaction of the discretized nonlinear equations within that time step.**

残差不是 CFD 解相对真实解的误差估计。很小的代数残差不能证明网格足够细、时间步足够小、边界条件符合真实生理，也不能声称解达到 1e-14 的物理精度。本次没有开展新的网格或时间步收敛研究。

## 源码依据与复现

- [Integrator.cpp](../../../vendor/svMultiPhysics_stage_q/Code/Source/solver/Integrator.cpp)：约 947–989 行，固定参考、本步参考与停止条件；约 151–162 行，求解、修正与日志输出的先后顺序。
- [petsc_impl.cpp](../../../vendor/svMultiPhysics_stage_q/Code/Source/solver/petsc_impl.cpp)：约 318–321 行，KSP 历史第一项赋给 `initNorm`。
- [output.cpp](../../../vendor/svMultiPhysics_stage_q/Code/Source/solver/output.cpp)：约 124–126 行，`Ri/R0`、`Ri/R1` 的打印定义。
- [utils.cpp](../../../vendor/svMultiPhysics_stage_q/Code/Source/solver/utils.cpp)：约 141–159 行，近零判定。

源码与日志 SHA-256 均保存在审计 JSON。新增绘图入口为 `scripts/flow_2mmps/residual_figures.py`，读取原日志后即可重绘，不需要 CFD、WSS 或轨迹重算。永久测试位于 `tests/flow_2mmps/test_residual_semantics.py`。
