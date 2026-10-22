# BraVa Domain A：BC-B 边界传递与收敛诊断

主诊断：**TECHNICALLY_INCONCLUSIVE_AFTER_BOUNDED_DIAGNOSIS**。结果状态：**BOUNDARY_SENSITIVITY_NOT_EVALUABLE**。

更换初值缓解了首档非线性停滞，但后续尝试仍受线性子问题未收敛限制。现有证据将范围缩小到 NS 的线性化/迭代求解路径，尚不能区分 Jacobian 病态、预条件器不足或稳态分支问题。

## 七个问题

1. **BC-B 为什么难收敛？** 更换初值缓解了首档非线性停滞，但后续尝试仍受线性子问题未收敛限制。现有证据将范围缩小到 NS 的线性化/迭代求解路径，尚不能区分 Jacobian 病态、预条件器不足或稳态分支问题。

2. **入口 profile 有没有传错？** 没有发现入口 profile 传错。原方法与严格 P2 trace L2 投影的相对差为 4.06e-16，入口流量及壁缘零速度均保持。

3. **三个 cut 的源流量是否自己守恒？** 三个截面不是严格零缺口：入口 17.524174 mL/min，两个出口合计 17.547795 mL/min，相对差 0.1348%。缺口与源 P1 场的积分散度一致，未发现法向或提取错误。

4. **单侧 P1 traction 是否明显不稳定？** 出口单侧牵引确实不唯一，切向分量相对差约 24.3% 和 33.0%；在固定求解压力基准下，总牵引差约 1.02% 和 0.40%。它是需要保留的边界不确定性，但本轮没有证据证明其极度不唯一或是失败主因。

5. **主要是边界数据、传递、初值还是求解器？** 已排除本轮检查范围内明显的入口传错、牵引装配错误和线性 Stokes 不相容。主要困难落在 NS 迭代求解阶段，具体根因仍需区分。

6. **是否值得再做正式 boundary sensitivity？** 目前不值得直接再启动一次正式边界敏感性：尚无合格 BC-B 字段，门槛应标为不可评估。

7. **下一步最小需要改什么？** 最小下一步是对失败那一档的 NS 线性化系统与原预条件器做有界诊断，保持边界与网格固定；本轮只提出建议，不自动改容差、重构牵引或尝试新 CFD 方法。

## 关键证据

| 试验 | 状态 | 最终残差 | 延拓档 | 最后 KSP 原因/次数 |
|---|---|---:|---:|---|
| 原 BC-B（直接复用） | 非线性迭代耗尽 | 0.04123 | 0.25 | 未记录完整值 |
| Stokes（仅诊断） | CONVERGED | 2.344e-12 | 1 | 2 / 157 |
| ns_full_projection | NOT_CONVERGED | 1.982 | 0.5 | -3 / 3000 |
| ns_stokes | NOT_CONVERGED | 1.982 | 0.5 | -3 / 3000 |


ns_full_projection 在 0.25 档最后记录残差 3.40364e-10；整个尝试耗时 9.27 分钟。

退出档位 0.5，最后线性迭代 3000 次、残差 0.0587749，分类 LINEAR_SOLVER_LIMITATION。SNES 原因值若为 0，表示 KSP 抛错时尚未设置最终 SNES 原因，不能解释为收敛。

ns_stokes 在 0.25 档最后记录残差 1.04107e-12；整个尝试耗时 9.43 分钟。

退出档位 0.5，最后线性迭代 3000 次、残差 0.0587749，分类 LINEAR_SOLVER_LIMITATION。SNES 原因值若为 0，表示 KSP 抛错时尚未设置最终 SNES 原因，不能解释为收敛。

旧日志保留了 41 个 SNES 状态（0–40）和每隔 50 次一次的 KSP 样本，不能补造阻尼因子、准确线性迭代次数或 reason。日志显示多次 Newton 更新完成后才在 SNES 层退出，且配置启用了线性失败报错，因此首先分类为 **NONLINEAR_CONVERGENCE_STALL**。这是结合日志过程与配置的诊断，不是仅根据最终残差猜测；也不等于证明线性预条件器最优。[PETSc 的线性失败报错说明](https://petsc.org/release/manualpages/KSP/KSPSetErrorIfNotConverged/)支持该配置解释；[收敛原因定义](https://petsc.org/release/manualpages/SNES/SNESConvergedReason/)用于解释本轮明确记录的原因值。在线文档现为 3.26，本轮执行环境始终为 3.25.5。

流量闭合差 0.1348% 与源场积分散度吻合到 2.81e-23 m³/s。1% 是本轮预先声明的“明显失配”诊断线，不是用户给出的精度标准，也不是物理精度证明；BC-B 只在入口强制速度，出口源流量仅作参考，不能把三个参考 Q 当作同时强制的通量条件。

入口 trace 的严格 L2 投影使用完整 P2 边界质量矩阵。原传递与投影最大矢量差 6.35e-16 m/s；未为匹配 18 mL/min 重新归一化。出口原 DG1 应力载荷与源 P1 牵引在 quadrature 上直接积分，向量相对差最大 2.38e-15，合力及相对 Domain A 中心的力矩一致。两侧牵引跳跃来自速度梯度的单侧不连续，P1 压力在同一面上连续；报告分别保留绝对压力、求解压力基准及切向比较，避免压力常数掩盖切向差异。

同一 BC-A fine 状态下，R_A=1.873e-10，R_IN=9.731，R_OUT=26.27，R_FULL=28.02。出口替换引起更大的代数扰动；这些范数包含与原 SNES 相同的 Dirichlet lifting，**不是物理误差**，除以接近零的 R_A 会产生很大的比例，不能当作场误差百分比。

全域初值经过实际体积质量矩阵 L2 投影，速度/压力相对 L2 误差分别为 1.29e-13/4.86e-14；220 点的源场细网格重建误差见 INITIAL_STATE_PROJECTION.json（历史引用 `INITIAL_STATE_PROJECTION.json`；结论见本报告和最终摘要，旧路径已清理）。原首档初值残差 6.968 降至 1.07，改善 6.51 倍，**未达到预设十倍“强相关”判据**；后续判断依据实际受控 NS 尝试，而非事后修改该判据。

Stokes 已收敛、场值有限、流量缺口 6.51e-13。它只证明这组数据在线性离散问题中可解，是 DIAGNOSTIC_ONLY；不作为 BraVa 科学流场，也不进入 Boundary Gate。它不能证明原非线性稳态解存在、唯一或具有足够的局部物理精度。

## 原冻结 Boundary Gate

**BOUNDARY_SENSITIVITY_NOT_EVALUABLE**。不生成虚假的六量比较，不把未取得收敛场写成 Boundary Gate FAIL。

## 范围、资源与停止

仅复用 Domain A：101,401 tet、447,679 mixed DOF、8 MPI、每进程 1 线程。原 BC-A 未重跑；原 BC-B 失败日志、边界数据、QA 和执行记录全部保留。只读阶段在约 15 分钟内完成。首次装配脚本中的 UFL 常数乘法 API 用法错误已在 CFD 前修正，2 秒退出记录保留；它不是一次 NS 尝试。

| 本轮求解 | 实际墙钟（分钟） | 进程树峰值（GiB） | cgroup 峰值（GiB） |
|---|---:|---:|---:|
| stokes | 0.84 | 8.040 | 15.009 |
| ns_full_projection | 9.27 | 10.487 | 17.674 |
| ns_stokes | 9.43 | 10.423 | 18.965 |


生成报告时本轮累计约 41.6 分钟（含诊断实现、等待、计算和结果整理）。资源值为实际采样峰值，不是事先外推。运行结束由监督程序和 MPI 进程退出检查确认。

Stokes 之后首次 NS 启动被额外 1.5 GiB 预测预留挡住，实际内存约 7.8 GiB，未启动 MPI。核对历史完整 P2-fine 峰值后将额外预测预留设为 1 GiB；**实际绝对 cgroup 20 GiB 停止线从未改变**，未改系统 cgroup 或全局缓存。每个实际 NS 尝试上限仍为 25 分钟，最多两种规定初值，rtol/atol、边界及 0.25/0.5/1 延拓未改。全部执行记录见 `monitoring/`。

本轮不编译 PETSc/DOLFINx，不重复 near-wall、parity、粒子或网格收敛验证。已有 `IMPLEMENTATION_MATCHES_DECLARED_MODEL` 和 `LOCAL_P2_NOT_YET_SPATIALLY_RESOLVED` 继续按原适用范围保留。没有开发 traction recovery、DG、Nektar++，没有启动 B/C、全域 CFD 或新轨迹。

因未取得收敛 NS 字段，本轮没有执行 GPU 字段比较；不为 GPU 利用率增加无必要计算。

四张图文入口（历史引用 `OPEN_RESULTS.html`；结论见本报告和最终摘要，旧路径已清理） · 原失败日志解析（历史引用 `SOLVER_FAILURE_DIAGNOSIS.json`；结论见本报告和最终摘要，旧路径已清理） · 质量检查（历史引用 `MASS_COMPATIBILITY.json`；结论见本报告和最终摘要，旧路径已清理） · 入口传递（历史引用 `INLET_TRANSFER.json`；结论见本报告和最终摘要，旧路径已清理） · 牵引两侧差（历史引用 `TRACTION_SIDE_SENSITIVITY.json`；结论见本报告和最终摘要，旧路径已清理） · 牵引载荷传递（历史引用 `TRACTION_TRANSFER.json`；结论见本报告和最终摘要，旧路径已清理） · 残差分解（历史引用 `RESIDUAL_DECOMPOSITION.json`；结论见本报告和最终摘要，旧路径已清理） · Stokes（历史引用 `STOKES_DIAGNOSIS.json`；结论见本报告和最终摘要，旧路径已清理） · NS 尝试（历史引用 `NS_DIAGNOSIS.json`；结论见本报告和最终摘要，旧路径已清理）

完成本报告后停止所有本轮任务，等待用户决定下一步。

## 退休状态

[最终摘要](FINAL_BC_SUMMARY.json)。原BC-B首次尝试在0.25档迭代40次、残差0.0412322，耗时15.65分钟；后续有界诊断如上。frozen Jacobian只是后续建议，现已退休、未运行，不将计划写成结果。
