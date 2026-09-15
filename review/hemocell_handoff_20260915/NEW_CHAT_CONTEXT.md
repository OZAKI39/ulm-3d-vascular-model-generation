# HemoCell 当前科研交接：先读本文件

当前任务链已到 **RBC Stage1 静态几何门禁失败**。`RBC_TIMESTEPS = 0`，没有 RBC runtime failure，也没有 IBM / MPI4 / RBC GPU 的运行正确性证据。下一阶段仅为 **RBC_GEOMETRY_AND_WALL_COMPATIBILITY_STAGE**。

本目录是审查用源码、合同、报告与状态交接。先读 [CURRENT_PROJECT_STATUS.md](CURRENT_PROJECT_STATUS.md) 和 [CURRENT_STATE.json](CURRENT_STATE.json)，再按下述证据定位；历史报告不回写新状态。

## 研究目的

最终对象是 **microbubble**：trajectory、velocity、margination、distance to wall、near-wall residence、adhesion。实验方向为 **PDMS microchannel**。mouse RBC 提供背景流体动力学环境（background hydrodynamic environment），不是最终研究对象。microbubble、adhesion、ultrasound 尚未进入当前计算。

## 冻结的源码与已完成基线

| upstream | commit |
|---|---|
| HemoCell | `5a410848bd5c57d5ae1c171112e78eab4a82e650` |
| original HemoCell Palabos | `05712164d940a42e06afdd705249912fa0c49f14` |
| GPU Palabos | `4127697e90169bbef982295f1d1c933cf6e90caa` |

STEP1 geometry、STEP2 voxelization / four-port mapping、pure-fluid BC、Step3C calibrated inlet、HUMAN_PARAVIEW_REVIEW、PURE_FLUID_BASELINE、GPU Stage4、新 RTX4090 恢复、PBS/BSA numerics smoke 均为 **PASS**。几何与 Step2/3 历史证据在仓库根目录相对路径 `review_bundle/step3_review/`，本次原样保留。

人工 ParaView 的当前 PASS 来自本次用户交接要求；历史 Step3/正式 Step3C 报告在出具时仍写 PENDING。该历史文字保持不变，不能据此替代当前状态，也不能将此 PASS 扩展到 **HUMAN_RBC_REVIEW（仍 PENDING）**。来源与范围见 [STATUS_RECONCILIATION.md](provenance/STATUS_RECONCILIATION.md)。

## Stage4 GPU 与 Step3C

生产 GPU candidate 为 **Stage4**。新 Vast RTX4090 correctness 及 5000-step sustained correctness 均 PASS；5000 步求解速率 **109.77304257495182 steps/s**。原 RTX5090 为约 **125.880863 steps/s**。GPU 源码数学未改变。RTX4090 每个 horizon 只有一次恢复验证，不能声称完成新的重复统计或 10000 步测试。见 [Stage4 原报告](gpu_stage4/original_rtx5090/GPU_STAGE4_REPORT.md) 和 [RTX4090 恢复报告](gpu_stage4/rtx4090_restore/VAST4090_STAGE4_RESTORE_REPORT.md)。

Step3C 的物理 `Qtarget = 2.7369132390905703e-15 m³/s`，固定数值入口倍率 `1.1197286861799598`。倍率补偿离散入口偏差，**不增加物理 Qtarget**。旧介质正式长期验证 PASS，375000 步 AUTO_CONVERGED，最终 `R_inlet = 1.2260984075843895e-07`；见 [正式长验证报告](pure_fluid/STEP3C_FORMAL_GPU_REPORT.md)。**此倍率尚未在新 PBS/BSA 介质下长期重新校准。**

## PBS/BSA 开发介质

1x PBS + 1% BSA，25°C；开发假设 `rho=1000 kg/m³`、`nu=1.0e-6 m²/s`、`tau=1`，冻结 `dx=1.9989918081065344e-7 m`。

`prepare_numerics.py` 是当前 **唯一 dt / outlet rho_LU 生成源**：从输入里的 nu、dx、tau 计算 dt，不读取既有 dt；先将 Pa 转换成 rho_LU，再让 solver 读取生成合同，不在 solver 内另算一套。代码在 [source/new_medium/scripts/prepare_numerics.py](source/new_medium/scripts/prepare_numerics.py)。

| 生成量 | 值 |
|---|---:|
| dt / s | 6.65994708146172e-9 |
| outlet_01 rho_LU | 1.0000484343922278 |
| outlet_02 rho_LU | 1.0004402376508774 |
| outlet_03 rho_LU | 0.9999543772756865 |

[数值合同](new_medium/NEW_MEDIUM_NUMERICS_CONTRACT.json) 与 [500/5000 步 smoke 报告](new_medium/PURE_FLUID_NEW_MEDIUM_SMOKE_REPORT.md) 已归档。smoke PASS；`EXPERIMENTALLY_MEASURED=NO`，`FINAL_EXPERIMENTAL_CONTRACT=PENDING`。smoke 约 163.82 steps/s 使用了不同介质和较低频率流量监测，**不能与上面的 109.77 steps/s 作为同口径 benchmark 比较**。

## MOUSE_RBC_DEV_V0 与失败证据

这是 mouse-like RBC 开发定义：HemoCell native biconcave / RBC_FROM_SPHERE，`RbcHighOrderModel`，`kLink=15, kArea=5, kBend=80, kVolume=20`，interior viscosity **OFF**。mouse 膜力学尚无实验校准。

固定名义目标 **50 µm³**、名义直径约 **6.43 µm**；实际离散体积为 **45.046330078320075 µm³**，不是 50。实际网格 **642 vertices / 1280 triangles**，三个轴向跨度约 **6.429 × 1.886 × 6.429 µm**。这约 −9.907% 的差异来自静态网格定义，不是运行体积漂移。future Hct 为 0%、5%、10% baseline、20%；均未开展 RBC 群体/Hct 验证。

**RBC_ONLY_VALIDATION_STAGE_1 = FAIL；RBC_STAGE1 = FAIL_GEOMETRY_GATE；RBC_GEOMETRIC_FIT = FAIL。** 原生 CPU/MPI 库和网格工具构建 PASS，但没有发布可执行耦合 timestep 的生产 case binary。MPI1 sanity、CASE A/C、MPI4 均未运行，**RBC_TIMESTEPS=0**。合同中的 monitor/output schedule 只是冻结计划，runtime monitors 未实现/未验证。

在冻结 spawn gate 下，**8,974,155 个确定性粗筛姿态，加预定细化**未找到安全放置。最佳失败候选仍有 **154/642 顶点位于 lumen 外**、最大穿透约 **0.597581 µm**、**186 个三角形接触对**。这是有界刚性放置搜索，**不是遍历所有连续姿态的数学不可嵌入证明，也不证明真实可变形 RBC 不能通过血管**。

请读 [RBC Stage1 报告](rbc_stage1/RBC_STAGE1_REPORT.md)、[几何适配报告](rbc_stage1/RBC_GEOMETRY_FIT_REPORT.md)、[冻结合同](rbc_stage1/RBC_STAGE1_CONTRACT.json)。`rbc_stage1/geometry/FROZEN_LUMEN_LU.vtp` 与 `BEST_RIGID_PLACEMENT.vtp` 可在同一 LU 坐标系共同打开；后者是 **INVALID diagnostic pose，绝非可运行 spawn**。

## 第二个 blocker：RBC-wall interaction

**WALL_INTERACTION_IMPLEMENTATION = ABSENT**，限定当前 Guo vascular production setup。HemoCell 库有 `enableBoundaryParticles`、`populateBoundaryParticles`、`applyBoundaryRepulsionForce`，但当前 Guo off-lattice 应用没有启用/提供兼容的壁面粒子集合。**流体 Guo no-slip 不等于 RBC-wall repulsion**。`CASE_B=BLOCKED_WALL_INTERACTION_IMPLEMENTATION`，同时几何门仍失败；见 [wall audit](rbc_stage1/provenance/WALL_INTERACTION_AUDIT.json)。

## NEXT_STAGE：RBC_GEOMETRY_AND_WALL_COMPATIBILITY_STAGE

1. 研究合理的 **upstream loading / inlet injection**，避免强行在狭窄 lumen 内放置未变形 RBC。
2. 为当前 **Guo STL wall** 设计并验证兼容的 RBC wall interaction / repulsion。
3. 重新审阅 **nominal 50 µm³ vs actual 45.046 µm³** 的定义。

完成上述工作并重新冻结有效初态与 Stage1 合同后，才能重新进行 Stage1。**不要进入 RBC Stage2、Hct 5/10/20% 或加入 microbubble。** 当前失败合同的 null positions 故意阻止运行。下一阶段不得通过悄悄缩小 RBC、改壁面、改物理参数或放宽门槛把失败改成 PASS。

源码定位、patch 应用范围和未包含依赖见 [source/README.md](source/README.md)；原始大型结果身份见 [OMITTED_LARGE_ARTIFACTS.tsv](OMITTED_LARGE_ARTIFACTS.tsv)。本次交接没有运行求解器、修改正式源工程、创建 PR 或合并分支。
