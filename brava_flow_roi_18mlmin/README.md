# BraVa CFD → microbubble trajectory correctness

从 [ACTIVE_DEVELOPMENT_SCOPE.md](ACTIVE_DEVELOPMENT_SCOPE.md) 开始，结果入口为 [当前报告与三张最终图](production_workflow/OPEN_RESULTS.html)。

已完成 continuous population simulator、GPU particle backend、冻结 Gmsh V2 manufactured reference 和 GPU svMultiPhysics steady field。当前科学状态：

**REFERENCE_PIPE_FIELD_FAIL_SCIENTIFIC + MIXED_FIELD_ERROR + P1_FIELD_TRAJECTORY_IMPACT_LARGE（DIAGNOSTIC_ONLY）**。

冻结 Gmsh V2 上的 GPU svMultiPhysics/PETSc 已达到 steady，正常保存第 71 步；这证明执行/稳态门槛通过，不代表局部场精度通过。0.377/149.195 µm 有 P1 representation limitation 证据，29.839/59.678 µm 更支持 CFD−P1 discrete-solution difference。误差向量有负 cross-term，不能直接相减百分数。

四组 exact-vs-same-mesh-P1 诊断中，A 的局部 10 µm transit 差为 5.148892%；四组均无径向分离、gap/contact 类别改变。它们不是正式 CFD-vs-oracle certification，也不是 dt convergence。

V3 = PROPOSAL_ONLY / NOT_RUN；正式 CFD-vs-oracle、0.5/0.25/0.125 ms dt convergence、新 BraVa CFD、新 8/46 human-track validation、新 population 均 NOT_RUN。

## 生产结果与协议

BraVa BG001 RMCA 四端口血管芯，入口18mL/min，rho=1056kg/m³、mu=0.00345312Pa·s；`cases/balanced_pressure_final/` 是保留的历史校准 pressure-outlet 解，三个出口约各1/3。candidate0/15 仅为不同打印姿态。以后新边界设计使用中心线阻力估算→出口压力。当前场的近壁/局部精度风险见 [CFD wall audit](population_v1/cfd_wall_audit/WALL_CFD_AUDIT_REPORT_ZH.md)，全局守恒不能证明局部正确。

正式旧群体结果 `population_v1/runs/realization_0000_contactfix/` 共1565条：458/527/437出流、141观察截止、2数值失败。当前固定灌注协议为1.240815660s灌注+10s观察，dt=0.5ms、浓度4e9m⁻³，尚未运行该协议新正式大批量 BraVa 结果。[协议报告](population_v1/fixed_window_protocol/FIXED_WINDOW_PROTOCOL_REPORT_ZH.md) · [旧群体报告](population_v1/POPULATION_VALIDATION_REPORT_ZH.md)。

[最新黑底蓝泡英文四视角](population_v1/candidate_0/views_0to1p2s_large_axes/OPEN_RESULTS.html) 沿用旧群体数据和逐帧低速显示筛选，不更改轨迹或终态。`microbubble/` 的旧独立1500条任务保持停止，不能与1565条群体结果混淆。

## 文件与复现边界

- `scripts/`、`vendor/`、`gpu_solver_fix/support/`：当前生产适配、解析/稳态/检查点、GPU wrapper与RA配置。没有更改 production physics。
- `production_workflow/`：冻结V2、GPU第71步、field/decomposition/sensitivity机器结果与当前validation源码。保留最终mesh、native checkpoint、FIELD_SAMPLES.npz和三个最终PNG。
- `inputs/`、`mesh/`、`cases/`、`population_v1/`：正式几何/旧场/群体记录。原始大数据仍保留本地和原服务器；具体同步范围见 [GitHub同步清单](GITHUB_SYNC_MANIFEST_ZH.md)。
- [退休路线索引](RETIRED_EXPERIMENTS_INDEX.md)：只用于按需历史审计；不恢复MMG/Omega/refine/MeshSim/全域加密/P2/BC-B。
- [清理计划](SYNC_CLEANUP_PLAN_ZH.md) · [大文件排除清单](EXCLUDED_LARGE_ARTIFACTS.json)。

服务器原位 `vast4090:/workspace/brava_flow_roi_18mlmin_20260928/`。冻结配置保留原绝对路径契约；此同步是可审查的源码/证据交付，并非空白机器一键部署。本轮没有安装环境或运行任何科学计算。
