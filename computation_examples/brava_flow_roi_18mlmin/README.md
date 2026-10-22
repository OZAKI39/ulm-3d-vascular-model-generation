# BraVa CFD → microbubble trajectory correctness

2026-10-06 重新执行结果：[中文主报告](population_v1/analysis/p2_linear_algebra_rerun_20261006/P2_LINEAR_ALGEBRA_AND_RESOURCE_REVIEW_ZH.md)与[三张300dpi图](population_v1/analysis/p2_linear_algebra_rerun_20261006/OPEN_RESULTS.html)。新候选保留GPU主系统，将Schur/GAMG预条件器放到CPU，P1同检查点回归通过（124.4秒、Δu约1.86e−16 m/s）。一次P2真实系统setup因任务RSS达到25.22 GiB触发预声明24 GiB保护停止；无主机OOM、无CUDA分配失败。瞬时保守余量9.90 GiB仍高于所需8.41 GiB，但内存持续增长且setup未完成，不能称安全门通过。最终 **NO / STOP_NATIVE_P2P1_ROUTE_RESOURCE_LIMIT**，只限定本轮资源保护范围；完整P2 CFD、BraVa、轨迹和群体均0。曲面压力回填gate仍FAIL。所有本轮任务退出，不push；下方上一轮结果保留为历史。

2026-10-06 最新线性代数审计：[主报告](population_v1/analysis/p2_linear_algebra_review_20261006/P2_LINEAR_ALGEBRA_AND_RESOURCE_REVIEW_ZH.md)与[三张300dpi核查图](population_v1/analysis/p2_linear_algebra_review_20261006/OPEN_RESULTS.html)。223,119个压力边占位已定位；原生P2实测111.79M结构项，ASM2子域NNZ重复1.85748倍。ASM0和Schur/GAMG均通过P1单步解等价回归；选定Schur在P2真实装配后、selfp CUDA临时工作区分配失败，尚无线性迭代和完整PC峰值。另独立复现原生边压力回填的曲面反解异常。最终 **P2_PIPE_LINEAR_ALGEBRA_FEASIBLE=NO / STOP_NATIVE_P2P1_ROUTE_RESOURCE_LIMIT**；一次非预期CUDA分配失败已保留并补充拒绝未知显存预算的入口保护，无重试，host OOM=0。146项工程检查通过；完整P2 CFD/BraVa/trajectory/population均0。所有自身任务退出，不自动继续、不push。下方旧轮次仅为历史记录。

2026-10-06 最新目标BC审计：[主报告](population_v1/analysis/target_bc_pipe_review_20261006/P2_TARGET_BC_AND_PIPE_CFD_REVIEW_ZH.md)与[两张300dpi核查图](population_v1/analysis/target_bc_pipe_review_20261006/OPEN_RESULTS.html)。实际XML经原生parser/setup和mixed重放确认：pipe与历史BraVa全部BC的weak装配调用0；上一轮体装配修复继续有效，无需补丁未使用的weak路径。唯一完整曲面elevation的174,785个TET10均通过正Jacobian证书，九项实现前置条件通过。固定8-rank GPU ASM/LU路线的资源预估超过当前cgroup，按本轮第9.3节运行前停止，不做OOM试错。最终 **P2_PIPE_NOT_EVALUABLE**，主建议 **P2P1_PIPE_RESULT_NOT_EVALUABLE**；新增CFD/trajectory/population均0，无actual P2误差改善结论、无自动后续。下方旧轮次停止规则仅为历史记录。

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

2026-10-06：[stabilization/局部守恒根因审查](population_v1/analysis/stabilization_local_conservation_review_20261006/STABILIZATION_AND_LOCAL_CONSERVATION_REVIEW_ZH.md)完成；[两张300 dpi核查图](population_v1/analysis/stabilization_local_conservation_review_20261006/OPEN_RESULTS.html)。实际生产源码身份已建立；1189 A时间份额0.594%，822 A/B为0.158/0.281%，触发Stop C，新增CFD/solver patch均为0。原468单元RT对照能强制守恒但大幅改变近壁速度，仅DIAGNOSTIC_ONLY，不替换production。现有Taylor–Hood流体链可用于受控原型，当前P1采样器仍不支持高阶DOF。本轮12项定向测试通过；所有自身作业收尾，不自动启动下一路线。

2026-10-06 原生 P2/P1 更新：[开发审查](population_v1/analysis/native_p2p1_review_20261006/NATIVE_P2P1_DEVELOPMENT_REVIEW_ZH.md)与[三张300dpi核查图](population_v1/analysis/native_p2p1_review_20261006/OPEN_RESULTS.html)已完成。A/B/C/D1通过：原生混合单元可执行，二次场与冻结V2的960原探针+4096随机点达到舍入误差表示精度；真实138584有1个自由边中点，但G仍受rank≤1限制。D2直接调用production construct_fluid，确认合法曲面单元的压力/速度映射及首个连续性积分Jacobian不一致，按Phase E 9.1运行前门槛停止。最终 **P2P1_NOT_EVALUABLE_DUE_TO_IMPLEMENTATION_BLOCKER**；58项工程检查通过不等于科学认证。新增CFD/粒子轨迹/完整曲面网格均0，生产源码、历史场与旧1565结果不变。仅保留validation-only CPU sampler，无GPU粒子迁移、无自动后续。

2026-10-06 common curved geometry 更新：[修复与验证总报告](population_v1/analysis/common_curved_geometry_review_20261006/COMMON_CURVED_GEOMETRY_MAPPING_REVIEW_ZH.md)及[三张300dpi图](population_v1/analysis/common_curved_geometry_review_20261006/OPEN_RESULTS.html)。独立、编译门控的 native mixed TET10 体装配补丁通过684项检查；四个P1/P1局部回归的R/K逐位不变。68项边界工程检查通过，但weak Dirichlet仍用四角点几何反解pressure：两个curved wall实际抛出原生异常，curved cap参考坐标不一致。按Phase F停止，最终 **COMMON_GEOMETRY_FIX_NOT_EVALUABLE**，不建议下一轮CFD；759项工程检查与10项边界语义/反解gate失败分别保存。初版Hessian失败及修正记录保留，容差未放宽。生产源码/binary、历史P1场和1565结果不变；新增CFD步、linear solve、trajectory、population均0。任务全部退出、无自动后续。
