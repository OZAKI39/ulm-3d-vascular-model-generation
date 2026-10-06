# BraVa 当前开发入口

2026-10-06 重新执行结果：[中文主报告](population_v1/analysis/p2_linear_algebra_rerun_20261006/P2_LINEAR_ALGEBRA_AND_RESOURCE_REVIEW_ZH.md)与[三张300dpi图](population_v1/analysis/p2_linear_algebra_rerun_20261006/OPEN_RESULTS.html)。新候选保留GPU主系统，将Schur/GAMG预条件器放到CPU，P1同检查点回归通过（124.4秒、Δu约1.86e−16 m/s）。一次P2真实系统setup因任务RSS达到25.22 GiB触发预声明24 GiB保护停止；无主机OOM、无CUDA分配失败。瞬时保守余量9.90 GiB仍高于所需8.41 GiB，但内存持续增长且setup未完成，不能称安全门通过。最终 **NO / STOP_NATIVE_P2P1_ROUTE_RESOURCE_LIMIT**，只限定本轮资源保护范围；完整P2 CFD、BraVa、轨迹和群体均0。曲面压力回填gate仍FAIL。所有本轮任务退出，不push；下方上一轮结果保留为历史。

2026-10-06 最新线性代数审计：[主报告](population_v1/analysis/p2_linear_algebra_review_20261006/P2_LINEAR_ALGEBRA_AND_RESOURCE_REVIEW_ZH.md)与[三张300dpi核查图](population_v1/analysis/p2_linear_algebra_review_20261006/OPEN_RESULTS.html)。223,119个压力边占位已定位；原生P2实测111.79M结构项，ASM2子域NNZ重复1.85748倍。ASM0和Schur/GAMG均通过P1单步解等价回归；选定Schur在P2真实装配后、selfp CUDA临时工作区分配失败，尚无线性迭代和完整PC峰值。另独立复现原生边压力回填的曲面反解异常。最终 **P2_PIPE_LINEAR_ALGEBRA_FEASIBLE=NO / STOP_NATIVE_P2P1_ROUTE_RESOURCE_LIMIT**；一次非预期CUDA分配失败已保留并补充拒绝未知显存预算的入口保护，无重试，host OOM=0。146项工程检查通过；完整P2 CFD/BraVa/trajectory/population均0。所有自身任务退出，不自动继续、不push。下方旧轮次仅为历史记录。

2026-10-06 最新目标BC审计：[主报告](population_v1/analysis/target_bc_pipe_review_20261006/P2_TARGET_BC_AND_PIPE_CFD_REVIEW_ZH.md)与[两张300dpi核查图](population_v1/analysis/target_bc_pipe_review_20261006/OPEN_RESULTS.html)。实际XML经原生parser/setup和mixed重放确认：pipe与历史BraVa全部BC的weak装配调用0；上一轮体装配修复继续有效，无需补丁未使用的weak路径。唯一完整曲面elevation的174,785个TET10均通过正Jacobian证书，九项实现前置条件通过。固定8-rank GPU ASM/LU路线的资源预估超过当前cgroup，按本轮第9.3节运行前停止，不做OOM试错。最终 **P2_PIPE_NOT_EVALUABLE**，主建议 **P2P1_PIPE_RESULT_NOT_EVALUABLE**；新增CFD/trajectory/population均0，无actual P2误差改善结论、无自动后续。下方旧轮次停止规则仅为历史记录。

2026-10-04：本轮仅完成清理、证据整理、smoke 与 GitHub 同步。停止科学开发，等待下一条明确指令。

**REFERENCE_PIPE_FIELD_FAIL_SCIENTIFIC + MIXED_FIELD_ERROR + P1_FIELD_TRAJECTORY_IMPACT_LARGE（DIAGNOSTIC_ONLY）**。

冻结 Gmsh V2 上的 GPU svMultiPhysics/PETSc 已达到 steady，正常保存第 71 步；这证明执行/稳态门槛通过，不代表局部场精度通过。0.377/149.195 µm 有 P1 representation limitation 证据，29.839/59.678 µm 更支持 CFD−P1 discrete-solution difference。误差向量有负 cross-term，不能直接相减百分数。

四组 exact-vs-same-mesh-P1 诊断中，A 的局部 10 µm transit 差为 5.148892%；四组均无径向分离、gap/contact 类别改变。它们不是正式 CFD-vs-oracle certification，也不是 dt convergence。

V3 = PROPOSAL_ONLY / NOT_RUN；正式 CFD-vs-oracle、0.5/0.25/0.125 ms dt convergence、新 BraVa CFD、新 8/46 human-track validation、新 population 均 NOT_RUN。

默认只按具体任务读取以下七条主线：

1. **Production microbubble simulator**：共享 `../_shared/particle_workflow/particle_3d/src/particle_3d/`、全部 active tests，本地 `population_v1/scripts/` 与 `population_v1/fixed_window_protocol/`。GitHub 对应 `../ulm_particle_formal_p9a5/particle_3d/`；科学源码逐字节保留。
2. **Frozen Gmsh V2 reference**：`production_workflow/gmsh_reference/reference_v2.py`、V2 plan、`candidate_v2/`、`GMSH_REFERENCE_MESH_V2_RESULT.json`、`CAP_TOPOLOGY_CORRECTION.json`。不继续研究 mesher。
3. **GPU steady field**：`production_workflow/manufactured_pipe/R1_GMSH_V2/cfd_gpu/`、`GPU_MANUFACTURED_RUN_IDENTITY.json`、`FIELD_VALIDATION.json`。当前 authority 是第71步，旧 CPU checkpoint 仅用于历史 restart assessment。保留现有 GPU wrapper、RA、solver parsing/steady/checkpoint/strict serialization 支持。
4. **V2 field-error decomposition**：`production_workflow/V2_P1_REPRESENTATION_DIAGNOSIS_ZH.md`、`V2_FIELD_ERROR_DECOMPOSITION.json`、`gmsh_reference/field_error_decomposition.py`、冻结 NPZ。
5. **P1 particle sensitivity**：`production_workflow/P1_FIELD_PARTICLE_SENSITIVITY.json`、`gmsh_reference/error_diagnosis/` 与 `exact_cylinder_particle.py`。仅 DIAGNOSTIC_ONLY。
6. **Conditional V3 proposal**：`production_workflow/PROPOSED_V3_MESH_POLICY.json`，仅提案，未生成/未运行。
7. **Future trajectory certification**：`production_workflow/PARTICLE_CFD_INPUT_CONTRACT.json`、`PIPE_PARTICLE_TRAJECTORY_VALIDATION_ZH.md`、`PIPE_PARTICLE_TIMESTEP_CONVERGENCE_ZH.md` 和当前快速接口测试。正式轨迹与收敛均未运行。

默认不要递归读取 retired experiments、旧报告/日志/mesher adapters 或清理清单；历史问题只按 [RETIRED_EXPERIMENTS_INDEX.md](RETIRED_EXPERIMENTS_INDEX.md) 定位。已删除 tracked 实现从各自 Git 历史恢复，不新建 archive，不恢复停止队列。

正式 BraVa 仍由 `reports/ACTIVE_FLOW.json` 指向 `cases/balanced_pressure_final/`，入口18mL/min、历史校准出口约各1/3；candidate0/15 是同一解的打印姿态。未来新边界设计使用中心线阻力估算→出口压力，不改写历史场身份。兼容符号链接保留。

旧正式 `population_v1/runs/realization_0000_contactfix/` 接受1565条：O1/O2/O3=458/527/437，141观察截止、2数值失败；不能称全部成功。当前固定协议为灌注1.240815660s+观察10s、dt0.5ms、浓度4e9m⁻³，尚无该协议新正式大批量结果。最新四视角在 `population_v1/candidate_0/views_0to1p2s_large_axes/`，显示筛选不改科学账本。

未创建 TASK.json、未启动任何科学作业。服务器操作先读 `/etc/vast-agents-guide.md`；不改变其它项目服务。

2026-10-05 独立事件诊断更新（上文为 2026-10-04 停止点记录）：已完成[原速度急变原因定位](population_v1/analysis/speed_event_cause_review/CAUSE_REPORT_ZH.md)及[事件区域背景流场核查](population_v1/analysis/field_reliability_event_review/FIELD_RELIABILITY_REPORT_ZH.md)，[五组可视化入口](population_v1/analysis/field_reliability_event_review/OPEN_RESULTS.html)。本轮复用原 M0 和经完整几何、历史 BC、稳态与真实载荷核验的固定 PLC 中/细档解，新增 CFD、网格和微泡积分均为 0。1189 的离壁环境仍有局部缺陷及敏感性；822 的原急增幅度/尺度未获得空间稳定性支持。三例近壁后段快流仅在已检查分辨率下得到特征支持，不是物理验证或全场精度通过。任务已收尾，无自动后续；原 Gmsh V2 科学 FAIL 和上文正式结果身份不变。

2026-10-06：[stabilization/局部守恒根因审查](population_v1/analysis/stabilization_local_conservation_review_20261006/STABILIZATION_AND_LOCAL_CONSERVATION_REVIEW_ZH.md)完成；[两张300 dpi核查图](population_v1/analysis/stabilization_local_conservation_review_20261006/OPEN_RESULTS.html)。实际生产源码身份已建立；1189 A时间份额0.594%，822 A/B为0.158/0.281%，触发Stop C，新增CFD/solver patch均为0。原468单元RT对照能强制守恒但大幅改变近壁速度，仅DIAGNOSTIC_ONLY，不替换production。现有Taylor–Hood流体链可用于受控原型，当前P1采样器仍不支持高阶DOF。本轮12项定向测试通过；所有自身作业收尾，不自动启动下一路线。

2026-10-06 原生 P2/P1 更新：[开发审查](population_v1/analysis/native_p2p1_review_20261006/NATIVE_P2P1_DEVELOPMENT_REVIEW_ZH.md)与[三张300dpi核查图](population_v1/analysis/native_p2p1_review_20261006/OPEN_RESULTS.html)已完成。A/B/C/D1通过：原生混合单元可执行，二次场与冻结V2的960原探针+4096随机点达到舍入误差表示精度；真实138584有1个自由边中点，但G仍受rank≤1限制。D2直接调用production construct_fluid，确认合法曲面单元的压力/速度映射及首个连续性积分Jacobian不一致，按Phase E 9.1运行前门槛停止。最终 **P2P1_NOT_EVALUABLE_DUE_TO_IMPLEMENTATION_BLOCKER**；58项工程检查通过不等于科学认证。新增CFD/粒子轨迹/完整曲面网格均0，生产源码、历史场与旧1565结果不变。仅保留validation-only CPU sampler，无GPU粒子迁移、无自动后续。

2026-10-06 common curved geometry 更新：[修复与验证总报告](population_v1/analysis/common_curved_geometry_review_20261006/COMMON_CURVED_GEOMETRY_MAPPING_REVIEW_ZH.md)及[三张300dpi图](population_v1/analysis/common_curved_geometry_review_20261006/OPEN_RESULTS.html)。独立、编译门控的 native mixed TET10 体装配补丁通过684项检查；四个P1/P1局部回归的R/K逐位不变。68项边界工程检查通过，但weak Dirichlet仍用四角点几何反解pressure：两个curved wall实际抛出原生异常，curved cap参考坐标不一致。按Phase F停止，最终 **COMMON_GEOMETRY_FIX_NOT_EVALUABLE**，不建议下一轮CFD；759项工程检查与10项边界语义/反解gate失败分别保存。初版Hessian失败及修正记录保留，容差未放宽。生产源码/binary、历史P1场和1565结果不变；新增CFD步、linear solve、trajectory、population均0。任务全部退出、无自动后续。
