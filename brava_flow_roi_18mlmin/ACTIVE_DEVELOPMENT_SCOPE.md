# BraVa 当前开发入口

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
