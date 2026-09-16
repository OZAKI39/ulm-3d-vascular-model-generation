# 当前科研状态

本表汇总到 2026-09-16。各阶段文件里的 NEXT_STAGE/PENDING 是当时快照；继续工作时以本表和 CURRENT_STATE.json 为入口，勿把旧计划当作新授权。

|阶段|状态|证据|实际含义|
|---|---|---|---|
|sonovue|PASS|[报告](sonovue/SONOVUE_SAMPLER_REPORT.md)|连续 empirical CDF / inverse-CDF；数目加权，未按 2.51 µm 重加权。|
|particle_engine|PASS|[报告](particle_engine/MINIMAL_LAMMPS_PARTICLE_ENGINE_REPORT.md)|最小技术引擎，技术接触不代表 SonoVue 接触材料参数已验证。|
|coupling|PASS|[报告](coupling/PALABOS_LAMMPS_COUPLING_V0_REPORT.md)|冻结流场单向读取；场来自 PBS/BSA 5000 步 smoke，仍是工程瞬态场。|
|passive_transport|PASS|[报告](passive_transport/PASSIVE_MICROBUBBLE_TRANSPORT_V0_REPORT.md)|overdamped RK2；真实 8 泡算例在粒子接近后触发 overlap safety stop。|
|migration_2026_rejected|FAIL|[报告](migration_2026_rejected/LAMMPS_2026_MIGRATION_AND_BUBBLE_BUBBLE_INTERACTION_V0_REPORT.md)|2026 release candidate 未通过源代码修正门槛；未提升为 active 安装。|
|rigid_pair|PASS_WITH_TWIST_PENDING|[报告](rigid_pair/STABLE_LAMMPS_RIGID_SPHERE_NEAR_FIELD_V0_REPORT.md)|稳定版刚性 no-slip sphere 近场；normal/shear/transverse rotation/TR 通过，twist 是 blocker。|
|wall_reference|PASS_WITH_LIMITATIONS|[报告](wall_reference/WALL_HYDRODYNAMICS_REFERENCE_AUDIT_REPORT.md)|选择 RMBW Lubrication；有限间隙 RR/TR 精度及远场分支限制保留。|
|wall_v0|BLOCKED_LOCAL_PLANE_VALIDITY|[报告](wall_v0/MICROBUBBLE_WALL_HYDRODYNAMICS_V0_REPORT.md)|核心实现 PASS_WITH_REFERENCE_LIMITATIONS；真实几何仅 4/10000 平面有效点且都在端盖，Case J 时间 0。|
|fixed_multiblob|FAIL_WALL_REPRESENTATION|[报告](fixed_multiblob/FIXED_MULTIBLOB_WALL_FEASIBILITY_AUDIT_REPORT.md)|仅指已测 Pecnut 固定硬球、名义平面及资源范围；不否定所有 regularized blobs。|

未启动曲壁、真实 STL fixed-sphere pilot、adhesion 或新的 dynamics。当前建议仅为 BEM 曲壁可行性审查，等待用户审阅。历史 Stage4 GPU 和 RBC Stage1 状态仍在 [上次交接](../hemocell_handoff_20260915/CURRENT_PROJECT_STATUS.md)；本次不更改这些基线。
