# Phase 1：现有模块与新 population 架构

本次权威流场：`../cases/balanced_pressure_final/frozen_flow/steady_flow.vtu`，SHA256：`73228d094c253cec75dbb0030476b687171d0d0fef7816a845392dad220285dc`。BraVa 已验收四端口血管芯，入口 18 mL/min；不修改 CFD。

## 执行链

事件寻址 Poisson source → 原始 SonoVue D≤4 μm 条件分布与精确 P1 正向通量位置采样 → 单次真实 WALL/active-pair 准入 → SoA active pool → 唯一 physical clock（birth 精确切分）→ FEM 采样 → P6.5 法向阻力 + P9-A.1 平面壁修正 + sphere-pair 阻力/接触 → 有限三角形连续路径证书 → 首次中心穿越出口移除 → 绝对时间轨迹、完整源账本与守恒账本。

GPU 顺序：先逐项 CPU reference parity，再 tetra walking/FEM、spatial hash pair、真实三角形 BVH、FP64 dynamics/contact，最后多独立 realization batching。CPU 仍负责来源事件、账本和必要的稀有几何/数值回退；回退必须记录，不能宣称全 GPU。

## 必须保留的语义

- 不调用旧 `attempt()` / `retry()`：一个 event 只抽一次 D、一次 x、检查一次，不产生 rejected particle_id。
- 旧 source 的 Network-H0 哈希限制不被 synthetic 标志绕过，新环境分别校验真实网格、流场和边界。
- 原先 SonoVue 使用 D≤4 μm 条件分布，不是把大泡裁到 4 μm。保持原始分布文件、归一化及取样实现。
- 沿用已冻结 wall/pair continuum handoff：a_ref = wall radius 或 pair min radius，h_lower=max(2 nm, 0.001 a_ref)，χ≤0.01 全开启、χ≥0.05 关闭、中间 C1 cubic。pair 阻力半径为 a_i a_j/(a_i+a_j)。该尺度是旧契约的模型交接，不是新增物理涂层或任意 gap floor。
- MB-MB 必须另立正式 contract，并验证后才能启用；不得把 P5 VALIDATION_ONLY 直接当生产。
- 固定球壁面支持静止仅按原 numerical/stationary gate 判断，终态移除，不能宣称生理堵塞/溶解。
- LAMMPS 保留为可选存储/邻居列表对照，正式积分由 particle_3d 负责。
- 旧 BraVa pilot 中出现长期近壁滞留（12 s 未出口），因此不能预先保证每个 source 都在有限观察窗完成。

## 已读源码位置

| 模块 | 关键入口及行号 |
|---|---|
| continuous_infusion.py | `canonical_bytes` L21, `sha256` L25, `positive_integer` L30, `counter_rng` L36, `poisson_interval` L45, `ContinuousInfusionSource` L60, `PopulationLedger` L125 |
| population_inlet_p9a4.py | `load_new_environment` L22, `make_contract` L46, `make_source` L85, `_proposal_shard` L109, `generate` L113 |
| inlet_flux.py | `positive_pieces` L6, `InletFluxSampler` L26, `frozen_boundary_flux` L66 |
| inlet_size_capacity.py | `InletClearanceTree` L15, `compute_inlet_global_size_capacity` L132 |
| injection_admission.py | `event_shape` L16, `FiniteSizeAdmission` L21, `plug_side_wall_margin` L91, `PlugPathAdmission` L97 |
| particle7_lifecycle.py | `PopulationLifecycle` L15 |
| particle7_bridge.py | `DynamicParticleBridge` L7 |
| field.py | `FlowSample` L8, `FlowBatch` L19, `FrozenFEMField` L30 |
| wall_geometry.py | `WallGeometry` L11 |
| wall_gap.py | `WallGapResult` L9, `wall_gap` L24, `touching_contacts` L62 |
| pair_broadphase.py | `all_pairs` L7, `broadphase` L11 |
| pair_geometry.py | `PairGeometryError` L19, `ParticlePairGapResult` L24, `_simplex` L43, `_gjk` L64, `_components` L83, `_stationary` L92, `_witness` L135, `_ordered_gap` L155, `pair_gap` L200 |
| hydrodynamic_resistance.py | `NonsphericalLubricationNotFrozen` L15, `require_sphere` L23, `viscosity_from_frozen` L28, `sphere_self_diagonal` L47, `PhysicalNearField` L56, `evaluate_block` L70 |
| resistance_assembly.py | `ResistanceSystem` L9, `assemble_resistance_system` L28 |
| resistance_solver.py | `ResistanceSystemIllConditioned` L12, `ResistanceSolution` L19, `contact_jacobian` L25, `solve_resistance` L43 |
| particle6_stepper.py | `bind_query_dependency` L21, `Particle6Stepper` L28 |
| particle65_motion.py | `assemble_v1` L25, `v1_trial` L45, `Particle65Stepper` L112 |
| particle9a_motion.py | `PlanarResistanceSystem` L14, `augment_planar_system` L30, `Particle9AStepper` L54 |
| formal_dynamics_p9a5.py | `safe` L20, `nonfinite` L26, `leaves` L32, `AuditSink` L43, `array_sha` L74, `audit_result` L79, `job` L124 |
| lammps_bridge.py | `new_lammps` L9, `validate_command` L21, `LammpsParticleBridge` L33 |
| lammps_neighbors.py | `ValidationNeighborPolicy` L11, `canonical_pairs` L27, `standalone_candidates` L31, `exact_physics_pairs` L38 |
| lammps_state.py | `BridgeParticle` L18, `schema_contract` L72 |
| nearfield_regularization.py | `NearFieldRegularizationV1` L21, `contract` L92 |
| nearfield_handoff.py | `InitialBelowContinuumHandoff` L16, `exact_interactions` L23, `require_admissible_initial` L36, `handoff_constraints` L44, `pair_handoff_certificate` L71, `wall_handoff_certificate` L80, `first_handoff_event` L104, `partitioned_wall_handoff_certificate` L174 |
| planar_wall_hydrodynamics.py | `stable_wall_frame` L17, `local_shear_vector` L27, `wall_weight` L33, `legacy_planar_coefficients_si` L41, `local_mobility_4d` L60, `planar_wall_affine_block` L72 |

完整路径、源哈希和符号范围见 `PHASE1_SOURCE_REVIEW.json`。这份文档是代码审查结果，不代表后续 parity 或生产计算已经通过。
