# 正文主图与建议附录图索引

正文共 16 张主图：14 张历史原图逐字节复制，2 张本轮整理图。原图数据不变；本轮没有启动新的轨迹计算。英文历史图题和阶段标识保留，正文使用中文学位论文式图注。原图记号与论文换记关系见符号表。

## 正文主图

| 图 / 节 | 图注标题与阅读重点 | 证据类别 | 历史来源 / 整理方式 |
|---|---|---|---|
| [图 1](figures/fig_01.png) / 第 1 节 | 从冻结血流到完整微泡轨迹累积的研究路线。先建立局部运动信息和有限尺寸约束，再处理近场阻力与连续入口；红细胞通行及多出口覆盖作为仍需独立验证的问题保留。箭头表示研究依赖，不表示所有物理层已在同一真实悬浮液中完成耦合。 | 方法路线示意 | 本次按既有研究逻辑绘制；非数值结果；本轮绘制，非新数值实验 |
| [图 2](figures/fig_02.png) / 第 2 节 | 冻结三维血管及入口、三个出口和壁面。读图时应先辨认开放端盖与实体壁面的区别，再观察分叉走向；整棵几何的可见性仅说明显示范围完整，并不代表粒子轨迹已经覆盖所有分支。 | 真实几何与边界 | [00_frozen_input_overview.png](../particle_3d/reports/particle0/figures/00_frozen_input_overview.png)；原样复制 |
| [图 3](figures/fig_03.png) / 第 3 节 | 解析线性场的四面体插值验证。解析值与查询值应落在一致性对角线附近，误差图用于观察浮点舍入量级；这一人工验证算例覆盖多种几何位置，不能替代真实血管的流场实验验证。 | 解析线性场验证算例 | [01_affine_field_validation.png](../particle_3d/reports/particle0/figures/01_affine_field_validation.png)；原样复制 |
| [图 4](figures/fig_04.png) / 第 4 节 | 真实冻结血管中单微泡从入口邻近位置至出口 02 的验证轨迹。应观察曲线是否沿管腔推进及末端是否到达正式出口；图中完整路径验证了这一具体初始化条件下的输运，尚不构成全尺寸分布或全部分支的覆盖证明。 | 真实几何单微泡验证轨迹 | [05_real_vessel_single_mb_trajectory.png](../particle_3d/reports/particle1/figures/05_real_vessel_single_mb_trajectory.png)；原样复制 |
| [图 5](figures/fig_05.png) / 第 5 节 | 红细胞直径和体积的模型分布及接受样本。应区分潜在分布参数、几何筛选后的样本分布与文献统计锚点；这里展示的是几何验证群体，并非某条血管内已经成功流动的红细胞人口。 | 统计几何模型验证样本 | [01_rbc_diameter_volume_distribution.png](../particle_3d/reports/particle2/figures/01_rbc_diameter_volume_distribution.png)；原样复制 |
| [图 6](figures/fig_06.png) / 第 6 节 | 轴对称扁球体在解析简单剪切流中的姿态验证。应观察不同形状比对应的方向演化与解析周期是否一致；该图属于人工流场验证，说明姿态模型的数学响应，不是红细胞在真实微血管中的变形影像。 | 解析剪切验证算例 | [06_simple_shear_jeffery_orbits.png](../particle_3d/reports/particle2/figures/06_simple_shear_jeffery_orbits.png)；原样复制 |
| [图 7](figures/fig_07.png) / 第 7 节 | 有限尺寸颗粒在人工壁面上的无摩擦接触验证。读者应同时观察法向继续靠近是否被限制、切向位移是否仍然保留；该算例验证几何约束机制，不表示真实血管壁具有刚性且无摩擦的全部生物学性质。 | 人工壁面验证算例 | [04_hard_contact_sliding_validation.png](../particle_3d/reports/particle3/figures/04_hard_contact_sliding_validation.png)；原样复制 |
| [图 8](figures/fig_08.png) / 第 8 节 | 双球迎面接近的有限尺寸接触验证。应观察球面间隙接近零后的相对法向运动，以及接受状态是否保持无穿透；这是构造的接触算例，不能据此推断真实血液中的碰撞频率或细胞相互作用强度。 | 双球接触验证算例 | [02_mb_mb_head_on_contact.png](../particle_3d/reports/particle4/figures/02_mb_mb_head_on_contact.png)；原样复制 |
| [图 9](figures/fig_09.png) / 第 9 节 | 球—壁法向近场阻力及正则化响应。应观察近场项如何随间隙缩小增强、在外区退出并在下限处封顶；下限以下的静态扫描仅用于比较公式，不属于允许接受的动态状态。 | 近场解析响应验证算例 | [03_sphere_wall_regularized_resistance.png](../particle_3d/reports/particle6_5/figures/03_sphere_wall_regularized_resistance.png)；原样复制 |
| [图 10](figures/fig_10.png) / 第 10 节 | 无量纲间隙与近场激活权重。随间隙由 0.05 减至 0.01，权重平滑升至一；端点零斜率用于减少人为开关效应，并不消除有限元梯度跨单元面的变化。 | 项目激活函数 | [01_nearfield_activation_weight.png](../particle_3d/reports/particle6_5/figures/01_nearfield_activation_weight.png)；原样复制 |
| [图 11](figures/fig_11.png) / 第 10 节 | 真实冻结场中同初态单球的旧模型与连续介质交接模型重放。应观察旧模型进入更小间隙而新模型保留下限附近的间隙；图中对照不能与更早双微泡历史最小值混为同一轨迹。该窗口检验近壁局部行为，没有出口事件。 | 真实场同初态局部重放 | [07_real_fem_subnanometer_replay.png](../particle_3d/reports/particle6_5/figures/07_real_fem_subnanometer_replay.png)；原样复制 |
| [图 12](figures/fig_12.png) / 第 12 节 | 入口流量加权采样验证。应比较采样密度与法向流量权重，而非仅比较点云是否均匀；图中人工解析对照用于核查概率机制，真实入口样本用于核查实际几何上的抽样。候选分布与接纳后的条件分布必须分别理解。 | 人工与真实入口采样验证 | [02_flux_weighted_inlet_sampling.png](../particle_3d/reports/particle7/figures/02_flux_weighted_inlet_sampling.png)；原样复制 |
| [图 13](figures/fig_13.png) / 第 12 节 | 人工开放控制区中连续计算与中断恢复的生命周期一致性。应观察出生、等待、活动和退出记录是否衔接，避免将恢复后的粒子误当作重新抽样；这一控制算例验证持续人口管理，不代表真实微血管通行。 | 人工生命周期重启验证 | [particle8_07_restart_parity.png](../particle_3d/reports/particle8/figures/particle8_07_restart_parity.png)；原样复制 |
| [图 14](figures/fig_14.png) / 第 14 节 | 完整三维血管及轨迹库中的代表性成功路径。应观察每条展示路径是否连接其入口出生与正式出口事件；背景显示全部血管，轨迹覆盖集中于第二出口相关路径。图中只展示代表性路径，完整样本数量由队列统计给出。 | 真实几何完整轨迹库 | [particle8_1_fig_01_full_geometry_trajectories.png](../particle_3d/outputs/particle8_1/figures/particle8_1_fig_01_full_geometry_trajectories.png)；原样复制 |
| [图 15](figures/fig_15.png) / 第 15 节 | 完成微泡轨迹在三个正交投影上的驻留时间累积。亮区表示在相应统计格中累计的驻留时间较大，同时受到路径数量和局部速度影响；未出现亮线的分支不能由背景血管形状补齐。图中密度属于计算轨迹占据量，并非探测到的超声回波强度。 | 真实轨迹驻留时间累积 | [particle8_1_fig_04_cumulative_localization_views.png](../particle_3d/outputs/particle8_1/figures/particle8_1_fig_04_cumulative_localization_views.png)；原样复制 |
| [图 16](figures/fig_16.png) / 第 16 节 | 冻结血流分配与完成微泡轨迹出口比例。出口 01 和 03 具有非零流量，却未出现完成轨迹；出口 02 的成功样本占比为 100%。两组柱的分母不同，图示揭示覆盖不足，不直接给出造成差异的唯一原因。 | 既有真实结果重绘；无新轨迹计算 | [PARTICLE8_1_VALIDATION.json](../particle_3d/outputs/particle8_1/PARTICLE8_1_VALIDATION.json)；本轮绘制，非新数值实验 |

图 1 另有 [矢量 PDF](figures/fig_01.pdf)，图 16 另有 [矢量 PDF](figures/fig_16.pdf)。这两张单图 PDF 由绘图工具导出，不是全文 LaTeX 编译预览。

## 建议附录图

以下列出未放入正文的历史静态图，按原阶段目录便于追溯。建议围绕具体审核问题择用，而非全部插入论文；含 scope、environment、audit、limitations 的图更适合技术补充材料。真实 FEM 与 synthetic 的判别须结合对应原 REVIEW，不凭文件名推断生理证据。

### particle8_1

- [particle8_1_fig_00_overview](../particle_3d/outputs/particle8_1/figures/particle8_1_fig_00_overview.png)
- [particle8_1_fig_02_inlet_sampling_audit](../particle_3d/outputs/particle8_1/figures/particle8_1_fig_02_inlet_sampling_audit.png)
- [particle8_1_fig_03_outlet_colored_ensemble](../particle_3d/outputs/particle8_1/figures/particle8_1_fig_03_outlet_colored_ensemble.png)
- [particle8_1_fig_05_path_statistics](../particle_3d/outputs/particle8_1/figures/particle8_1_fig_05_path_statistics.png)
- [particle8_1_fig_06_representative_journeys](../particle_3d/outputs/particle8_1/figures/particle8_1_fig_06_representative_journeys.png)
- [particle8_1_fig_07_saved_scene_audit](../particle_3d/outputs/particle8_1/figures/particle8_1_fig_07_saved_scene_audit.png)
- [particle8_1_fig_08_exported_animation_storyboard](../particle_3d/outputs/particle8_1/figures/particle8_1_fig_08_exported_animation_storyboard.png)
- [particle8_1_fig_09_scientific_limitations](../particle_3d/outputs/particle8_1/figures/particle8_1_fig_09_scientific_limitations.png)

### particle0

- [02_node_back_sampling_error](../particle_3d/reports/particle0/figures/02_node_back_sampling_error.png)
- [03_shared_face_continuity](../particle_3d/reports/particle0/figures/03_shared_face_continuity.png)
- [04_inside_outside_classification](../particle_3d/reports/particle0/figures/04_inside_outside_classification.png)
- [05_real_flow_velocity_vectors](../particle_3d/reports/particle0/figures/05_real_flow_velocity_vectors.png)
- [06_real_flow_scalar_diagnostics](../particle_3d/reports/particle0/figures/06_real_flow_scalar_diagnostics.png)
- [07_sampling_runtime_smoke](../particle_3d/reports/particle0/figures/07_sampling_runtime_smoke.png)

### particle1

- [00_particle1_scope_and_provenance](../particle_3d/reports/particle1/figures/00_particle1_scope_and_provenance.png)
- [01_single_mb_size_provenance](../particle_3d/reports/particle1/figures/01_single_mb_size_provenance.png)
- [02_uniform_flow_trajectory](../particle_3d/reports/particle1/figures/02_uniform_flow_trajectory.png)
- [03_pure_rotation_validation](../particle_3d/reports/particle1/figures/03_pure_rotation_validation.png)
- [04_single_step_dataflow](../particle_3d/reports/particle1/figures/04_single_step_dataflow.png)
- [06_real_vessel_single_mb_diagnostics](../particle_3d/reports/particle1/figures/06_real_vessel_single_mb_diagnostics.png)
- [07_validation_timestep_comparison](../particle_3d/reports/particle1/figures/07_validation_timestep_comparison.png)

### particle2

- [00_particle2_scope_and_literature](../particle_3d/reports/particle2/figures/00_particle2_scope_and_literature.png)
- [02_rbc_derived_geometry_distribution](../particle_3d/reports/particle2/figures/02_rbc_derived_geometry_distribution.png)
- [03_rbc_joint_geometry_scatter](../particle_3d/reports/particle2/figures/03_rbc_joint_geometry_scatter.png)
- [04_static_flow_orientation_validation](../particle_3d/reports/particle2/figures/04_static_flow_orientation_validation.png)
- [05_rigid_rotation_orientation_validation](../particle_3d/reports/particle2/figures/05_rigid_rotation_orientation_validation.png)
- [07_distribution_wide_jeffery_validation](../particle_3d/reports/particle2/figures/07_distribution_wide_jeffery_validation.png)
- [08_real_fem_rbc_orientation_trajectory](../particle_3d/reports/particle2/figures/08_real_fem_rbc_orientation_trajectory.png)
- [09_real_fem_orientation_diagnostics](../particle_3d/reports/particle2/figures/09_real_fem_orientation_diagnostics.png)
- [10_particle2_validation_timestep_comparison](../particle_3d/reports/particle2/figures/10_particle2_validation_timestep_comparison.png)

### particle3

- [00_particle3_scope_and_wall_contract](../particle_3d/reports/particle3/figures/00_particle3_scope_and_wall_contract.png)
- [01_wall_triangle_normals](../particle_3d/reports/particle3/figures/01_wall_triangle_normals.png)
- [02_sphere_wall_gap_validation](../particle_3d/reports/particle3/figures/02_sphere_wall_gap_validation.png)
- [03_ellipsoid_wall_gap_validation](../particle_3d/reports/particle3/figures/03_ellipsoid_wall_gap_validation.png)
- [05_physical_time_subdivision](../particle_3d/reports/particle3/figures/05_physical_time_subdivision.png)
- [06_capillary_deformation_surrogate](../particle_3d/reports/particle3/figures/06_capillary_deformation_surrogate.png)
- [07_deformation_feasibility_map](../particle_3d/reports/particle3/figures/07_deformation_feasibility_map.png)
- [08_real_wall_gap_samples](../particle_3d/reports/particle3/figures/08_real_wall_gap_samples.png)
- [09_real_wall_controlled_contact](../particle_3d/reports/particle3/figures/09_real_wall_controlled_contact.png)
- [10_real_fem_mb_wall_clearance](../particle_3d/reports/particle3/figures/10_real_fem_mb_wall_clearance.png)
- [11_real_fem_rbc_wall_clearance](../particle_3d/reports/particle3/figures/11_real_fem_rbc_wall_clearance.png)
- [12_particle3_timestep_comparison](../particle_3d/reports/particle3/figures/12_particle3_timestep_comparison.png)

### particle4

- [00_particle4_scope_and_contact_contract](../particle_3d/reports/particle4/figures/00_particle4_scope_and_contact_contract.png)
- [01_pair_gap_geometry_matrix](../particle_3d/reports/particle4/figures/01_pair_gap_geometry_matrix.png)
- [03_mb_mb_glancing_contact](../particle_3d/reports/particle4/figures/03_mb_mb_glancing_contact.png)
- [04_rbc_mb_contact_geometry](../particle_3d/reports/particle4/figures/04_rbc_mb_contact_geometry.png)
- [05_rbc_rbc_orientation_contact](../particle_3d/reports/particle4/figures/05_rbc_rbc_orientation_contact.png)
- [06_off_center_angular_contact](../particle_3d/reports/particle4/figures/06_off_center_angular_contact.png)
- [07_capsule_pair_contact](../particle_3d/reports/particle4/figures/07_capsule_pair_contact.png)
- [08_pair_physical_time_subdivision](../particle_3d/reports/particle4/figures/08_pair_physical_time_subdivision.png)
- [09_three_particle_simultaneous_contact](../particle_3d/reports/particle4/figures/09_three_particle_simultaneous_contact.png)
- [10_contact_order_invariance](../particle_3d/reports/particle4/figures/10_contact_order_invariance.png)
- [11_wall_plus_particle_multicontact](../particle_3d/reports/particle4/figures/11_wall_plus_particle_multicontact.png)
- [12_all_pairs_vs_broadphase](../particle_3d/reports/particle4/figures/12_all_pairs_vs_broadphase.png)
- [13_mixed_four_particle_contact_trajectory](../particle_3d/reports/particle4/figures/13_mixed_four_particle_contact_trajectory.png)
- [14_real_fem_two_mb_contact_smoke](../particle_3d/reports/particle4/figures/14_real_fem_two_mb_contact_smoke.png)
- [15_particle4_timestep_comparison](../particle_3d/reports/particle4/figures/15_particle4_timestep_comparison.png)

### particle5

- [00_particle5_scope_and_resistance_contract](../particle_3d/reports/particle5/figures/00_particle5_scope_and_resistance_contract.png)
- [01_isolated_sphere_stokes_resistance](../particle_3d/reports/particle5/figures/01_isolated_sphere_stokes_resistance.png)
- [02_sphere_wall_lubrication](../particle_3d/reports/particle5/figures/02_sphere_wall_lubrication.png)
- [03_sphere_sphere_lubrication](../particle_3d/reports/particle5/figures/03_sphere_sphere_lubrication.png)
- [04_rbc_mb_resistance_symmetry](../particle_3d/reports/particle5/figures/04_rbc_mb_resistance_symmetry.png)
- [05_sparse_resistance_matrix](../particle_3d/reports/particle5/figures/05_sparse_resistance_matrix.png)
- [06_resistance_dissipation_validation](../particle_3d/reports/particle5/figures/06_resistance_dissipation_validation.png)
- [07_resistance_metric_contact](../particle_3d/reports/particle5/figures/07_resistance_metric_contact.png)
- [08_lubrication_transit_delay](../particle_3d/reports/particle5/figures/08_lubrication_transit_delay.png)
- [09_real_fem_near_wall_mb_lubrication](../particle_3d/reports/particle5/figures/09_real_fem_near_wall_mb_lubrication.png)
- [10_real_fem_two_mb_resistance_smoke](../particle_3d/reports/particle5/figures/10_real_fem_two_mb_resistance_smoke.png)
- [11_particle5_timestep_comparison](../particle_3d/reports/particle5/figures/11_particle5_timestep_comparison.png)

### particle6

- [00_particle6_scope_and_lammps_environment](../particle_3d/reports/particle6/figures/00_particle6_scope_and_lammps_environment.png)
- [01_lammps_state_roundtrip](../particle_3d/reports/particle6/figures/01_lammps_state_roundtrip.png)
- [02_lammps_neighbor_equivalence](../particle_3d/reports/particle6/figures/02_lammps_neighbor_equivalence.png)
- [03_one_step_standalone_vs_lammps](../particle_3d/reports/particle6/figures/03_one_step_standalone_vs_lammps.png)
- [04_multistep_bridge_parity](../particle_3d/reports/particle6/figures/04_multistep_bridge_parity.png)
- [05_resistance_system_parity](../particle_3d/reports/particle6/figures/05_resistance_system_parity.png)
- [06_lammps_zero_force_audit](../particle_3d/reports/particle6/figures/06_lammps_zero_force_audit.png)
- [07_checkpoint_restart_parity](../particle_3d/reports/particle6/figures/07_checkpoint_restart_parity.png)
- [08_shape_metadata_restart](../particle_3d/reports/particle6/figures/08_shape_metadata_restart.png)
- [09_real_fem_standalone_vs_lammps](../particle_3d/reports/particle6/figures/09_real_fem_standalone_vs_lammps.png)
- [10_neighbor_rebuild_stress](../particle_3d/reports/particle6/figures/10_neighbor_rebuild_stress.png)

### particle6_5

- [00_particle6_5_scope_and_literature](../particle_3d/reports/particle6_5/figures/00_particle6_5_scope_and_literature.png)
- [02_h_lower_vs_particle_radius](../particle_3d/reports/particle6_5/figures/02_h_lower_vs_particle_radius.png)
- [04_sphere_sphere_regularized_resistance](../particle_3d/reports/particle6_5/figures/04_sphere_sphere_regularized_resistance.png)
- [05_sphere_wall_continuum_handoff](../particle_3d/reports/particle6_5/figures/05_sphere_wall_continuum_handoff.png)
- [06_sphere_sphere_continuum_handoff](../particle_3d/reports/particle6_5/figures/06_sphere_sphere_continuum_handoff.png)
- [08_handoff_scale_sensitivity](../particle_3d/reports/particle6_5/figures/08_handoff_scale_sensitivity.png)
- [09_handoff_physical_time_refinement](../particle_3d/reports/particle6_5/figures/09_handoff_physical_time_refinement.png)
- [10_particle6_5_lammps_bridge_parity](../particle_3d/reports/particle6_5/figures/10_particle6_5_lammps_bridge_parity.png)
- [11_particle6_5_timestep_comparison](../particle_3d/reports/particle6_5/figures/11_particle6_5_timestep_comparison.png)

### particle7

- [00_particle7_scope_and_literature](../particle_3d/reports/particle7/figures/00_particle7_scope_and_literature.png)
- [01_frozen_inlet_flux_map](../particle_3d/reports/particle7/figures/01_frozen_inlet_flux_map.png)
- [03_mb_cumulative_flux_scheduler](../particle_3d/reports/particle7/figures/03_mb_cumulative_flux_scheduler.png)
- [04_rbc_volume_flux_scheduler](../particle_3d/reports/particle7/figures/04_rbc_volume_flux_scheduler.png)
- [05_rbc_isotropic_orientation_validation](../particle_3d/reports/particle7/figures/05_rbc_isotropic_orientation_validation.png)
- [06_finite_size_inlet_admission](../particle_3d/reports/particle7/figures/06_finite_size_inlet_admission.png)
- [07_pending_injection_queue](../particle_3d/reports/particle7/figures/07_pending_injection_queue.png)
- [08_injection_timestep_independence](../particle_3d/reports/particle7/figures/08_injection_timestep_independence.png)
- [09_particle_lifecycle_and_outlet_removal](../particle_3d/reports/particle7/figures/09_particle_lifecycle_and_outlet_removal.png)
- [10_injection_checkpoint_restart](../particle_3d/reports/particle7/figures/10_injection_checkpoint_restart.png)
- [11_long_injection_population_balance](../particle_3d/reports/particle7/figures/11_long_injection_population_balance.png)
- [12_scheduled_vs_admitted_distributions](../particle_3d/reports/particle7/figures/12_scheduled_vs_admitted_distributions.png)
- [13_real_frozen_inlet_injection_smoke](../particle_3d/reports/particle7/figures/13_real_frozen_inlet_injection_smoke.png)
- [14_local_tube_hematocrit_diagnostic](../particle_3d/reports/particle7/figures/14_local_tube_hematocrit_diagnostic.png)

### particle8

- [particle8_00_overview](../particle_3d/reports/particle8/figures/particle8_00_overview.png)
- [particle8_01_real_inlet_flux](../particle_3d/reports/particle8/figures/particle8_01_real_inlet_flux.png)
- [particle8_02_rbc_volume_scheduler](../particle_3d/reports/particle8/figures/particle8_02_rbc_volume_scheduler.png)
- [particle8_03_mb_scheduler](../particle_3d/reports/particle8/figures/particle8_03_mb_scheduler.png)
- [particle8_04_real_mixed_smoke](../particle_3d/reports/particle8/figures/particle8_04_real_mixed_smoke.png)
- [particle8_05_synthetic_lifecycle](../particle_3d/reports/particle8/figures/particle8_05_synthetic_lifecycle.png)
- [particle8_06_pending_identity](../particle_3d/reports/particle8/figures/particle8_06_pending_identity.png)
- [particle8_08_animation_storyboard](../particle_3d/reports/particle8/figures/particle8_08_animation_storyboard.png)
- [particle8_09_limitations](../particle_3d/reports/particle8/figures/particle8_09_limitations.png)

## 动画补充材料

动画不计入正文 16 张静态图，也不因显示帧数增加轨迹样本数。真实场景、人工控制及真实未通行的混合接纳场景应保持原分类。
- [完整轨迹、累积和出口着色动画总览](../particle_3d/outputs/particle8_1/index.html)
- [此前完整几何与三维旋转回放总览](../particle_3d/outputs/particle8_full3d/index.html)
