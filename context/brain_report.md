<!-- PRODUCTION_WORKFLOW_CURRENT_STATUS_START -->
**2026-10-04 当前：REFERENCE_PIPE_FIELD_FAIL（CFD_RUNTIME_OR_EXECUTION / TIME_LIMIT），已停止。VALIDATION_ONLY=true。**

单次端面core径向3→6修正PASS（速度0.210666%、Q 0.165549%）；唯一V2全部mesh gates PASS，永久冻结。wall altitude24.9983 µm，中央bulk215.4203 µm/P0=1.173873。唯一CFD在1800s硬限内完成15步，只有第10步快照，没有稳态field；局部u/G/s/curl未评估，oracle/trajectory/dt均未运行。

修复了validation reader的float32序列化严格对应与监督器pidfd退出异常。CFD没有重启，输入和结果未改。监督记录显式恢复，原生exit code与完整内存峰值未知。所有本次native进程已停止，TASK.json移除。不要续算旧CFD、改reference mesh、生成V3或恢复其它路线/BraVa。

入口：[ACTIVE_DEVELOPMENT_SCOPE.md](computation_examples/brava_flow_roi_18mlmin/ACTIVE_DEVELOPMENT_SCOPE.md)；[总报告](computation_examples/brava_flow_roi_18mlmin/production_workflow/PRODUCTION_WORKFLOW_REVIEW_ZH.md)；[图表](computation_examples/brava_flow_roi_18mlmin/production_workflow/OPEN_RESULTS.html)。下文仅历史；原production physics、固定协议与旧1565（141截止、2失败）身份不变。0.5 ms仍无canonical convergence结论。
本轮提交：mesh `9c8a7551a8c9`；CFD停止与接口修复 `82dd025593b0`。
<!-- PRODUCTION_WORKFLOW_CURRENT_STATUS_END -->

# BraVa 流场求解与微泡轨迹：代码、数据和报告路径指南

初始索引核对日期：2026-09-29，当时仅整理文档，没有启动数值计算。2026-09-30 补充第 15 节：正在进行的网格、流场与固定队列微泡验证，该阶段已在服务器执行新 CFD 和验证重放。表中本地链接均使用绝对路径；相对名称以各节注明的根目录为准。当前执行状态以第 15 节所列报告与监控记录为准。

## 1. 先分清当前代码与已有结果

当前对象是 **BraVa BG001 RMCA 的 BALANCED 四端口血管芯**，包含人工接管、去除五面盒体。candidate 0 与 candidate 15 是同一几何的两种刚体打印姿态，复用同一个无重力、刚壁、牛顿流体物理解，不能按两套独立 CFD 计数。

| 对象 | 当前状态 | 正确入口 |
|---|---|---|
| 已求解正式流场 | 全压力出口复算已完成；入口 18 mL/min，三出口约各 1/3 | `cases/balanced_pressure_final/`，由 `reports/ACTIVE_FLOW.json` 指定 |
| 已有正式微泡群体 | `realization_0000_contactfix`，1565 条接受轨迹；旧规则为各出口达到至少 20 条后停灌，随后观察 12 s | `population_v1/runs/realization_0000_contactfix/` |
| 当前微泡生产代码 | 已改为固定时间协议：灌注 1.240815660 s，停灌后观察 10 s，结束于 11.240815660 s | `population_v1.run_fixed_window_population` 模块及固定协议合同 |
| 固定时间协议的结果 | 已完成小规模、源事件重放和真实环境零源验证；**尚无该协议的正式大批量 BraVa 轨迹** | `population_v1/fixed_window_protocol/validation/` |
| 最新微泡展示 | candidate 0、物理时间 0–1.2 s、四视角、放大坐标文字；仍读取上述旧协议已算轨迹 | `population_v1/candidate_0/views_0to1p2s_large_axes/` |
| 网格与输运收敛验证 | 另行冻结 46 个出生参数，关闭泡间作用；正在比较新网格与 CFD，尚未选定新正式流场或认证 0.5 ms | `mesh_particle_validation_v2/`，详见第 15 节 |
| 最早独立 1500 条微泡方案 | 已按用户要求停止，只保留局部记录，不能当完成批次 | `microbubble/data/USER_STOP.json` |

已有群体结果的 1565 条不等于 1565 条成功出流：O1/O2/O3 分别为 **458/527/437**，141 条观察截止、2 条数值失败、0 条支持静止。`scientific_integrity_pass=false` 原样保留。新固定窗协议没有回填或改写这批历史结果。

## 2. 常用总入口

- **BraVa 算例根目录**：[/home/lzy/projects/computation_examples/brava_flow_roi_18mlmin](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin>)

- **BraVa 算例身份说明**：[/home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/CASE_README.md](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/CASE_README.md>)

- **正式流场及数据身份**：[/home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/reports/ACTIVE_FLOW.json](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/reports/ACTIVE_FLOW.json>)

- **candidate 0 流场动画**：[/home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/visualization/candidate_0/OPEN_RESULTS.html](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/visualization/candidate_0/OPEN_RESULTS.html>)

- **candidate 15 流场动画**：[/home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/visualization/candidate_15/OPEN_RESULTS.html](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/visualization/candidate_15/OPEN_RESULTS.html>)

- **微泡群体总预览**：[/home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/OPEN_RESULTS.html](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/OPEN_RESULTS.html>)

- **candidate 0 当前微泡预览**：[/home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/candidate_0/OPEN_RESULTS.html](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/candidate_0/OPEN_RESULTS.html>)

- **最新四视角大字号微泡动画**：[/home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/candidate_0/views_0to1p2s_large_axes/OPEN_RESULTS.html](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/candidate_0/views_0to1p2s_large_axes/OPEN_RESULTS.html>)

- **微泡生产验证主报告**：[/home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/POPULATION_VALIDATION_REPORT_ZH.md](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/POPULATION_VALIDATION_REPORT_ZH.md>)

- **当前固定时间协议主报告**：[/home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/fixed_window_protocol/FIXED_WINDOW_PROTOCOL_REPORT_ZH.md](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/fixed_window_protocol/FIXED_WINDOW_PROTOCOL_REPORT_ZH.md>)

- **近壁 CFD 审计主报告**：[/home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/cfd_wall_audit/WALL_CFD_AUDIT_REPORT_ZH.md](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/cfd_wall_audit/WALL_CFD_AUDIT_REPORT_ZH.md>)

- **微泡突然加速原因说明**：[/home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/analysis/acceleration_diagnosis/ACCELERATION_EXPLANATION_ZH.md](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/analysis/acceleration_diagnosis/ACCELERATION_EXPLANATION_ZH.md>)

- **新网格、CFD 与微泡输运验证**：[/home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/mesh_particle_validation_v2/OPEN_RESULTS.html](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/mesh_particle_validation_v2/OPEN_RESULTS.html>)

## 3. 几何、ROI、端口与网格输入

### 3.1 几何上游项目

独立几何/打印工程仍在 [/home/lzy/projects/vascular_printing](</home/lzy/projects/vascular_printing>)，不在 BraVa CFD 目录中重复保存整套原始数据。

| 路径入口 | 功能、输入输出或使用范围 |
|---|---|
| [开发与使用指南.md](</home/lzy/projects/vascular_printing/开发与使用指南.md>) | BraVa 原始数据、语义 ROI、BALANCED 候选、四端口血管芯及 candidate 0/15 的完整来源说明。 |
| [vessel_model/T - Brava/swc_files](</home/lzy/projects/vascular_printing/vessel_model/T - Brava/swc_files>) | BraVa 原始 SWC 和 ColorCoded 输入；当前 CFD 采用工程处理后的 ROI，并非直接求解全脑原始树。 |
| [tools/compact_bg001_rmca.py](</home/lzy/projects/vascular_printing/tools/compact_bg001_rmca.py>) | 从 refined ROI 生成 MINI/BALANCED/RICH 制造候选的代码入口。 |
| [s1-4_sacrificial_box_and_ports.py](</home/lzy/projects/vascular_printing/s1-4_sacrificial_box_and_ports.py>) | 四端口血管芯、人工接管及对齐处理；当前 CFD 几何的上游。 |
| [s1-6_abs_casting_mold.py](</home/lzy/projects/vascular_printing/s1-6_abs_casting_mold.py>) | 完整倒模盒与打印姿态生成；盒体不进入当前 CFD 流体域。 |
| [outputs/topbrain_brava_transfer/nn_production/BG001/RMCA/compact_manufacturing_roi/print_fixture_design_all_ports_aligned/core/BG001_RMCA_BALANCED_core_with_ports_all_aligned.stl](</home/lzy/projects/vascular_printing/outputs/topbrain_brava_transfer/nn_production/BG001/RMCA/compact_manufacturing_roi/print_fixture_design_all_ports_aligned/core/BG001_RMCA_BALANCED_core_with_ports_all_aligned.stl>) | 实际采用的已验收四端口芯，原始单位 mm。 |
| [outputs/topbrain_brava_transfer/nn_production/BG001/RMCA/compact_manufacturing_roi/final_abs_casting_mold/candidate_0_and_15](</home/lzy/projects/vascular_printing/outputs/topbrain_brava_transfer/nn_production/BG001/RMCA/compact_manufacturing_roi/final_abs_casting_mold/candidate_0_and_15>) | 两种最终打印姿态、模型和制造记录；用于核对姿态，不是两份独立流场。 |


### 3.2 BraVa CFD 的冻结输入和求解网格

本节根目录：[/home/lzy/projects/computation_examples/brava_flow_roi_18mlmin](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin>)

| 路径入口 | 功能、输入输出或使用范围 |
|---|---|
| [inputs/geometry.json](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/inputs/geometry.json>) | 源 STL 路径/哈希、mm→m 单位、平移原点、端口身份、中心坐标、法向、材料、流量和两姿态变换；端口命名以此为准。 |
| [inputs/source_core_mm.stl](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/inputs/source_core_mm.stl>) | 源血管芯的本算例副本，单位 mm。 |
| [inputs/labelled_core_m.vtp](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/inputs/labelled_core_m.vtp>) | 缩放到 m 并赋予边界身份的表面。 |
| [inputs/fitted_centerline.swc](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/inputs/fitted_centerline.swc>) | 实际拟合中心线，用于几何与阻力设计辅助分析。 |
| [inputs/original_radius.swc](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/inputs/original_radius.swc>) | 原始半径参考。 |
| [inputs/compensated.swc](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/inputs/compensated.swc>) | 制造补偿半径参考。 |
| [inputs/port_attachment_summary.json](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/inputs/port_attachment_summary.json>) | 人工接管与原血管连接关系。 |
| [inputs/boundary_design.json](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/inputs/boundary_design.json>) | 中心线阻力/管长近似得到的初步边界设计；不是最终出口压力的唯一依据。 |
| [mesh/generation](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/mesh/generation>) | 官方 SimVascular/TetGen 原始生成产物及生成记录。 |
| [mesh/SV_MESH/mesh-complete.mesh.vtu](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/mesh/SV_MESH/mesh-complete.mesh.vtu>) | 权威三维四面体网格：135852 节点、692313 四面体，坐标 m。 |
| [mesh/SV_MESH/mesh-complete.exterior.vtp](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/mesh/SV_MESH/mesh-complete.exterior.vtp>) | 体网格对应外表面；包含壁面和人工端口截面。 |
| [mesh/SV_MESH/mesh-surfaces](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/mesh/SV_MESH/mesh-surfaces>) | INLET.vtp、OUTLET_01.vtp、OUTLET_02.vtp、OUTLET_03.vtp、WALL.vtp；WALL 与开口截面分开。 |
| [mesh/SV_MESH/mesh_arrays.npz](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/mesh/SV_MESH/mesh_arrays.npz>) | 规范四面体连接、坐标、边界面和标签；梯度计算使用规范连接顺序。 |
| [mesh/SV_MESH/wall_distances.npz](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/mesh/SV_MESH/wall_distances.npz>) | 网格相关壁面距离辅助数据。 |
| [reports/geometry_qc.json](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/reports/geometry_qc.json>) | 源表面、端口和坐标准备检查。 |
| [reports/mesh_validity.json](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/reports/mesh_validity.json>) | 实际网格有效性检查。 |
| [reports/mesh_quality.json](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/reports/mesh_quality.json>) | 体网格质量统计。 |


`inputs/candidate_0_*` 和 `inputs/candidate_15_*` 保存各姿态的 manifest、打印变换和还原变换；全部在 [/home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/inputs](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/inputs>)。局部端口身份不能按 face_id 的数值大小重新猜测。

## 4. 流场求解代码入口与执行链

实际顺序为：准备几何 → 生成/检查网格 → 初步阻力设计 → 等流量出口校准 → 测量并统一平移出口参考压力 → 全压力出口正式 CFD → 独立后处理与冻结结果。求解程序为 svMultiPhysics，CPU/MPI 负责装配及部分预条件计算，PETSc CUDA 参与线性代数；不是所有工作都在 GPU。

| 路径入口 | 功能、输入输出或使用范围 |
|---|---|
| [scripts/prepare_geometry.py](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/scripts/prepare_geometry.py>) | 读取已验收血管芯，复制来源、识别端口、建立 SI 单位和两种打印姿态映射。 |
| [scripts/generate_mesh.py](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/scripts/generate_mesh.py>) | 调用官方 SimVascular/TetGen 生成体网格。 |
| [scripts/audit_mesh.py](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/scripts/audit_mesh.py>) | 检查连接、闭合性、正体积、质量和端口，导出实际 SV_MESH。 |
| [scripts/design_boundary.py](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/scripts/design_boundary.py>) | 用中心线阻力近似作设计辅助，输出 boundary_design.json / network_edges.csv。 |
| [scripts/prepare_case.py](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/scripts/prepare_case.py>) | 创建独立 calibration 或 pressure 算例、solver.xml、policy 和输入身份；不是求解器本体。 |
| [scripts/start_calibration_remote.py](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/scripts/start_calibration_remote.py>) | 早期校准启动封装；用于理解执行来源，已有结果不应从该入口无差别重启。 |
| [scripts/run_flow_sequence.py](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/scripts/run_flow_sequence.py>) | 校准结果后处理 → 读取设计压力 → 创建 balanced_pressure_final → 正式求解 → 独立后处理；会启动 CFD，且拒绝已有正式目录覆盖。 |
| [scripts/solve_remote.py](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/scripts/solve_remote.py>) | 实际服务器运行器：MPI 启动、残差/稳态/守恒监视、停止请求、检查点与执行记录。 |
| [scripts/postprocess_flow.py](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/scripts/postprocess_flow.py>) | 读取原生速度/压力、核对 mesh/flow、检查相邻稳态结果、恢复 WSS、积分端口流量，生成 frozen_flow 和 INDEPENDENT_CHECK。 |
| [vendor/flow_solver_support/wss.py](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/vendor/flow_solver_support/wss.py>) | WSS 数学核：p1_gradients → boundary_owners → wall_geometry → tangential_traction；nodal_average 仅用于节点显示。 |
| [vendor/flow_solver_support/wss_case.py](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/vendor/flow_solver_support/wss_case.py>) | 材料读取、坐标身份及逐算例 WSS 恢复封装 recover。 |
| [vendor/flow_solver_support](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/vendor/flow_solver_support>) | 残差解析、稳态监视、检查点、并行和数值质量检查的本算例代码副本。 |
| [vendor/vascular_validation](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/vendor/vascular_validation>) | 几何、网格诊断、结果积分、来源身份和验证的本算例代码副本。 |
| [scripts/build_gpu_sync_fix.py](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/scripts/build_gpu_sync_fix.py>) | 独立构建 GPU 解向 CPU 显式同步的修复版可执行程序；不覆盖原始 solver。 |
| [gpu_solver_fix/gpu_host_sync.patch](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/gpu_solver_fix/gpu_host_sync.patch>) | 上述 PETSc 同步修复的补丁。 |
| [gpu_solver_fix/build.json](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/gpu_solver_fix/build.json>) | 构建命令、原/新 executable 路径及源码哈希；其 BUILT_NOT_YET_VALIDATED 是构建时状态，后续验证见下项。 |
| [scripts/validate_gpu_backend.py](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/scripts/validate_gpu_backend.py>) | 比较已保存 CPU8/GPU8 的同条件短步流场及 WSS，不是新 CFD 启动器。 |
| [gpu_solver_fix/backend_equivalence.json](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/gpu_solver_fix/backend_equivalence.json>) | CPU8 原始后端与 GPU8 修复后端两步比较结果，passed=true；不能由此宣称网格无关。 |
| [scripts/check_internal_sections.py](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/scripts/check_internal_sections.py>) | 在人工直管段内部截面上积分真实速度通量。 |
| [scripts/check_internal_divergence.py](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/scripts/check_internal_divergence.py>) | 独立散度定理核对，解释内部截面流量差。 |
| [scripts/plot_internal_flux.py](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/scripts/plot_internal_flux.py>) | 把内部通量和散度检查绘成图。 |


正式计算的实际命令及 PETSc 选项以 [execution.json](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/cases/balanced_pressure_final/reports/execution.json>) 和 [PETSC_OPTIONS.txt](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/cases/balanced_pressure_final/run/PETSC_OPTIONS.txt>) 为准：8 MPI ranks、OMP 1、FGMRES、CUDA 向量/稀疏矩阵，局部 ASM/LU 在 CPU。实际 CFD 时间步为 0.01 s，与微泡名义 0.0005 s 不同。

WSS 从冻结 P1 FEM 速度梯度恢复，每个 WALL 三角形使用相邻四面体的梯度：

$$G=\nabla u,\qquad t_\mu=\mu(G+G^T)n,\qquad \tau_w=t_\mu-(t_\mu\cdot n)n,\qquad WSS_{raw}=\|\tau_w\|.$$

壁面面片原值与面积加权节点显示值分别保存；`WSS_display_Pa` 的平滑外观不表示原始面片值相同。

## 5. 正式流场配置、结果、报告与日志

正式算例根目录：[/home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/cases/balanced_pressure_final](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/cases/balanced_pressure_final>)

| 路径入口 | 功能、输入输出或使用范围 |
|---|---|
| [purpose.json](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/cases/balanced_pressure_final/purpose.json>) | 实际正式出口压力及其校准来源，明确三个出口流量是求解后测量值。 |
| [policy.json](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/cases/balanced_pressure_final/policy.json>) | 执行与验收策略；不要用旧模板默认值替代它。 |
| [input_hashes.json](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/cases/balanced_pressure_final/input_hashes.json>) | 此次计算输入哈希。 |
| [run/solver.xml](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/cases/balanced_pressure_final/run/solver.xml>) | 真正传入 svMultiPhysics 的材料、步长、边界和求解设置。 |
| [run/PETSC_OPTIONS.txt](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/cases/balanced_pressure_final/run/PETSC_OPTIONS.txt>) | 实际 PETSc 配置。 |
| [run/solver.log](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/cases/balanced_pressure_final/run/solver.log>) | 原始残差和线性/非线性求解日志。 |
| [run/petsc_profile.txt](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/cases/balanced_pressure_final/run/petsc_profile.txt>) | CPU/GPU PETSc 性能记录。 |
| [run/8-procs](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/cases/balanced_pressure_final/run/8-procs>) | 原始并行结果、checkpoint 和时间步输出；日常后处理优先使用冻结场。 |
| [reports/execution.json](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/cases/balanced_pressure_final/reports/execution.json>) | 执行命令、进程配置、退出码、残差解析、守恒、稳态、GPU 采样及输入/输出哈希。 |
| [reports/executed_runner.py](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/cases/balanced_pressure_final/reports/executed_runner.py>) | 本次实际运行器副本，追踪历史执行版本。 |
| [reports/executed_postprocess.py](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/cases/balanced_pressure_final/reports/executed_postprocess.py>) | 本次实际 WSS/流场后处理脚本副本。 |
| [frozen_flow/steady_flow.vtu](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/cases/balanced_pressure_final/frozen_flow/steady_flow.vtu>) | 正式 SI 冻结解：节点 Velocity（m/s）、Pressure（Pa），坐标 m；微泡读取的就是此场。 |
| [frozen_flow/flow_arrays_si.npz](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/cases/balanced_pressure_final/frozen_flow/flow_arrays_si.npz>) | 同一场的 NumPy 数组，包括规范 tetra 连接、边界及速度/压力。 |
| [frozen_flow/candidate_0](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/cases/balanced_pressure_final/frozen_flow/candidate_0>) | candidate 0 姿态下的场副本，仅刚体变换。 |
| [frozen_flow/candidate_15](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/cases/balanced_pressure_final/frozen_flow/candidate_15>) | candidate 15 姿态下的场副本，仅刚体变换。 |
| [postprocess/wall_wss_si.vtp](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/cases/balanced_pressure_final/postprocess/wall_wss_si.vtp>) | 原始壁面三角形 WSS_raw_Pa、切向矢量、法向、父 tetra，以及节点显示 WSS_display_Pa。 |
| [postprocess/pressure_surface_si.vtp](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/cases/balanced_pressure_final/postprocess/pressure_surface_si.vtp>) | 表面压力/速度展示输入。 |
| [postprocess/boundary_summary.csv](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/cases/balanced_pressure_final/postprocess/boundary_summary.csv>) | 每个端口实测流量、入口占比和截面平均压力。注意平均压力不等同于 prescribed traction BC 数值。 |
| [postprocess/INDEPENDENT_CHECK.json](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/cases/balanced_pressure_final/postprocess/INDEPENDENT_CHECK.json>) | 正式场独立检查、场身份、WSS 和稳态数值摘要。 |


### 5.1 当前数值身份与边界

流场 SHA256：`73228d094c253cec75dbb0030476b687171d0d0fef7816a845392dad220285dc`。

权威体网格 SHA256：`e755a6ce04eca518fc02df5cd6139c98887e94617bdf9410d72f55b99f443c94`。

材料为 ρ=1056 kg/m³、μ=0.00345312 Pa·s；刚性无滑移壁面、无重力、牛顿流体。这是当前数值模型参数，不代表实际打印灌注液已测物性。

| 边界 | 规定压力/牵引 BC（Pa） | 求解后流量（μL/min） | 相对入口流量（%） |
|---|---:|---:|---:|
| INLET | 不规定；给定流量 18000 μL/min | 18000.000000 | 100.000000 |
| OUTLET_01 | 556.970655549 | 6000.016253 | 33.333424 |
| OUTLET_02 | 0.000000000 | 5998.129842 | 33.322944 |
| OUTLET_03 | 17.137427154 | 6001.853905 | 33.343633 |

### 5.2 校准与总体检查

| 路径入口 | 功能、输入输出或使用范围 |
|---|---|
| [brava_calibration_equal_split/flow](</home/lzy/projects/computation_examples/brava_calibration_equal_split/flow>) | 已迁出的等流量校准实体目录；BraVa 的 cases/calibration_gpu8_production 是兼容链接。校准强制流量不是最终自然分流结论。 |
| [brava_calibration_equal_split/flow/postprocess/designed_outlet_pressures_Pa.json](</home/lzy/projects/computation_examples/brava_calibration_equal_split/flow/postprocess/designed_outlet_pressures_Pa.json>) | 从校准截面平均压力提取、减去共同最小值后的正式压力输入。 |
| [brava_calibration_equal_split/flow/postprocess/pressure_design_note.json](</home/lzy/projects/computation_examples/brava_calibration_equal_split/flow/postprocess/pressure_design_note.json>) | 压力设计方法及统一参考平移说明。 |


| 路径入口 | 功能、输入输出或使用范围 |
|---|---|
| [reports/FLOW_SEQUENCE_COMPLETE.json](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/reports/FLOW_SEQUENCE_COMPLETE.json>) | 完整校准→正式场顺序完成记录。 |
| [reports/LOCAL_FLOW_VERIFICATION.json](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/reports/LOCAL_FLOW_VERIFICATION.json>) | 既有传回核验记录；是记录时点的收据，不应假定其旧媒体哈希描述所有后续展示改版。 |
| [reports/internal_sections/INTERNAL_SECTION_CHECK.json](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/reports/internal_sections/INTERNAL_SECTION_CHECK.json>) | 内部截面流量检查汇总。 |
| [reports/internal_sections/DIVERGENCE_THEOREM_CHECK.json](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/reports/internal_sections/DIVERGENCE_THEOREM_CHECK.json>) | 体积分/面通量独立对应。 |
| [reports/internal_sections/section_flux.csv](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/reports/internal_sections/section_flux.csv>) | 各端口向内 2/5/10 mm 截面通量。 |
| [reports/internal_sections/divergence_theorem.csv](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/reports/internal_sections/divergence_theorem.csv>) | 内部通量偏差和散度积分对照。 |
| [reports/internal_sections/internal_flux_and_divergence.png](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/reports/internal_sections/internal_flux_and_divergence.png>) | 内部通量与散度结果图，旁有 PDF。 |
| [logs](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/logs>) | 网格、各次 CFD、顺序控制、流场绘图和回传日志。 |


正式场全局边界守恒通过，但已有内部截面偏差可达入口流量约 2.49%；仍为单档网格。使用局部近壁速度或轨迹解释时，应同时阅读第 10 节近壁审计，不能把总守恒或 CPU/GPU 一致当成局部 CFD 精度已证明。

## 6. 流场动画与后处理入口

| 路径入口 | 功能、输入输出或使用范围 |
|---|---|
| [scripts/export_print_pose_fields.py](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/scripts/export_print_pose_fields.py>) | 把同一合格流场转换到 candidate 0/15 的打印坐标姿态；速度等矢量同时旋转。 |
| [scripts/prepare_flow_visualization.py](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/scripts/prepare_flow_visualization.py>) | 从真实 P1 FEM 场生成流线与局部速度 glyph 采样；流线不是有限尺寸微泡动力学轨迹。 |
| [scripts/render_flow_parallel.py](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/scripts/render_flow_parallel.py>) | 两姿态流场动画/静态图的渲染调度器；会生成展示文件，不启动 CFD。 |
| [visualization/render_flow.py](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/visualization/render_flow.py>) | 速度流线、局部矢量、压力、显示 WSS、原始 WSS 的主体渲染。 |
| [visualization/style.py](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/visualization/style.py>) | 共享配色、文字、坐标网格、相机、旁置标签等风格实现；最新微泡坐标轴也复用其布局方法。 |
| [visualization/style_provenance.json](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/visualization/style_provenance.json>) | 来自既有 rotate_visualization 的风格来源、调整与保留项。 |
| [visualization/render_streamlines_fitted.py](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/visualization/render_streamlines_fitted.py>) | candidate 0 流线的 36 s 慢旋转、整周最大安全取景改版；publish 才更新页面引用文件。 |
| [visualization/data](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/visualization/data>) | 真实流线、局部速度采样及准备清单。 |
| [visualization/candidate_0/animations](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/visualization/candidate_0/animations>) | candidate 0 的 streamlines_rotation、velocity_vectors_rotation、pressure_rotation、wss_rotation、wss_raw_rotation 五个 MP4。 |
| [visualization/candidate_15/animations](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/visualization/candidate_15/animations>) | candidate 15 对应五个 MP4。 |
| [visualization/candidate_0/figures](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/visualization/candidate_0/figures>) | candidate 0 高清图片。 |
| [visualization/candidate_15/figures](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/visualization/candidate_15/figures>) | candidate 15 高清图片。 |
| [visualization/candidate_0/MEDIA_VALIDATION.json](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/visualization/candidate_0/MEDIA_VALIDATION.json>) | 当前各视频各自的帧数、时长、来源及媒体校验。 |
| [visualization/candidate_15/MEDIA_VALIDATION.json](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/visualization/candidate_15/MEDIA_VALIDATION.json>) | candidate 15 媒体记录。 |
| [visualization/candidate_0/streamlines_adjustment_20260929](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/visualization/candidate_0/streamlines_adjustment_20260929>) | 流线慢旋转/最大取景改版的报告、CAMERA_FIT、验证、改版前备份及渲染输出。 |
| [scripts/verify_flow_delivery.py](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/scripts/verify_flow_delivery.py>) | 既有传回流场及媒体核验入口；会写验证收据，使用时要注意媒体清单所属版本。 |


candidate 0 当前流线为 **36 s/周、864 帧、24 fps**；该姿态其他四段仍为原 18 s/432 帧。流线最终取景为保证整周完整入镜而最大化，不能继续称固定放大 30%。`wss_rotation.mp4` 为节点显示值，`wss_raw_rotation.mp4` 为逐三角形原值。

## 7. 当前微泡计算代码、合同及 GPU/CPU 分工

### 7.1 真正的代码与正式入口

本地当前实现：[/home/lzy/projects/computation_examples/_shared/particle_workflow/particle_3d/src/particle_3d/population_v1](</home/lzy/projects/computation_examples/_shared/particle_workflow/particle_3d/src/particle_3d/population_v1>)。BraVa 算例中的 [/home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/code/particle_core](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/code/particle_core>) 指向同一实现。BraVa 的 `population_v1/` 主要保存算例配置、运行数据、分析和展示，不包含第二套物理核心。

| 路径入口 | 功能、输入输出或使用范围 |
|---|---|
| [run_fixed_window_population.py](</home/lzy/projects/computation_examples/_shared/particle_workflow/particle_3d/src/particle_3d/population_v1/run_fixed_window_population.py>) | 当前唯一正式 CLI：main 读取固定合同、校验流场/源身份、选择 CPU/CUDA、拒绝覆盖已有 run，并记录实际参数与代码身份。模块名 particle_3d.population_v1.run_fixed_window_population。 |
| [run.py](</home/lzy/projects/computation_examples/_shared/particle_workflow/particle_3d/src/particle_3d/population_v1/run.py>) | 兼容入口，仅转发到固定时间 CLI；旧 --minimum-exits、--post-infusion、--fixed-infusion 不再接受。 |
| [fixed_window.py](</home/lzy/projects/computation_examples/_shared/particle_workflow/particle_3d/src/particle_3d/population_v1/fixed_window.py>) | run_fixed_window 负责固定源时间窗、出生切步、精确终点；censor_remaining 记录观察截止末态；FixedWindowOutput 分开保存源提议与准入账本。 |
| [environment.py](</home/lzy/projects/computation_examples/_shared/particle_workflow/particle_3d/src/particle_3d/population_v1/environment.py>) | load_environment：验证 flow/mesh/五个边界哈希，建立 P1 场、有限三角形 WALL、出口分类器、真实通量入口及原 SonoVue 分布。 |
| [source.py](</home/lzy/projects/computation_examples/_shared/particle_workflow/particle_3d/src/particle_3d/population_v1/source.py>) | 事件寻址 Poisson 源、粒径/入口位置、一次准入检查；拒绝后不重试、不按出口挑选微泡。 |
| [engine.py](</home/lzy/projects/computation_examples/_shared/particle_workflow/particle_3d/src/particle_3d/population_v1/engine.py>) | CPU 共享物理时钟参考引擎；所有活跃泡同一接受时间区间推进。 |
| [gpu_engine.py](</home/lzy/projects/computation_examples/_shared/particle_workflow/particle_3d/src/particle_3d/population_v1/gpu_engine.py>) | CUDA 群体引擎、事件步长、有限尺寸壁面/球对检查及困难小子集 CPU 回退。 |
| [gpu_field.py](</home/lzy/projects/computation_examples/_shared/particle_workflow/particle_3d/src/particle_3d/population_v1/gpu_field.py>) | 规范 P1 四面体 GPU 定位、邻接 walking、BVH 恢复及速度/压力/梯度插值。 |
| [gpu_wall.py](</home/lzy/projects/computation_examples/_shared/particle_workflow/particle_3d/src/particle_3d/population_v1/gpu_wall.py>) | 有限三角形壁面最近点、间隙、法向和连续路径证明。 |
| [gpu_neighbors.py](</home/lzy/projects/computation_examples/_shared/particle_workflow/particle_3d/src/particle_3d/population_v1/gpu_neighbors.py>) | GPU 空间散列、邻居候选、精确球间隙及 realization 隔离。 |
| [gpu_resistance.py](</home/lzy/projects/computation_examples/_shared/particle_workflow/particle_3d/src/particle_3d/population_v1/gpu_resistance.py>) | FP64 批量近壁阻力、梯度 forcing 和接触约束解。 |
| [gpu_state.py](</home/lzy/projects/computation_examples/_shared/particle_workflow/particle_3d/src/particle_3d/population_v1/gpu_state.py>) | GPU SoA 状态、稳定 ID、出生/移除及输出数据转移。 |
| [gpu_safety.py](</home/lzy/projects/computation_examples/_shared/particle_workflow/particle_3d/src/particle_3d/population_v1/gpu_safety.py>) | GPU 时间进展/终态候选证据，不随意把低速轨迹当失败。 |
| [gpu_batch.py](</home/lzy/projects/computation_examples/_shared/particle_workflow/particle_3d/src/particle_3d/population_v1/gpu_batch.py>) | 多个独立 realization 批处理；来源与时钟独立，禁止跨 realization 相互碰撞。 |
| [gpu_extension.py](</home/lzy/projects/computation_examples/_shared/particle_workflow/particle_3d/src/particle_3d/population_v1/gpu_extension.py>) | 本项目 CUDA 扩展的加载/编译，float64 核。 |
| [cuda](</home/lzy/projects/computation_examples/_shared/particle_workflow/particle_3d/src/particle_3d/population_v1/cuda>) | 实际 CUDA/C++ 查询及计算核源码。 |
| [bvh.py](</home/lzy/projects/computation_examples/_shared/particle_workflow/particle_3d/src/particle_3d/population_v1/bvh.py>) | CPU 预处理静态精确几何 BVH。 |
| [components.py](</home/lzy/projects/computation_examples/_shared/particle_workflow/particle_3d/src/particle_3d/population_v1/components.py>) | 接触连通组件阻力装配与耦合约束求解；不创建全群体 dense 6N×6N 系统。 |
| [pairs.py](</home/lzy/projects/computation_examples/_shared/particle_workflow/particle_3d/src/particle_3d/population_v1/pairs.py>) | 精确球几何、近场作用与连续扫掠候选。 |
| [reference.py](</home/lzy/projects/computation_examples/_shared/particle_workflow/particle_3d/src/particle_3d/population_v1/reference.py>) | 调用原始单泡物理的 CPU 权威参考。 |
| [safety.py](</home/lzy/projects/computation_examples/_shared/particle_workflow/particle_3d/src/particle_3d/population_v1/safety.py>) | 继承的时间进展/支持静止证据规则。 |
| [state.py](</home/lzy/projects/computation_examples/_shared/particle_workflow/particle_3d/src/particle_3d/population_v1/state.py>) | CPU 活跃状态、身份和终态移除。 |
| [output.py](</home/lzy/projects/computation_examples/_shared/particle_workflow/particle_3d/src/particle_3d/population_v1/output.py>) | 接受时间端点状态分块输出，源账本、终态和交互记录。 |
| [provenance.py](</home/lzy/projects/computation_examples/_shared/particle_workflow/particle_3d/src/particle_3d/population_v1/provenance.py>) | 输入、源码、合同、种子和输出的可复查身份。 |


### 7.2 当前合同、配置、启动和测试

| 路径入口 | 功能、输入输出或使用范围 |
|---|---|
| [fixed_window_protocol/FIXED_WINDOW_POPULATION_EXPERIMENT_V1.json](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/fixed_window_protocol/FIXED_WINDOW_POPULATION_EXPERIMENT_V1.json>) | 当前生产固定合同：C=4×10⁹ m⁻³、Q=3×10⁻⁷ m³/s、名义 dt=0.0005 s、T_infusion=1.240815660 s、停灌后 10 s、T_end=11.240815660 s。 |
| [fixed_window_protocol/CODE_VERSION.json](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/fixed_window_protocol/CODE_VERSION.json>) | 固定协议开发版本记录；该记录不替代具体 run 的逐源码 SHA。 |
| [data/environment.local.json](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/data/environment.local.json>) | WSL 流场/网格/边界及粒径数据路径；仍含旧绝对路径，但由迁移兼容链接定位新实体文件，身份哈希未变。 |
| [data/environment.json](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/data/environment.json>) | 服务器环境路径与同一套输入哈希。 |
| [contracts/MB_MB_INTERACTION_V1.json](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/contracts/MB_MB_INTERACTION_V1.json>) | 微泡间 V1 近场/接触模型契约。 |
| [contracts/POPULATION_RUN_REQUEST.json](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/contracts/POPULATION_RUN_REQUEST.json>) | 旧覆盖停止批次的请求记录；不能替代当前固定时间合同。 |
| [scripts/launch_production.py](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/scripts/launch_production.py>) | 服务器新批次启动封装；当前使用固定合同、CUDA 后端、新目录，记录 PID 和 CPU/GPU 遥测；执行会开始轨迹积分。 |
| [scripts/analyze_population.py](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/scripts/analyze_population.py>) | 读取已有 run，独立检查出口/壁面与账本，导出逐条轨迹和物理时间展示数据。 |
| [scripts/collect_results.sh](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/scripts/collect_results.sh>) | 服务器结果回传到本地；是文件传输入口。 |
| [scripts/final_qc.py](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/scripts/final_qc.py>) | 半径、身份、来源哈希、出口端点与编码视频终检。 |
| [scripts/inventory_delivery.py](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/scripts/inventory_delivery.py>) | 生成代码/结果/审阅附件清单。 |
| [scripts/build_report.py](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/scripts/build_report.py>) | 根据指定已有 run 和独立检查生成中文报告；会写派生报告。 |


测试入口：[/home/lzy/projects/computation_examples/_shared/particle_workflow/particle_3d/tests/population_v1](</home/lzy/projects/computation_examples/_shared/particle_workflow/particle_3d/tests/population_v1>)；当前协议重点为 [test_fixed_window.py](</home/lzy/projects/computation_examples/_shared/particle_workflow/particle_3d/tests/population_v1/test_fixed_window.py>)。真实场/FEM/壁面与 CPU/GPU 比较脚本在 [/home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/scripts](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/scripts>)：`validate_gpu_field.py`、`validate_gpu_wall.py`、`validate_real_reference.py`、`validate_real_gpu_dynamics.py`、`validate_real_source.py`、`validate_contact_fallback.py`。其中有些会执行小规模粒子验证，不应当成只读文件查看器。

### 7.3 原单泡物理依赖和粒径数据

以下仍被当前 population 模块复用，不是仅凭旧阶段名称可以删除的文件。根目录：[/home/lzy/projects/computation_examples/_shared/particle_workflow/particle_3d/src/particle_3d](</home/lzy/projects/computation_examples/_shared/particle_workflow/particle_3d/src/particle_3d>)

| 路径入口 | 功能、输入输出或使用范围 |
|---|---|
| [field.py](</home/lzy/projects/computation_examples/_shared/particle_workflow/particle_3d/src/particle_3d/field.py>) | 冻结 FEM 场、规范网格对应、P1 插值和逐四面体速度梯度。 |
| [wall_geometry.py](</home/lzy/projects/computation_examples/_shared/particle_workflow/particle_3d/src/particle_3d/wall_geometry.py>) | 真实有限三角形壁面几何。 |
| [validation_boundary.py](</home/lzy/projects/computation_examples/_shared/particle_workflow/particle_3d/src/particle_3d/validation_boundary.py>) | 正式出口身份与穿越分类。 |
| [inlet_flux.py](</home/lzy/projects/computation_examples/_shared/particle_workflow/particle_3d/src/particle_3d/inlet_flux.py>) | 真实入口正通量积分/采样。 |
| [injection_method_c.py](</home/lzy/projects/computation_examples/_shared/particle_workflow/particle_3d/src/particle_3d/injection_method_c.py>) | TruncatedSonoVue 等原条件粒径与入口方法。 |
| [injection_admission.py](</home/lzy/projects/computation_examples/_shared/particle_workflow/particle_3d/src/particle_3d/injection_admission.py>) | 有限尺寸微泡在 WALL/流体域附近的准入判定。 |
| [open_boundary_rim.py](</home/lzy/projects/computation_examples/_shared/particle_workflow/particle_3d/src/particle_3d/open_boundary_rim.py>) | 开放端口 rim 拓扑及壁面/出口边缘处理。 |
| [particle9a_motion.py](</home/lzy/projects/computation_examples/_shared/particle_workflow/particle_3d/src/particle_3d/particle9a_motion.py>) | 原单泡近壁动力学权威实现。 |
| [planar_wall_hydrodynamics.py](</home/lzy/projects/computation_examples/_shared/particle_workflow/particle_3d/src/particle_3d/planar_wall_hydrodynamics.py>) | 平面近壁平动/转动水动力修正与作用范围。 |
| [hydrodynamic_resistance.py](</home/lzy/projects/computation_examples/_shared/particle_workflow/particle_3d/src/particle_3d/hydrodynamic_resistance.py>) | 阻力模型基础。 |
| [resistance_assembly.py](</home/lzy/projects/computation_examples/_shared/particle_workflow/particle_3d/src/particle_3d/resistance_assembly.py>) | 阻力装配。 |
| [resistance_solver.py](</home/lzy/projects/computation_examples/_shared/particle_workflow/particle_3d/src/particle_3d/resistance_solver.py>) | 阻力/约束系统求解基础。 |
| [kinematic_contact.py](</home/lzy/projects/computation_examples/_shared/particle_workflow/particle_3d/src/particle_3d/kinematic_contact.py>) | 不可穿透接触运动学约束。 |


粒径数据入口：[/home/lzy/projects/computation_examples/_shared/sonovue_size_distribution_v0](</home/lzy/projects/computation_examples/_shared/sonovue_size_distribution_v0>)。生产采用原 SonoVue **D≤4 μm 条件分布**，不是裁掉大泡后篡改轨迹。当前 BraVa population 模型没有加入 RBC、声学振荡、泡壳破裂、聚并或双向 CFD 反馈；不能将其他项目的 RBC 参考文件当作此批结果。

GPU 实际负责 FEM/壁面查询、邻居搜索、批量阻力/接触及状态推进；CPU 负责源事件、准入、严格出口判定、困难球对/有限边缘小组件回退和输出。旧正式 run 内部计数约 99.83% 接受动力学更新由 GPU 求解；这是该次运行记录，不是对所有工况的性能保证。动画采用 NVIDIA EGL 渲染、CPU libx264 编码。

## 8. 已有正式微泡数据：应该读取哪些文件

已算群体的唯一正式 run：[/home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/runs/realization_0000_contactfix](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/runs/realization_0000_contactfix>)。生产身份对应历史 commit `98d705f13babdfa0cba6c12331623292a588cb77`；当前固定协议开发记录为 `9836c264863eca7f34806c22e020b4c174ef3a7c`。两者不是同一次计算。

| 路径入口 | 功能、输入输出或使用范围 |
|---|---|
| [SUMMARY.json](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/runs/realization_0000_contactfix/SUMMARY.json>) | 实际源事件/接受/出流/失败/观察截止计数、停灌和结束时刻、GPU/CPU 更新计数、完整性标志。 |
| [IDENTITY.json](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/runs/realization_0000_contactfix/IDENTITY.json>) | 真实生产参数、种子、源/流场/逐源码哈希、设备和运行身份。 |
| [POPULATION_CONTRACT.json](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/runs/realization_0000_contactfix/POPULATION_CONTRACT.json>) | 该历史 run 真正使用的合同；与新固定窗合同分开。 |
| [MANIFEST.json](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/runs/realization_0000_contactfix/MANIFEST.json>) | 本次 run 输出文件清单/身份。 |
| [PROGRESS.json](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/runs/realization_0000_contactfix/PROGRESS.json>) | 最后保存的进度；不是自动继续任务的授权。 |
| [states](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/runs/realization_0000_contactfix/states>) | chunk_*.npz 原始接受子步群体状态；包括 stable ID、出生时间、x/v/ω/q、半径、tetra、wall gap 等。state_time_s 与 state_offsets 表示同一物理时刻的群体行段。 |
| [source.jsonl.gz](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/runs/realization_0000_contactfix/source.jsonl.gz>) | 真实源事件及一次准入接受/拒绝账本。 |
| [terminals.jsonl.gz](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/runs/realization_0000_contactfix/terminals.jsonl.gz>) | 正式出流、数值失败、观察截止等终态和真实时刻。 |
| [interaction_audit.jsonl.gz](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/runs/realization_0000_contactfix/interaction_audit.jsonl.gz>) | 球对 IDs、近场/接触、handoff、拒绝试步及 CPU 回退证据。 |


后处理结果根目录：[/home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/analysis/realization_0000_contactfix](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/analysis/realization_0000_contactfix>)

| 路径入口 | 功能、输入输出或使用范围 |
|---|---|
| [tracks](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/analysis/realization_0000_contactfix/tracks>) | 逐条轨迹 mb_*.npz，SI 单位；适合查看单泡三维位置、速度、壁面间隙和单元身份。 |
| [trajectory_catalog.csv](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/analysis/realization_0000_contactfix/trajectory_catalog.csv>) | 逐微泡编号、粒径、出生、终态、出口等目录。 |
| [population_timeline.csv](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/analysis/realization_0000_contactfix/population_timeline.csv>) | 源、活跃群体和各终态的物理时间演化。 |
| [visualization_data.npz](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/analysis/realization_0000_contactfix/visualization_data.npz>) | 用于原总览/动画的派生采样，不能替代所有接受子步的原始数据。 |
| [pair_interactions.csv](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/analysis/realization_0000_contactfix/pair_interactions.csv>) | 球对近场/接触摘要。 |
| [independent_outlet_checks.csv](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/analysis/realization_0000_contactfix/independent_outlet_checks.csv>) | 真实出口穿越独立核查。 |
| [independent_wall_checks.csv](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/analysis/realization_0000_contactfix/independent_wall_checks.csv>) | 轨迹与真实壁面独立核查。 |
| [INDEPENDENT_AUDIT.json](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/analysis/realization_0000_contactfix/INDEPENDENT_AUDIT.json>) | 账本、轨迹、几何检查汇总。 |
| [failure_diagnostics.json](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/analysis/realization_0000_contactfix/failure_diagnostics.json>) | 两条数值失败的定位及原因证据。 |
| [contact_fix_comparison.json](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/analysis/realization_0000_contactfix/contact_fix_comparison.json>) | 历史接触修复前后对照记录。 |
| [resource_statistics.json](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/analysis/realization_0000_contactfix/resource_statistics.json>) | 该次执行的资源遥测统计。 |
| [performance_summary.csv](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/analysis/realization_0000_contactfix/performance_summary.csv>) | 吞吐及阶段性能摘要。 |


轨迹使用 SI：位置 m、时间 s、速度 m/s；CSV 中显式命名 `diameter_um` 等列采用对应标注单位。显示姿态可转换到 mm。应按出生≤t<终态显示微泡，不能把不同出生时刻按 age=0 对齐冒充持续灌注，也不能把终态后的点继续停留在画面中。

## 9. 微泡动画、静态图与筛选来源

### 9.1 当前推荐展示

[四视角总预览](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/candidate_0/views_0to1p2s_large_axes/OPEN_RESULTS.html>) · [四视角末帧对照](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/candidate_0/views_0to1p2s_large_axes/CONTACT_SHEET.png>) · [字号/角度/时间窗与复现说明](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/candidate_0/views_0to1p2s_large_axes/README_ZH.md>)

| 相对基准方位角 / 仰角 | MP4 绝对路径入口 |
|---|---|
| 240° / 25.0° | [/home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/candidate_0/views_0to1p2s_large_axes/az240_el25p0/animations/microbubble_az240_el25p0_0to1p2s.mp4](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/candidate_0/views_0to1p2s_large_axes/az240_el25p0/animations/microbubble_az240_el25p0_0to1p2s.mp4>) |
| 330° / 25.0° | [/home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/candidate_0/views_0to1p2s_large_axes/az330_el25p0/animations/microbubble_az330_el25p0_0to1p2s.mp4](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/candidate_0/views_0to1p2s_large_axes/az330_el25p0/animations/microbubble_az330_el25p0_0to1p2s.mp4>) |
| 60° / 25.0° | [/home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/candidate_0/views_0to1p2s_large_axes/az060_el25p0/animations/microbubble_az060_el25p0_0to1p2s.mp4](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/candidate_0/views_0to1p2s_large_axes/az060_el25p0/animations/microbubble_az060_el25p0_0to1p2s.mp4>) |
| 150° / 25.0° | [/home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/candidate_0/views_0to1p2s_large_axes/az150_el25p0/animations/microbubble_az150_el25p0_0to1p2s.mp4](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/candidate_0/views_0to1p2s_large_axes/az150_el25p0/animations/microbubble_az150_el25p0_0to1p2s.mp4>) |

每段 1920×1080、24 fps、262 帧、约 10.916667 s，显示物理时间 0–1.2 s。黑底、蓝泡、英文文字、流场风格坐标网格、无外框，主图保持已校验取景；右上角显示角度。坐标轴标题 40 px、刻度 30 px。角度为相对基准方位角，不应当成世界坐标绝对方位角。

该展示继续隐藏 2 条失败轨迹和 58 条加速窗口经过高散度单元的数值风险轨迹，并逐帧隐藏当前时刻前后各 0.05 s 内三维平均速度低于 2 mm/s 的微泡。58 条是风险关联，不表示逐条物理错误已确认；隐藏仅影响显示，不修改原数据及总账。

| 路径入口 | 功能、输入输出或使用范围 |
|---|---|
| [data/QC.json](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/candidate_0/views_0to1p2s_large_axes/data/QC.json>) | 四段视频解码、标签布局、取景及与前版逐帧数据一致性检查。 |
| [data/views.json](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/candidate_0/views_0to1p2s_large_axes/data/views.json>) | 视角、输出文件、帧数/时长列表。 |
| [data/view_summary.csv](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/candidate_0/views_0to1p2s_large_axes/data/view_summary.csv>) | 各视角展示摘要。 |
| [az240_el25p0/data/render_manifest.json](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/candidate_0/views_0to1p2s_large_axes/az240_el25p0/data/render_manifest.json>) | 该视角完整渲染参数、数据身份与筛选语义；其他角度保存相同结构。 |
| [az240_el25p0/data/animation_samples.npz](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/candidate_0/views_0to1p2s_large_axes/az240_el25p0/data/animation_samples.npz>) | 真实物理时刻、位置、平均速率和 displayed 掩码；可核查为何某泡某帧未显示。 |
| [az240_el25p0/data/frame_summary.csv](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/candidate_0/views_0to1p2s_large_axes/az240_el25p0/data/frame_summary.csv>) | 每帧 alive/prior_hidden/slow_hidden/displayed 数量。 |
| [az240_el25p0/data/coordinate_annotations.json](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/candidate_0/views_0to1p2s_large_axes/az240_el25p0/data/coordinate_annotations.json>) | 坐标名称、数字和实际屏幕边界框。 |


### 9.2 当前与既有绘图代码

| 路径入口 | 功能、输入输出或使用范围 |
|---|---|
| [scripts/render_population_four_views.py](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/scripts/render_population_four_views.py>) | 当前四角度批量渲染入口；调用 render_view05_flow_axes。默认小字号，复现最新版本须带 40/30/38 参数。 |
| [scripts/render_view05_flow_axes.py](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/scripts/render_view05_flow_axes.py>) | 读取既有逐条轨迹、逐帧低速筛选、固定视角、流场式坐标轴和角度文字；不积分新轨迹。 |
| [scripts/check_population_four_views.py](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/scripts/check_population_four_views.py>) | 独立核查四角度视频、数值采样、角度、边界、文字；会输出 QC/抽帧等检查文件。 |
| [scripts/render_static_t6s.py](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/scripts/render_static_t6s.py>) | 物理 t=6 s 的多角度静态图，保留白框阶段的显示版本；不是当前无框动画。 |
| [scripts/render_population.py](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/scripts/render_population.py>) | 两姿态原群体动画及早期 candidate 0 显示调整入口；重现其默认输出不会自动得到最新四角度版本。 |
| [scripts/style_reference.py](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/scripts/style_reference.py>) | 保留原透明血管、相机和标签布局参考。 |
| [scripts/spatial_axes_reference.py](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/scripts/spatial_axes_reference.py>) | 此前按人脑 ROI 右侧子图提取的空间坐标风格；当前 flow-axes 版本使用第 6 节共享流场布局。 |
| [scripts/plot_summary.py](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/scripts/plot_summary.py>) | 群体账本、轨迹和性能静态统计图。 |
| [candidate_0/view05_flow_axes](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/candidate_0/view05_flow_axes>) | 旧 0–6 s、view05 无外框流场坐标轴动画与记录。 |
| [candidate_0/views_0to1p2s](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/candidate_0/views_0to1p2s>) | 四视角 0–1.2 s 的小字号前版。 |
| [candidate_0/static_t6s_white_frame](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/candidate_0/static_t6s_white_frame>) | t=6 s 白色外框多角度静态图及低速筛选记录。 |
| [candidate_0/animations](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/candidate_0/animations>) | 早期 candidate 0 群体动画及展示改版；最新四视角 MP4 不在此目录。 |
| [candidate_15/animations/microbubble_physical_time.mp4](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/candidate_15/animations/microbubble_physical_time.mp4>) | candidate 15 的既有物理时间动画；未应用后续仅 candidate 0 的全部展示调整。 |
| [figures](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/figures>) | 群体数量、性能等汇总图。 |


最新四视角复现参数为 `--output-dir <独立新目录> --axis-title-size 40 --axis-tick-size 30 --axis-tick-offset 38`；`render_population_four_views.py` 固定使用四角度和 1.2 s 终点。渲染运行环境与服务器路径见第 12 节。本次仅记录入口，没有执行这些命令。

## 10. CFD 近壁审计、加速诊断与验证报告

### 10.1 整体说明及架构

| 路径入口 | 功能、输入输出或使用范围 |
|---|---|
| [POPULATION_VALIDATION_REPORT_ZH.md](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/POPULATION_VALIDATION_REPORT_ZH.md>) | 旧正式群体模型、算法、CPU/GPU 验证、生产结果和限制的主报告；其中旧展示描述/旧停灌协议不替代当前显示或新固定协议。 |
| [reports/ARCHITECTURE_ZH.md](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/reports/ARCHITECTURE_ZH.md>) | 原物理模块复用关系、新架构、实际调用链和代码审查。 |
| [reports/PHASE1_SOURCE_REVIEW.json](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/reports/PHASE1_SOURCE_REVIEW.json>) | 逐模块来源、符号及审查记录。 |
| [reports/DEVELOPMENT_FAILURES.md](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/reports/DEVELOPMENT_FAILURES.md>) | 开发失败及修复经过，保留未通过证据。 |
| [validation/server](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/validation/server>) | 实际服务器 CPU/GPU、入口源、接触、壁面、batch 和性能验证 JSON。 |
| [validation/server/final_delivery_qc.json](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/validation/server/final_delivery_qc.json>) | 原正式批次的数据一致性和媒体交付记录；不把两条失败改成通过，且媒体部分属于当时版本。 |
| [FILE_INVENTORY.csv](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/FILE_INVENTORY.csv>) | 该阶段交付文件清单。 |
| [DELIVERY_MANIFEST.json](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/DELIVERY_MANIFEST.json>) | 该阶段交付哈希/清单，不代替之后显示改版各自的 manifest。 |


### 10.2 当前固定时间协议（没有正式大批量轨迹）

| 路径入口 | 功能、输入输出或使用范围 |
|---|---|
| [FIXED_WINDOW_PROTOCOL_REPORT_ZH.md](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/fixed_window_protocol/FIXED_WINDOW_PROTOCOL_REPORT_ZH.md>) | 固定源窗/终点、右删失、旧 CLI 禁用、CPU/GPU 和源统计验证，及实际可用新生产入口。 |
| [OPEN_RESULTS.html](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/fixed_window_protocol/OPEN_RESULTS.html>) | 固定时间线和源计数图展示。 |
| [validation/server](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/fixed_window_protocol/validation/server>) | 固定窗测试 A–H 等服务器记录；部分为合成调度/少量粒子，不是新 BraVa 生产批次。 |
| [validation/real_source_only](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/fixed_window_protocol/validation/real_source_only>) | 真实 BraVa 来源重放，1565 源提议；没有运行准入和轨迹积分。 |
| [validation/empty_real_environment_fixedwindow_current](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/fixed_window_protocol/validation/empty_real_environment_fixedwindow_current>) | 真实 BraVa 环境 C=0 的 CLI 检查，精确到达固定终点；不是有泡正式结果。 |
| [validation/history_immutability.json](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/fixed_window_protocol/validation/history_immutability.json>) | 历史群体数据未被新协议改写的哈希证据。 |
| [migration](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/fixed_window_protocol/migration>) | 旧覆盖停灌入口的只读 .py.txt 存档和保护记录；不是当前生产入口。 |
| [logs](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/fixed_window_protocol/logs>) | 固定协议开发、测试、源重放、空源 CLI 和报告日志。 |


### 10.3 近壁 CFD 审计

审计根目录：[/home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/cfd_wall_audit](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/cfd_wall_audit>)

| 路径入口 | 功能、输入输出或使用范围 |
|---|---|
| [WALL_CFD_AUDIT_REPORT_ZH.md](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/cfd_wall_audit/WALL_CFD_AUDIT_REPORT_ZH.md>) | 近壁速度、完整梯度、散度、法向分量、WSS、网格结构与轨迹访问对应的中文审计。 |
| [OPEN_RESULTS.html](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/cfd_wall_audit/OPEN_RESULTS.html>) | 审计图和附件入口。 |
| [AUDIT_CONFIG.json](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/cfd_wall_audit/AUDIT_CONFIG.json>) | 各阈值、区域/层定义、统计单位和分母。 |
| [AUDIT_IDENTITY.json](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/cfd_wall_audit/AUDIT_IDENTITY.json>) | 原流场、网格、边界和审计输入身份。 |
| [SUMMARY.json](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/cfd_wall_audit/SUMMARY.json>) | 主数值摘要：壁面节点速度、近壁层、WSS 重算、体积和散度。 |
| [SUPPLEMENT_SUMMARY.json](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/cfd_wall_audit/SUPPLEMENT_SUMMARY.json>) | 补充结构/分区/邻域检查摘要。 |
| [all_cells.csv](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/cfd_wall_audit/all_cells.csv>) | 全部 692313 个四面体的诊断指标。 |
| [wall_adjacent_cells.csv](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/cfd_wall_audit/wall_adjacent_cells.csv>) | W0/W1/W2 邻壁层相关单元指标及映射。 |
| [wall_nodes.csv](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/cfd_wall_audit/wall_nodes.csv>) | 壁面节点与速度检查。 |
| [wall_face_structure.csv](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/cfd_wall_audit/wall_face_structure.csv>) | 壁面面片及其父单元结构。 |
| [wall_sample_metrics.csv](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/cfd_wall_audit/wall_sample_metrics.csv>) | 近壁取样、法向/切向速度等指标。 |
| [wall_audit.vtp](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/cfd_wall_audit/wall_audit.vtp>) | 可在 ParaView/PyVista 打开的壁面诊断数据。 |
| [wall_layers_audit.vtu](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/cfd_wall_audit/wall_layers_audit.vtu>) | 带层标签和诊断数组的体数据。 |
| [recomputed_arrays.npz](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/cfd_wall_audit/recomputed_arrays.npz>) | 从原坐标/节点速度独立重算的诊断数组。 |
| [particle_exposure.csv](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/cfd_wall_audit/particle_exposure.csv>) | 1565 条微泡与风险流场区域的访问/暴露汇总。 |
| [acceleration_overlap.csv](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/cfd_wall_audit/acceleration_overlap.csv>) | 加速候选与高散度区域的对应；当前展示隐藏 58 条风险轨迹的来源。 |
| [anomaly_clusters.csv](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/cfd_wall_audit/anomaly_clusters.csv>) | 高风险单元空间簇摘要；相应单元/壁面映射在同目录 cluster_*.csv。 |
| [figures](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/cfd_wall_audit/figures>) | 12 张关键图及 PDF、FIGURE_INDEX.json。 |
| [scripts](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/cfd_wall_audit/scripts>) | audit_field、validate_recomputation、supplement、particle_overlap、plot_audit、build_report、verify_delivery 等检查实现；只后处理原有场/轨迹。 |
| [validation/FINAL_DELIVERY_QC.json](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/cfd_wall_audit/validation/FINAL_DELIVERY_QC.json>) | 文件、表格、图和身份核查。 |


该审计确认壁面节点速度为零、所有体单元有正体积，独立 WSS 重算与原始值接近机器精度；同时保留局部高散度、法向分量及近壁精度风险。前几项一致不能抵消后者，也不能证明微泡加速是精确生理行为。

### 10.4 “先慢后快”专门诊断

| 路径入口 | 功能、输入输出或使用范围 |
|---|---|
| [ACCELERATION_EXPLANATION_ZH.md](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/analysis/acceleration_diagnosis/ACCELERATION_EXPLANATION_ZH.md>) | 两个旧视频播放 3 s 以后明显加速的解释、证据和限制。 |
| [OPEN_RESULTS.html](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/analysis/acceleration_diagnosis/OPEN_RESULTS.html>) | 加速案例图与说明预览。 |
| [speedup_diagnosis_all_candidates.csv](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/analysis/acceleration_diagnosis/speedup_diagnosis_all_candidates.csv>) | 201 个按明确速度阈值筛出的加速事件及 FEM/间隙/交互信息。 |
| [visible_speedup_candidates.csv](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/analysis/acceleration_diagnosis/visible_speedup_candidates.csv>) | 可见加速候选清单。 |
| [CPU_STATE_REPLAY.json](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/analysis/acceleration_diagnosis/CPU_STATE_REPLAY.json>) | 在已保存状态上用原 CPU 方程复核，不是重新积分完整轨迹。 |
| [WALL_RELEASE_AND_FEM.json](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/analysis/acceleration_diagnosis/WALL_RELEASE_AND_FEM.json>) | 有限三角形 FACE→EDGE、离壁过程与同位置 FEM 诊断。 |
| [speedup_examples.png](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/analysis/acceleration_diagnosis/speedup_examples.png>) | 四个具体 ID 的速度/间隙时间序列图，旁有 PDF。 |
| [scan_acceleration.py](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/analysis/acceleration_diagnosis/scan_acceleration.py>) | 原轨迹速度变化扫描。 |
| [diagnose_acceleration.py](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/analysis/acceleration_diagnosis/diagnose_acceleration.py>) | FEM、壁面和 CPU 状态复核。 |
| [check_release_geometry.py](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/analysis/acceleration_diagnosis/check_release_geometry.py>) | 几何边缘/离壁路径检查。 |


## 11. 日志、性能、历史批次与版本对应

| 路径入口 | 功能、输入输出或使用范围 |
|---|---|
| [logs/server/realization_0000_contactfix.log](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/logs/server/realization_0000_contactfix.log>) | 正式群体计算主日志。 |
| [logs/server/realization_0000_contactfix_PROCESS.json](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/logs/server/realization_0000_contactfix_PROCESS.json>) | 当时启动命令、PID 与运行身份。历史 PID 不用于当前停止任务。 |
| [logs/server/realization_0000_contactfix_telemetry.csv](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/logs/server/realization_0000_contactfix_telemetry.csv>) | 正式运行 CPU/GPU 时间序列；GPU 是设备级采样，可能含同期进程。 |
| [logs/server/tests_contact_fallback.log](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/logs/server/tests_contact_fallback.log>) | 原正式版本接触回退测试记录。 |
| [logs/server/four_views_large_axes_render.log](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/logs/server/four_views_large_axes_render.log>) | 最新四视角大字号渲染记录。 |
| [logs/server/four_views_large_axes_qc.log](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/logs/server/four_views_large_axes_qc.log>) | 相应视频及显示核查。 |
| [scripts/benchmark.py](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/scripts/benchmark.py>) | 真实 BraVa 上已认证单步、CPU/GPU 分阶段基准；执行会进行数值更新，不是只读报告查看。 |
| [scripts/benchmark_pair_components.py](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/scripts/benchmark_pair_components.py>) | 球对组件 CPU 回退阶段性能。 |
| [scripts/benchmark_sparse_batch.py](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/scripts/benchmark_sparse_batch.py>) | 多独立 realization 稀疏 batch 性能验证。 |
| [runs/realization_0000](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/runs/realization_0000>) | 首次中止开发批次，保留错误和部分状态。 |
| [runs/realization_0000_guardfix](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/runs/realization_0000_guardfix>) | 中止的时间进展保护修复批次，不能作为正式结果。 |
| [runs/realization_0000_final](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/population_v1/runs/realization_0000_final>) | 中止的近共线接触容差诊断批次；名字带 final 仍不是有效生产批次。 |


另外，[/home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/microbubble/data/USER_STOP.json](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/microbubble/data/USER_STOP.json>) 明确标记旧独立微泡批次被用户停止。[/home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/scripts/run_microbubble_sequence.py](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/scripts/run_microbubble_sequence.py>) 和 `microbubble/scripts/` 属于该旧方案；当前持续灌注任务从第 7 节入口开始，不能恢复旧自动衔接链。[/home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/cases](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/cases>) 中 calibration_fgmres、calibration_gpu8_lu、diagnostic_* 等是校准/后端尝试或诊断目录，不能仅凭存在 solver.log 当作正式流场。

Git 开发工作树保留于 [/home/lzy/projects/computation_examples/_archives/github_sync_20260928_latest](</home/lzy/projects/computation_examples/_archives/github_sync_20260928_latest>)；目录名中的日期并不等于每个文件或 run 的最新科学版本。当前计算代码从本指南第 7 节共享源码入口读取；追溯具体已有 run 以其 IDENTITY.json 和逐源码哈希为准。迁移记录见 [/home/lzy/projects/computation_examples/MIGRATION_REPORT.md](</home/lzy/projects/computation_examples/MIGRATION_REPORT.md>) 和 [/home/lzy/projects/computation_examples/PATH_MIGRATION.csv](</home/lzy/projects/computation_examples/PATH_MIGRATION.csv>)。

## 12. 服务器位置、环境与本地对应

SSH 别名：`vast4090`。服务器算例根：`/workspace/brava_flow_roi_18mlmin_20260928`。2026-09-29 的初始索引仅核对下表涉及的 17 项服务器入口存在；该次检查没有启动、停止或变更服务器任务。随后运行的独立验证工程见第 15 节。服务器路径不是本机文件链接，须先通过 SSH 登录。

| 服务器绝对路径 | 功能 |
|---|---|
| `/workspace/brava_flow_roi_18mlmin_20260928` | BraVa 服务器总目录；对应本地 computation_examples/brava_flow_roi_18mlmin。 |
| `/workspace/brava_flow_roi_18mlmin_20260928/scripts` | CFD 几何适配、启动、后处理和流场渲染调度。 |
| `/workspace/brava_flow_roi_18mlmin_20260928/gpu_solver_fix/svmultiphysics` | 本次真正执行的独立 GPU 修复版 CFD 二进制；本地主要保留源码补丁、构建身份及日志。 |
| `/workspace/brava_flow_roi_18mlmin_20260928/cases/balanced_pressure_final` | 正式 CFD 配置、结果、checkpoint 和日志。 |
| `/workspace/brava_flow_roi_18mlmin_20260928/mesh/SV_MESH` | 真实求解网格与五类边界。 |
| `/workspace/brava_flow_roi_18mlmin_20260928/visualization/candidate_0/OPEN_RESULTS.html` | 服务器流场预览入口。 |
| `/workspace/brava_flow_roi_18mlmin_20260928/population_v1/source/particle_3d/src/particle_3d/population_v1/run_fixed_window_population.py` | 服务器当前固定窗生产模块；其父目录是实际 population 核心。 |
| `/workspace/brava_flow_roi_18mlmin_20260928/population_v1/source/particle_3d/tests/population_v1` | 服务器 CPU/GPU 测试。 |
| `/workspace/brava_flow_roi_18mlmin_20260928/population_v1/env/bin/python` | population 独立运行环境，包含服务器 CUDA Torch 依赖。 |
| `/workspace/brava_flow_roi_18mlmin_20260928/population_v1/runs/realization_0000_contactfix` | 已有 1565 条群体原始生产结果。 |
| `/workspace/brava_flow_roi_18mlmin_20260928/population_v1/fixed_window_protocol` | 新协议、测试、合同和报告；没有新正式大批量结果。 |
| `/workspace/brava_flow_roi_18mlmin_20260928/population_v1/cfd_wall_audit` | 近壁场审计和轨迹风险关联。 |
| `/workspace/brava_flow_roi_18mlmin_20260928/population_v1/candidate_0/views_0to1p2s_large_axes` | 最新四视角大字号动画。 |
| `/workspace/brava_flow_roi_18mlmin_20260928/supervisord.conf` | 本算例既有 supervisor 配置，区别于实例系统服务配置；本次仅核实存在。 |
| `/root/particle8_2_runs/env/bin/python` | 既有 CFD 数据处理和流场渲染科学 Python。 |
| `/venv/main/bin/python` | 服务器预装 Torch/CUDA 环境，population 隔离环境继承相关依赖。 |
| `/workspace/particle9a5_formal_trajectories_20260925T080115Z_dt1ms` | 旧受保护单泡物理及 SonoVue 参考源；不是 BraVa 当前群体结果目录。 |

服务器校准仍位于 `/workspace/brava_flow_roi_18mlmin_20260928/cases/calibration_gpu8_production/`，与本地迁入 `brava_calibration_equal_split/flow/` 对应。源代码的模块搜索路径应指向 `/workspace/brava_flow_roi_18mlmin_20260928/population_v1/source/particle_3d/src`。

本地科学 Python：[/home/lzy/projects/computation_examples/_shared/python_environment/bin/python](</home/lzy/projects/computation_examples/_shared/python_environment/bin/python>)，实际目标为 [/home/lzy/projects/temp_storage/ulm_particle_3d_particle0/.venv/bin/python](</home/lzy/projects/temp_storage/ulm_particle_3d_particle0/.venv/bin/python>)。本地核心模块搜索路径是 [/home/lzy/projects/computation_examples/_shared/particle_workflow/particle_3d/src](</home/lzy/projects/computation_examples/_shared/particle_workflow/particle_3d/src>)。GPU 运行和 EGL 渲染使用已验证服务器环境；本地可用的 NumPy/PyVista 不等于有可用 CUDA Torch。

## 13. 阅读与使用顺序

- **只找网格/速度/压力/WSS**：第 3.2 节网格 → 第 5 节冻结场/后处理 → 第 6 节动画。
- **只找轨迹及统计**：第 8 节正式 run 和 analysis；展示用第 9 节。筛选后画面不是完整群体计数。
- **继续开发微泡模型**：第 7 节共享核心、固定合同和测试；先读固定协议报告，再参考旧生产验证报告。
- **分析近壁慢速/突然加速**：第 10.3、10.4 节；原始运行身份、单元/三角形 ID 和时间语义需一起核对。
- **继续新网格/输运收敛验证**：第 15 节；先看开发状态、收敛评估和暂停恢复记录，避免重复启动活动队列。
- **追踪 CPU/GPU 执行证据**：CFD 的 execution/PETSc profile、微泡的 SUMMARY/backend_counts 和 telemetry；不要用 GPU 利用率单独证明某项算法在 GPU。

旧 `README.md`、`REPRODUCE.md` 和原验证报告为阶段记录，其中覆盖停止参数、旧动画风格或旧路径可能已过时。尤其 `population_v1/REPRODUCE.md` 中旧 `--minimum-exits/--post-infusion` 命令不适用于当前代码；新生产入口以固定时间协议报告为准。上述历史文件保留原字节，本文负责提供当前索引，不批量改写冻结证据。

旧 `/home/lzy/projects/brava_flow_roi_18mlmin`、`/home/lzy/projects/ulm_particle_formal_p9a5` 等路径目前通过兼容链接可达。新工作优先采用本文 `computation_examples` 路径；不要删除这些兼容链接而不先处理既有配置依赖。

## 14. 本次文档核对范围

2026-09-29 初始索引逐项检查了当时 **276 个不同本地文件/目录**，并通过 SSH 检查上述 **17 项服务器入口**存在。参数和状态来自实际 ACTIVE_FLOW、purpose、boundary_summary、SUMMARY、固定协议报告及渲染清单；该次检查没有重新进行数值求解，也不将文件存在性检查称为数值验证。第 15 节为后续开发新增索引，不包含在上述历史计数内。

## 15. 新网格、CFD 与微泡输运验证（2026-09-30 已恢复开发）

独立工程根目录：[/home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/mesh_particle_validation_v2](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/mesh_particle_validation_v2>)。服务器对应 `/workspace/brava_flow_roi_18mlmin_20260928/mesh_particle_validation_v2`，使用该目录自己的 `supervisord.conf` 和 `supervisor.sock`，不能与父目录或实例系统服务混用。

冻结的46个出生参数用于验证，不是新的正式大批量群体。出生时刻、位置、半径、取向和ID相同；初速度及角速度仍按原规则从各候选场采样。截至2026-09-30 22:44 UTC，已完成10套新CFD及11批46条单泡重放（含同协议M0），结果均已回收并核验。25 µm coplanar细档已经完成：41出流、5观察截止、0数值失败；50/25 µm完整空间比较FAIL，0/46通过全部逐轨迹门槛，12条两档均出流者改变出口。20 µm候选M1_coplanar_fine_h020网格QA已通过，CFD从零场运行中，尚无其合格流场或轨迹。活动状态以[最近监控](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/mesh_particle_validation_v2/monitoring/LATEST.json>)中的UTC时间为准。旧暂停与恢复均保留在后文时间线；不能重启已关闭队列。不得用零数值失败、全局守恒或CPU/GPU重现替代局部场与输运收敛，ACTIVE_FLOW未变，dt/P2阶段未启动。

| 文件/目录入口 | 内容与适用范围 |
|---|---|
| [OPEN_RESULTS.html](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/mesh_particle_validation_v2/OPEN_RESULTS.html>) | 本轮中文报告、300 dpi PNG/PDF 与候选状态入口。 |
| [FINAL_VALIDATION_REPORT_ZH.md](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/mesh_particle_validation_v2/FINAL_VALIDATION_REPORT_ZH.md>) | 随已完成数据更新的科学报告；文件名中的 FINAL 不表示活动计算已完成或 gate 已通过。 |
| [DEVELOPMENT_STATUS_ZH.md](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/mesh_particle_validation_v2/DEVELOPMENT_STATUS_ZH.md>) | 各阶段实际状态；正式时间步与新 FEM 阶段受前序 gate 限制。 |
| [MESH_EXPERIMENT_PLAN_ZH.md](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/mesh_particle_validation_v2/MESH_EXPERIMENT_PLAN_ZH.md>) | M1/M2/M3 设计、资源上限、执行次序与验收规则。 |
| [00_baseline](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/mesh_particle_validation_v2/00_baseline>) | 原网格/场/源码身份、两个异常目标、固定空间采样与 46 条队列；不覆盖旧生产目录。 |
| [01_mesh_capability](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/mesh_particle_validation_v2/01_mesh_capability>) | 服务器实际 mesher API、尺度与边界层审计。 |
| [02_mesh_candidates](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/mesh_particle_validation_v2/02_mesh_candidates>) | 候选网格与各自 QA；失败尝试保留，存在目录不等于通过。 |
| [03_cfd_runs](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/mesh_particle_validation_v2/03_cfd_runs>) | 每候选求解输入、日志、检查点、执行身份和独立后处理。 |
| [04_field_audit](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/mesh_particle_validation_v2/04_field_audit>) | 原始逐单元 G/div/η/rQ/χ、固定区域与取样比较，另有明确标记的只读诊断。 |
| [05_particle_replay](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/mesh_particle_validation_v2/05_particle_replay>) | MB-MB 关闭的确定性重放、全部接受状态、原始绝对时间比较与事件记录。 |
| [MESH_CONVERGENCE_ASSESSMENT.json](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/mesh_particle_validation_v2/05_particle_replay/MESH_CONVERGENCE_ASSESSMENT.json>) | 完整中细档收敛状态；待完成项与失败项明确保留，不自动选生产场。 |
| [M2_MEDIUM_FINE_FIELD_REVIEW_ZH.md](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/mesh_particle_validation_v2/M2_MEDIUM_FINE_FIELD_REVIEW_ZH.md>) | M2 固定位置 u/G/WSS 比较及入口离散剖面、入口邻壁峰值分析。 |
| [M2_MEDIUM_FINE_TRAJECTORY_REVIEW_ZH.md](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/mesh_particle_validation_v2/M2_MEDIUM_FINE_TRAJECTORY_REVIEW_ZH.md>) | 两档 46 条均出流但仅 1 条通过完整门槛、6 条改变出口；不把零失败改称收敛。 |
| [M2_PATH_FIELD_DECOMPOSITION_ZH.md](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/mesh_particle_validation_v2/M2_PATH_FIELD_DECOMPOSITION_ZH.md>) | 在相同记录坐标上交叉采样，区分同位置场差与已经分离的路径所采样的场差；代数分解不是因果比例。 |
| [M2_PATH_GRADIENT_REVIEW_ZH.md](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/mesh_particle_validation_v2/M2_PATH_GRADIENT_REVIEW_ZH.md>) | 46 条轨迹 × 两条基准路径的速度/梯度 RMS；沿途统计不替换固定区域 gate。 |
| [M2_EVENT_INPUT_TIMING_ZH.md](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/mesh_particle_validation_v2/M2_EVENT_INPUT_TIMING_ZH.md>) | 微泡记录速度与实际区间起点输入的时序核对，区分筛查触发、此前减速与后续恢复。 |
| [M2_VELOCITY_COMPONENT_REVIEW_ZH.md](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/mesh_particle_validation_v2/M2_VELOCITY_COMPONENT_REVIEW_ZH.md>) | CUDA 重现 ID301/636 的 11165 个原接受区间并分解近壁求速；不重新积分或宣称物理准确。 |
| [M3_MEDIUM_FINE_FIELD_REVIEW_ZH.md](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/mesh_particle_validation_v2/M3_MEDIUM_FINE_FIELD_REVIEW_ZH.md>) | M3 中细档固定点 u/G/WSS、域外记录与原异常中心审查；场部分 FAIL，完整输运 gate 随重放更新。 |
| [M3_MEDIUM_FINE_TRAJECTORY_REVIEW_ZH.md](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/mesh_particle_validation_v2/M3_MEDIUM_FINE_TRAJECTORY_REVIEW_ZH.md>) | 中细档0/46通过；两档各45出流但删失对象不同，10条两档均出流的轨迹改出口。 |
| [FIXED_PLC_FIELD_REVIEW_ZH.md](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/mesh_particle_validation_v2/FIXED_PLC_FIELD_REVIEW_ZH.md>) | 原边界精确保留后的同坐标 CFD 对照：零速结构消失，但散度与梯度仍需中细档及输运收敛检查。 |
| [FIXED_PLC_TRAJECTORY_REVIEW_ZH.md](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/mesh_particle_validation_v2/FIXED_PLC_TRAJECTORY_REVIEW_ZH.md>) | 已关闭轨迹原始末状态核对，分别报告观察截止、出流与失败；含 ID440/1505、301/636 速度和壁间隙，细档完成后更新中细档门槛。 |
| [FIXED_PLC_WSS_REVIEW_ZH.md](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/mesh_particle_validation_v2/FIXED_PLC_WSS_REVIEW_ZH.md>) | 固定几何下两处相邻壁面 WSS、对壁高度和切向速率交叉核对；峰值在原高 η 区域，不在入口。 |
| [WALL_REPRESENTATION_REVIEW_ZH.md](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/mesh_particle_validation_v2/WALL_REPRESENTATION_REVIEW_ZH.md>) | 全部 417292 个 M0 记录状态的等价壁面查询；ID393 双壁近等距法向切换及同输入瞬时求速诊断，含 gap/normal 独立交换。 |
| [NATIVE_OUTPUT_IDENTITY_REVIEW_ZH.md](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/mesh_particle_validation_v2/NATIVE_OUTPUT_IDENTITY_REVIEW_ZH.md>) | 新coplanar网格的Float32原生输出坐标身份审查；精确编码/拓扑对应、反例拒绝、检查点连续性与派生场身份。 |
| [INTERNAL_P1_ZERO_FACES_ZH.md](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/mesh_particle_validation_v2/INTERNAL_P1_ZERO_FACES_ZH.md>) | 排查内部三WALL点插值零速面；M0为2，全部10个通过QA候选为0，不作为物理精度认证。 |
| [COPLANAR_MESH_REVIEW_ZH.md](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/mesh_particle_validation_v2/COPLANAR_MESH_REVIEW_ZH.md>) | 同原分片平面域的面内与体网格细化实际尺度；两档 QA 通过仅允许进入 CFD，未认证收敛。 |
| [COPLANAR_VALIDATION_REVIEW_ZH.md](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/mesh_particle_validation_v2/COPLANAR_VALIDATION_REVIEW_ZH.md>) | 原分片平面保持后的同位置场、冻结微泡终态与全部接受状态速度/间隙曲线；已关闭50/25 µm完整比较FAIL；20 µm后续比较另行生成。 |
| [COPLANAR_PLC_PROBE_PLAN_ZH.md](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/mesh_particle_validation_v2/COPLANAR_PLC_PROBE_PLAN_ZH.md>) | Phase 2G 后续面内/体网格细化能力试验，保持输入 PLC 与端口，原生输出须独立 QA；不自动 CFD。 |
| [SOLVER_RESOURCE_REVIEW_ZH.md](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/mesh_particle_validation_v2/SOLVER_RESOURCE_REVIEW_ZH.md>) | 已关闭 PETSc 执行段的 CPU/GPU 事件耗时与 LU rank 不均衡；暂停、重启验证和续算分别保留，嵌套事件不相加，不作为物理精度证明。 |
| [FIXED_PLC_REFINEMENT_RESULTS_ZH.md](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/mesh_particle_validation_v2/FIXED_PLC_REFINEMENT_RESULTS_ZH.md>) | 原边界三角面完全保留的局部体网格设计及 QA；网格通过不等于 CFD/轨迹通过。 |
| [MOUSE_HUMAN_MESH_SCALE_COMPARISON.md](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/mesh_particle_validation_v2/MOUSE_HUMAN_MESH_SCALE_COMPARISON.md>) | 历史小鼠与 BraVa 的相对尺度对照，不恢复旧小鼠计算。 |
| [RESUME_REPORT_ZH.md](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/mesh_particle_validation_v2/RESUME_REPORT_ZH.md>) | 安全暂停、原生重启连续性验证及后续开发时间线。 |
| [PAUSE_RESUME_ZH.md](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/mesh_particle_validation_v2/PAUSE_RESUME_ZH.md>) | 已有安全暂停/恢复过程；停止队列与活动 CFD 必须区分。 |
| [monitoring/LATEST.json](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/mesh_particle_validation_v2/monitoring/LATEST.json>) | 最近同步的服务器监控快照；本地副本不是实时查询，须看其 UTC 时间。 |
| [scripts](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/mesh_particle_validation_v2/scripts>) | 本轮开发脚本；`executed_versions/` 保存带 SHA256 的版本归档，实际执行身份由各产物记录。 |
| [CURRENT_STUDY_SCRIPT_DEPLOYMENT.json](</home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/mesh_particle_validation_v2/CURRENT_STUDY_SCRIPT_DEPLOYMENT.json>) | 当前脚本与服务器归档校验记录；不替代各科学结果的执行源码身份。 |

继续工作时先完整读取 [/home/lzy/projects/BRAVA_CONTEXT_HANDOFF.md](</home/lzy/projects/BRAVA_CONTEXT_HANDOFF.md>)，再核对服务器当前队列与结果身份。阶段未通过时保留失败证据，不能自动跳到时间步认证、P2/P1 或 exact-divfree 后端。


### 2026-09-30 14:45 UTC BraVa暂停更新

BraVa网格/粒子验证已按用户指令安全暂停：M1_coplanar_fine原生第37步检查点已验证并备份，本研究服务器科学任务全部停止。已完成9套新CFD、10批冻结46初态回放；完整中细档收敛尚未通过。最新完成结果/未完成任务及恢复边界见 [PAUSE_RESUME_ZH.md](computation_examples/brava_flow_roi_18mlmin/mesh_particle_validation_v2/PAUSE_RESUME_ZH.md)，机器记录PAUSED_STATE_20260930T144554Z.json。恢复必须等待新指令，旧队列不可直接重启。


### 2026-09-30 16:05 UTC：恢复细档验证，保留重启门槛FAIL并改用更早检查点

用户已明确恢复开发。原细档30→37重放native exit0、u/p差约4e-15，但A相对L2约1e-7，未通过1e-8门槛；不是新的微泡失败。原FAIL和24个闭合执行文件已保存在本地，完整诊断见 `mesh_particle_validation_v2/FINE_RESTART_CONDITIONING_ZH.md` 与 `figures/restart_A_conditioning.{png,pdf}`。

当前新任务 `check_coplanar_fine_restart_30` 正从原20步验证至原30步，保持原u/p与逐rank Y/A门槛不变；只有PASS才能从原30步独立续算，不放行原37步恢复。技术审查 `RECOVERY_REVIEW_FINE30.json` 记录预算：扣除原段7422.862194秒和第一次验证1017.505819秒，第二次验证开始前剩13159.631987秒，仍共用原6小时。新队列 `fine_user_resume_followup` 尚未启动，等待该门槛。

截至16:13本地报告为16:09:36UTC版：167链接、96图文件SHA核验PASS，报告入口 `mesh_particle_validation_v2/OPEN_RESULTS.html`。80脚本快照b761d50efd3a385e7540dd57f3aebeb6b621478f49fe01af6f9d343e9ac2aece已部署且当前/归档源核验PASS。M3旧等待回放文字已修正为实际完成FAIL。当前仍9套新CFD/10批冻结46单泡回放，无新生产网格或正式大批量；时间步、P2/P1、严格无散度阶段均未启动。实时任务状态读服务器monitoring/LATEST.json，本地副本不是自动推送。


## 2026-09-30 16:55 UTC：最新用户指令为暂停，现已安全停止并回收

本条取代前述活动恢复状态。用户16:48:55 UTC要求在下一个可停止点暂停；原生续算段user_pause_resume_from30_0001已完成31至37步并退出0，health_failures=[]。fine_user_resume_followup和report_watch均STOPPED，全部研究队列STOPPED/EXITED，GPU计算进程为空；实例保留。第37步Y/A完整检查点SHA 276f72b5a23ec2203aa6829026c0322017cb929b3959f38e7251d2c1b1c79fe7，26个封闭文件本地回收哈希和独立checkpoint_audit均PASS。原37步文件、首次30→37 FAIL、20→30 PASS全部保留。

本续算耗时1042.891990184784秒，原21600秒总预算现剩10567.13324713707秒，休眠等待不计入。当前仍未稳态合格；原execution.json保留FAIL，不能改称成功或微泡失败。完成量仍9套新CFD、10批46单泡回放，无细档正式场或细档微泡结果，无新生产网格，dt/P2/无散度阶段未开始。

最新暂停记录：mesh_particle_validation_v2/PAUSED_STATE_20260930T165531Z.json；可读总结mesh_particle_validation_v2/PAUSE_RESUME_ZH.md；图形入口OPEN_RESULTS.html。CURRENT_EXECUTION_STATE为PAUSED_BY_USER，等待新的明确恢复指令，禁止自动重启任何队列。恢复前需审查本次37检查点的连续性路径和剩余预算，不直接复用现有段或旧启动器。


## 2026-09-30T17:56:09.883833+00:00：用户再次恢复，复用已验证第30步的新续算已运行

最新指令“恢复开发，在开发停止之前不要停止 thinking 的信息显示”解除16:55暂停。服务器fine_resume_followup_0002已于17:54:47 UTC启动，report_watch运行；旧fine_user_resume_followup及其他暂停队列保持STOPPED。新段为M1_coplanar_fine/segments/user_pause_resume_from30_0002，从字节未改动且20→30原门槛已PASS的原第30步检查点恢复。没有重跑验证，没有放行新旧第37步检查点；两个原生暂停记录及30→37 FAIL全部保留。

新启动预算10567.13324713707秒，原6小时扣除原执行、两次检查和第一续算1042.891990秒。31→37重复计算仍扣预算。技术依据RECOVERY_REVIEW_FINE30_0002.json、授权RESUME_AUTHORIZATION_FINE30_0002.json，正向预检及放宽门槛/漏算已执行段/重置预算三个负向控制均PASS，276个本地保护输入SHA一致。实际XML只改restart=true、仍500步、8MPI+CUDA，PETSc/物理边界/科学模型不改。83脚本快照0414b516cb4341fb3c2b5d93c2374b64a351f7f12670bcfbc070b1fa5685664e服务器当前及归档SHA/AST全PASS。

新入口scripts/resume_validated_fine.py、fine_resume_followup.py串行推进稳态、后处理、精确原生/规范场身份、全场审计、46初态单泡重放和分析。只有逐项成功才前进；不自动发布生产网格或进入dt/P2。report_fine_restart_recovery与render_report已支持多次暂停续算的预算和历史耗时，不隐藏上一续算成本。当前仍9个合格新CFD、10批回放；完整空间gate尚未完成。持续监控脚本monitor_progress.py；暂停时应先冻结新producer，再向当前新段run/STOP_SIM请求原生边界停机，并停止report_watch及后续队列。不要把旧段当成活动段。


## 2026-09-30T19:30:10.428855+00:00：细档稳态与全场审计关闭，46条回放正在运行

M1_coplanar_fine的segments/user_pause_resume_from30_0002已在第81步原生exit0，health_failures=[]，execution PASS；40/50/60/70/80五连续区间满足原稳态门槛。续算耗时4899.31183719635秒，原6小时预算扣除全部原执行、两次验证、两次续算后剩5667.82140994072秒。后处理与独立规范/原生输出身份核验PASS，规范flow SHA d049122a0e77c001a1d70e0b1ef09a88faae7b634b5b0bcb99e461cc55bf168a。原暂停execution FAIL与30→37 FAIL证明仍原样保留，不追溯改写。

全场审计已完成：1064977 tetra，严格零速单元0，1059固定点均在域内；CPU/GPU G最大绝对差9.095e-13 s⁻¹。此数值一致性不等于局部物理精度。中细档固定点u相对L2为1.9433%，G为19.3892%，WSS为20.9880%；部分局部区域u也超过2%，空间场门槛未过。原388769中心细档|div|约86.119 s⁻¹，η约0.037544，相比原场改善但区域梯度/WSS仍不稳定。

181个已关闭CFD/审计文件已本地回收并SHA PASS，见CLOSED_FINE_0002_LOCAL_COLLECTION_VERIFIED.json。当前合格新CFD计10套，已完成分析回放仍10批；第11批M1_coplanar_fine_dt0.00050000正运行，固定46初态/0.5ms/关闭泡间作用，不能提前称其完成或无失败。fine_resume_followup_0002与report_watch继续运行。原正式生产流场、1565条旧协议population和显示参数均未更改；dt/P2阶段未启动。

若用户再次暂停，原生CFD已经完成，不要重复从30步启动，也不要改写已PASS执行；先核对当前回放/分析安全边界。下一步关闭回放和分析、收集轨迹、复核完整空间gate，并按实际资源评估进一步mesh refinement。新增只读资源测算位于01_mesh_capability/COPLANAR_NEXT_REFINEMENT_FEASIBILITY，不生成候选/不运行CFD；25µm表面数组复现完全一致，20/18.75/12.5µm仅计算候选表面规模与非上界资源估计。


## 2026-09-30T19:50:01.761401+00:00：25 µm细档全部关闭，20 µm有界网格实验已启动

M1_coplanar_fine第81步稳态CFD、规范输出身份、全场审计及冻结46单泡回放/分析均完成。细档回放41出流、5观察截止、0数值失败、0有支撑静止；O1/O2/O3=16/17/8，截止ID5/393/534/1329/1377。20次固定阈值slow–fast事件（中档25次）。ID440 O2停留0.963458s，ID1505 O1停留0.845126s，ID636 O3停留1.800663s；ID1189 O3停留9.937649s，ID393/1329观察截止。原旧1565结果和协议身份均不改。

完整50/25µm比较FAIL：全部逐条门槛通过0/46；两档均出流但出口改变12条（123,168,175,341,487,542,650,708,1077,1208,1240,1505），终止类别改变6条（4,5,534,1329,1449,1522）。速度历史门槛46条均未通过；39条终点差、39条持续时间差、27条最小间隙差、11条路径长度差超门槛。固定点u/G/WSS相对L2=1.9433%/19.3892%/20.9880%，无正式网格选择，dt/P2仍未启动。不能把无数值失败称为空间收敛。

376个细档回放/分析文件回收SHA PASS，全部ID及终止记录交叉核验PASS，见CLOSED_FINE_REPLAY_0002_LOCAL_COLLECTION_VERIFIED.json。当前完成量为10套新CFD、11批46轨迹回放（包括同协议M0）。总报告19:39:59UTC版本地171链接/98图SHA PASS，REPORT_MANIFEST SHA d6f5173056b550fd793de532e026d99c9c6c8f4baa16542d7cba22d1754af116。新细档曲线已人工查看，无平滑/时间对齐；完整解释见COPLANAR_VALIDATION_REVIEW_ZH.md。fine_resume_followup_0002与report_watch均已正常EXITED，不能再次启动这些独占输出入口。

依据用户Phase2G继续可行空间细化，19:48:41UTC启动独立coplanar_h020_mesh_only。候选M1_coplanar_fine_h020，对照M1_coplanar_fine；只生成及QA，不自动CFD。表面272904面/136454点，原顶点和端口保留；显式表面准备资源上限30万，旧默认25万不变，所有科学门槛不变。原生TetGen已输出295691点/1412899 tetra，尚待完整拓扑与严格域对应QA和实际CFD资源复核；不可先计作新合格网格或CFD。

依据COPLANAR_H020_MESH_REVIEW.json（SHA49dad3a905716e457dc8c05c27fec4fa4418f9183abcc235b8506fd174982ea3）、COPLANAR_H020_MESH_PLAN_ZH.md；资源预测及25µm输入逐数组复现见COPLANAR_NEXT_REFINEMENT_RESOURCES_ZH.md。prepare_coplanar_plc_probe只新增显式表面资源上限参数，audit_coplanar_transport_mesh新增独立plan入口及源记录；新coplanar_mesh_probe_followup按源/证据锁定串行执行并禁止自动CFD。84脚本快照fd04695f582befc62430fe18de1932bb32e477abc74b69a005d59d6bff50ee3a当前/归档AST及SHA全PASS。旧83快照及旧执行源码保留，粒子/CFD物理未修改。

当前活动任务只有coplanar_h020_mesh_only。如收到新的暂停指令，先核对该网格队列原生生成/QA边界，不要向已关闭CFD写停止标记或从30步重复恢复。网格阶段结束后检查实际质量、规模和内存再决定下一步；20µm未运行CFD，不存在它的流场或轨迹结果。


## 2026-09-30T20:31:04.066971+00:00：20 µm网格通过，新的原始零场CFD正在运行

最新用户指令仍为恢复开发。20 µm候选M1_coplanar_fine_h020已完成网格及严格域对应QA：296162节点、1414875四面体、全WALL零速风险单元0，负/零体积和非流形面0；原分片平面壁面和4端口严格保留。原生Float32 probe的FAIL证据保留，输入精度恢复及471次完整内部边扇分割后最终GATE_MESH PASS。minSICN=0.01165949，比25 µm的0.01811更低；不能只因更细就声称更准。46个封闭网格文件本地SHA核验PASS，见CLOSED_H020_MESH_LOCAL_COLLECTION_VERIFIED.json。

资源复核、85源文件和30证据哈希预检及预算/漏源/错网格/放宽门槛四项负向控制均PASS后，coplanar_h020_cfd_followup于20:18:56.261662 UTC启动。实际case为03_cfd_runs/M1_coplanar_fine_h020（根case，从零场开始，不是旧25 µm续算段）。实际XML/PETSc与正式基线逐字节相同：8MPI+CUDA、P1/P1、dt=0.01、500步、21600秒、5连续稳态区间；压力及微泡模型未变。当前仍首次CPU ASM/LU分解，尚无合格h020流场或轨迹。20:29实际求解器RSS约34.44 GiB，cgroup总42.05 GiB含约7.61 GiB可回收文件缓存，OOM/oom_kill均0，继续监测实际峰值而不把预测当上界。

当前85脚本快照SHA ca522c5f7e9b654306c4d5b49302c85ac2283e41b82ea9ca9dbeb7120d25da33；CFD review SHA f248f5c8092ddcaae934791a489e6f4ae151a2ac5819853ec3bf95020475efb2。见COPLANAR_H020_CFD_REVIEW.json/.md（可读文件为COPLANAR_H020_CFD_REVIEW_ZH.md）、H020_CFD_PREFLIGHT_CONTROLS.json、H020_CFD_STARTUP_LOCAL_VERIFIED.json。启动配置10文件及12个输入哈希已收回本地核验。当前运行producer引用的85脚本不得运行中修改。

队列仅在CFD稳态PASS后顺序执行后处理、独立规范/原生输出身份核验、全场审计、冻结46初态0.5ms回放、分析及空间比较。当前完成科学量仍10套新CFD/11批回放，25 µm已全部结束，50/25 µm完整空间gate仍FAIL（逐轨迹全门槛0/46）。无新生产网格、无dt/P2/无散度试验。20:16:49报告（启动前静态快照）已本地验证171链接/98图SHA PASS，manifest SHA dbd00dfdc05c7ed80afb52cf56ee704df323fa5fa3056f7cdb065ec34e30fb06；其active_CFD=[]只代表生成时刻，实际运行状态见monitoring/LATEST.json与CURRENT_EXECUTION_STATE.json。

只有coplanar_h020_cfd_followup活动，fine_resume_followup_0002/report_watch/coplanar_h020_mesh_only已EXITED，不能重启独占输出旧入口。若收到新的暂停指令，应核对当前阶段并冻结后续推进，向活动h020根case/run/STOP_SIM请求原生安全边界退出，保留检查点与已耗预算；不能用旧25 µm暂停辅助原样操作，不能恢复已关闭25 µm或重置6小时。新的h020检查点须另行审查连续性，旧25 µm证明不可套用。


### 2026-09-30 20:53 UTC 活动监控补充

20 µm CFD已完成第1时间步（14次非线性迭代，末Ri/R0=4.530e-11）并进入第2步；首步观察原文与哈希见H020_FIRST_STEP_OBSERVATION.json。两次CPU LU构建已通过，第二次由ITERATION_GROWTH触发；近时cgroup总内存约38.68GiB，solver进程树RSS约37.42GiB（共享页可能重复计数），OOM/oom_kill=0。未有稳态资格、未有h020流场或轨迹，继续coplanar_h020_cfd_followup，不重启旧任务。

只读原生日志复核：SV13Q age在同一时间步的多次logical求解间保持不变，进入下一时间步才增长；例如25 µm原段logical3至14均step1/age0，logical15为step2/age1。因此不能把max_age=5口头写成“每5次线性求解必重建”；还存在ITERATION_GROWTH等提前触发条件。以上是实际日志观察，不修改任何策略。

待当前85脚本producer完全退出后再部署的本地草稿在pending_development/：report_coplanar_h020_validation.py生成独立25/20 µm比较及6张PNG/PDF；report_fine_restart_recovery.py修正已合格25 µm仍显示“尚需完成”的旧文字并列实际预算。均只做报告，不改科学门槛或结果，尚未部署/运行。两个brava_collect_closed_h020_*.py是对应封闭阶段的独占回收脚本，已做语法检查，未执行；执行前必须核对阶段和manifest不存在。旧25 µm回收脚本不能直接用于新case。


### 2026-09-30 21:22 UTC 首小时资源图已核验

20 µm CFD已完成前3时间步，当前第4步迭代，3次预条件器构建已完成，约38.8GiB总内存且OOM/求解健康告警均0。首小时20:18:56至21:18:56的60份监控快照已封闭保存为monitoring/H020_FIRST_HOUR_STATUS.jsonl（SHA9fcf480688ed8bdf285e9aaba811d30beeca4d6a169e241935a39dcb5559b936），本地只读生成资源PNG/PDF并人工查看，7个文件本地/服务器SHA PASS见monitoring/H020_FIRST_HOUR_RESOURCE_ARTIFACTS_VERIFIED.json。资源脚本位于monitoring/report_h020_first_hour.py，不在运行中85脚本producer目录，不改变其身份。

图文入口monitoring/H020_FIRST_HOUR_RESOURCE_REVIEW_ZH.md。CPU区间加权7.680012核（配额7.68），采样cgroup峰值42.055626GiB（限制42.055664），OOM/kill增量0；包含文件缓存/共享内存，不等于solver私有RSS。GPU是稀疏瞬时样本，不推算平均利用率、时间占比或完整瓶颈。此图只涵盖首小时，不声称完整CFD内存上界、稳态资格或输运改善；主报告98图manifest未因此修改。


## 2026-09-30 22:11 UTC：20 µm首个保存状态通过完整性核验，CFD继续

第10步已完成并进入第11步。原生stFile_010.bin的8-rank Y/A检查及result_010.vtu有限性、壁面无滑移核验PASS，两个原文件及证明已回收本地SHA PASS，见03_cfd_runs/M1_coplanar_fine_h020/reports/SAVED_STATE_010_VERIFIED.json与SAVED_STATE_010_LOCAL_COLLECTION.json。检查点SHA e1e56de0118089ab487ab5dfdb70edde22f9832ce18f78c404e68c5046c2edf8。这只是保存状态完整性，绝不声称重启连续性或稳态已通过。

第10步time=0.1s：入口18mL/min，出口归一化分流0.3332182392/0.3331620949/0.3336196659；epsilon_mass=0.0009765625，尚大于1e-6稳态物理门槛。压力0.439712…2166.73231Pa，umax0.54997159m/s，壁面最大速度0。尚无可比较的保存间隔，连续稳态区间0/5。输出适配器读取无重试或失败；实际后处理/全场审计/微泡回放尚未开始。

活动任务仍仅coplanar_h020_cfd_followup，实际case仍根case，从零场连续推进。总内存约39GiB、OOM=0；4次预条件器构建均完成。实际源85快照不变；未修改物理、容差、预算或生产场。完成量仍10套新CFD/11批46回放。

后续保存状态回收辅助为/tmp/brava_archive_h020_saved_step.py（副本pending_development/brava_archive_h020_saved_step.py），参数10至500中的10倍数；等待原生进入下一步才读，首个状态interval=null合法。10已执行且输出独占，不要重跑。以后20/30等可独立回收，但这些证明都不是重启连续性认证。


## 2026-09-30 23:19 UTC：短暂停止、汇总、自动恢复均已完成；当前继续开发

用户23:00要求的“下个安全点短暂停止并汇总，然后恢复开发”已全部执行。h020第20步stFile_020.bin完整8-rank Y/A及原生流场核验PASS并回收本地；检查点SHA cdbc69781439fc20fe0c1d1c1303ab5796049a3909ca8914ec95d8ed03d9f16a，证明reports/SAVED_STATE_020_VERIFIED.json与SAVED_STATE_020_LOCAL_COLLECTION.json。step20质量误差9.536743159e−7、wall no-slip PASS；10→20区间Eu=0.0006024017628、EQ=0.0005338033864，仍高于1e−5/1e−6，连续稳态0/5。10/20已回收且独占，不重复执行collector。

在完整step20后挂起同一组12进程45.202254秒，23:18:02.622881 UTC挂起，23:18:47.825130 UTC恢复；当前MPI及8rank PID/start_ticks未变，日志已继续增长。没有native退出，没有检查点重启，没有STOP_SIM，没有改写budget或科学配置，短暂停顿计入原21600秒。审计monitoring/H020_BRIEF_STATUS_HOLD_20260930.json与H020_BRIEF_HOLD_RESUMPTION_VERIFIED.json，恢复及3个报告文件本地/服务器SHA一致。后者status=PASS_RESUMED_SAME_PIDS_LOG_ADVANCED。当前状态恢复RUNNING_REVIEWED_COPLANAR_CFD；不是PAUSED，也无需再等恢复指令。

用户要求的完成量、结果、后续和时间估计已通过commentary答复并写入PROGRESS_REVIEW_AND_NEXT_STEPS_ZH.md/.json（9个链接存在，SHA验证），已同步服务器。MD SHA1101e49cdab0f5284761a1c8520ac6b35ca196356e8ea8d1930d627c204ca919。当前阶段估计4–6h，空间通过后全剩余工作约6–10h；若预算触发并需重启验证，额外1–2h或更多；若需P2代表域则规划1–3工作日以上，非保证时限。ETA实测证据留在monitoring/H020_ETA_EVIDENCE_20260930T2304Z/，预算可能不足仍未修改。报告生成器首次使用相对--hold路径报ValueError，随后改用绝对路径运行成功；科学结果不受影响，实际使用时--hold传绝对路径。

23:19监控第21步已完成2次非线性迭代，CPU区间平均包含45秒挂起故短时4.54核，不能误判性能退化；恢复后的独立区间待继续采样。内存39.18GiB、health/read retry=0。继续原producer；完成量仍10套新CFD/11批46回放；h020未稳态、未回放。开发和持续进度指令仍生效。


## 2026-10-01 00:13 UTC：第30步检查点已收回；启动收敛图完成

h020第30步已完成、当前第32步NL2，health/read retry/OOM均0，内存39.18GiB，唯一活动producer仍coplanar_h020_cfd_followup。第30步stFile SHA cced807e052650ca93e87d35c1306ae6b07346aa18c3f99a566af54ac44d42b2，size20209024；8rank Y/A、流场有限性、无滑移和本地文件SHA均PASS，local collection 2026-10-01T00:10:20.304084UTC。服务器证明SHA5cd6c5d2122f5117d1049f7c3143c9e68dcab400f3a1f7a721d8bc9a53fb0016。10/20/30均已独占回收，不重跑对应collector；下一个40。此为保存完整性，不是重启连续性或合格流场认证。

20→30：Eu=7.0319738117395925e−6低于1e−5，EQ=1.1857101942204126e−6仍高于1e−6，连续稳态0/5，未改门槛。mass=9.313226075571546e−10、WALLu=0，umax=0.5499245861335055m/s。第30步pressure_range_pa=[0.4384423082164724, 2725.8372827409903]；outlet_fractions={'OUTLET_01': 0.33308154594874734, 'OUTLET_02': 0.33337173066055575, 'OUTLET_03': 0.3335467233906969}。

新增只读启动收敛图：monitoring/H020_STARTUP_STEP30/REVIEW_ZH.md、REVIEW.json、startup_convergence.png/.pdf。本地标准matplotlib生成、已目视检查，绑定三个完整保存状态及其本地回收证明；第10步缺相邻区间，保留缺失不补零，图不是合格稳态场。5文件本地/服务器SHA PASS见monitoring/H020_STARTUP_STEP30_ARTIFACTS_VERIFIED.json，源monitoring/report_h020_startup.py SHA580453e79dc8f7ef76ec5dd1402d06972623abd2734b47a872954cca263af490；PNG SHA272dbe1d5f137625ef3632aaaf0d9c004964e903a07fde4df47553a67155b93c。图在monitoring，不修改主报告98图manifest或活动85源快照。

用户23:00要求的短暂停止/汇总/恢复已完成，当前继续开发，无待确认事项。完成科学量仍10套新CFD/11批回放，h020未稳态；原21600秒预算未更改。当前进度报告PROGRESS_REVIEW_AND_NEXT_STEPS_ZH.md保持23:18阶段快照，最新进度读monitoring/LATEST.json。


## 2026-10-01 00:41 UTC — 20 µm继续运行，预算恢复仅为已验证保护的备用草稿

- 最新用户的“短暂停止→输出完成项/结果/待办/ETA→自动恢复”已于23:18完成（原PID暂停45.202254秒后恢复）。当前开发继续，不再等用户恢复；若收到新的暂停指令则按新指令执行。
- 当前仍只有 coplanar_h020_cfd_followup 运行，原6小时/500步时限和85个源文件不变。00:41监控：原生第38步预条件器第10次重建；37步完成，累计15730秒，CPU约7.679核，内存约37.06 GiB，健康标志全0，OOM/oom_kill=0。memory.max累计事件从746607增到747602，不能把累计历史事件全部归于本次；本次观察间隔新增995次，pressure累计新增约0.017秒。仍0/5稳态区间；10/20/30已收集，下一次40。
- 新增pending_development三个备用草稿（未调度、未启动、未修改活动代码）：recover_h020_wall_budget.py (91b8532ebc826527cded97932e88ba952704e7584e29bad99621d2f4f56f63cc)、h020_wall_budget_followup.py (dd513c9458aa790d10491ff27c02ae85139e0b199d8498fb42f3b18d4ed1cfc8)、review_h020_closed_budget.py (685cae6c74638bef9a7d8d04866aa3c87127034258b57d6201cacf6f0553fff3)。均AST PASS，在真实活动状态下执行保护检查均拒绝：原execution.json未关闭，未创建审查或新段、未启动CFD。
- 方案见pending_development/H020_WALL_BUDGET_RECOVERY_PLAN.json及同名_ZH.md；本地/远端8个草稿/证据SHA一致，证明H020_BUDGET_DRAFT_FILES_VERIFIED.json SHA916da85b22ce2dab74d34bbf19ac26b88d4dba5eb653d846065e42f0cb4076c0。正向预检和原生连续性尚未执行。
- 原用户附件Phase2明确允许随DOF增加记录后提高wall-clock allowance。备用拟增加6小时使h020累计原生额度12小时；这不是用户指定数值，也不是ETA/收敛保证。必须先等原producer完全退出，确认原本段确为预算FAIL、exit0/health[]，保存关闭原始资料、核对原生最终checkpoint，才运行review_h020_closed_budget.py --study <server study>生成独立H020_WALL_BUDGET_REVIEW.json；现在该文件不存在。严禁修改活动policy或绕过其STOP_SIM。
- 具体恢复方案用h020自身20→30重放，1小时上限、原u/p/各rank Y/A相对L2门槛1e-8；仅PASS后从原h020第30步续算，并重新满足5区间稳态。全部原段/验证/续算耗时累计扣除12小时，原FAIL不改写。25 µm的连续性证明不转用。失败不自动重试。
- 闭合审查创建后还须先实际正向preflight并核对资源，再调度pending_development/h020_wall_budget_followup.py --study <study> --review <absolute review>；它依次检查/续算/后处理/独立身份/场审计/冻结46条回放/分析/原有报告。未创建新的supervisor job。若初始段本身PASS，则完全不启用此备用路径。
- 已准备的/tmp/brava_collect_closed_h020_cfd.py仍只接受fresh root PASS；若使用续算路径必须先适配segment和本地SV_MESH兼容链接，不能直接运行旧草稿。报告整合草稿仍要等producer闭合后部署；85脚本当前禁止变更。


## 2026-10-01 00:58 UTC — h020第40步首次通过区间门槛，继续运行

40步已完成独立核验/本地回收（不再重复执行collector40）。原生cp SHA a0bc9fad418af43e2177c38a32775de7ec5dda876352a3911a5495a7b34f5735。30→40 E_u=1.4440512382013205e-7、E_Q=3.6555320171826055e-9；mass=9.093251619169604e-13，壁面u=0，u_max=.549924585365713 m/s。连续通过1/5，仍无合格h020场/轨迹。00:57:57原生进入43步重建第11个PC，原段elapsed16737秒，CPU7.677核/RAM35.05GiB，health全0/OOM0，memory.max累计747602未再增长。原6h政策不变，备用追加预算仍未激活。

新科学图表 monitoring/H020_STARTUP_STEP40/REVIEW_ZH.md + REVIEW.json + startup_convergence.png/.pdf 已运行、看图并双端SHA验证。证明 monitoring/H020_STARTUP_STEP40_ARTIFACTS_VERIFIED.json SHA5a9a35c16ac86b0413c87d302fdb921d6db3d78e2744536d30744c92dd89f795。新通用绘图程序 monitoring/report_h020_startup_series.py SHAf88b1c2e7ab38ee3d4b266e16a2eb6ac40e95e69ed89168556650189b401e5bf，可在各保存状态收集后用 --through-step 50/60/...生成新独立目录（已存在则拒绝覆盖）。原生比较为<=，图标明连续区间数，不单独认证稳态或局部物理精度；不修改原30步图。

时间估计已随实测更新：monitoring/H020_ETA_STEP40/REVIEW_ZH.md及ESTIMATE.json/monitor_snapshot.json，双端证明SHA cfe1f951967eced22a3448ead78087217f2917ab99594538f2f2ed7658f0625b。00:53快照原6h剩5156秒；若假定50..80均通过，到81步不中断尾段约7472秒；因此很可能预算停机。备用h02020→30验证约3090秒、原30→81续算约10677秒，加原段剩余后原生总剩约18923秒；含审计/46条回放/报告，本轮预计仍需5–7小时。假定未来NL次数像已关闭25µm、普通迭代26秒/重建额外661秒，非收敛/时限保证。若仍需P2/P1，全项目预留1–3工作日或更长。当前用户未再暂停。

新的通用关闭CFD回收草稿 pending_development/brava_collect_closed_h020_cfd_any_segment.py，/tmp同名，SHA3aeefd66b49d89f2423780d606be00cca7792188941fe899c85cbcc51e205e11；外层及嵌入远端probe均AST PASS、已同步但未执行。兼容fresh PASS或经过独立验证的预算续算PASS，校验原FAIL保留，回收全部原生段；只把已知远端绝对SV_MESH链接转为本地相对链接。优先使用此新collector；旧root-only版本保留。下一保存状态50，预算到点前继续原producer；不要把备用草稿当已运行或已正向预检。


## 2026-10-01 01:10 UTC — 备用恢复草稿输入保护补充（仍未启用）

pending_development/recover_h020_wall_budget.py更新为SHA6d66846b4504e10ab5eccb3541d4ee39125412caf3f85e2f9a08f287a62c7c42：在每个验证/续算阶段重新核验BUDGET_STOP_PRESERVED内所有原段文件（包括原生分区文件），并拒绝未计入账目的额外分段。其余followup/review_creator源码不变。旧9项草稿/保护证明完整保存于pending_development/previous_drafts/h020_budget_20261001T0108Z，两端均已归档。更新后的三个真实活动状态负向保护均PASS，8项草稿双端SHA一致，新的H020_BUDGET_DRAFT_FILES_VERIFIED.json SHA3adbc16d02651ca5bdbab32bc634e9158bf9ebed6950cfe07316d5a01bec1fd0。85个活动源文件、原6小时政策未变；正向闭合审查、重启连续性、续算均未执行。01:09:34主计算到46步、原生elapsed17431秒、health全0、仍1/5。若以后进入分段，/tmp/brava_h020_memory.py仍读取根case的progress，需要先适配RESUME_INDEX再用于新段资源监控；普通brava_monitor_live.py已经能识别分段。


## 2026-10-01 01:26 UTC — h020第50步已归档，连续稳态2/5

50步collector已执行成功，不重复：cp SHA0d8ef7f82dfcfb7770f8d8fadc218793cf0df98a0d7093e4720752ae36218bab；40→50 E_u=3.256860679840286e-9，E_Q=1.3692905220556472e-11，mass=7.058607893785836e-16，wallu=0，u_max=.5499245853645829。累计2/5，仍未正式合格/未回放。monitoring/H020_STARTUP_STEP50下MD/JSON/PNG/PDF已生成、开图检查并双端SHA验证，证明monitoring/H020_STARTUP_STEP50_ARTIFACTS_VERIFIED.json SHAef871879fa0ab24a854880d8b11ef5427f8237a690dc68b3f26a20728552e2ca。01:25:31监控进入52步，elapsed18391.804秒、CPU7.6805核、RAM39.164GiB、health全0，原6h预算仍有效。

下一保存节点60；10/20/30/40/50都已独立收集。当前用户要求的短暂停止/总结/自动恢复早已完成，开发仍ACTIVE。备用预算恢复代码版本以01:10记录为准，仍未激活；正向review_creator须用/root/particle8_2_runs/env/bin/python（负向测试在缺少closed execution时已提前拒绝，未验证科学依赖），先等原producer EXITED、原预算FAIL关闭且HEAVY锁可用。原runner只写STOP_SIM并等待原生步完成；若6小时点落在耗时约11分钟的PC重建中，实际关闭可能超出时限数分钟，全部真实耗时仍计入累计12小时。不能因deadline到达就假定已保存退出。


## 2026-10-01 01:58 UTC — h020第60步归档，连续稳态3/5，原预算将到期

60步collector已成功执行，不重复：cp SHA9f8b36202e705b29d2b4f5242fc06e5f09c911121508fa09cc5e03e21988207c；50→60 E_u=6.920730656668296e-11，E_Q=5.461597857816791e-14，mass=1.764651973446459e-16，wallu=0，u_max=.5499245853642368。累计3/5，不是正式合格场。10/20/30/40/50/60均已独立回收，下次70（也可能先触发6h预算停止）。01:57:37已进入63步重建第15个PC，elapsed20318.287秒，health全0，原6h剩约1282秒；不能在deadline到达时假定原生已经退出。

monitoring/H020_STARTUP_STEP60的MD/JSON/PNG/PDF已生成、看图并同步验证。最初排版的边缘标签靠近y轴刻度，初版完整保存至H020_STARTUP_STEP60_INITIAL_LAYOUT；新monitoring/report_h020_startup_series_v2.py SHA4042a972e8f8060ac8d1283234c7c5a7a2a73b976f753cd412de50427d596911只调整首末标签对齐，输入哈希及所有数值与初版一致。40/50图和旧程序f88原样保留。9文件（当前5+初版4）双端证明monitoring/H020_STARTUP_STEP60_ARTIFACTS_VERIFIED.json SHAac58814618cd631e1aa2fe87defb1c895c03e94b94fec9683c1afa229b8c082b。后续绘图优先v2。

最近memory.max累计748084，较747602增加482次，发生于58步重建期间；pressure累计增加约6.5ms，OOM/oom_kill仍0。峰值不是局部物理精度证据，CPU/GPU资源策略和源文件保持原样。当前用户没有再次要求暂停；开发ACTIVE。备用预算方案及三个程序仍未正向执行，使用01:10版本driver6d66846b...和followup dd513c...，创建具体闭合审查必须先见原producer EXITED、原执行budget FAIL/exit0/health[]及完整final checkpoint。建议闭合审查生成后先回收原rootcase/审查并逐项校验BUDGET_STOP_PRESERVED文件哈希，再正向preflight/调度新job。不要直接运行只接受CFD PASS的collectors来回收这段预算FAIL；本地原始元数据需要保留，然后再启动恢复。


## 2026-10-01 02:32 UTC — 原h020预算段已关闭，独立20→30连续性验证正在运行

原coplanar_h020_cfd_followup已EXITED。原生最终68步、exit0、health[]、checkpoint PASS；status FAIL仅因为wall time budget exceeded，原6小时停止请求02:19:02.040367 UTC，实际原生wall_time21918.14679002762秒（6h05m18s），等待当前PC/步完成的额外318秒如实计费。原execution SHA8f9d30ca32a002e0d00b217e19e5e9e026228e0ea7e6d6ce6ba8537af0f8e970；68步cp SHAf5d67eaea2ad635088ff18ca615b21e1bf06c9b087d055c8aaacc968879a5bfa。原最后完整稳态保存60步为3/5；未资格化、未回放。原20/30/40/50/60和68数据不覆盖。

已用solverpy运行review_h020_closed_budget.py生成H020_WALL_BUDGET_REVIEW.json SHA362c0507e7fe7a8ac58b4b84653d91000e64ad8555dbe99b878912fd0f6daa67。原47项关闭文件已回收并逐项SHA通过：H020_INITIAL_BUDGET_LOCAL_COLLECTION_VERIFIED.json；原元数据单独归档monitoring/H020_INITIAL_BUDGET_CLOSED_METADATA，包含原身份/policy/input_hashes/execution/progress/保存清单和恢复review。保存清单BUDGET_STOP_PRESERVED SHA ddc37d844e6ee89c315d8f75c988dd9dff751dac85e0ffc15e2891c935ffcb12。原生/源文件/网格均未修改。

正向preflight PASS：H020_WALL_BUDGET_POSITIVE_PREFLIGHT.json SHA7746651eaca38fb4961217b8a855f3ca9ac143dd951d4c3b6e60bf761163a90d。新额度累计43200秒，追加21600，扣原实际21918.14679后在验证前剩21281.85320997238秒。所有验证/续算实际耗时继续累计，原FAIL不改。中文review H020_WALL_BUDGET_REVIEW_ZH.md SHA220676f3f6e837583e9995442d7fc9f5139f27ebce5bfc733f66a37808e92ece。

已通过原scripts/add_job.py调度唯一新任务coplanar_h020_budget_followup_0001，producer PID776396；H020_WALL_BUDGET_FOLLOWUP_JOB.json SHA3fcb2ef12f7ccd6fd74ed921659ce595e97561d8becfd3e05308d656a6df48a2；启动证明H020_WALL_BUDGET_LAUNCH.json SHAdfc37c22009e6102f784f01813c2ec0b8daf699053e1d2a4c25ce4dec757303c。CURRENT=VALIDATING_H020_BUDGET_RESTART，RESUME_INDEX指segments/wall_budget_check_20_to_30_0001，MPI PID776416。先20→30原生重放，上限3600秒，u/p与每rank Y/A原1e-8门槛。仅PASS才由followup自动从原30步继续；失败直接停止，不能自动重试/放宽/转用25µm证明。30→最终续算、审计、46回放均尚未运行。初步资源约7.68核，验证段第一PC构建中，health全0。

代码继续使用driver6d66846b4504e10ab5eccb3541d4ee39125412caf3f85e2f9a08f287a62c7c42、followup dd513c9458aa790d10491ff27c02ae85139e0b199d8498fb42f3b18d4ed1cfc8；活动85脚本保持原快照ca522...，直到新producer也完全结束前仍不可改动。新分段所有文件独立保存，根case仍原预算FAIL。检查段即便原生raw execution FAIL（30步到顶而非稳态停止）也不能单凭该字段认定连续性失败，要读独立RESTART_CONTINUITY_BUDGET_20_TO_30.json；validate会严格校验最终步、全部实际步、health及原1e-8误差。

资源工具/tmp/brava_h020_memory.py已经按RESUME_INDEX选段，当前SHA93ac06843284a86fb1b6bfce4648091ea5808e1a796595edc86120a8ae7d0cc4，原版本和一次首键解析失败版均归档monitoring/local_helpers。一次远端查询成功但本地原解析器假定utc首键导致StopIteration，已改为与JSON键顺序无关并验证9个新原生进程；无CFD中断。记录monitoring/MEMORY_MONITOR_SEGMENT_UPDATE.json SHAfe995fbeb94be1eebc64b73ddf5d6327834ba593b7bde9bcca4dc1d5669fca03。普通monitor_live已经选新段，初次低CPU均值包含前几分钟回收/预检空档，之后7.68正常。不要继续用旧根progress判断新段耗时/进程。

后续：监控新验证与自动续算；验证关闭后独立回收校验，别重跑root已归档10..60 collector来收新段。任何新用户暂停优先按新段实际状态处理。待最终producer闭合后再运行wall query、部署报告整合并做25/20空间审查。pending报告report_coplanar_h020_validation.py最后关于“没有从检查点重启”的句子原指45秒短暂停顿；正式报告必须缩窄该句并补充此后的预算恢复执行链、原FAIL、连续性和累计耗时，避免误称整个h020从未重启。collector使用any_segment版本3aeefd...，原root-only已不适用。


## 2026-10-01 02:39 UTC — 最终h020报告草稿补充预算重启执行链

pending_development/report_coplanar_h020_validation.py当前SHA710cc43f541badecab70616a980a6130ad3f75b75c48903ae1243da8b1f88e17，旧545459...保存在previous_drafts同哈希名源码中；H020_REPORT_INTEGRATION_DRAFT.json已更新。新增原预算FAIL保留、同网格20→30原1e-8连续性、最终合格段与canonical identity、所有执行SHA/耗时账目的交叉核验；缩窄45秒短停未重启的措辞为“仅就那次短停而言”，另列后续预算检查点重启表。尚未部署/完整运行，85活动源码未改。真实远端VALIDATING_RESTART状态下，单独执行报告身份保护函数会拒绝宣称合格且不写报告；证明H020_REPORT_LINEAGE_DRAFT_GUARD.json SHA8e65afb909fc8bbad8260f2e22c310b9ccc061e4861e552e075212bcaa18aa2c。正向真实闭合场报告验证仍待最终CFD/回放结束。02:38:28验证段elapsed586.635秒，仍第21步FIRST_SOLVE PC构建，health全0、CPU7.6806核、RAM35.93GiB。新段自己的稳态计数0不覆盖原段3/5，也不等于连续性失败。


## 2026-10-01 02:58 UTC — 连续性验证第26步，补充闭合原段资源证据

开发仍ACTIVE；23:18短停/进度汇总/自动恢复已履行，无新暂停待处理。唯一活动job仍coplanar_h020_budget_followup_0001，MPI776416，20→30验证段。02:57:38累计1737秒，25步7次非线性迭代完成，26步按MAX_AGE重建第二个PC，尚无新残差；CPU7.682核，health全0，OOM/oom_kill=0，memory.max累计749293未新增。36e2秒验证上限、原1e-8连续性门槛、累计12h账目均未改；原FAIL保留。原运行85源码和活动pending恢复三程序均不修改。

新增只读资源报告01_mesh_capability/COPLANAR_H020_CLOSED_INITIAL_RESOURCES/REVIEW_ZH.md、REVIEW.json、resource_planning.png/.pdf已实际生成、看图并本地/远端SHA核验；证明monitoring/H020_CLOSED_INITIAL_RESOURCE_ARTIFACTS_VERIFIED.json SHA64ea9b5f282998e4637c2b578ee14ef711f5b37b7bb064efd9f7584f64f1f9c0（5文件含源码、5文档链接）。来源monitoring/review_h020_closed_initial_resources.py dbd30c1bb49679b3620407bc7fcf22cc29aa4f0cbfdd47d9882c6ab77dd080c0。20µm原闭合段进程树RSS采样峰值41682.796875 MiB=40.7058553467 GiB（包含runner和后代，5秒采样，共享页可能重复计）；之前37–38GiB是常态而非峰值。25µm全部闭合段同口径最大23.3538131714GiB。节点增长1.383倍、RSS增长1.743倍；原20µm含10%余量审查估计35.53008279GiB，被实测超出14.5673%。独立首小时cgroup采样峰42.055625916/42.055664063GiB，OOM增量0；cgroup与RSS口径不可等同。

18.75µm新两档PLC面数/体网格割线外推326135节点、1542728tetra，两种RSS情景44.8255/47.0463GiB；12.5µm52.1783/58.3627GiB。这些不是实测、上界或所有P1细化不可行证明。未生成新网格/启动新CFD/批准P2或dt；25/20空间判定仍待闭合场+同46回放。若仍失败，下一档需要以当前资源证据重新审查，不能沿用旧线性预测直接批准。

新增pending_development/collect_closed_h020_budget_continuity.py SHAe88755f9fdda5f3923c556187c175298c5faa1628b7755059d4952736ea693d5，外层及远端probe AST PASS、已同步、尚未执行。仅在独立RESTART_CONTINUITY_BUDGET_20_TO_30.json和检查段execution/ledger均闭合后运行；支持真实PASS/FAIL结论且保留raw execution状态，核对原47文件、20/30检查点、完整步骤/8rank误差、原1e-8门槛与账目，创建monitoring/H020_BUDGET_CONTINUITY_CLOSED独立快照，归档整个检查段，不碰后续续算。只补对应segment的本地../../SV_MESH链接。当前不要执行缺闭合数据的collector。

新增pending_development/report_h020_budget_continuity.py SHAcebb86396519c262aa4aec91b9ba1fbd774949b1258e4e3ca39c27efd669a498，AST PASS、已同步、尚未执行。等collector成功后在本地scientific venv用--study运行；生成该独立目录REPORT_ZH.md/REPORT.json/continuity_errors.png/.pdf，不覆盖原始证据，分别列u/p和8rank Y/A相对L2、raw执行状态、独立连续性结论，不称稳态/空间精度通过。生成后需看图、文档链接/来源哈希和双端同步核验。原主报告整合草稿仍710cc43...，等新producer退出再部署。


## 2026-10-01 03:24 UTC — h020原生连续性PASS，已自动进入原30步续算

20→30检查段已关闭，原生exit0/health[]，wall_time3125.663186311722秒（52分05.663秒），在3600秒限额内，checkpoint30完整。独立连续性PASS：u相对L2=3.7038325221178695e-15，p=2.9178908968086167e-15，8rank A最大9.96650696968471e-10，所有u/p/Y/A均<=原1e-8门槛。原生检查段raw status仍FAIL（指定30步结束未稳态资格化），没有改写。proof reports/RESTART_CONTINUITY_BUDGET_20_TO_30.json SHAadb7a7a89e79b57000fd750259a230830ea8c2301d9472ecb0a84b0d656d761a；检查段execution SHAe6e787d3ba1c3d626b7936143641593aab07ec3e7e0cc781417463d03e793afa。原30cp cced807e...，重算30cp34960bc48dd08aa9105cb51659810a63009e303ae986d07f9654fc5c62b7a1dd；检查点非位级相同，但通过原状态误差门槛。

followup已自动从原30步进入segments/wall_budget_resume_from30_0001。producer仍776396、job coplanar_h020_budget_followup_0001；新的MPI PID778658，验证MPI776416已退出。CURRENT RUNNING_H020_REVIEWED_BUDGET_CONTINUATION，RESUME_INDEX status RUNNING，checkpoint30/validation_only=false。03:24:15新段elapsed197.32秒，31步FIRST_SOLVE PC构建，CPU7.680核/RAM30.95GiB，health全0。仍无合格h020场、无h020回放，稳态必须在新段再满足5区间。

原段+验证已实际累计25043.80997633934秒，43200秒总额度在续算启动前剩18156.19002366066秒（约5.0434h）。账目reports/H020_WALL_BUDGET_LEDGER.json当前SHA439aed9e4373821113ce6e70a20a2227d64ec49cc99ea3ec0cabf8a02299408e，后续续算关闭才追加其实际耗时。不得把这个额度当剩余耗时保证，也不得预算重置。

已执行pending_development/collect_closed_h020_budget_continuity.py，30个检查段/原始引用/快照文件本地SHA核验PASS：monitoring/H020_BUDGET_CONTINUITY_CLOSED/MANIFEST.json SHA9488ca510b4cfbc5bbc10b2c47ca6e8f1dbfa197dbd0d2ebdb33a0a4aec40ea9，LOCAL_COLLECTION_VERIFIED SHA 4fb48cdd28b92f0e15b004bae5514a668e9496c6571beb1da2375c52aafb642c。完整验证case已在本地segments/wall_budget_check_20_to_30_0001，SV_MESH指本地../../SV_MESH；不要重跑此collector覆盖独立快照。root原预算FAIL不覆盖。

已执行report_h020_budget_continuity.py生成上述目录REPORT_ZH.md/REPORT.json/continuity_errors.png/.pdf，18个相对L2数据点及门槛使用实值，无平滑。图已目视检查，9文件含source/原manifest/快照/报告双端SHA PASS：monitoring/H020_BUDGET_CONTINUITY_ARTIFACTS_VERIFIED.json SHAb06bd8805f3f55f969e3679f837fe7378575663007144c4afb9cd862476ca36a。图只证明重启连续性，不证明稳态或局部精度。

新增pending_development/archive_h020_continuation_saved_step.py SHAb52724659ed7b91bf1ec8051272dd06d2c18fc3514790ecdb2cb01169ea0f1f7，外层/嵌入probe AST PASS，尚未执行。用于续算40/50/60/70/80完整保存并进入下一步后归档，固定选择新continuation目录且要求index/proof/cp身份；校验checkpoint、复制原生VTU/bin，单独保存段内SAVED_STATE_XXX证明和本地回收证明。不得使用旧root collector归档这些续算状态。完整CFD闭合后仍用any_segment collector 3aeefd...，回放后用原replaycollector。

最终报告草稿report_coplanar_h020_validation.py已加入实测资源报告入口和输入哈希检查，当前SHA9e3580cf3252afc1ee9e50f13f1913577a55bc40179e166e27580fe6e3444279；旧710cc43...已保存previous_drafts。H020_REPORT_INTEGRATION_DRAFT.json已更新，AST PASS、尚未部署/正式执行，旧负向lineage guard记录对应710cc43版（核心身份保护未改变）。活动85源和恢复driver/followup全部不动，直到新producer完全结束。原用户短停/总结/自动恢复已完成，开发ACTIVE。


## 2026-10-01 04:04 UTC — h020续算40步归档，稳态区间1/5

当前仍coplanar_h020_budget_followup_0001 / MPI778658，在segments/wall_budget_resume_from30_0001运行；04:03:34已进入41步MAX_AGE第3次PC构建，elapsed2554.94秒、health全0。续算30→40 E_u=1.444051238026335e-7、E_Q=3.6555320171826055e-9，质量误差9.09501627114305e-13，壁面u=0、u_max=.5499245853657136。新段连续区间1/5，不能与原预算段3/5相加，不是合格稳态场。原生31..35均各5NL、36..40各4NL，与原闭合记录一致；之前总结中35步4NL应以原execution实际5NL为准。

已实际运行pending_development/archive_h020_continuation_saved_step.py 40：段内reports/SAVED_STATE_040_VERIFIED.json SHA1863e9f8934c5945b6f47d6385b2cad19a95bea8f0284be4d7700681e1b48ca7；40cp SHA56c9e018c191e74387914aca0c881153cc5c6a5621a1583f8bd3684df385418d；reports/SAVED_STATE_040_LOCAL_COLLECTION.json SHA450d03f37d8ceaabcf5d4679e9bb36c1e45fafb68bdc0e9c0b4609e724b07a2e。两个原生文件本地SHA已通过，server同步本地核验记录。不重跑collector40。下一保存状态50，使用新continuation脚本而非旧root工具。新段本地目录目前仅40步资料，尚未完整同步活动case，也未创建SV_MESH链接；最终any_segment collector会收齐并处理已知链接。

最近内存约38.2GiB，memory.max累计749293未新增、OOM/oom_kill=0，CPU约7.68核。验证报告及资源报告均已完成并双端核验，最终报告草稿仍9e3580cf...，活动85代码未变。当前用户短停/汇总/自动恢复已完成，开发ACTIVE；下一步继续稳态区间、合格后审计+46回放，再闭合后部署报告/壁面查询和25/20空间审查。


## 2026-10-01 04:40 UTC — h020续算50步归档，稳态2/5，更新进度/ETA报告

同一续算segment和MPI778658继续运行。50步已完整保存并进入51步第5个PC构建；04:40:02新段elapsed4745.62秒，CPU7.679核/RAM34.1GiB，health全0。40→50 E_u=3.2568607014569315e-9、E_Q=1.3692861104257136e-11；mass5.293955920339377e-16，wallu0，u_max=.5499245853645819；新段连续稳态2/5，仍未正式资格化或回放。

已运行archive_h020_continuation_saved_step.py 50，原生文件和证明本地/服务器已回收核验。段内reports/SAVED_STATE_050_VERIFIED.json SHA25253c6dfd4e4c6394c37790c370bf3d7de4b52b444a9e99c4298d59123bc5c4；50cp SHA96e75eecf6e4be749bfe098c4165427210aa2e124e5e9f21cbebeec29144745c；SAVED_STATE_050_LOCAL_COLLECTION.json SHA10f4362dd155cb1cb41ea8c332a430cf78737d03709c49104187d41689fdcc8c。新续算40/50均已收集，不重复collector；下一次60。

新增monitoring/H020_CONTINUATION_STEP50_PROGRESS/PROGRESS_ZH.md、STATUS.json、monitor_snapshot.json，04:38:06快照。实查10套合格新CFD/11批分析完成回放，20µm未计完成。该快照沿用先前原30→81条件模型10677秒，减当前实际续算4629.166秒，原生估计剩6047.834秒约1.68h；含后处理、场审计、46回放、比较/回收/报告预计本轮20µm还需2–3h（从04:38算，不每次重置）。前提60/70/80继续通过并正常停止/审计/回放顺利，非收敛或时限保证。全项目若需P2/P1代表区域验证仍预留1–3工作日或更长，未自动启动P2/dt。原生histor.dat的reported_elapsed可能包含检查点继承的旧计时；预算/ETA使用runner真实wall elapsed与closed ledger，不用native打印绝对elapsed冒充本段耗时。

报告源monitoring/report_h020_continuation_step50_progress.py 1b9726f8358280cc1634f79a0b53b3cbc7b5eac7a029b858a73c11ba749ba1d0；已执行、6链接核验、4文件含源双端SHA PASS，证明monitoring/H020_CONTINUATION_STEP50_PROGRESS_VERIFIED.json SHA09e3c90bd59862fa45317a5509a921d7a80a6f7a7a584b1f35d8d1fc73b7b5d0。不覆盖23:18原用户暂停汇总快照。当前开发ACTIVE；85活动源及恢复三程序原样，报告整合草稿9e3580cf...仍未部署。最近memory.max仍749293、OOM0，无新增事件。


## 2026-10-01 05:08 UTC — h020续算60步归档，稳态3/5

同一segments/wall_budget_resume_from30_0001、MPI778658持续运行；05:07:24进入61步第7次MAX_AGE PC构建，elapsed6388.53秒、health全0。60步已完整保存并独立回收SHA核验，50→60 E_u=6.92073671421056e-11、E_Q=5.43953970814871e-14，mass1.764651973446459e-16，wallu0，u_max=.549924585364235。新段连续稳态3/5，不能合并原预算段的3/5，仍未资格化/回放。

archive_h020_continuation_saved_step.py 60已实际完成，不重复运行：段内SAVED_STATE_060_VERIFIED.json SHA0f32984c64fe337f445462df823d4f10a47b356fd63fef51b5a2400db0929a29；原生60cp SHA44b3e13841779c6dfa54ae46646cb63e56ac58e71da423e8d2ca0447e3325c68；SAVED_STATE_060_LOCAL_COLLECTION.json SHA07c3c966daa9e93dc9e662696b775b592b3f619be159c38dd4bd2c9f1920c386。新续算40/50/60均已收集，下一70；旧root10..60文件独立保留。

最近普通迭代内存约38.5GiB，CPU约7.68核，memory.max仍749293，OOM/oom_kill为0；无新增数值失败。04:38进度/ETA快照保持历史，不每次重置其2–3h估计。全部正式微泡1565及ACTIVE_FLOW未更改，h020仍是独立实验，主报告整合待producer完全结束。开发ACTIVE，无新暂停请求。


## 2026-10-01 05:39 UTC — h020续算70步归档，稳态4/5

同一续算segment、MPI778658持续运行；05:38:51进入71步第9次MAX_AGE PC构建，elapsed8273.42秒、health全0。新段60→70 E_u=1.3919122895068985e-12、E_Q=1.0587911840678754e-15，mass3.529303946892918e-16，wallu0，u_max=.5499245853642302。连续稳态4/5，仍须80步第五区间以及原生正常退出，当前不是正式合格场/没有h020回放。

archive_h020_continuation_saved_step.py 70已实际执行并本地SHA核验：段内SAVED_STATE_070_VERIFIED.json SHAc5daa2f577fdb6ccfd9f18667e69c88bfa30de551acdeb55b789b4576dbb31a9；70cp SHA9be7af45bdca09053506392e63ea92de15cfb5b713b65bf66d342a8fa419f2f7；SAVED_STATE_070_LOCAL_COLLECTION.json SHA9625a0fd9d40a2422c8378b9cfdee01304208c074a6990f25f0515db51bc604a。新续算40/50/60/70已收集，不重复。下一80保存且进入81后收集80；若80触发稳态STOP_SIM仍要等待原生当前步完成，不能先假定退出。

目前CPU约7.68核，普通阶段内存约38.7GiB，memory.max仍749293、OOM0，无新增数值失败。活动源85及恢复程序未改，报告整合草稿9e3580cf...仍待整个producer退出后部署。后续流水线将独立身份核验/场审计/同46条回放/分析/原报告；闭合后再执行实际路径wall query与新报告整合。用户短停/汇总/自动恢复已完成，开发ACTIVE，无待暂停指令。

### 2026-10-01 06:09 UTC：20 µm续算已完成五次稳态检查，等待原生正常退出

实际续算第80步已保存；40/50/60/70/80连续5/5通过，程序于06:08:57 UTC发出five consecutive steady intervals and physical gates正常停止请求。原生求解器已进入第81步预条件器重建，仍须完成该步才能正常退出；此时尚无最终合格场、场审计或20 µm轨迹结果，不能宣称整轮完成。唯一活动生产器仍为coplanar_h020_budget_followup_0001，未启动dt/P2或正式群体批次。

第80步原生检查点/VTU已独立核对并收回本地：E_u=3.0838435797285595e-14，E_Q=2.205814966808074e-16，质量误差3.529303946892918e-16，无滑移壁点最大速度0；O1/O2/O3分流0.3330815453214332/0.3333717340193113/0.3335467206592554。仅表示保存状态和全局稳态检查，不认证局部场精度。

位于03_cfd_runs/M1_coplanar_fine_h020/segments/wall_budget_resume_from30_0001/：检查点080 SHA74e177a3e2f464de57960bbb0626b5125bb1779d8874474948567934bffed770；reports/SAVED_STATE_080_VERIFIED.json SHA df9bdfbf620f874d592ca73e6f013dff56769685fe14d9861c3a2e34e8b1696d；reports/SAVED_STATE_080_LOCAL_COLLECTION.json SHA0e3baded88518cbc41ecdf3b592ba591b8991d86a55db7449ca4b39a98a2905a。收集脚本正常exit0，原数据不变。

用户要求的短暂停止/进展报告/自动恢复早已完成，当前继续开发。待原生PASS退出后，生产器自动执行后处理、canonical身份检查、场审计、冻结46初态回放及比较；仅当整个生产器EXITED后才能部署pending_development中的新报告版本或运行额外重任务。最新报告草稿9e3580cf...已同步但未部署，85项运行源码仍冻结。

### 2026-10-01 06:39 UTC：20 µmCFD与场审计已完成；场门槛FAIL，46条回放进行中

20 µm原生续算最终81步、exit0、health_failures=[]、status PASS，wall_time_s=10816.527707576752；execution SHA47619e36478eec152efe1273976d8f7435511855597f35455af302b02795f465。全部原生段累计35860.33768391609s，12h总预算剩7339.662316083908s；ledger SHAfa458c8b1ef54245e51f785f98fe63822cfe99c945440b955666f616a2702053。原6h失败段仍保留，未改写。最终续算段仍segments/wall_budget_resume_from30_0001，RESUME_INDEX已QUALIFIED_STUDY_FLOW（SHA26f7d7deadfe9af26f6576940554e57c542fabd5f0da54e12014ba39f9891566）。

最终flow SHA cb36025fece24f087f550f4e8753c8c868a1894f67ec444fb498af6b42220925；native SHA be4b2bc24f13140d91863b8fd9d30eaceea45c7d0a79cf9313a69ceedb0870c2；mesh SHA e6fc8966bf25d8d2b2190265ed98dd53cef90c055115bad24b8e317a564f56e3。STUDY_FLOW报告SHA15fad6ae17622ed54bb1083ea7c7f42f25456f345ea3b8738a2f12de6e345ce4，独立canonical身份PASS，证明SHA9eb672b99f6b360918628a67c15aee72e9279e24c1bcc15d6c841f5995c99477。最终质量误差1.764651973446459e-16，最大壁点u0，WSS均值3.0751888565863297Pa、P955.717266737924176Pa、最大27.55816780237174Pa。正式ACTIVE_FLOW未替换。

场审计634.8188254833221s完成，1059固定点全部在域内，zero_tetra=[]、eta_zero_gradient_count0、centroid_zero_speed_count0。CPU/GPU G最大差1.5916157281026244e-12仅为计算一致性。全场max|div|4584.603041364238s^-1。SUMMARY SHA b16a7a22a4aae72961e0c4ae6c798c6218f7689d74c09c6ce73fecf40e5903d2。已用any_segment收集器完成141文件逐一SHA验证：CLOSED_H020_CFD_OUTPUT_MANIFEST.json SHA dcb44eb10747ea6fa93243c3ce77906ec2584c516f716fbecc92121103fc3cb9；CLOSED_H020_CFD_LOCAL_COLLECTION_VERIFIED.json SHA318e298012723275b5f9a8fb54e0e3e99e3be22fa09c98b42f6f8e28fcab152e。无需重跑该排他收集器。

新增独立场比较monitoring/H020_CLOSED_FIXED_POINT_COMPARISON/REVIEW_ZH.md、COMPARISON.json、field_convergence.png/pdf，已执行、查看并双端验证5文件/3链接。源码pending_development/report_h020_closed_fixed_points.py SHA7152095053674db90a1aedee37661fb42a027c4ccaca0cd767584f6a19c37b1b；比较JSON SHA5e117d6cd5b09857b762cd7f0478f8ace4bbb29a66c2357d8cdb2c6a510764b3；双端证明monitoring/H020_CLOSED_FIXED_POINT_ARTIFACTS_VERIFIED.json SHA39241997be3fea0ae7374254830d60708e5b22108c47a07c9671b9404e50713e。此为回放关闭前的场比较快照，不是完整轨迹判定。25/20 ALL：u1.680135%、G16.588714%、WSS22.112642%；六区域u均<2%，G/WSS均超过各5%门槛，场FAIL。high_eta区G25.796676%、WSS37.520462%；zero区G20.226279%、WSS15.439108%。不得选production或直接做正式dt/P2。

原388769中心20 µm：u0.078366828544123m/s，div-64.80922289030491s^-1，G1986.5017187410656s^-1，eta0.032624800813853626；原138584中心u0.02347642529326879m/s，div6.089754709468679s^-1，G914.156253017135s^-1，eta0.006661612486234861。单点改善不代表区域/轨迹收敛。

06:39:01监控：唯一生产器coplanar_h020_budget_followup_0001仍RUNNING，CURRENT_EXECUTION_STATE=REPLAYING_REVIEWED_COPLANAR_COHORT；20 µm冻结46初态回放elapsed376.93s、5526迭代、24/46终止，其余继续，未有整批终止汇总或分析，不得把24条计成出流。GPU瞬时36%，RAM7.56GiB。旧85项源码仍冻结ca522c...。

新增待部署辅助程序pending_development/deploy_closed_h020_reports.py SHA00fcda40b4516de44fc1ac535c1974a869aae188c6e449a1f057f6437f10516b，语法通过并已同步。只允许整个生产器EXITED、所有活跃作业为空、46分析完成、壁面查询完成后部署原三个报告草稿，并先备份旧报告/figures/源码快照到monitoring/PRE_H020_REPORT_INTEGRATION（尚不存在）。实际负向验证H020_REPORT_DEPLOYMENT_ACTIVE_GUARD.json：活动原生持锁时拒绝预检查，85源码/SNAPSHOT不变、无发布输出，不能把它说成部署成功。正向部署尚未执行。

接下来：等待46回放/分析/自动比较/旧报告生产器整体EXITED；收集完整replay；再独占运行原版本audit_coplanar_wall_queries.py（25/20、reference25、输出COPLANAR_FINE_H020_PATH_WALL_QUERIES），然后运行报告部署辅助程序及新版render_report，检查图片/链接/双端哈希。完成25/20完整评估和资源可行性复核后，按原Phase2G决定进一步P1细化/有条件代表域P2，不跳过前置门槛。当前未再暂停，用户要求的短暂停止/报告/自动恢复已经履行。

### 2026-10-01 07:03 UTC：20 µm验证与报告全部闭合；18.75 µm仅网格检查已启动

20 µm冻结46初态回放734.7184898853302s完成，44出流、2观察截止、0数值失败、0有支撑静止；O1/O2/O3=16/15/13，截止ID393、1377，冻结slow-fast筛查24事件（25 µm为20）。完整25/20比较FAIL：0/46全部指标通过；共同出流但出口变化12条[2,4,9,168,440,487,650,1041,1077,1189,1208,1505]；终止类别变化3条[5,534,1329]。速度史46条均未过5%门槛，端点与停留时间各40、最小间隙23、路径长度12未过。历史失败ID440现在O3、residence1.353973870579249s；1505现在O2、1.7897540451226592s，不能把出流当成收敛。旧1565身份不变，也没有新正式群体批次。

376回放文件已收回并逐SHA验证，不能再运行排他收集器：CLOSED_H020_REPLAY_OUTPUT_MANIFEST.json SHAa14672d44b5e96e068ba659ece21d09013f240d7e9d5355db3b7f021bb6e75a7；IDENTITY SHA6c379e3cbc4f98349f15802170e9d0083749eae5fd5ca48e85613a84c99e76a7；analysis/SUMMARY SHA0a17a96bf0e6035565a06fe1764923d244b97e6c0efca20e60ddb86c613e8f18。生产器coplanar_h020_budget_followup_0001于06:47已EXITED，不得重启。

追加GPU只读壁面查询已完成：01_mesh_capability/COPLANAR_FINE_H020_PATH_WALL_QUERIES/REVIEW.json SHA6e819fe82d92ef8f311afe32368e514b61f2608c3efcd0fb90b35742af4cd903。在25 µm全部407106记录位置（242336近壁查询）上，25/20最大间隙差8.51098705401121e-18m、法向向量差0、接触数量变化0、feature/rim变化0。耗时12.764s，无重新积分；只覆盖这些坐标，不证明完整接触求解或所有20 µm路径等价。

已部署三个报告源码，旧85快照与报告/figures保存在monitoring/PRE_H020_REPORT_INTEGRATION；现在CURRENT_STUDY_SCRIPT_SNAPSHOT为86项、SHAdb2266fe8a4033612fe15dd317e41783afe33a81b09857b91df681b1746bfae2。render_report SHA799f09b3f3cd7e175defe898ebe70555c8bd9da6818c5a65f9f2ef2ec0331a91（增加独立场图链接），report_coplanar_h020_validation SHA9e3580cf3252afc1ee9e50f13f1913577a55bc40179e166e27580fe6e3444279，report_fine_restart_recovery SHA44f922737f688296d6f2783ce663ee13ad99c1b5d7a6e1da08790c9e2425f3d6。部署助手已正向执行，勿重跑排他归档。

H020_INTEGRATED_REPORT_EXECUTION.json记录新版报告实际PASS、exit0、58.674601793289185s，ACTIVE_FLOW SHAde168191873967374caea00083ea6c34baec164dd873795a5fef1f6e3c19cb6d不变。报告时刻06:49:55.968724UTC，REPORT_MANIFEST SHA6c755ccf0f3ab033b0d845f02f3a2f3467bd92ff72a279b29da3ead8f92c20f2。pending_development/collect_integrated_h020_report.py SHA4302cacbe24b3d76f3bc9e6b9221cd750e43d91dca0d9a9212d474df0343b8db已实际完成：86源码、106图文件、HTML+h020MD共201链接通过。证明H020_INTEGRATED_REPORT_LOCAL_VERIFIED.json及更新的REPORT_LINKS_AND_HASHES_VERIFIED.json；旧链接证明另存monitoring/H020_REPORT_LINK_PROOF_PREVIOUS.json。三张新h020图已目视核对，显示原始接受状态、无平滑/事件对齐。OPEN_RESULTS.html和COPLANAR_H020_VALIDATION_REVIEW_ZH.md已发给用户。该静态报告的active=[]属于其生成时刻，当前新网格任务见下。

遵照Phase2G，在确认完整25/20 FAIL后继续下一档有界网格可行性检查，尚未开发P2或正式dt。新计划COPLANAR_H01875_MESH_PLAN_ZH.md SHA c307e529580b0294a436eda18e2fd9e831365964d8c7a80023279888272546d9；审查COPLANAR_H01875_MESH_REVIEW.json SHAe4055ebee1ca4a1d0e61b29b2425c4fec0e6f0feedfd07db3dc984af11e5792c；pending_development/h01875_mesh_only_followup.py SHA2710c49491d3ac179e60fc18a114bba1c44c4c6081b2f14be272892502638558。首次未执行草稿5f309...及review ab726...保留于pending_development/previous_drafts/h01875_before_state_guard，当前版本增加启动状态及步骤间PAUSE检查。

新job coplanar_h01875_mesh_only_0001于07:02:41UTC启动，唯一活动job。预检查证明SHA72e99ccaef6e3d92ca6b226e5b95b0bf2e0301b3ee7f1184b32d314d7deb5310、job SHA1b721e403dac7863f2763843d33c20496bca5b251b0e73ba59fe7936762b01e4。CURRENT_EXECUTION_STATE=RUNNING_H01875_MESH_ONLY_RESOURCE_REVIEW。目标probe M1_coplanar_h01875_probe、candidate M1_coplanar_fine_h01875，comparison_medium=h020。表面上限300000→330000、tet上限1900000、native网格900s、prelaunch RAM<30GiB，科学/物理门槛不变。预计mesh+QA约2–20min，非保证时限。CFD_automatic_launch/particle_automatic_launch均false。当前86项源码与新pending驱动冻结，任务结束前勿修改。

07:03观测：输入表面实际154739点/309474面，与已测精确面数一致；原顶点/端口三角形不变，原面面积误差6.83e-15、平面误差2.71e-20m。原生网格已生成337592点/1613894 tetra，原始QA/随后共形拓扑处理仍未闭合，不能算合格候选。旧资源外推18.75为44.83/47.05GiB，超过37GiB CFD规划与42.06GiB容器额度；即使GATE_MESH PASS也不授权CFD，尤其不能用旧线性RSS QA结果绕过最新资源审查。下一步待实际网格QA后重新评估内存，可独立研究更多MPI分区，但目前没有改变MPI/求解器配置或声称内存改善。

原生base runner已核对可从policy.MPI_ranks读取进程数，默认启动command为wrapper+launcher -n ranks+wrapper+exe；而study prepare_candidate_cfd.py的身份metadata固定CPU_MPI_ranks=8，coplanar_output_runner.py候选白名单也仅含50/25/20。未来更多MPI/新候选需要独立准备与身份检查，不能只改数字后直接投入长CFD。目前未启动此类实验。用户没有新的暂停请求，最早短暂停止/报告/自动恢复已完成，开发持续。


### 2026-10-01 07:53 UTC：18.75 µm 网格 QA 闭合；16 MPI 资源短试验已启动，开发 ACTIVE

18.75 µm mesh-only job coplanar_h01875_mesh_only_0001 已于07:40:48.961752UTC结束，H01875_MESH_ONLY_CLOSED=CLOSED_MESH_QA_PASS，实际2288.6561284065247s（38.1min），原2–20min规划偏短。原生337592节点/1613894tet有453 all-WALL且float32墙面精度FAIL，原记录保留；恢复声明PLC输入精度并完成405内部边扇共形细分后，337997节点/1615458tet、all-WALL0、负/零体积及非流形0、1063固定几何探针全部在域内、体积差0、五角色细分PLC三角形坐标集合精确一致、四端口原三角形不变。严格域QA PASS。minSICN最低.020611120735138218，中位.8134400522315466，低质量607。不能将mesh PASS当CFD/空间收敛通过。

72文件（含原始失败、共形候选、输入PLC/ancestry、MPI启动能力记录）已本地/服务器逐SHA回收。collector pending_development/collect_closed_h01875_mesh.py SHA368a21231bfe76ec9e7c29fe3b18113a8792d1ade126c4556203ca164c5f7ebb；CLOSED_H01875_MESH_OUTPUT_MANIFEST SHAa4ac9a4573ebcfb164607299ba254e9ce6fb9ce2e487c1b9ee7d74533b7a5804；CLOSED_H01875_MESH_LOCAL_COLLECTION_VERIFIED SHA406f7b7d1866297b7c75b9e50bf062da31c1341593fb9f9341269879b7db8a63；GATE_MESH SHAd8f1e40cae163f7520daeae77eac81f1ffbf8d07173614a3ca125c7b69a59ac9。收集器已完成，勿重跑排他proof。

实际节点数下8rank外推RSS46.45584956775984/49.55544275313739GiB，高于37GiB规划及42.0556640625GiB容器限额；是规划情景，不是实测或上界。新增01_mesh_capability/COPLANAR_H01875_MESH_RESOURCE_REVIEW/REVIEW_ZH.md、REVIEW.json、mesh_and_resources.png/pdf，图已目视，8文件9链接双端核验，证明H01875_MESH_RESOURCE_REPORT_VERIFIED SHAe1bcf279f80650c2face1cc6acdf154b5c7727553e337a89635b6f469c393c46。报告JSON SHAc9537f94e76b2d6c2c4b80eb53ed0cc65f5664b2dbcb692b1103667119940b70。

原生启动器8rank默认PASS，16默认slot不足FAIL已保留；显式--use-hwthread-cpus可启动16个进程（每个仍affinity0..15，CPU总额度7.68核），只是启动能力，原记录位于01_mesh_capability/MPI_LAUNCH_CAPABILITY_8_16。没有修改全局MPI配置。

为实测更多分区资源和数值一致性，新独立case位于01_mesh_capability/MPI_PARTITION_RESOURCE_PROBE_H020_16/case。使用现有20µm mesh（mesh_arrays f93bff504ef4bd035b58157f9e39e5d73205c082a0cf9f3c934670fcda9c3f2d），零初态到第10步；MPI8→16，启动命令显式硬件线程槽位，总步500→10，其余XML/物性/BC/dt.01/PETSc及非线性容差不变。预算9000s、正常停止宽限900s、采样cgroup39.5GiB软停/41GiB硬停，仅终止本试验进程组。原先未执行的3600s草稿及初次guard保存在pending_development/previous_drafts/mpi16_before_startup_cost_review；修改依据原8rank第10步reported elapsed6660s/137次线性求解，不预设16rank加速。

审查H020_MPI16_RESOURCE_PROBE_REVIEW SHA0c21fb9ec6357cf851fca9ccc2a691a7b8e09169668e5710a87c6e02096acb45；计划H01875_RESOURCE_AND_MPI16_PLAN_ZH.md。驱动pending_development/probe_h020_mpi16_resources.py SHA422e3a9047aff913272065ac9ac75139d2290f20ca33962939612c84bc7f8efb；语法/5处adapter完全逆向还原base源检查通过；实际并发锁负向guard V2 PASS（此时mesh已结束，由test_parent持同一锁，不称仍由mesh持锁）。H020_MPI16_RESOURCE_PREFLIGHT_VERIFIED SHA f81ae286737abaa0ed1eaba81dd6388d31614a77aec08489c5f19780c1ba5fbc 已正向PASS。执行adapter SHA62c844eecf9a5d2b0eb7ebff29c6e28f7563d69073d011255d13c7f0a1678053，11输入/源码本地proof SHA760bd892744e539617239b2789679ed2c15e30733ee345fc3cef59b825028267（输入证明不等于CFD完成）。

job h020_mpi16_resource_probe_0001 于07:48:47.982028UTC启动，PID788441，原生MPI进程组788449。仅此昂贵job活动。CURRENT_EXECUTION_STATE清理旧暂停/续算字段，状态RUNNING_H020_MPI16_RESOURCE_PROBE。运行期间冻结86活动脚本快照db2266fe...及新pending驱动422e3a.../review0c21...，不得改活动源/配置。监控使用monitoring/monitor_mpi16_resource_probe.py，写MPI16_LATEST.json与MPI16_HISTORY.jsonl；旧/tmp/brava_monitor_live.py虽active正确，其CFD摘要仍硬编码h020已完成PASS，不能拿它判断新短试验完成。

原始固定10步execution预计仍FAIL（不满足稳态停止），必须保留。独立RESULT只有在原生正常step10/残差health/无停止或OOM、同canonical u/p/G相对L2与Linf均≤1e-8、采样cgroup<39.5GiB时才PASS_RESOURCE_AND_PARTITION_PROBE。G L2体积权重，u/p同节点原值，不平滑/插值。8rank参考原始step10 VTK SHA1f1acf461a968c7011e500b90b478c1def255c9e4fdc5582d19a347b411c731e；cp10 e1e56de0118089ab487ab5dfdb70edde22f9832ce18f78c404e68c5046c2edf8。试验PASS也不自动授权18.75µm全CFD；仍须资源审查，不能声称稳态/物理/空间精度通过。

07:50首次监控elapsed87.6s，首次PC分解，cgroup32.36GiB/treeRSS26.26GiB，GPU6420MiB，OOM0、无STOP、无保存步；仅早期样本，不能当峰值。预计本短试验加核验2–3h（从07:48启动计，不每次重置）；若后续可行，18.75µm完整CFD/场审计/46回放另6–12h。整个任务仍依赖空间gate，不作全部完成保证；正式dt/P2/P1及生产群体尚未开始。此前用户短暂停止/进度说明/自动恢复已履行，没有待处理暂停请求。


### 2026-10-01 08:12 UTC：用户新指令暂停，等待明确恢复；取消先前自动恢复语义

最新用户要求“下一个可停止点暂时停止，输出已完成/结果/后续/预计时间；完成后我会给指令是否恢复开发”。已执行，**现在 PAUSED，禁止自动恢复或继续开发计算**。上文旧的 ACTIVE/自动恢复描述已被本条覆盖。

16 MPI 首次资源试验于08:03:30UTC结束，881.8395447731018s，cgroup采样41.0255012512207GiB触发预设41GiB硬停止；原生exit1、final_step=null、没有完成第1步/没有可比较输出。无OOM/max事件增量；treeRSS采样35.01164245605469GiB，仅是中止前样本，不是完整峰值上界。RESULT=FAIL_RESOURCE_OR_PARTITION_PROBE；equivalence_metrics=null，不能声称16MPI不一致或一致。execution SHA09e95578b643b20204b58193634c0f444645cc9f8fe8d37c92fc449d3c6631bf，RESULT SHA f4fc74fda621f86b9059abc1b9757cc7a125179e4d8adf108bbf6f559710b4b1。29文件已双端SHA回收：CLOSED_MPI16_RESOURCE_PROBE_MANIFEST SHA3d6bee7b8c76745ed285d32d97a061c9c1d95f60caf8a856ebe49d2353ba4343；CLOSED_MPI16_RESOURCE_PROBE_LOCAL_VERIFIED SHAc7235ef7c661e3553c0f03d20ae215d1f124b056e17351f2c0bcfa31c8aa5ff4。收藏脚本pending_development/collect_closed_mpi16_resource_probe.py SHA3fca03580624dbd98a1d7e89fc6b76a32bd23723915530b3b99697d7610ee0a8已完成，勿重跑排他proof。

期间完成只读缓存初查：pending_development/audit_closed_evidence_page_cache.py，输出01_mesh_capability/MPI16_PAGE_CACHE_AUDIT/RESIDENCY.json。583个已归档证据文件有1.601959228515625GiB驻留页；mincore不读取大文件内容、不请求回收缓存。停前MEMORY_STAT_AFTER_SOFT_STOP记录anon约30.50GiB/file约8.30GiB；这些类别可能交叠，不能任意求和或把1.60GiB直接从cgroup中减掉。文件缓存调查尚不完整，**没有执行posix_fadvise/清缓存，没有资源重试或18.75µm CFD**。只读查询MPI安装得默认mpi_yield_when_idle=false，记录pending_development/MPI_IDLE_PROGRESS_CAPABILITY.json（服务器）；未改变当前MPI配置，未验证该调度参数的数值或性能影响。

新增固定点分组补充：monitoring/COPLANAR_FIXED_POINT_STRATIFICATION/REVIEW_ZH.md、REVIEW.json、fixed_point_strata.png/pdf；源pending_development/report_fixed_sample_strata.py。5文件5链接双端哈希核验、图目视通过，证明monitoring/COPLANAR_FIXED_POINT_STRATIFICATION_VERIFIED SHA faa18074d06d02cf1ed717f161e8a9f5b521901ad54a971d94fed1bd0d4da82f；JSON SHAd0b8961512e254355124f30b70926938e5daadd3d6b8d9bf74e5ad882168681d。25/20最靠壁两个130点组u相对L2均14.223548%（ALL1.680135%）；20场中位壁距.135981/1.359812µm，中位u.100970/1.009701mm/s，最大速度向量差.113110/1.131103mm/s。50/25这两组u16.333698%。各网格130对最近射线点全部同单元；ray fraction是原壁面面心→对顶点比例，不一定法向。球域并集内345点u/G/WSS1.333372/21.975626/27.715555%；并集外714点2.259107/10.974946/6.395081%。这是事后分层诊断，不改变冻结gate、不把采样误差份额解释为因果或全域概率。

实际暂停08:12:32.539850UTC，CURRENT_EXECUTION_STATE=PAUSED_BY_USER_AFTER_MPI16_RESOURCE_PROBE、active_jobs=[]、resume_requires_explicit_user_instruction=true。全部study supervisor无活动任务，重计算锁可取得，旧MPI788449不再存在。暂停快照pause_snapshots/20261001T081232Z/，PAUSE_VERIFIED SHA7666cec50f106c017ad41012169be365158fe818f85d63cb09527b35557dc61a；保存PRE_PAUSE_EXECUTION_STATE、PROGRESS_ZH.md、LOCAL_PAUSE_REPORT_VERIFIED，并追加最终已知原生CFD可执行文件/proc/exe无残留进程和报告远端SHA核验。实例/服务保留。正式ACTIVE_FLOW SHAde168191...、86源码快照db2266fe...、旧1565及动画均未更改。

完整已完成/结果/后续/ETA见pause_snapshots/20261001T081232Z/PROGRESS_ZH.md，包含现有综合图文链接。本轮用户要求已履行后应结束，**等待下一条明确恢复指令**。恢复后先读这份暂停报告及16MPI失败/缓存证据，再审查资源方案；不能直接重启旧失败job或自动启动18.75 CFD。条件规划：资源复核/有界验证2–4h；若可行下档CFD+场审计+46回放6–12h；空间通过后dt约1–2h。顺利情况下剩余约10–18h；若需局部P2/P1另1–3工作日或更久，不保证收敛或绝对完成时间。


<!-- LOCAL_P2_CURRENT_STATUS_START -->
## 2026-10-01 局部 P2 / 近壁审查已按 A 区门槛结束（2026-10-01T13:35:42.769575+00:00）

用户的新指令已替代旧全域 P1 路线：不恢复 18.75 µm CFD、17.5/15 µm 网格、旧微泡批次或旧资源探测。新目录 `computation_examples/brava_flow_roi_18mlmin/p2p1_nearwall_cross_validation/`；服务器对应 `/workspace/brava_flow_roi_18mlmin_20260928/p2p1_nearwall_cross_validation/`。新边界使用原“中心线阻力网络 → 出口压力”，未使用校准 CFD；局部入口按网络流量、下游按沿程网络压力。

结论 **LOCAL_P2_NOT_YET_SPATIALLY_RESOLVED**。A 的局部 P1、P2-medium、P2-fine 已完成；P2 两档为 28,268/101,401 tet、130,794/447,679 DOF。内部速度差 1.30%，但 0.1/1/10 µm 近壁速度差 21.76/21.62/20.22%，G 差 21.75/21.47/19.31%，WSS 差均 15.95%，未通过冻结门槛。fine 正常退出，最终非线性残差 1.87e-10、全局流量缺口 9.36e-10；全局守恒不能代替局部精度。异常单元中心细分后有 6 个单侧梯度，η 范围 0.092–0.409，不能择一作改善证明。

近壁五类夹具及 ID1189/637 共 33,198 个保存状态重算完成：**IMPLEMENTATION_MATCHES_DECLARED_MODEL**，最大速度复现误差 3.06e-18 m/s；没有新积分、近壁源码或参数修改。硬接触启用的速度非光滑属于声明模型局限。按 A 区停止规则，B/C 只提取未求解，7 条新轨迹未执行；本轮科学作业已退出，没有自动下一阶段。

[中文报告](computation_examples/brava_flow_roi_18mlmin/p2p1_nearwall_cross_validation/FINAL_P2_NEARWALL_REVIEW_ZH.md) · [六张图文入口](computation_examples/brava_flow_roi_18mlmin/p2p1_nearwall_cross_validation/OPEN_RESULTS.html) · [机器状态](computation_examples/brava_flow_roi_18mlmin/p2p1_nearwall_cross_validation/CURRENT_EXECUTION_STATE.json)。原正式 ACTIVE_FLOW、旧 1565 条（1422 出流/141 截止/2 失败）和动画身份保持；不能改称新固定灌注协议结果。旧任务暂停记录只属于被替代路线，不是恢复授权。
<!-- LOCAL_P2_CURRENT_STATUS_END -->

<!-- LOCAL_FORMULATION_BOUNDARY_STATUS_START -->
## 2026-10-01 最新任务：Domain A 边界 benchmark 已技术性停止

本段是最新执行状态；旧全域加密、旧资源探测和上一轮局部 P2 任务均不恢复。本轮请求存档为 `computation_examples/brava_flow_roi_18mlmin/local_formulation_benchmark/provenance/USER_REQUEST.txt`。新目录 `local_formulation_benchmark/`，服务器对应 `/workspace/brava_flow_roi_18mlmin_20260928/local_formulation_benchmark/`。

状态 **TECHNICAL_STOP_BC_B_SOLVER_FAILED**；边界门槛 **NOT_EVALUATED_BC_B_SOLVER_FAILED**，科学主结论 `main_case=null`。仅新跑一次 Domain A 的 BC-B P2-fine，复用原 fine 网格、P2 空间、求解器参数、core 和全部固定探针，BC-A 未重跑。BC-B 按本轮最新提示使用正式全域三维截面数据：入口三分量速度、出口完整应力牵引；BC-A 仍为中心线阻力边界。正式 CFD 未重新校准或修改。

BC-B 在首档延拓 0.25、第 40 次 SNES 迭代后未收敛退出，残差 0.04123218，目标约 6.97e-8。耗时 15.645 分钟；8 MPI，实测使用约 7.62/7.68 CPU 核配额；进程树峰值 10.102 GiB、绝对 cgroup 峰值 16.786 GiB。无收敛 FIELDS，不能计算有效的 BC-A/B core 差异，不能宣称 BOUNDARY_DOMINATED，也不能用此失败证明 DC 优劣。未调参、重试或增加迭代上限。

边界提取/施加 QA 通过；原始 P1 应力内侧单侧导数仍存在离散相容性和物理精度限制。BC-B 入口原始截面流量 17.524174 mL/min 未强制归一化到 18。PETSc 实测不支持 CUDA/HIP/Kokkos Vec/Mat；因无合格 BC-B 字段，GPU 字段比较未运行。不要把预备脚本或历史 GPU 结果写成本轮新计算。

服务器监督程序 `brava_boundary_benchmark` 已 EXITED、autostart=false、autorestart=false，MPI 主进程 796084 已消失。未启动 DC、Nektar++、新轨迹或任何 B/C 计算。已有 `LOCAL_P2_NOT_YET_SPATIALLY_RESOLVED` 与 `IMPLEMENTATION_MATCHES_DECLARED_MODEL` 冻结不变；141 个原 particle 源文件身份检查通过，未重复 near-wall 验证。

[最终判断](computation_examples/brava_flow_roi_18mlmin/local_formulation_benchmark/FINAL_DECISION_ZH.md) · [中文边界报告](computation_examples/brava_flow_roi_18mlmin/local_formulation_benchmark/FINAL_BOUNDARY_REVIEW_ZH.md) · [两张图文入口](computation_examples/brava_flow_roi_18mlmin/local_formulation_benchmark/OPEN_RESULTS.html) · [停止依据](computation_examples/brava_flow_roi_18mlmin/local_formulation_benchmark/WHY_CONTINUE_OR_STOP_ZH.md)。建议下一阶段先做有界的人工截面牵引传递/非线性收敛相容性诊断，仅为建议，未经新指令不得自动执行或恢复旧路线。
<!-- LOCAL_FORMULATION_BOUNDARY_STATUS_END -->


<!-- BOUNDARY_TRANSFER_DIAGNOSIS_STATUS_START -->
## 2026-10-01 最新：Domain A BC-B 有界诊断完成，全部停止

完成时间 2026-10-01T20:23:13.213210+00:00。本段覆盖上文旧任务的 ACTIVE/恢复表述；未经新指令不得自动继续。最新请求存档于 `computation_examples/brava_flow_roi_18mlmin/boundary_transfer_diagnosis/provenance/USER_REQUEST.txt`。本轮目录 `boundary_transfer_diagnosis/`，服务器 `/workspace/brava_flow_roi_18mlmin_20260928/boundary_transfer_diagnosis/`。

主诊断 **TECHNICALLY_INCONCLUSIVE_AFTER_BOUNDED_DIAGNOSIS**；结果状态 **BOUNDARY_SENSITIVITY_NOT_EVALUABLE**。旧 `local_formulation_benchmark/` 失败案例 48 项冻结文件原样保留。BC-A 未重跑；网格、447679 DOF、8 MPI/1 线程、物理、NS 容差及 0.25/0.5/1 延拓均保持。

源 P1 三截面流量：入口 17.524174 mL/min、出口合计 17.547795 mL/min，缺口 0.134792%，与原 P1 积分散度一致。入口原方法与严格 P2 trace L2 投影相对差 4.06e-16；出口原载荷与 quadrature 直接积分相对差约 2.38e-15。两侧总牵引跳跃在固定压力基准下为 1.02%/0.40%，但切向分量为 24.26%/33.02%；保留不确定性，未证明它是失败主因。固定 BC-A 状态的出口替换代数扰动大于入口替换，不能将残差当物理误差。

全域投影初值残差改善 6.51 倍，未达到预设十倍强相关判据。Stokes 48.25 秒收敛、残差 2.34e-12，仅 DIAGNOSTIC_ONLY。ns_full_projection：NOT_CONVERGED，最后档位 0.5，KSP -3/3000，9.27 分钟；ns_stokes：NOT_CONVERGED，最后档位 0.5，KSP -3/3000，9.43 分钟。旧失败首因分类 NONLINEAR_CONVERGENCE_STALL；本轮详细证据和线性退出原因见 NS_DIAGNOSIS.json，不把失败字段当有效 Boundary Gate 比较。

原正式 CFD、旧 1565 条（1422 出流、141 截止、2 数值失败）、动画及固定灌注协议的结果身份不变。未重跑 near-wall/parity/网格收敛/粒子；`LOCAL_P2_NOT_YET_SPATIALLY_RESOLVED` 与 `IMPLEMENTATION_MATCHES_DECLARED_MODEL` 仍按历史范围保留。没有启动 DC、Nektar++、B/C、extra-fine 或全域加密。Supervisor `brava_boundary_diagnosis` 已退出，autostart=false、autorestart=false；双端交付校验文件位于本轮目录。

[中文诊断报告](computation_examples/brava_flow_roi_18mlmin/boundary_transfer_diagnosis/FINAL_BOUNDARY_DIAGNOSIS_ZH.md) · [四张图文入口](computation_examples/brava_flow_roi_18mlmin/boundary_transfer_diagnosis/OPEN_RESULTS.html) · [当前状态](computation_examples/brava_flow_roi_18mlmin/boundary_transfer_diagnosis/CURRENT_EXECUTION_STATE.json)。下一步仅建议保持边界/网格固定，针对失败档位的线性化系统与预条件器做最小有界审查；不自动执行。
<!-- BOUNDARY_TRANSFER_DIAGNOSIS_STATUS_END -->
