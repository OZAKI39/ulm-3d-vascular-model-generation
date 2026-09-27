# 当前 Network-H0 3D FEM 流场

当前算例：`formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-H0-pressure-v1/`。

- `SV_MESH/`：四面体网格和边界表面。
- `run/`：当前 solver.xml、PETSc 参数、计算日志、边界历史及求解输出。
- `frozen_flow/steady_flow_mean_2p0_mmps_A_H0.vtu`、`flow_arrays_si.npz`：当前新边界条件稳态场。
- `reports/physics_validation_H0.json`：物理验收。
- `field_diagnostics/excel_export_20260925/`：当前残差/收敛表格导出。

1D/0D 网络与算例生成、远程求解、验收代码位于 [`../ulm_3D_vascular`](../ulm_3D_vascular/README.md)。本工作树保留当前求解器源码 `vendor/svMultiPhysics_stage_q/`、共享验证模块和网格准备代码。

`mean-2p0-mmps/` 仅保留生成 H0 算例、输入哈希验证及现有 RBC 展示复现实际需要的文件；旧日志、诊断和结果已删除。`frozen_reference` 为构建与校验依赖，不是当前新边界条件场。

[当前工作流](../ACTIVE_VASCULAR_WORKFLOW.md) · [旋转可视化](../formal_3D_flow_solver/FEM_SimVascular/rotate_visualization/OPEN_RESULTS.html)
