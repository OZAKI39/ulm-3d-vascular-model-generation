# 血管 FEM 网格、求解器支持与可视化

完整用途与当前工作区分工见 [formal_3D_flow_solver 的用途](../README.md)。

- [mesh_generate](mesh_generate/README.md)：表面导入、TetGen 四面体生成、网格拓扑/质量检查和网格图片。
- [solver_support](solver_support/README.md)：日志解析、稳态/检查点/并行验收及 GPU/PETSc 支持。
- `inputs/fem_reference/`：准备好的表面、参考网格及端口约定。
- `configs/`：`mesh_policy.json`、`face_map.json`、`reference_physics.yaml`、`flow_solver.xml` 等基础配置。
- [outputs/mesh_and_flow](outputs/mesh_and_flow/README.md)：当前导入表面、TetGen 原始输出和 `solver_mesh/`；旧未验收流场与早期示例测试数据已删除。
- [reports](reports/README.md)：当前网格质量、几何、端口和工具 API 报告，以及对应网格图片。
- [rotate_visualization/OPEN_RESULTS.html](rotate_visualization/OPEN_RESULTS.html)：当前 H0 新边界条件流场的流线、速度矢量、压力、WSS 视频与 4K 图片。
- `external/SimVascularDistribution/`：官方网格工具；`external/flow_solver_source/`：求解器源码。
- `.venv/`：本地依赖和共享 Python 环境的依赖来源。

当前 H0 正式求解配置、冻结流场和求解日志位于 [独立流场工作区](/home/lzy/projects/ulm_flow_mean_2p0_mmps/README.md)。本目录的 `reports/` 仅保留当前网格与工具检查文件。

Python 包名为 `vascular_validation` 与 `flow_solver_support`，安装项目名为 `vascular-workflow-tools`。官方 `import sv` 仍由 SimVascular 内嵌 Python 提供。自有路径命名及校验清单已更新；旧阶段目录、包名与兼容链接均未保留。
