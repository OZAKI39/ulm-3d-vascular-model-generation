# 当前网格工具、共享环境和流场可视化

- [旋转可视化](rotate_visualization/OPEN_RESULTS.html)：Network-H0 新边界条件流线、速度矢量、压力和 WSS，保留原有渲染风格。
- `rotate_visualization/input_data/`：当前可视化输入；`results/`：视频、4K 图片和矢量数据。
- [mesh_generate](mesh_generate/README.md)：表面准备/导入、TetGen 四面体生成、网格拓扑与质量检查、网格图片；代码已从原 `scripts/`、`src/` 迁入。
- `inputs/`、`configs/`、`outputs/sv1/`、`reports/sv1/`：原有输入、配置、网格和检查结果。迁移代码时保持这些数据路径与内容不变。
- [solver_support](solver_support/README.md)：求解器日志、检查点、稳态检查与 GPU/PETSc 支持，按功能命名。
- `external/SimVascularDistribution/`：官方内嵌 Python / TetGen 网格工具。
- `external/svMultiPhysics-reuse/`：保留原 Git 元数据和本地修改的求解器源码。
- `.venv/`：被当前共享 Python 环境引用的依赖，不能作为旧阶段缓存删除。

当前 FEM 算例位于 [ulm_flow_mean_2p0_mmps](/home/lzy/projects/ulm_flow_mean_2p0_mmps/README.md)。旧求解阶段、构建环境、性能试验、日志和旧渲染副本已删除。

2026-09-27 的代码迁移清单及验证记录：[mesh_code_migration_20260927](/home/lzy/projects/temp_storage/mesh_code_migration_20260927/README.md)。原 `scripts/sv13n`、`scripts/sv13q`、`configs/sv1_3q`、`external/sv13q` 容器目录已移除。独立 H0 算例中的冻结代码副本和服务器上的实际路径不在本次迁移范围内。
