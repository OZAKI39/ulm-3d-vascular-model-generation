# 表面导入、四面体网格生成与检查

本目录集中保存正式 SimVascular / TetGen 网格工具。所有脚本按自身位置解析上一级 `FEM_SimVascular` 数据根目录，可从任意工作目录调用。

| 功能 | 入口 |
|---|---|
| 将已准备的 ROI 表面转为带端口标签的 VTP | `scripts/prepare_surface.py` |
| 用官方 SimVascular 导入表面、核对端口 | `scripts/import_surface_model.py` |
| 冻结网格尺寸策略 | `scripts/freeze_mesh_policy.py` |
| 官方 TetGen 四面体网格生成 | `scripts/generate_tetra_mesh.py` |
| 拓扑、正体积、端口、几何误差和 minSICN 检查；导出求解网格 | `scripts/audit_mesh.py` |
| 网格及质量图片 | `scripts/mesh_figures.py` |
| 启动官方内嵌 Python | `scripts/run_meshing_python.py` |
| 官方 API 可用性与接口说明 | `scripts/probe_meshing_api.py`、`scripts/inspect_meshing_api.py` |
| 几何、网格诊断、通用校验及测量模块 | `src/vascular_validation/` |

输入与结果仍使用以下相对 `FEM_SimVascular/` 的路径：

- 输入：`inputs/fem_reference/`。
- 配置：`configs/face_map.json`、`configs/mesh_policy.json`。
- 表面：`outputs/mesh_and_flow/model/`。
- TetGen 原始输出：`outputs/mesh_and_flow/mesh_generation/primary/`，以及仅在允许条件下使用的 `fallback/`。
- 检查后的求解网格：`outputs/mesh_and_flow/solver_mesh/`。
- 质量、几何检查与图片：`reports/mesh_and_flow/`；运行日志：`logs/mesh_and_flow/`。

表面准备和官方导入均要求 `reports/mesh_and_flow/simvascular_api.json` 检查通过且网格 API 可用。旧示例流体测试报告已移除；当前报告清单见 [reports/README.md](../reports/README.md)。

普通 Python 脚本使用 `FEM_SimVascular/.venv/bin/python`。包含 `import sv` 的脚本必须通过官方内嵌 Python 启动，例如从 `FEM_SimVascular/` 运行：

```bash
.venv/bin/python mesh_generate/scripts/run_meshing_python.py probe_meshing_api.py
.venv/bin/python mesh_generate/scripts/run_meshing_python.py import_surface_model.py --label geometry_import
.venv/bin/python mesh_generate/scripts/run_meshing_python.py generate_tetra_mesh.py --label mesh_generation
.venv/bin/python mesh_generate/scripts/audit_mesh.py --attempt primary
```

这些命令是工具入口说明，不是对现有成果重复执行的批处理。网格生成保留原有“已存在体网格则停止”保护，网格策略保留冻结保护；本次迁移没有重新生成正式网格。当前 H0 新边界条件流场仍由独立 [H0 算例](/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-H0-pressure-v1/) 提供。
