# 表面导入、四面体网格生成与检查

本目录集中保存正式 SimVascular / TetGen 网格工具。所有脚本按自身位置解析上一级 `FEM_SimVascular` 数据根目录，可从任意工作目录调用。

| 功能 | 入口 |
|---|---|
| 将已准备的 ROI 表面转为带端口标签的 VTP | `scripts/prepare_surface.py` |
| 用官方 SimVascular 导入表面、核对端口 | `scripts/sv_import_model.py` |
| 冻结网格尺寸策略 | `scripts/freeze_mesh_policy.py` |
| 官方 TetGen 四面体网格生成 | `scripts/sv_generate_mesh.py` |
| 拓扑、正体积、端口、几何误差和 minSICN 检查；导出求解网格 | `scripts/audit_mesh.py` |
| 网格及质量图片 | `scripts/mesh_figures.py` |
| 启动官方内嵌 Python | `scripts/run_sv_python.py` |
| 官方 API 可用性与接口说明 | `scripts/sv_api_probe.py`、`scripts/sv_inspect_api.py` |
| 几何、网格诊断、通用校验及测量模块 | `src/sv_validation/` |

输入与结果仍使用以下相对 `FEM_SimVascular/` 的路径：

- 输入：`inputs/fem_reference/`。
- 配置：`configs/face_map.json`、`configs/mesh_policy.json`。
- 表面：`outputs/sv1/model/`。
- TetGen 原始输出：`outputs/sv1/mesh_generation/primary/`，以及仅在允许条件下使用的 `fallback/`。
- 检查后的求解网格：`outputs/sv1/SV_MESH/`。
- 质量、几何检查与图片：`reports/sv1/`；运行日志：`logs/sv1/`。

普通 Python 脚本使用 `FEM_SimVascular/.venv/bin/python`。包含 `import sv` 的脚本必须通过官方内嵌 Python 启动，例如从 `FEM_SimVascular/` 运行：

```bash
.venv/bin/python mesh_generate/scripts/run_sv_python.py sv_api_probe.py
.venv/bin/python mesh_generate/scripts/run_sv_python.py sv_import_model.py --label geometry_import
.venv/bin/python mesh_generate/scripts/run_sv_python.py sv_generate_mesh.py --label mesh_generation
.venv/bin/python mesh_generate/scripts/audit_mesh.py --attempt primary
```

这些命令是工具入口说明，不是对现有成果重复执行的批处理。网格生成保留原有“已存在体网格则停止”保护，网格策略保留冻结保护；本次迁移没有重新生成正式网格。当前 H0 新边界条件流场仍由独立 [H0 算例](/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-H0-pressure-v1/) 提供。
