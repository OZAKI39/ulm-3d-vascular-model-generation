# formal_3D_flow_solver 的用途

这个文件夹保存当前血管工作流使用的 **3D FEM 网格工具、求解器支持工具、共享依赖，以及新边界条件流场的旋转可视化**。实际内容集中在 `FEM_SimVascular/`。

| 目录或入口 | 作用 |
|---|---|
| [mesh_generate](FEM_SimVascular/mesh_generate/README.md) | 导入已准备的 ROI 表面，生成 TetGen 四面体网格，检查拓扑、单元质量、端口和几何误差，输出网格检查图 |
| [solver_support](FEM_SimVascular/solver_support/README.md) | 解析求解器日志，检查稳态、原生检查点和并行一致性，保存 GPU/PETSc 环境包装与预条件器复用支持代码 |
| [inputs/fem_reference](FEM_SimVascular/inputs/fem_reference/) | 网格工具所需的表面、参考网格、端口约定及输入来源信息 |
| [configs](FEM_SimVascular/configs/) | 网格尺寸策略、面标签映射和基础参考求解配置；`flow_solver.xml` 是保留的基础参考算例配置 |
| [outputs/mesh_and_flow](FEM_SimVascular/outputs/mesh_and_flow/README.md) | 仅保留当前表面 `model/`、原始四面体 `mesh_generation/primary/` 和检查后的 `solver_mesh/`；旧流场及示例测试输出已删除 |
| [reports](FEM_SimVascular/reports/README.md) | 当前网格质量、端口、几何一致性、工具 API 检查及网格图片；旧流场、构建和测试报告已清理 |
| [旋转可视化总入口](FEM_SimVascular/rotate_visualization/OPEN_RESULTS.html) | 当前 Network-H0 新边界条件的流线、局部速度矢量、压力和 WSS 视频及 4K 图片；配套输入在 `input_data/` |
| `external/SimVascularDistribution/` | 官方 SimVascular 内嵌 Python 和 TetGen 网格工具 |
| `external/flow_solver_source/` | 带预条件器复用修改的 svMultiPhysics 求解器源码，保留原 Git 元数据 |
| `.venv/` | 本目录 Python 依赖，也被共享 Python 环境引用 |

表中未写完整前缀的路径均相对 `FEM_SimVascular/`。

当前正式计算分布在以下工作区：

- 血管 A、ROI、1D/0D 管网：[ulm_3D_vascular](/home/lzy/projects/ulm_3D_vascular/README.md)。
- 新边界条件 3D FEM 正式算例：[mean-2p0-mmps-A-H0-pressure-v1](/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-H0-pressure-v1/)。
- 当前微泡生成与轨迹：[ulm_particle_formal_p9a5](/home/lzy/projects/ulm_particle_formal_p9a5/README.md)，时间步长为 1 ms。
- RBC 代码、轨迹和展示：[当前结果索引](/home/lzy/projects/ulm_particle_formal_p9a5/CURRENT_RESULTS.md)。
- 服务器算例、二进制与运行环境：[CURRENT_SERVER_PATHS.md](/home/lzy/projects/CURRENT_SERVER_PATHS.md)。

2026-09-27 已将自有代码、配置与数据目录改为功能命名：`vascular_validation`、`flow_solver_support`、`import_surface_model.py`、`generate_tetra_mesh.py`、`mesh_and_flow/`、`solver_mesh/` 等。当前目录没有保留旧名字的兼容文件或软链接。

官方依赖内部的 `sv` Python API、XML 格式、C++ 文件名和实际服务器路径保留原定义；历史运行记录保留当时的原始路径。它们与本项目自有文件的命名区分管理。命名整理没有改变血管几何、网格数值、物理参数、正式轨迹或可视化风格。

命名映射、原始哈希和复核结果见 [整理记录](/home/lzy/projects/temp_storage/descriptive_names_20260927/README.md)。
