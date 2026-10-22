# 当前有效的表面与网格数据

本目录只保留当前网格工具链使用的数据：

| 目录 | 保留用途 |
|---|---|
| `model/` | 带端口标签的源表面和官方 SimVascular 导入表面；供表面导入与 TetGen 生成使用 |
| `mesh_generation/primary/` | 已成功生成的原始四面体网格、分面及生成记录；`audit_mesh.py` 直接读取它们进行完整复核 |
| `solver_mesh/` | 检查后的体网格、外表面、入口/出口/壁面、数组及壁面距离；供网格检查和绘图使用 |

`solver_mesh/` 的 9 个文件已与当前 [Network-H0 正式算例](/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-H0-pressure-v1/SV_MESH/) 核对，逐字节一致：70,363 个节点、371,402 个四面体。

正式新边界条件流场位于该算例的 [frozen_flow/](/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-H0-pressure-v1/frozen_flow/)。

2026-09-27 已删除旧未验收血管流场、早期官方示例测试输出、旧绘图缓存和 API 说明快照，共 21 个文件路径。剩余 19 个数据文件哈希保持不变。绘图缓存和 API 说明可由对应工具按需重新生成。

删除清单和哈希记录：[清理回执](/home/lzy/projects/temp_storage/mesh_outputs_cleanup_20260927/result.json)。
