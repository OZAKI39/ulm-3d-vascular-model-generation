# 当前血管、流场、微泡与 RBC 工作流

2026-09-27 已执行 WSL 历史内容清理。本文件是当前路径入口；旧同步快照和旧目录总表已删除。服务器内容和独立的血管打印项目保留。

| 当前内容 | 本地路径 |
|---|---|
| 血管 A、ROI、表面与 1D/0D Network H0 | [ulm_3D_vascular](ulm_3D_vascular/README.md) |
| 原始小鼠几何与共享预处理代码 | `vascular_printing/`；`ulm_3D_vascular` 中保留兼容链接 |
| 当前 FEM：ROI-only-balanced-pressure-v1 | [新流场算例](temp_storage/github_sync_20260927/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-ROI-only-balanced-pressure-v1/) |
| 保留 H0 FEM 基线 | [H0 算例](ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-H0-pressure-v1/) |
| 表面导入、TetGen 四面体生成及网格检查 | [mesh_generate](formal_3D_flow_solver/FEM_SimVascular/mesh_generate/README.md) |
| formal_3D_flow_solver 目录用途与结构 | [目录说明](formal_3D_flow_solver/README.md) |
| 求解器日志、检查点、GPU/PETSc 支持 | [solver_support](formal_3D_flow_solver/FEM_SimVascular/solver_support/README.md) |
| 当前微泡生成、计算和轨迹，1500 条、dt=0.5 ms | [ulm_particle_formal_p9a5](ulm_particle_formal_p9a5/README.md) |
| RBC 代码、轨迹和展示 | [RBC/微泡结果索引](ulm_particle_formal_p9a5/CURRENT_RESULTS.md) |
| 流线、速度、压力、WSS 旋转展示 | [OPEN_RESULTS.html](formal_3D_flow_solver/FEM_SimVascular/rotate_visualization/OPEN_RESULTS.html) |
| 保留 H0 残差表格 | [A_Global_Transient_Evolution_A_H0.xlsx](A_Global_Transient_Evolution_A_H0.xlsx) |
| 当前服务器路径 | [CURRENT_SERVER_PATHS.md](CURRENT_SERVER_PATHS.md) |
| 本次清理回执 | [清理回执](temp_storage/cleanup_20260927/result.json) |

2026-09-28 最新流场 SHA256：`fb3c6ad0815d156b09b289fc347766486f24ddb521b8e51e8bbc4cf1fd2ace12`，对应 `mean-2p0-mmps-A-ROI-only-balanced-pressure-v1/frozen_flow/steady_flow_mean_2p0_mmps_A_ROI_only_balanced_pressure.vtu`。新微泡 1500 条及动画已完成，具体数据与核验见 [最新结果](ulm_particle_formal_p9a5/CURRENT_RESULTS.md)。H0 的残差表格和 CORE500 是单独保留的基线。

`temp_storage/ulm_particle_3d_particle0/` 现在仅保留共享 Python 环境与 Git 标记；旧工程内容已删除。该环境还被独立血管打印项目使用，并引用 `formal_3D_flow_solver/FEM_SimVascular/.venv/`。`temp_storage/ulm-3d-vascular-model-generation/.git/` 是当前工作树共用的 Git 元数据库，不能随旧主检出内容一起删除。

仍带旧阶段名的少量文件是当前代码直接读取的依赖：v1 网络图与端口数据、P9A4 冻结出生样本、NEW 流场适配器、网格/黏度契约、RBC 原始场以及部分轨迹时长验证输入。它们不代表仍保留旧开发流程。

辅助目录现统一存放于 [temp_storage](temp_storage/README.md)：Git 元数据、共享 Python、旧 FEM 溯源契约和清理回执。原目录已移除，当前代码直接使用新路径。
