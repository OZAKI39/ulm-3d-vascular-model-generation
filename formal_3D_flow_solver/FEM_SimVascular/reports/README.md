# 当前有效报告

本目录保留当前表面、网格和官方工具依赖对应的 12 个文件，均在 `mesh_and_flow/`。

| 文件 | 用途 |
|---|---|
| `surface_adapter.json` | 输入表面转换记录，绑定原始表面及当前 `model/source.vtp` 的哈希 |
| `geometry_import.json` | 官方表面导入结果、面标签及点/三角形数量；网格生成的前置检查 |
| `geometry_qc.json` | 当前网格与输入表面的几何误差、端口面积及法向检查 |
| `mesh_policy_provenance.json` | 当前网格尺寸策略的路径与哈希；完整网格检查直接读取 |
| `mesh_validity.json` | 体网格有效性、连通性、体积及原始/导出网格哈希 |
| `mesh_quality.json` | 当前四面体 minSICN 质量分布 |
| `simvascular_distribution.json` | 当前官方发行包安装位置；内嵌 Python 启动器直接读取 |
| `simvascular_api.json` | 当前官方网格 API 检查结果；表面准备与导入的前置检查 |
| `source_geometry.png` | 当前输入表面展示 |
| `simvascular_surface.png` | 当前检查后的血管表面展示 |
| `mesh_cutaway.png` | 当前四面体网格剖切展示 |
| `mesh_quality.png` | 当前网格质量分布图 |

这些报告已与当前表面、网格策略及网格哈希核对，网格为 70,363 个节点、371,402 个四面体。保留依据是当前有效性和工具依赖，文件创建日期不作为单独的删除依据。

2026-09-27 删除了旧未验收流场的速度/压力/收敛图、旧运行及边界条件报告、过时构建和测试记录、旧 FEM 对比快照及旧总报告，共 39 个文件。旧示例流体测试报告也已删除，表面入口改为读取当前 `simvascular_api.json`，几何处理逻辑保持不变。

新边界条件 H0 流场的正式物理验收见 [H0 流场工作区](/home/lzy/projects/ulm_flow_mean_2p0_mmps/README.md) 和 [H0 验收说明](/home/lzy/projects/ulm_3D_vascular/reports/a_network_1d0d_boundary_v2_idealized/ROI_3D_A_H0_PRESSURE_VALIDATION_ZH.md)。当前流线、速度、压力和 WSS 展示见 [可视化总入口](../rotate_visualization/OPEN_RESULTS.html)。

本次删除与验证记录：[清理回执](/home/lzy/projects/temp_storage/mesh_reports_cleanup_20260927/result.json)。
