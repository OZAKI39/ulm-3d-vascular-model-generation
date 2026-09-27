# 网格报告目录清理回执

范围：`formal_3D_flow_solver/FEM_SimVascular/reports/`。

已删除 39 个旧文件，包括失败流场的验收报告及图表、旧构建/测试日志、旧 FEM 对比快照、过时总报告和示例测试记录。保留 8 个当前网格/工具 JSON 与 4 张当前几何/网格图片，保留文件字节不变。

- `before.json`：删除/保留清单、文件哈希、大小，以及受保护输入和网格数据哈希；没有复制被删除的旧报告。
- `before_prepare_surface.py`、`before_import_surface_model.py`：两处前置条件更新前的少量源码。
- `verify_surface_entrypoints.py` 与对应日志：在隔离目录实际执行表面准备和官方导入；不用旧示例报告，输出几何、标签、端口映射与现有数据一致；API 状态失败或网格 API 不可用时停止。
- `retained_reports_verification.json`：保留报告与实际表面、网格质量、网格文件、策略和官方启动器的一致性检查。
- `result.json`：删除及验证汇总。

前置条件从旧流体示例的 `official_smoke.json` 改为当前 `simvascular_api.json`。几何转换、网格生成、流场和轨迹算法未修改。隔离验证产生的副本在检查后清理。
