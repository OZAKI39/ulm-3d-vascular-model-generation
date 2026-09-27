# 2026-09-27 网格代码迁移回执

当前工具入口：[mesh_generate](/home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/mesh_generate/README.md)、[solver_support](/home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/solver_support/README.md)。

- `before.json`：迁移映射、代码原始 SHA256、数据/展示/求解器源树原始 SHA256 和 Git 身份。
- `before_code/`、`before_metadata/`：本次涉及的少量原代码/路径配置备份，不在活跃源码搜索范围。
- `migrate.py`：本次已执行的迁移脚本，拒绝重复执行。
- `reference_recipes/`：原 `sv13q` 的补丁生成与远程构建脚本原文。它们依赖已删除的历史输入及部署模块，仅供溯源，不是当前可运行入口。
- `verification.json`：迁移后的路径、数据完整性与代码核对结果。
- `result.json`：最终汇总，含官方 Python/API、现有网格、包路径、49 项管网测试和 7,760 个当前科学文件的完整性检查。
- `mesh_check_reports/`、`mesh_check.log`：在隔离目录中执行完整网格检查所得报告。节点、单元、面标签、邻接和质量数组与原网格完全一致，测试生成的大文件已清理。

正式输入、网格、流场和展示数据保持原位置；没有删除其内容，也没有重新计算或替换当前场。服务器路径和独立 H0 工作树中的冻结代码未迁移。
