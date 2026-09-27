# 求解器支持工具

原阶段名称已改为功能名称；这里的工具不参与表面导入或 TetGen 网格生成。

| 当前路径 | 用途 |
|---|---|
| `src/sv_solver_support/flow_parser.py` | 求解器日志解析 |
| `src/sv_solver_support/solver_checks.py` | 线性/非线性求解记录与验收检查 |
| `src/sv_solver_support/transient_checks.py` | 瞬态、重启与相邻场变化检查 |
| `src/sv_solver_support/steady_monitor.py` | 稳态停止条件与场等价性检查 |
| `src/sv_solver_support/parallel_checks.py` | 并行数值一致性检查 |
| `src/sv_solver_support/checkpoint_checks.py` | 原生检查点及 PETSc 运行语义检查 |
| `use_gpu_env.sh` | 既有服务器 CUDA/MPI 环境包装脚本 |
| `include/preconditioner_reuse.h` | PETSc 预条件器复用实现 |
| `configs/preconditioner_reuse_policy.json` | 保留的复用策略与测量记录 |
| `scripts/flow_figures.py` | 原始流场检查图脚本 |
| `../external/svMultiPhysics-reuse/` | 求解器源码与原 Git 元数据 |

Python 包从 `sv_solver_support` 导入；通用几何、后处理和文件校验仍从 `sv_validation` 导入。项目 `pyproject.toml` 和本地可编辑安装已指向新的源码位置。

本次仅迁移、重命名并修正路径/导入。检查模块中读取历史阶段报告的辅助函数仍需要调用者提供当时的文件；本目录不构成新的远程一键部署入口。原补丁生成与远程安装脚本已经缺少历史源树及 `runner_remote`、`environment_remote` 等依赖，原文迁至 [迁移回执的 reference_recipes](/home/lzy/projects/temp_storage/mesh_code_migration_20260927/reference_recipes/)，未作为可运行入口保留。

服务器 `/workspace/...` 路径沿用服务器实际位置，见 [CURRENT_SERVER_PATHS.md](/home/lzy/projects/CURRENT_SERVER_PATHS.md)。C++ 内部类型、函数及停止回执协议名称保留，以保持与当前二进制/日志的对应关系；头文件本身已更名为 `preconditioner_reuse.h`，包含路径同步更新。本次没有重编译或重新部署远程求解器。

当前 H0 求解任务仍使用 [独立流场工作树](/home/lzy/projects/ulm_flow_mean_2p0_mmps/README.md) 中的冻结代码。本地源码迁移不修改其输入哈希、服务器二进制或当前结果。当前流线、速度、压力、WSS 展示入口是 [rotate_visualization](../rotate_visualization/OPEN_RESULTS.html)。
