# 工作流辅助存储

这里存放当前项目所需的共享环境、Git 元数据、溯源契约和操作回执。业务代码仍在各自的当前项目目录中。

| 目录 | 用途 |
|---|---|
| `ulm-3d-vascular-model-generation/` | 当前流场、P9A5 等工作树共用的 Git 元数据库 |
| `ulm_particle_3d_particle0/` | 共享 Python 环境 `.venv/` 及其 Git 工作树标记 |
| `formal_3D_flow_solver_FEM/` | 旧 FEM 元数据和当前网络仍读取的几何溯源契约 |
| `cleanup_20260927/` | 上一次删除操作的清单、哈希和检查回执 |
| `migration_20260927/` | 本次目录迁移、路径修复和验证回执 |

共享 Python：`/home/lzy/projects/temp_storage/ulm_particle_3d_particle0/.venv/bin/python`。

本次采用同文件系统目录重命名，完整迁移原有内容；修复 Git 工作树指针、虚拟环境启动器、激活脚本、依赖链接及当前入口说明。原位置没有留下兼容目录或软链接。服务器路径与服务器文件没有变更。

[当前业务工作流](../ACTIVE_VASCULAR_WORKFLOW.md) · [迁移验证](migration_20260927/verification.json)
