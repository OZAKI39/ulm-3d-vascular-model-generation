> 历史记录：2026-10-22。当前删除结果和待确认项目请看 [2026-10-23 清理报告](CLEANUP_REPORT_ZH.md)。本页不代表当前磁盘状态。

# 旧目录待清理清单

本轮没有删除旧工程、科研数据、兼容链接或服务器内容。下表是待用户确认的清单，不是已经执行的删除记录。文件级迁移证明见 `audit/data_manifest.json`、`audit/source_map.json` 和 `audit/original_files_after.json`。

| 原位置 | 已迁入内容 | 建议与前置条件 |
|---|---|---|
| `/home/lzy/projects/computation_examples/mean-2p0-mmps-A-ROI-only-balanced-pressure-v1` | 冻结场、原生状态、配置日志、1500 条原始轨迹、队列、验收、图表/动画 | 本轮保留；新目录完成备份且确认所有外部引用后，才评估撤下旧副本 |
| `/home/lzy/projects/computation_examples/mean-2p0-mmps-A-best-feasible-balance-v1` | 独立流场及另 1500 条原始轨迹、全部纳入的正式证据 | 同上，不能与 ROI-only 合并成一个算例 |
| `/home/lzy/projects/computation_examples/mean-2p0-mmps-A-H0-pressure-v1` | H0 场、500 条正式轨迹及现有报告 | 同上；未完成的历史媒体交付状态保留 |
| `/home/lzy/projects/computation_examples/mean-2p0-mmps`、`fem_reference_stage_q` | 原始基线和 Stage-Q 冻结参考、checkpoint 与已有验证 | 科研溯源数据，先保留原位置；不得称为可随意删除的缓存 |
| `/home/lzy/projects/computation_examples/_shared/particle_workflow/particle_3d/src/particle_3d` | 当前依赖闭包、轻量验证模块；114 受保护原源码的完整字节证据 | 仅从“小鼠当前运行”角度已解除依赖。共享代码可能仍服务 BraVa/其他项目，不能直接删整个 `_shared/particle_workflow` |
| `/home/lzy/projects/computation_examples/_shared/active_network_development/ulm_3D_vascular/network_1d0d` | 当前 21 个网络模块及对应输入/设计报告 | 比旧 `geometry_network/network_1d0d` 更完整的权威实现已迁入；删除旧 worktree 前必须检查 Git 元数据和其他消费者 |
| `/home/lzy/projects/computation_examples/_shared/geometry_network` | 必要几何/网络输入已选择迁入；重复旧网络实现未再复制为运行代码 | 该目录存在打印项目链接，不得递归删除链接目标；原始 Git 历史仍有价值 |
| `/home/lzy/projects/computation_examples/_shared/formal_3D_flow_solver` | 自有网格/FEM/WSS 支持代码、有效 native 源、补丁和必要网格 | 旧 SimVascular GUI、已安装原生库和共享环境未复制。新原生环境尚未验证，不能直接删除整个旧目录 |
| `/home/lzy/projects/computation_examples/_shared/legacy_flow_workflow` | 冻结数据与实际执行方法已在新工程归一；原 native runner 逻辑迁入 | 原构建资料与独立实验仍可能依赖此目录，暂保留 |
| `_shared/wss_audit`、`_shared/wss_validation_v2`、`validation_cases` | 已执行证据与取消状态迁入 `validation/` | 保留旧记录；不能恢复细档/压力敏感性。未迁入的失效中间帧详见排除清单 |
| `/home/lzy/projects/temp_storage/mouse_flow_github_sync_20261022` | 发布清单和必要身份证据已复制 | 是已发布快照及本地 Git worktree，删除必须按 worktree 管理处理；本轮保留 |
| 各旧工作区 `__pycache__`、`.pytest_cache`、中间渲染帧、可再生 GPU 缓存 | 不作为当前运行依赖 | 可列为低风险清理候选，但本轮没有执行删除；不要把原始 `.npz` 轨迹或 checkpoint 当作缓存 |
| 共享 Python 环境、`_shared/python_environment` 与旧顶层兼容链接 | 新工程不调用它们 | 其他项目可能仍依赖，未确认前保留 |
| `/home/lzy/projects/vascular_printing`、`new_vessel_printing`、BraVa、声学微泡目录 | 只提取确需的小鼠源样本及几何函数，没有整理这些工程 | 不在本次清理范围内 |

## 已确认的重复与未知项

`audit/duplicate_data_review.json` 按内容哈希列出新项目中的相同数据。部分相同网格位于不同原生算例包或成品媒体输入包中，原始清单使用包内相对路径校验；当前保留这些必要副本，不通过软/硬链接维持对旧目录的依赖。不同身份契约原文全部保留。

`audit/published_EXCLUDED.json` 保留此前未同步文件的来源、原因与已记录哈希；`audit/selection_decisions.json` 记录本次未作为运行数据复制的项目。未执行目录的细档大网格、可重建缓存、上游测试 LFS 指针、安装环境等仍在旧位置或原记录中，不应自动推断它们已无用途。

历史 `particle6_5/data/07_real_v1_dt1.json` 在本轮开始前就不在当前本地/快照中；三个相关旧测试明确 BLOCKED。不能将其“补造”后替代原始验证证据。

需要继续保留/单独配置的外部依赖：原 svMultiPhysics/PETSc/MPI/CUDA 栈、SimVascular TetGen API、Ultraliser/VMTK；可选 LAMMPS bridge。它们与已经完成独立安装的 Python 环境、CPU C++ 几何内核区分记录。
