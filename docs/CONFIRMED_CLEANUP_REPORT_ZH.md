# BraVa 待确认清单处理报告

日期：2026-10-23。根据用户“确认处理【待确认清单】中的内容”，**9 项已全部处理，待处理项为 0**。

旧工程和旧备份已实际删除。正式数据、必要失败记录和代码历史保存在唯一工程 `/home/lzy/projects/new_human_vessel_flow`。这次没有建立另一套备份或同步工程。

## 1. 实际删除了什么

以下路径均相对于 `/home/lzy/projects`，9 个目标互不包含；删除后均不存在。

| 目标 | 原分配空间 GB | 处理 |
|---|---:|---|
| `computation_examples/brava_flow_roi_18mlmin` | 25.566 | 旧工程删除；必要内容和代码历史先保留 |
| `computation_examples/brava_calibration_equal_split` | 0.277 | 旧校准目录删除；该工况独立冻结结果保留 |
| `server_migration/brava_new_server_20261009` | 0.025 | 旧迁移工作目录删除；连接、环境说明及必要验证记录保留 |
| `server_migration/brava_vast4090_20261009T094459Z` | 46.374 | 旧服务器完整本地备份删除；必要科学证据先保留 |
| `computation_examples/_archives/github_sync_20260928_latest/brava_flow_roi_18mlmin` | 1.978 | 只删除 BraVa 子目录，保留混合工作区及共享 Git |
| `github_sync_audits/flow_projects_20261022/BRAVA_INCLUDED.json` | 0.0045 | 旧发布清单删除 |
| 同目录 `BRAVA_EXCLUDED.json` | 0.0024 | 旧发布清单删除 |
| 同目录 `BRAVA_PUBLISH_PLAN.json` | 0.0024 | 旧发布清单删除 |
| 同目录 `BRAVA_PUBLISH_PLAN_BEFORE_LFS_EXCLUSION.json` | 0.0027 | 旧发布清单删除 |

完整路径与完成时间见[删除清单](cleanup_confirmed_20261023/DELETED_FILES.csv)及[逐项执行凭据](cleanup_confirmed_20261023/deletion_receipt.json)。[原待确认清单](cleanup_20261023/PENDING_REVIEW.csv)已全部标为 `DELETED_AFTER_USER_CONFIRMATION`。

## 2. 释放了多少空间

GB 为十进制。按实际分配空间测量，没有把逻辑大小、父子目录或硬链接重复相加。完整扫描了 428,795 个文件系统条目，没有发现目标之外的硬链接，也没有跨项目删除共享对象。

| 不重叠的测量范围 | 清理前 GB | 清理后 GB |
|---|---:|---:|
| 本轮 9 个旧目标 | 74.232 | 0 |
| 新 Human 工程，含环境、保留数据和 Git | 11.026 | 12.284 |
| 合计 | 85.258 | 12.284 |

**上述范围净减少 72.974 GB，约 73.0 GB。** 新工程增加的空间主要是导入的原仓库历史、补齐的必要科学资料及审计记录。当前 `.git` 约 1.043 GB，保留历史不代表又建了一个运行工程。

整块 WSL 文件系统可用空间从 684.306 GB 变为 759.689 GB；它同时受其他活动影响，因此不能把其全部增量归给本任务。本轮没有压缩或收缩 Windows 的 WSL 虚拟磁盘文件。数据来自[清理前](cleanup_confirmed_20261023/disk_before.json)、[清理后](cleanup_confirmed_20261023/disk_after.json)。测量在最终文档和 Git 提交前，后续少量元数据会使字节数略有变化。

## 3. 删除前专门保留了什么

在已有正式资料之外，补齐 **628 个必要文件、238,992,404 字节（约 239 MB）**；相同内容只引用已有一份，没有把旧工作区整套复制回来。来源和逐文件 SHA 见[保留凭据](cleanup_confirmed_20261023/preservation_receipt.json)。

- 两份事件区中、细档对照流场、原体网格、最后保存步 81 的 checkpoint 和必要执行记录。以前仅保留了部分分析输出，本轮补齐了原始对照输入，详见[输入索引](../diagnostics/field_events/reference_fields/README.md)。它们仍是历史研究解，不是新的正式生产结果。
- 三个已停止微泡批次的契约、停止/失败记录与关键日志。它们的大量完整旧状态及中间计算结果已删除；正式 `contactfix` 批次完整原始状态全部保留。
- P2/P1 研究实际编译的源码变体、必要边界输入、解一致性检查输出及失败记录。原来的失败和资源停止结论不变。
- 历史正式求解器二进制约 6.3 MB，仅保留其身份。当前求解软件及环境仍在 `new_vessel_flow_solvers`，没有删除或换成这个旧程序。
- 旧服务器的必要连接入口、环境配方、版本信息和验证证据。旧完整备份已退出使用，不再声称可以直接恢复整套旧服务器工作区；本轮没有登录或修改服务器。

原仓库的 4 个分支、全部 **37 个历史/reflog 提交**及 37 个原先未挂接对象已经保存到新工程本地 `.git`；连通性检查通过。旧代码可用 `git show` 查询，无需保留旧运行目录。详见[历史凭据](cleanup_confirmed_20261023/history_preservation_receipt.json)和[Git 使用方式](GIT_WORKFLOW_ZH.md)。这些历史引用仅在本地，普通开发分支 push 不会上传它们；已有远端历史分支未改写。

四份旧发布 JSON 的原文件已删除，删除凭据保留各自 SHA。它们不被冒称为新清单的同一版本。

## 4. 删除后的检查

| 检查 | 实际结果 |
|---|---|
| 当前保护清单 | **3960 个文件全部 SHA 一致**，共 4,560,205,398 字节 |
| 删除前与原件再次比较 | 2759 个来自目标目录的保留文件一致 |
| 原有源码及原生文件 | 1207 个文件没有改变 |
| 小规模功能测试 | **74 通过、0 失败、0 跳过**，45.09 s；含已有 CPU/GPU 和固定时间窗小测试 |
| 独立性 | 127 个模块及真实数据读取通过，旧工作区访问尝试为 0 |
| 两套冻结流场、网格、WSS、checkpoint | 只读核对通过 |
| 正式旧微泡状态和账本 | 1565 接受，458/527/437 出流，141 观察截止、2 数值失败 |
| 当前原生求解程序 | 二进制身份检查通过；没有执行新的 CFD |
| 小鼠工程 | HEAD 和工作树状态与清理前一致 |

原正式场 SHA 为 `73228d094c253cec75dbb0030476b687171d0d0fef7816a845392dad220285dc`；校准场 SHA 为 `f5d7a2ccf7cdcc08cbfebd725d0a17e23e5780df23e4ef6488a26ffa5983d4d2`。

固定协议仍为灌注 1.240815660 s、随后观察 10 s、dt=0.5 ms、浓度 4×10⁹ m⁻³，尚无新正式批量结果。旧群体的 `scientific_integrity_pass=false` 保留。

这里的通过表示文件完整、程序能独立工作，**没有解决或取消原科学问题**：`REFERENCE_PIPE_FIELD_FAIL_SCIENTIFIC`、`MIXED_FIELD_ERROR`、`P1_FIELD_TRAJECTORY_IMPACT_LARGE`、`STOP_NATIVE_P2P1_ROUTE_RESOURCE_LIMIT` 继续保留，V3 仍仅提案。

证据：[验证汇总](cleanup_confirmed_20261023/validation_summary.json)、[文件核对](cleanup_confirmed_20261023/integrity_after_summary.json)、[独立性](../validation/confirmed_cleanup_independence/INDEPENDENCE.json)、[只读核查](../validation/confirmed_cleanup_readonly/SUMMARY.json)、[小测试日志](../validation/confirmed_cleanup_small_tests/pytest.log)。

## 5. 路径、其他项目与 Git

顶层 `/home/lzy/projects/brava_flow_roi_18mlmin` 快捷入口保留，现指向 `new_human_vessel_flow`。它只提供新工程入口，旧的内部目录结构不再兼容。当前开发应直接使用新目录 README 中的路径。

更新了共享算例目录中仅属 BraVa 的两行索引及图表链接，并更新 `server_migration/CURRENT_BRAVA_SERVER.json` 的本地位置。其它算例条目未修改。原始历史报告里的旧路径保留原文，不作为当前运行入口。差异见[导航更新](cleanup_confirmed_20261023/navigation_update_receipt.json)和[补丁](cleanup_confirmed_20261023/catalog_updates.patch)。

删除前全项目符号链接与 `.pth` 检查未发现其他运行依赖指向目标；72158 个范围内代码/入口文件检索只发现退休迁移工具和声学项目的历史**远端**路径，后者未修改。本次可访问的 22 个进程未使用删除目标，另有 28 个进程因权限/退出无法检查，没有把它们计为已检查。

混合历史 worktree 保留，其 Git 显示 2015 个 BraVa 路径删除，这是本轮授权删除的结果；该工作区 HEAD 和其他文件未变，没有顺手提交、重置或删除其他项目。Mouse、Clot、打印、声学、RBC 及共享 Git/环境继续保留。

源码与当前方法没有改动。修改内容是必要证据补齐、清单状态、忽略规则、路径导航和报告。临时整套工程、clone、worktree、ZIP 均未创建；详细本地审计和小测试证据保留，没有留下待运行任务。

本轮直接在 `dev/human-vessel-flow` 检查、提交和推送轻量文件及 SHA 索引；大场、原始微泡状态、历史二进制和旧仓库本地历史不随普通 push 上传。最终提交、远端一致性和推送状态写入本地 `docs/cleanup_confirmed_20261023/GIT_STATE.json`，不把自身 commit ID 写回同一个提交。

新增的当前说明和清单通过空白格式检查。原样保留的历史数据、日志、源码片段含原有行末空格/换行格式提示，逐项核对其 SHA 后保留原字节，没有为消除格式提示改写科学证据。详见[检查记录](cleanup_confirmed_20261023/whitespace_summary.json)。

未执行：服务器操作、新 CFD、新网格、新正式微泡群体、新轨迹认证，均为 **NOT_RUN**。
