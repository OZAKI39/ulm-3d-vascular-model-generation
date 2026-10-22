# BraVa 工作区清理与直接 Git 开发报告

日期：2026-10-23。唯一日常开发目录：`/home/lzy/projects/new_human_vessel_flow`。

**已实际删除两处重复同步工作区，并把 Git 直接建在新 Human 工程里。** 被删除内容原占用约 **7.871 GB**；把新 Git、补齐证据和验证文件计入后，三个明确变化目录实测净减少约 **7.793 GB**。没有建立新的同步工程或完整备份副本。

## 1. 真实空间变化

GB 均为十进制；按 `du -B1` 的分配空间计算，不把逻辑文件大小当作释放空间。测量时间和精确字节见 [清理前](cleanup_20261023/disk_before.json)、[清理后](cleanup_20261023/disk_after.json)。

| 不重叠的变化目录 | 清理前 GB | 清理后 GB |
|---|---:|---:|
| 新 Human 工程，含独立环境、数据和 Git | 10.950 | 11.026 |
| `temp_storage/brava_flow_github_sync_20261022` | 5.927 | 0 |
| `_archives/github_sync_brava_20261006` | 1.941 | 0 |
| 三处合计 | 18.819 | 11.026 |

原 BraVa 工程仍占约 25.566 GB，独立校准工况约 0.277 GB，未计入已清理空间。旧共享 Git 根目录约 14.065 GB、共享环境约 3.627 GB，继续供其他工作区使用，未删除。新 Human `.git` 在本次测量时约 0.055 GB。

整块磁盘可用空间从 676.517 GB 到 684.315 GB。磁盘还可能有其他活动，所以本任务净变化按上面三个不重叠目录计算。库存中的父子目录不能再次相加；后续提交会使 Git 元数据有少量变化。

完整根目录盘点见 [DISK_INVENTORY.csv](cleanup_20261023/DISK_INVENTORY.csv)；71 个原始子目录的占用和用途保存在本地 `docs/cleanup_20261023/original_directory_inventory.csv`。

## 2. 实际删除与判断依据

- **2026-10-22 发布工作区**：8,749 个发布来源重新做 SHA 比较，8,745 个仍与原来源一致；另外两项是发布属性差别，两项是原位置不再存在的历史说明。历史说明和属性有已验证 Git 恢复途径。工作区另有 22 份发布元数据/入口文件，已核对对应 Git 对象。
- **2026-10-06 早期工作区**：2,880 个普通文件逐一对照。57 项没有直接对应的当前同内容文件；其中 39 项必要验证资料先迁入新工程，其余 18 项是旧发布资料或继承的历史表格/说明，Git 可追溯。两份小鼠相关表格另核对了仍保留的原件，没有混入 Human 数据或改变 Mouse 工程。
- 两个工作区删除前均无未提交、未追踪或忽略文件；HEAD 与实际远端分支一致。用 **`git worktree remove` 正常注销**，没有强制删除、reset、clean 或清理共享 Git 对象。
- 新工程内 8 个已完成的一次性迁移工具已移除，明确的 Python/pytest 缓存已清理；复测时重建的少量缓存继续由 Git 忽略。工具原文保存在初始提交 `4a2d4e7`，不再保留可重放迁移的活动入口。临时 Git 暂存列表也已删除。

全部目标、原占用、保留位置与证据见[实际删除清单](cleanup_20261023/DELETED_FILES.csv)和[执行凭据](cleanup_20261023/deletion_receipt.json)。删除前没有发现可访问进程使用两处同步目录；27 个受权限限制的进程无法检查，未计为已检查。没有发现指向这两处目录的兼容链接或当前运行源码依赖。

## 3. 补齐了哪些必要资料

共补齐 **43 个文件、7,920,074 字节**，每个均在复制前后核对 SHA。其中 39 个来自早期同步工作区：参考管第 10 步场与 checkpoint、P2 资源审查的边界输入/求解配置、曲面失败检查的实际执行记录等；另 4 个是原 Gmsh V2 生成器及其少量冻结前置输入。

这些文件保存到已有的 `diagnostics/reference_pipe`、`diagnostics/p2p1`、`diagnostics/nearwall` 和 `results/reports` 中，没有复制整套历史工程。补齐的第 10 步是历史验证材料，GPU V2 的正式参考仍为原第 71 步。

原 `reference_v2.py` 保持原 SHA，作为冻结来源证据；当前读取和分析仍通过 `src/human_vessel_flow/reference_pipe/` 的普通模块。没有安装新 mesher 或执行生成器。完整来源和 SHA 见[保留凭据](cleanup_20261023/preservation_receipt.json)。

第一次通用审计把迁移清单已经声明移除的 22 个旧控制器误判为缺失；按原清单的 `PRUNED_UNUSED_CONTROLLER_IN_NEW_PROJECT_ONLY` 分类核对后排除该误判。另有旧根目录 handoff 文件已不存在，但新工程保留了原文。**原有 3,622 个交付文件在本轮开始时全部与交付 SHA 一致**，没有通过修改科研数据来消除报错。见[分类处理记录](cleanup_20261023/identity_resolution.json)。

## 4. 当前源码与长期使用方式

权威运行实现仍在 `src/`，原科学实现没有改动。本轮只改进了独立性检查脚本：输出到新目录、拒绝覆盖旧记录，并按真实路径检查文件读取、目录列举和扫描。已创建简洁的 [AGENTS.md](../AGENTS.md)，明确禁止普通 push 创建同步副本。

[当前方法清单](CURRENT_METHODS_ZH.md)逐项列出几何、TetGen 网格、中心线压力、原生 FEM、冻结场、WSS、微泡、固定时间窗、P1/P2 研究和可视化的实际源码位置。[目录树](cleanup_20261023/DIRECTORY_TREE.txt)来自实际文件系统。

生产原生源码在 `native/svmultiphysics`，GPU 同步补丁在 `native/production_patch`；当前独立安装由 `configs/cfd/native.json` 指定。P2 实验补丁保持诊断身份，没有并入生产求解器。共享 native 安装、其他项目环境和原服务器目录均未清理。

## 5. 删除后的数据与功能检查

- **74 项现有小测试通过，0 失败、0 跳过**，实际包含原 CPU/GPU 人工场景、固定时间窗和接口测试。
- 当前保留清单 **3,332 个文件、4,321,212,994 字节**，全部 SHA 一致。清单包括正式输入/结果、必要验证证据及受保护原生源码等不同类别，不能把所有文件都称为轨迹文件。
- **1,201 个原有源码及原生文件**与前次交付身份一致；计算方法和冻结参数没有修改。
- 两处同步目录实际删除后，再阻断旧工作区访问，**127 个模块导入和真实数据读取通过**，旧目录访问尝试为 0。
- 两套冻结场、原网格和端口、速度压力、原始面片 WSS、流量记录与两套必要最终 checkpoint 的只读复核通过。
- 旧群体完整状态及账本复核通过：**1565 接受，458 / 527 / 437 出流，141 观察截止，2 数值失败**，原 `scientific_integrity_pass=false` 保留。
- 固定协议仍为灌注 **1.240815660 s**、随后观察 **10 s**、dt **0.5 ms**、浓度 **4×10⁹ m⁻³**。只有协议和小测试，没有新正式大批量结果。
- 原生可执行文件、PETSc 库和 wrapper 的身份检查通过；本轮没有发起原生 CFD。以前的构建/加载和 2,146 点真实场 CPU/GPU 比较记录保留，没有将其冒称为本轮重跑。

正式场 SHA：`73228d094c253cec75dbb0030476b687171d0d0fef7816a845392dad220285dc`。

校准场 SHA：`f5d7a2ccf7cdcc08cbfebd725d0a17e23e5780df23e4ef6488a26ffa5983d4d2`。

科学状态不变：`REFERENCE_PIPE_FIELD_FAIL_SCIENTIFIC`、`MIXED_FIELD_ERROR`、`P1_FIELD_TRAJECTORY_IMPACT_LARGE`，以及 `STOP_NATIVE_P2P1_ROUTE_RESOURCE_LIMIT`。V3 仍仅提案，曲面压力/高阶输出问题未被宣布解决。没有新增 CFD、网格、正式微泡群体、轨迹重放或时间步收敛验证。

证据：[验证汇总](cleanup_20261023/validation_summary.json)、[独立性](../validation/cleanup_independence_after/INDEPENDENCE.json)、[只读核查](../validation/cleanup_readonly_after/SUMMARY.json)、[小测试日志](../validation/cleanup_small_tests/pytest.log)。

## 6. Git 与旧同步入口

Git 根目录就是新 Human 工程，分支为 **`dev/human-vessel-flow`**，远端为 `OZAKI39/ulm-3d-vascular-model-generation`，upstream 使用 `origin/dev/human-vessel-flow`。本地 `.git` 独立，没有外部 alternates、commondir 或第二工作树。

初始独立根提交 `4a2d4e7` 保存已重组的当前工程；旧 `sync/brava-flow-project-20261022` 和 `sync/brava-p2-algebra-20261006` 远端历史均保留，没有改写。此结构通过独立开发分支管理完整项目，不需要新建仓库或向 main 塞入多套工作区。

6 个旧同步入口已改成启动即说明退休并退出的小文件，不再复制、发布或调用旧服务器。逐个实际调用确认在产生副作用前退出；原差异见[入口处理记录](cleanup_20261023/retired_entrypoints.json)。其他项目共享的旧历史工具没有擅自执行。没有发现当前用户 crontab；这一检查不代表扫描了所有系统级服务。

Git 管理源码、配置、必要轻量资料、补丁与报告，大批网格/场/原始状态/动画及环境被 `.gitignore` 排除。本地大数据没有因 Git 排除而删除，当前 SHA 索引为 `configs/active_data_manifest.json`。没有修改 LFS 或拆分大文件。

初始导入含原样保存的上游源码、历史补丁和原始文本数据，其已有空白格式问题被单独记录，没有为格式检查改写这些身份文件。本轮新修改的当前文件单独执行 diff 检查。提交和推送的最终实测状态在本地 `docs/cleanup_20261023/GIT_STATE.json`，不纳入自身描述的提交。具体操作见 [Git 说明](GIT_WORKFLOW_ZH.md)。

## 7. 保留并待确认的范围

[最新待确认清单](cleanup_20261023/PENDING_REVIEW.csv)列出 **9 项**：原 BraVa 工程、独立校准工况、两处服务器迁移资料、混合历史工作区中的 BraVa 部分，以及四份独有旧发布清单。原始工程和备份按本轮指令等待确认，没有整目录删除。

其中旧工程约 25.566 GB，旧服务器迁移备份约 46.374 GB。旧 mesh/field 研究目录约 13.099 GB；停止的 `realization_0000`、`guardfix`、`final` 原状态分别约 0.756 / 0.898 / 2.126 GB。这些数字包含于父目录，**不能相加推算可释放空间**。它们含原始状态、旧试验或独有失败证据，未因“停止”或“大文件”自动删除。

Mouse、Clot、打印、RBC、声学独立工程未改动；Mouse HEAD 和工作树状态在前后核对一致。共享 Git、Python、native 依赖继续保留。没有连接服务器，也没有依靠远端副本作为唯一删除理由。

**在已删除的两处同步工作区不存在时，新 Human 工程仍可独立开发、读取冻结场、分析旧微泡数据并运行本轮通过的小型功能。** 新的原生 CFD 数值对照、完整 GPU 群体和科学精度认证仍是未执行或未通过的独立事项。
