# 小鼠项目清理报告

日期：2026-10-23。唯一日常开发目录：`/home/lzy/projects/new_mouse_vessel_flow`。

后续更新：用户已确认删除，原清单中的 15 项已实际删除。最新状态见[确认后删除报告](CONFIRMED_CLEANUP_REPORT_ZH.md)。下文空间统计和验证记录属于前一次清理，不能与后续删除重复计算。

**已实际删除重复同步工程，并把 Git 直接建在新小鼠工程内。** 已删除目标原占用约 **5.11 GB**；按本次发生变化的两个工作区计算，扣除补回正式数据、Git 和审计记录后，净减少约 **2.86 GB**。没有再次建立同步工程、备份镜像或压缩包。

## 1. 清理前后空间

以下 GB 为十进制，来自实际分配空间（`du -B1`），不是文件的逻辑大小。精确字节数和测量时间见 [disk_after.json](cleanup_20261023/disk_after.json)。

| 范围 | 清理前 GB | 清理后 GB | 变化 |
|---|---:|---:|---:|
| 新工程 + 本次小鼠同步副本，两个不重叠目录 | 9.391 | 6.527 | 减少 2.864 GB |
| `temp_storage/mouse_flow_github_sync_20261022`，也是本次冗余 Git 工作树 | 5.104 | 0 | 减少 5.104 GB |
| 本次清理的缓存 | 0.00313 | 0.00256 | 部分在复测时正常重建 |
| 新工程，含独立环境、补回的结果、报告和 Git | 4.287 | 6.527 | 增加 2.240 GB |
| 新工程自己的 `.git` | 0 | 0.017 | 新建独立 Git 管理 |
| 整块磁盘可用空间 | 661.206 | 672.031 | 观察到增加 10.825 GB，不能全归于本任务 |

表中明细互相包含，**不能相加**。整块磁盘在任务期间还有其他占用变化，因此本次净节省以两个明确工作区的前后实测为准。删除操作本身的短时间窗口，可用空间增加 5,115,146,240 字节；已核准目标分配空间为 5,107,482,624 字节，差额还涉及 Git 注销的管理文件及同期变化。

约 14.065 GB 的共享主 Git 库、约 2.675 GB 的旧几何 Git 库、约 3.627 GB 的旧共享环境仍保留。它们不是本次可以任意删除的小鼠独占副本。原先只测量了共享主 Git 库所在父目录，未单独测量其清理前对象库大小；不补造这个数值。

完整目录调查见 [磁盘清单](cleanup_20261023/DISK_INVENTORY.csv)和 [Git 对象库调查](cleanup_20261023/git_object_stores_summary.json)。这些原始/共享目录没有被算作“已清理”。

## 2. 实际删除了什么

- `/home/lzy/projects/temp_storage/mouse_flow_github_sync_20261022`：用 `git worktree remove` 正常注销并删除，未使用强制选项。已核对其远端提交 `987d0638944f09e1acc1f415b26bf5b0479c3f1d`，没有未提交、忽略但未保存的文件或未推送提交。
- 新工程及两个明确小鼠算例中的 Python、pytest、原生查询临时缓存。未清理共享 Python 环境。复测产生的少量缓存仍由 `.gitignore` 排除。
- 旧共享审计目录的 `MOUSE_LINKS.json`：与新工程保留的 `docs/audit/published_LINKS.json` 逐字节相同。
- 新工程内四个已完成任务的旧入口：`scripts/migrate_local.py`、`migrate_sources.py`、`finalize_audit.py`、`inventory_data_review.py`。原文可从首个提交 `e68fe08` 追溯。此次一次性调查/恢复/清理助手也已移除，只保留记录与当前可用检查脚本。

逐项路径、理由、原占用和保留位置见 [实际删除清单](cleanup_20261023/DELETED_FILES.csv)、[执行凭据](cleanup_20261023/deletion_receipt.json)。这些列表中的缓存目录项包含其内部文件，不额外重复计数。

对同步副本的 34,645 个纳入文件做了 SHA 核对：34,542 个原位置仍有相同内容；3 项是有记录的发布属性差别；100 项是打印项目以前已退休的代码副本，已逐一验证可以从保留的本地和远端历史取回。它们不是新小鼠运行依赖，也没有删除或修改独立打印工程。171 个已迁入源码的原位置没有出现遗漏的新修改。详见 [副本判断](cleanup_20261023/snapshot_resolution_summary.json)。

## 3. 为什么新工程还补回了数据

你在 2026-10-22 要求删除过 best-feasible、H0、原始 2 mm/s 和 Stage-Q 在新工程中的副本。这次附件明确要求保留它们的必要正式数据，因此本轮按最新指令执行。

从仍保留的原件补回 17,445 个必要数据/证据文件，逻辑大小 2,108,502,003 字节，源和目标 SHA 全部一致。没有恢复旧 Python、shell 执行工程，也没有重新计算结果。只有四份独立算例配置加入当前入口；计算方法仍共用一套源码。详见 [恢复凭据](cleanup_20261023/restoration_summary.json)。

## 4. 正式结果和实际验证

清理后重新检查 **31,400 个保留文件，4,863,923,028 字节，零哈希不一致**。五个保存流场、必要网格、原生 checkpoint、原始微泡记录和已保存验证证据均按现有清单核查。

| 流场 | 已核验 SHA256 |
|---|---|
| ROI-only | `fb3c6ad0815d156b09b289fc347766486f24ddb521b8e51e8bbc4cf1fd2ace12` |
| Best-feasible | `fcc692caa74c70d7ecbf9ae6662d29d45abbae926b04783aa52886d382a4bdb2` |
| Network-H0 | `064cbd28f3efa72f426fc946b2f29da21f056c596609095e7283d39070aa55f4` |
| 原始 2 mm/s | `129ebb77396550b300b20ce29cde7b03919e6037b7a81de60fe75d6c6efb616d` |
| Stage-Q | `373ff549430cab710d57eb4406f2e7588f7a5bb68ebaa68cf32f8da5662e26f1` |

| 微泡数据 | 数量 | 时间步 | O1 / O2 / O3 出流 | 接触支持静止 | 执行失败 |
|---|---:|---:|---:|---:|---:|
| ROI-only | 1500 | 0.5 ms | 787 / 61 / 467 | 185 | 0 |
| Best-feasible | 1500 | 0.5 ms | 0 / 133 / 1182 | 185 | 0 |
| Network-H0 | 500 | 1 ms | 0 / 139 / 291 | 70 | 0 |

上述全部 3500 条记录的队列顺序、粒径、出生位置、流场身份、时间步、原始轨迹、完成标记和实际分类都已检查。静止不算出流。没有把历史数据冒称为新计算。

- 清理前后均运行现有小测试；最终 **116 通过，0 失败，0 跳过**。
- 禁止访问旧 `computation_examples`、`temp_storage`、打印、BraVa 等目录后，**190 个模块导入通过**；旧目录访问尝试为 0。
- 在同一禁止旧路径的进程内，几何读取、小网络回归、五个场、三组轨迹和原 WSS 恢复核对通过。
- 原生 CFD/网格安装身份检查通过，读取以前实际小规模试算的证据；本轮没有重跑原生求解器。
- 新几何重建、网格生成、完整 CFD、正式微泡计算：**NOT_RUN**。细档血管与 O2 敏感性仍为 **SKIPPED_BY_USER**。旧 Ultraliser 特殊极小模型问题按你之前要求暂缓处理，未声称修复。

结论仅适用于这些实际检查，不增加完整血管局部物理精度或新算例成功的保证。参见 [验证汇总](../validation/cleanup_20261023/SUMMARY.json)、[独立性验证](../validation/cleanup_20261023/independence.json)、[最终数据核查](../validation/cleanup_20261023/integrity_after_removal.json)。

## 5. 原待确认事项的后续结果

见[已更新的原待确认清单](cleanup_20261023/PENDING_REVIEW.csv)。用户随后授权删除，处理结果如下：

- `computation_examples/mean-2p0-mmps*` 和 `fem_reference_stage_q` 五个原始工程：已删除旧目录，新工程所需数据保留且核验通过。
- `_shared/active_network_development`、`particle_workflow`、`geometry_network`、旧流场和 WSS 工作区：含共享关系、原始资料或未提交状态，不能整体视作纯临时副本。
- `temp_storage/ulm-3d-vascular-model-generation/.git`：仍支撑多个独立项目；`ulm_particle_3d_particle0/.venv` 仍有旧项目消费者。
- `_archives/github_sync_20260928_latest`：混有 BraVa 内容，且相对已缓存远端有 4 个本地提交，本轮不删除。
- 2026-09-27 的六个迁移/清理目录、四份旧发布清单：已按后续授权删除，包括不再需要的独有历史整理记录。

原清单标为“保留”的 18 项继续保留。14 个指向已删除小鼠目录的兼容链接改指向新工程中对应数据；共享环境、混合项目、正式保留数据和服务器文件没有删除。

## 6. Git 和以后怎样开发

- Git 根目录就是 `/home/lzy/projects/new_mouse_vessel_flow`。
- 当前分支：`dev/mouse-vessel-flow`；upstream：`origin/dev/mouse-vessel-flow`。
- 远端：`https://github.com/OZAKI39/ulm-3d-vascular-model-generation.git`。
- 已从本目录正常 commit 并 push 成功，旧 `main`、历史同步分支和其他项目分支未覆盖；没有 force push。
- `.git` 独立位于本目录，无外部对象库链接；这里只有一个工作树。旧来源以 Git 提交和文件身份对应表追溯。
- 根目录已建立 `AGENTS.md`；`.gitignore` 排除环境、缓存、构建产物和大批科研数据。代码、配置、测试、小输入、身份和精选报告进入 Git。没有改变已有 LFS 规则。
- 共享旧导出、发布、封装工具的**小鼠路径**已禁用；两项小鼠命令实测在复制/推送前退出。四个共享脚本的 BraVa 执行路径未改动；完整差异在本地审计记录中。没有发现小鼠同步定时任务。

最新提交与远端一致性实测记录在本地 `docs/cleanup_20261023/GIT_STATE.json`，该文件在最后推送后生成，不纳入自身描述的提交。使用 `git log -1` 也能直接查看。

当前新增文件中有原封保留的上游源码、原始 CSV 与历史补丁，所以整树 `git diff --cached --check` 曾报告已有空白格式问题；没有为通过格式检查改变科学/上游证据。700 余份原生文件已核对原 SHA；本轮修改的当前文件单独通过检查。详见 [格式检查说明](cleanup_20261023/inherited_whitespace_summary.json)。

日常操作见 [Git 使用说明](GIT_WORKFLOW_ZH.md)。**clone 这个分支只获取所管理的代码和轻量材料，不自动下载本地全部科研数据**；完整结果的迁移以 `configs/active_data_manifest.json` 为准。

## 7. 最终目录与独立性结论

[真实目录树](cleanup_20261023/DIRECTORY_TREE.txt)；[当前方法清单](CURRENT_METHODS_ZH.md)；[完整函数索引](CURRENT_SOURCE_FUNCTIONS.csv)。主要入口保持为：

```text
new_mouse_vessel_flow/
  .git/                 直接管理当前文件
  AGENTS.md             后续任务规则
  mouse_vessel_flow/    唯一活跃 Python 实现
  cases/                五个独立算例配置
  configs/              当前环境、源码和数据身份
  inputs/               正式输入与受保护证据
  results/              各场、正式轨迹、原图表
  external/             原生来源和安装说明
  scripts/              当前小检查、独立性检查
  tests/                当前可执行测试
  validation/           历史证据与本轮实测记录
  docs/                 来源、清理与 Git 说明
```

**已删除旧小鼠同步目录后，新工程仍可独立修改代码、使用 Git、加载保存结果及执行本轮已验证功能。** 这已经通过实际删除后的阻断旧路径检查。未来完整新计算仍需要已安装的 `/home/lzy/projects/new_vessel_flow_solvers` 和正常系统库；这套明确的软件安装不是旧同步工程。没有声称未运行的大计算已通过。
