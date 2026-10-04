# CLEANUP_RECOVERY_AUDIT_V1

审计对象：`/home/lzy/projects/temp_storage/clot_clearance_github_sync_20260930`。

**结论：安全恢复与依赖分类已完成，但恢复验收未通过；停在 `CRITICAL_UNRECOVERED`，等待人工审核。** 原 123 项中，51 项已按原 SHA256 恢复，72 项仍缺原字节。12 份科学字段等价的 SUMMARY 候选仅保存在暂存区，不覆盖历史路径，也不替代旧哈希断言。

用户已明确回复 **“没有”清理前备份**。本轮不继续 repository cleanup、activation audit 或科学开发；没有执行文件系统级恢复。

## 1. 最终必报指标

| 指标 | 结果 |
|---|---|
| missing original | **123** |
| exact-byte recovered | **51**：本轮开始前 48，本轮新增 3 |
| semantic-exact recovered | **0** |
| scientifically-equivalent regenerated | **12**：从原历史记录重组科学字段，仅暂存 |
| manifest-only no-recovery-needed | **0** |
| temporary/cache no-recovery-needed | **0** |
| critical unrecovered | **72**：真实的旧测试字节依赖 |
| unknown | **0** |
| tests：本轮 before / after | **407 passed / 8 failed / 0 skipped → 407 / 8 / 0** |
| tests：清理前基线 | **412 passed / 3 failed / 0 skipped** |
| C5 state loadable | **yes** |
| C7 state loadable | **yes** |
| old manifests modified | **NO**：141 份原清单哈希未变 |
| production modified | **NO** |
| old test assertions modified | **NO** |
| scientific baseline intact | **YES，限当前生产源码与接受的 C5/C7 科学基线** |
| full pre-cleanup byte integrity restored | **NO** |
| pre-deletion snapshot found | **NO**，限本机可见且获授权的检查范围 |
| post-deletion read-only backup created | **YES** |
| 是否需要用户提供 backup 才能继续精确恢复 | **YES；用户已确认没有，不重复索要** |
| filesystem recovery | **FILESYSTEM_RECOVERY_MAY_BE_REQUIRED** |

计数口径：原文件的互斥处置为 **51 已精确恢复 + 72 未精确恢复 = 123**。12 个科学等价候选属于这 72 项的子集，因此不可再加到 123 的分区总数中。它们缺失原运行 telemetry，且旧测试比较整文件 SHA；既不是全容器 `SEMANTIC_EXACT`，也没有消除测试依赖。另有 15 个不完整 SUMMARY 候选判为 `NOT_EQUIVALENT`，不计作恢复成功。

上次已恢复的 5,871 项不在 123 项统计分母内。加上本次分母内的 51 项，误删测试依赖累计精确恢复 **5,922 项**。

## 2. 冻结状态与备份

- Git HEAD：`f0e891686b9b107510ab3292aa27ce83638b3339`。
- 分支：`clot-clearance/energy-regularization-20260930`。
- 初始 HEAD/status：[FROZEN_GIT_STATE.txt](FROZEN_GIT_STATE.txt)；结束状态另存 `FINAL_GIT_STATE.txt`。
- 副本：`/home/lzy/projects/temp_storage/clot_cleanup_recovery_snapshot_20261004T203008Z`。
- 24,740 个普通文件、143 个符号链接、3,906,563,979 bytes；`cp -a --reflink=auto --no-dereference`，无与原仓库共享的硬链接，普通文件和目录去掉写权限。

备份工作区文件逐项 SHA 校验无差异。复制与校验之间，原仓库 Git 的 loose objects/refs 转为 packed 存储，导致原始逐文件检查列出 985 个 `.git` 存储差异；差异没有涉及工作区科学文件。未查明并发 packing 的触发者，本轮没有运行 `git gc/repack`。两份仓库 HEAD、status 相同，`git fsck --full --unreachable` 均退出 0。原始失败检查和后续逻辑验证同时保留在 [BACKUP_VERIFICATION.json](BACKUP_VERIFICATION.json)，不把底层存储相等与 Git 对象完整性混为一谈。

该副本是**误删后**的冻结备份，只能防止后续变化，无法提供当时已缺失的原文件；权限只读也不等于文件系统不可变快照。

## 3. 123 项依赖判定

[RECOVERY_MISSING_123.csv](RECOVERY_MISSING_123.csv) 包含全部要求字段，并增加 `dependency_class`、`secondary_class`、`git_ever_tracked`、`semantic_status`。原分母取自只读保留的 `../repository_cleanup_v1/RECOVERY_UNAVAILABLE.json`，未改写为 72 项。

123 项均为 `DIRECT_TEST_DEPENDENCY`：测试读取 PRE manifest 后，逐路径执行 `read_bytes()` 并验证 SHA256。间接经过 manifest 并不等于 `MANIFEST_ONLY`。

| 次级文件性质 | 原 123 项 | 现仍缺失 |
|---|---:|---:|
| REGENERABLE_CACHE：Numba / Matplotlib | 45 | 45 |
| TEMPORARY_ARTIFACT：历史 pytest 输出 | 78 | 27 |
| 合计 | 123 | 72 |

仍缺失文件为 **41 `.nbc`、3 `.nbi`、27 `SUMMARY.json`、1 font cache JSON**，合计 **8,644,268 bytes**。逐文件列表及预期哈希见 [CRITICAL_UNRECOVERED.csv](CRITICAL_UNRECOVERED.csv)。

未发现这些路径被正式科学 REVIEW/STATUS/SUMMARY 或运行时代码按字面路径读取；Numba 对可选缓存的自动读取单独记录。运行时可以再生缓存，但无法由此证明再生字节与历史编译缓存相同。没有任何一项被判为当前 C5/C7 接受态或必需的生产运行时输入。

证据：[DEPENDENCY_EVIDENCE.json](DEPENDENCY_EVIDENCE.json) 记录扫描范围、逐测试读取链、manifest 哈希和语义重建来源。字面搜索不能排除一切动态访问；本次没有以“搜不到文件名”作为删除许可。

关键读取链：

| 旧测试 | 实际保护范围 |
|---|---|
| `test_spatial_discretization_v2.py::test_legacy_files_are_unchanged` | 整个旧 branch/spatial V1 树，包括测试临时输出和缓存 |
| `test_spatial_discretization_v2.py::test_historical_sources_and_tests_unchanged` | `pd_clot/` 前缀下的编译缓存也被读取 |
| `test_pfpd_spectral_domain.py::test_all_preexisting_production_and_candidate_sources_unchanged` | 同上 |
| `test_c5_stability_bracket.py::test_saved_roots_production_and_existing_tests_unchanged` | 同上 |
| `test_c7_early_branch.py::test_production_and_historical_geometry_unchanged` | 同上 |

## 4. 安全恢复与搜索范围

按原 SHA 已恢复的 51 项包括 27 个 gzip 序列化文件、21 个恢复历史 fixture 路径前缀的 JSON，以及本轮新增 3 个 `IDENTITY.json`。新增 3 项用精确保留的 CONFIG 与 source_snapshot 恢复配置哈希和源码哈希，写入 staging、复验原 SHA 和长度，再复制原路径并复验。没有用新运行结果或随意数值填补旧哈希。

本轮新增恢复日志：[EXACT_RECOVERY_LOG.jsonl](EXACT_RECOVERY_LOG.jsonl)。此前 48 项的恢复证据仍在 `../repository_cleanup_v1/SERIALIZED_ARTIFACT_RESTORATION_LOG.csv`。

只读检索结果：

- 遍历 `/home/lzy` 与 `/tmp`，覆盖 `/home/lzy/projects`、`temp_storage` 及实际 sibling 目录；元数据检查 782,544 个文件。
- 对现仓库及 `temp_storage` 实体 sibling 目录建立全文件 SHA 索引；其余 home/tmp 文件按预期大小筛选后计算 SHA。共哈希 **55,626 个文件、13,949,878,489 bytes**，不依赖文件名相同。[FILE_SHA256_INDEX.csv](FILE_SHA256_INDEX.csv)
- 检查 **43 个 Git 对象库**；现仓库 branches/tags/reflog/all commits、1,050 个可达对象记录、1,052 个本地对象元数据及 fsck unreachable 均已检查。两个 unreachable blobs 也包含在 all-objects 扫描内。
- 123 项当前均未 tracked；现有所有分支及 reflog 可见历史中也未 tracked。4 个 partial clone 对象库有 promisor 警告，本轮未联网拉取缺失远端对象；本地 Git 检索不代表证明所有远端对象都不存在。[GIT_OBJECT_SEARCH.json](GIT_OBJECT_SEARCH.json)
- 只读逐成员检查 **166 个 ZIP/TAR 归档**，未解包到原目录，未找到缺失 SHA。[ARCHIVE_SEARCH.json](ARCHIVE_SEARCH.json)
- 对 3 个 Numba index 检查 112 个保留索引的文件编号排列候选，无原 SHA 匹配；没有编译新内核。[CACHE_INDEX_RECONSTRUCTION.json](CACHE_INDEX_RECONSTRUCTION.json)

搜索限制：不跟随目录符号链接、不访问 `.ssh`/`.gnupg`、跳过本轮输出和新冻结副本的重复扫描；6 个 `/tmp` 系统私有目录无读取权限，未绕过。冻结副本另有完整哈希验证。没有发现更早备份；不能由此断言未挂载磁盘或宿主系统不存在备份。

## 5. 容器语义恢复

本 123 项没有缺失 `.npz/.npy` 文件；数组容器名称、shape、dtype 的检查用于 C5/C7 加载 smoke，未启动模型求解。

[SEMANTIC_RECOVERY_COMPARISON.json](SEMANTIC_RECOVERY_COMPARISON.json) 保存每个 JSON 候选的输入文件原 SHA、字段比较、丢失字段及输出位置。

- **12 个 `SCIENTIFIC_EQUIVALENT`**：ADR SUMMARY 的科学字段全部由同次历史 `history.json`、`damage_jumps.json`、`equilibrium_*/SUMMARY.json`、CONFIG 与原 writer 的纯字典表达式重组。输入均与清理前 SHA 一致；另与已存在的测试输出科学字段比较相同。未调用 mechanics、damage 或 solver。原 `wall_time_s`、`cpu_time_s`、`peak_process_rss_bytes` 无法还原，候选明确省略，不填写假值。
- **15 个 `NOT_EQUIVALENT`**：旧 streaming SUMMARY 可恢复 final/history 部分，但全程 `minimum_J_or_intrinsic_stretch`、`maximum_acceleration_m_s2` 每个机械步更新，现存端点和 settling samples 不覆盖所有步骤。无法验证原全程极值，因此不把相似新运行认作历史等价。
- **0 个全容器 `SEMANTIC_EXACT`**。缓存二进制未获得可验证原语义对象；Matplotlib 缓存也没有足够原证据。

所有非精确候选仅在 `semantic_recovery_staging/`。没有因为 SHA 不同而直接判科学字段失败，也没有因为科学字段相同就宣称旧字节完整性测试通过。

## 6. 不可恢复项与异常记录

[RECOVERY_EXCEPTION_REPORT.json](RECOVERY_EXCEPTION_REPORT.json) 记录损失原因、原路径/哈希、当前科学依赖、真实测试依赖、候选替代物和语义状态。其状态为 `NOT_A_WAIVER`，没有授权忽略任何现存测试依赖。

`CRITICAL_UNRECOVERED=72` 的含义是 **旧测试所需的历史字节仍未恢复**，不是 72 个当前科学接受态损坏。即使其中 12 项有科学字段等价的副本，仍无法通过原 SHA 断言，所以保守地保留在阻塞列表中。

旧科学数据没有被改写，但 15 个旧 pytest SUMMARY 的独有全程极值目前无法重建；不能将“当前 C5/C7 科学基线完整”扩大表述成“所有历史临时结果信息完整”。

## 7. 新保护基线提案

在 123 项全部分类、unknown=0 后生成 [SCIENTIFIC_PROTECTED_MANIFEST.json](SCIENTIFIC_PROTECTED_MANIFEST.json)，共保护 **7,970 文件、1,147,813,239 bytes**。包括源码、测试、配置、接受/诊断状态、最终科学数据、REVIEW/STATUS、图件和关键 provenance。

依据同一策略反查清理前库存，没有遗漏策略选中的旧科学文件，已存文件相对清理前无哈希变化。438 项原有 source/test/config 及 141 份旧 manifest 均通过哈希检查，详见 [PROTECTED_BASELINE_VERIFICATION.json](PROTECTED_BASELINE_VERIFICATION.json)。

此清单为 **`PROPOSED_BASELINE_REQUIRING_HUMAN_REVIEW`**，不替代旧 manifest，不消除 72 项旧依赖。临时/永久边界见 [TRANSIENT_ARTIFACT_POLICY.md](TRANSIENT_ARTIFACT_POLICY.md)。供现有测试和报告使用的 C5 Hessian `cache/` 继续保护；不凭目录名称判断可删除性。

## 8. 测试与接受态检查

本轮运行未修改的 `pytest tests -q`，附加 `-p no:cacheprovider` 和独立 `/tmp` basetemp；Python/Numba/Matplotlib 运行缓存重定向到本轮临时目录，不写入旧科学结果。未启动独立科学 simulation。

| 时点 | passed | failed | skipped |
|---|---:|---:|---:|
| 清理前原基线 | 412 | 3 | 0 |
| 本轮开始（复用上轮最后实测记录） | 407 | 8 | 0 |
| 本轮恢复后实测 | 407 | 8 | 0 |

三份记录收集到的 415 个测试 ID 完全相同。本轮 before/after 的 8 个失败 ID、消息、traceback **逐文本相同**；原有 staggered 3 个失败同样保持原消息与 traceback。**本轮没有新回归，但上次清理引入的 5 项失败仍未消除；“恢复至清理前无新增失败”的验收不通过。**

证据：[TEST_COMPARISON.json](TEST_COMPARISON.json)、[AFTER_TESTS.json](AFTER_TESTS.json)、[AFTER.log](AFTER.log)。

C5/C7 smoke 读取 **52 个状态 + 6 份报告**，shape/dtype/有限性检查通过，文件 SHA 与清理前一致。包含 C5 15 个稳定性审计态、C5 激活相关 6 态、C7 接受分支及激活相关 31 态。[SCIENTIFIC_SMOKE_AFTER.json](SCIENTIFIC_SMOKE_AFTER.json)

## 9. Snapshot 与文件系统恢复边界

只读 `findmnt`/`lsblk -f` 显示工作区所在卷是 `/dev/sdd` 上的 ext4，挂载含 `discard`。未发现可见 ZFS/Btrfs/LVM snapshot，也未发现可用 restic/borg/rsnapshot 工具或已检查路径中的备份。

检查路径：`/.snapshots`、`/snapshots`、`/mnt/snapshots`、`/backup`、`/backups`、`/home/.snapshot`、`/home/lzy/.snapshots`、`/home/lzy/backup`、`/home/lzy/backups`；均不存在。用户已确认没有备份。宿主机或未挂载介质未检查，无法对其下结论。

**FILESYSTEM_RECOVERY_MAY_BE_REQUIRED**。若仍要争取恢复旧测试所需原字节，建议停止向该 ext4 文件系统大量写入，再由人工决定后续磁盘恢复方案；现有信息不足以保证可恢复。

本轮未运行 extundelete、photorec、testdisk、debugfs 恢复/写操作，也未执行文件系统级 restore。没有修改 physics、production、旧测试、旧 manifest 或旧科学结果来制造通过。

## 10. 停止点

恢复报告可供审核，恢复验收仍被 72 项旧测试字节依赖阻塞。没有清理前备份可继续搜索。**停止本轮工作，等待人工审核；不继续 cleanup、activation audit、simulation 或 solver 开发。**
