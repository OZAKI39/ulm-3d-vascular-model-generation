# Clot 历史目录与 Git 审计

审计日期：2026-10-20；执行日期：**2026-10-21**。用户已明确批准 **A+B+C** 及固定分支方案。**A、B 已完成删除，旧项目目录清理前占用合计 8.594 GiB**；两份独有说明已保存并提交，C 也已按精确清单全部删除。 第一批未访问或清理服务器 `/workspace/y_clot/`，没有重跑正式 3.6 ms 轨迹。唯一日常开发工程是本目录；保留的旧目录仅用于科研数据和历史证据。

阅读说明：以下“没有 push / 没有 upstream”是第一批清理执行时的历史状态。首次正式提交 `982a0d1` 已正常推送至长期分支 `clot-simulation-v4` 并设置 upstream。随后用户明确批准第二批；执行结果见本报告末节，未清理其他项目或服务器。

## 空间与结论

按 `du -x -B1` 统计实际分配空间，不跟随符号链接。1 GiB = 2³⁰ bytes。开始审计时磁盘可用 **445.71 GiB**，`temp_storage` 占 **79.76 GiB**；生成临时只读审计索引后约 **445.56 GiB**。系统及其他项目活动可能使数字变化。

- A：已删除干净 worktree **4.880 GiB**，共同 Git 元数据也由 Git 正常注销。
- B：已保留两份独有说明并删除恢复快照 **3.714 GiB**。
- 其余需判断用途的环境、文献和审计目录 **3.979 GiB**，本轮不建议删除。
- 必须暂时整体保留的工程/科研数据目录 **58.864 GiB**，包括新工程及下表四个旧 Git 工程和本地 `clot_clearance`；这是目录占用，不能解释为全部都是唯一数据。
- A+B 删除前 `df` 可用 478,415,798,272 bytes，删除后立即为 487,647,141,888 bytes，增加 **9,231,343,616 bytes（8.597 GiB）**，包含额外注销的 worktree 元数据。本轮临时索引约 0.145 GiB 的回收单列，不计入旧项目成果。最终磁盘读数见本报告末尾；测试和其他进程可能同时产生文件。

## 已批准并执行的精确删除清单

| 编号 / 原始路径 | 实际占用 | 文件用途 | 判断依据 | 建议处理 |
|---|---:|---|---|---|
| A `/home/lzy/projects/temp_storage/clot_github_sync_20261020` | 5,239,787,520 bytes / 4.880 GiB | 已推送的精选源码、科研结果快照 worktree | `.git` 是文件；工作树干净；HEAD 与实时远端一致；未跟踪文件为零，仅两份忽略的 `.pyc`。全部其余普通文件的 Git blob 在核对的远端 Clot 分支中存在，符号链接受 Git 跟踪 | 已从共同主仓库执行 `git worktree remove`，未使用 force；原分支和共同 `.git` 保留，未 prune 对象 |
| B `/home/lzy/projects/temp_storage/clot_cleanup_recovery_snapshot_20261004T203008Z` | 3,987,644,416 bytes / 3.714 GiB | 10 月 4 日只读恢复副本，含独立 `.git` 和 11 个嵌套上游仓库 | 其中 `clot_clearance/` 的 **23,171 个普通文件**与保留的 9 月 30 日工程同相对路径、同 SHA-256；符号链接文字完全相同。12 个 Git 对象库的全部对象（含未被 refs 引用的对象）及全部 refs 在对应保留仓库中存在 | 已保存下列两份说明并验证 SHA、提交后删除整个 B；README 已在远端历史分支，未另复制工程 |

B 的必要前置步骤已经完成，保留提交为 `0aea532`：

| 独有文件 | 文件大小 | 现保留位置 |
|---|---:|---|
| `B/PROJECT_ACTIVE_CONTEXT.md` | 7,773 bytes | `archive/reports/cleanup_20261004_PROJECT_ACTIVE_CONTEXT.md` |
| `B/HISTORICAL_INDEX.md` | 5,825 bytes | `archive/reports/cleanup_20261004_HISTORICAL_INDEX.md` |

这两份说明记录当时清理恢复的未解决问题，不用当前说明覆盖它们。除这两份必要资料外，没有把 B 再复制到 archive。删除前重新核对了 A/B 共 44,758 个普通文件的 SHA、目录内新增文件、符号链接、实时远端 HEAD 及 B 所有嵌套仓库没有外接 worktree。当前用户可读的 `/proc` 中未发现 A/B 的打开文件或工作目录；26 个系统进程的文件描述符不可读取，不把这一检查称为全系统排他保证。

两份说明的 SHA-256 分别为 `e5fe3b2ecc34d41312fd01d08338827bec879fd595ba57d8bd46889887d32b55`、`ca4568ca8b1b29b09413b749eba96af29722babdc837fa402b7191ec1bdffec3`。

目录内容指纹：按相对路径排序的 `[路径, SHA-256]` 列表，以 UTF-8 紧凑 JSON 计算 SHA-256，包含普通文件、不包含 `.git` 目录；Git 对象另核对。

- A：`e4b39247282be09f17970df3353989a349ca67e069b56e8c5a66bd6b070fd4fd`
- B：`12a1ff904a8ca16d237a3270073309d34b5f8a677dec8038d1b7b6e53e0d67d8`

## 必须保留和无法确认的目录

以下路径前缀 `T` 为 `/home/lzy/projects/temp_storage`，`P` 为 `/home/lzy/projects`。均不属于本次 A/B 删除授权范围。

| 原始路径 | 占用 | 用途与依据 | 建议 |
|---|---:|---|---|
| `T/clot_clearance_github_sync_20260930` | 28.377 GiB | 共同 Git 主仓库、旧实验全量数据、原生求解器；1 个 tracked 修改、15,364 个 untracked 文件。其 `.git` 供四个外部 worktree 使用 | 必须保留，尤其不能单删 `.git` |
| `T/clot_github_sync_20261015` | 21.013 GiB | V4 原开发 worktree；1 个 tracked 修改、2,646 个 untracked 文件；大量完整原生科研数据 | 保留科研档案，停止在此日常开发 |
| `T/clot_github_sync_20261011` | 5.509 GiB | worktree；4 个 tracked 源码修改、3,609 个 untracked 文件；1,333 个文件在其他核对的保留工程中未找到同字节副本，约 3.16 GiB | 保留 |
| `T/y_clot_github_sync_20261010` | 1.706 GiB | worktree；1 个 tracked 修改、1,111 个 untracked 文件；174 个文件在其他核对的保留工程中未找到同字节副本 | 保留 |
| `P/clot_clearance` | 1.080 GiB | 无 `.git` 的早期工程；3,791 个文件中 1,876 个未在核对的保留工程中找到同字节副本，约 580 MiB，含早期 VTK 结果 | 保留，不能整目录判成重复 |
| `P/new_clot_simulation` | 1.179 GiB | 唯一开发工程、已验证本地动态库、输入与测试、独立虚拟环境 | 必须保留 |
| `T/curved_clot_discrepancy_gpu_env_20261006` | 2.813 GiB | 历史 GPU 环境，含原生第三方库；新 V4 不使用，历史复现需求未确认 | 暂留，需用户判断历史用途 |
| `T/clot_literature_fem_env_20261007` | 0.145 GiB | 历史文献/FEM 环境 | 暂留 |
| `T/nonlocal_boundary_literature_20261005` | 0.058 GiB | 论文和研究材料 | 保留文献，不按临时目录删除 |
| `P/.venvs/clot_analysis` | 0.312 GiB | 旧分析环境 | 暂留，未作为本轮可删除缓存 |
| `P/.task_envs` | 0.327 GiB | `traction_failure_v1`、`material_audit` 两个旧环境 | 暂留 |
| `P/.task_work` | 0.170 GiB | 一次性脚本、PDF 环境及服务器导出数据；`clot_github_sync_20261011/server_selected` 有 13 个压缩文件未找到同字节副本 | 暂留，未把 gzip 字节不同推断为内容可删除 |
| `P/github_sync_audits` | 0.031 GiB | 历史传输、排除与验证记录、一次性同步脚本 | 保留证据，不再执行旧同步流程 |
| `P/server_migrations/y_clot_20261009` | 0.122 GiB | 迁移清单和验证记录，无须再次生成备份 | 保留 |

**其他项目，排除删除：** `T/ulm-3d-vascular-model-generation` 8.121 GiB 是其他项目多个 worktree 的共同仓库，不能删除其 `.git`；`T/ulm_particle_3d_particle0` 3.378 GiB 是对应 worktree，存在大量 staged/unstaged 变更。它们的非 Clot 提交未逐一与实时远端历史比对，不声称已全部推送。

`T/github_sync_20260927` 和 `T/github_sync_20260928_latest` 是指向其他项目的符号链接，不是独立可回收的完整目录。`T/formal_3D_flow_solver_FEM` 也是链接。其余 `descriptive_names_20260927`（37.65 MiB）、`cleanup_20260927`（2.97 MiB）、`mesh_code_migration_20260927`（0.79 MiB）及更小的 mesh/migration 目录用途涉及其他项目，全部保留。`server_migration` 的 72 GiB 主要属于 BRAVA、打印和声场等项目，不在本次范围。根目录论文和 `clot_rupture.md` 不删除。

## 科研数据与缓存的边界

共对 **164,825 个文件**计算 SHA-256 与 Git blob 身份，读取约 65.94 GB；对实时远端的 **7 个 Clot 分支 tip**核对其 Git 文件树。以下是“在这些分支 tip 中没有同字节原文件”的占用，包含可能有无损 gzip 表示或精简 NPZ 表示的情况，**不能直接称作 GitHub 全历史完全没有任何表示**：

| 工程内科研目录 | 原字节未直接出现在上述分支 tip 的数据 |
|---|---:|
| 9 月 30 日工程 `verification/results/inputs/runs` | 约 21.28 GiB |
| 10 月 11 日工程相应目录 | 约 3.17 GiB |
| 10 月 15 日工程相应目录 | 约 17.59 GiB |

须保留的具体例子：

- `T/clot_clearance_github_sync_20260930/clot_clearance/verification/y_clot_point_cloud_and_early_flow_coupling_v1/cases/C1_GAUSSIAN_A200/saved_frame_force_audit.bin`：约 912.5 MiB。
- 同工程 `verification/lammps_contact_transfer_and_fragment_size_v1/PER_FRAME_DETACHED_COMPONENTS.csv`：约 376.4 MiB。
- `T/clot_github_sync_20261011/clot_clearance/verification/lps_variational_softening_v3_integration_and_single_pilot_v1/20261014T232612Z/full/native_frames.bin`：约 168.3 MiB。
- `T/clot_github_sync_20261015/clot_clearance/verification/clot_fragment_surface_coupling_v3/20261018T170400Z/cases/biphasic/CONNECTED_COMPONENT_HISTORY.jsonl.gz`：约 1.739 GiB。
- `T/clot_github_sync_20261015/clot_clearance/verification/clot_biphasic_heterogeneous_shear_fragmentation_v2/20261017T231900Z/cases/biphasic/CONNECTED_COMPONENT_HISTORY.jsonl.gz`：约 1.122 GiB。

10 月 20 日导出记录明确含 **391 个本地排除项、75 个服务器排除项、2,663 个字段精简状态**。精简状态不是完整原生检查点。旧排除清单曾把 restart 归为可重建二进制；本次不采纳它作为删除理由。服务器完整检查点仍在服务器保留，本轮未核验服务器是否有新增数据。

发现旧非环境目录的 `__pycache__`、`.pytest_cache`、`CMakeFiles` 共约 **137.7 MiB**。三个旧 LAMMPS `CMakeFiles` 合计约 120 MiB，其 `.a/.so/lmp` 均在别处保留，但历史 PRE manifest/测试可能把缓存视为固定输入。本次未执行全部旧历史测试，因此**不批准式地将它们列为可删缓存**。10 月 4 日历史报告记录过缓存误删后仍有 72 个依赖未恢复；本轮不改写该报告或原 SHA。

当前新库 `native/build/liblammps.so` SHA-256：`49154b2b1825b50925b314ba39e72d502c0b1c2ec7f5943c627958923b7fe831`，以及 `libmapped_geometry.so`、7 份自定义 C++/头文件、固定核心源码包、输入和回归检查点均保留。旧 `tools/lammps-contact-onset-v1` 是另一个 worktree，含 3 个未提交源码/文档修改，不能当普通下载副本删除。

## Git 状态和后续方案

新工程是普通独立 Git 仓库，`.git` 为目录；分支 **clot-simulation-v4**，origin 为原 GitHub 仓库。初始重构提交 **29202fe**，规则及审计提交 **88bf0a3**，独有历史说明保留提交 **0aea532**；本轮维护与清理结果也在本地正常提交。**没有 push、没有 upstream；新工程这些本地提交尚未保存在 GitHub。**

删除 A 前，旧 Clot 四个外部 worktree 均归属于 `T/clot_clearance_github_sync_20260930/.git`；相应已提交 HEAD 均被本次核对的实时远端分支历史覆盖。但这不包含各工作树的未提交内容。10 月 15 日开发 HEAD `5036268` 已被 `project-sync-20261020` 的历史覆盖，即使没有同名远端开发分支，也不应误称该提交未推送。

A 的已核对远端：`https://github.com/OZAKI39/ulm-3d-vascular-model-generation.git`，分支 `clot-clearance/project-sync-20261020`，提交 `032cb75ecf2b99290fceb0f4ef672f91eb3df981`。

**已按用户确认**在新工程配置同一仓库的 origin，并将本地分支改名为 `clot-simulation-v4`；剩余三个旧外部 worktree 保持原状。本工程 root 布局与旧仓库的 `clot_clearance/` 布局不同，且是独立根提交；不合并无关历史、不覆盖旧 main、不 force push。只在用户明确要求推送时首次 `git push -u origin clot-simulation-v4`，随后一直用普通 `git push`。正常 commit/push 命令见 README。

新 `.gitignore` 不再整体忽略 `results/`，正式结果、日志、检查点不会被宽泛扩展名规则隐藏。当前生成目录被忽略不代表可删除。`AGENTS.md` 已禁止自动同步副本、推送、创建普通提交专用 worktree、无授权的大型仿真和不明用途删除。

旧模式入口已定位：

- `T/clot_github_sync_20261020/clot_clearance/AGENTS.md:3` 指向快照导航；10 月 15 日版本第 3 行也要求先读 `github_sync_20261015/README.md`。它们**没有明确要求每次新建同步目录**，不能把历史目录名误称成强制规则。
- `P/github_sync_audits/clot_20261015/export_payload.py`、`write_navigation.py`、`verify_remote.py` 写死同步目录并执行导出/导航改写。
- A 内 `clot_clearance/github_sync_20261020/export_files.py`、`finalize_snapshot.py` 形成旧复制和筛选流程。
- 新 `clot_sim/`、`native/build.py`、`native/env.sh`、配置均没有继承上述同步执行逻辑；旧路径只在来源报告和测试禁止访问的断言中出现。

## 验证与限制

- 删除前小型回归 **9 passed，17.72 s**；删除后 `source native/env.sh; .venv/bin/python -m pytest -q --junitxml=results/validation/junit.xml` 再次 **9 passed，27.40 s**。114 条 VTK/NumPy 废弃提示，无测试失败。
- 覆盖节点/键与原 ID、积分体积与材料身份、材料能量/力、表面与接触力、损伤/kappa、短步推进、连通性、能量账本、保存和恢复。测试中的文件审计钩子禁止读取旧 `temp_storage` 或服务器历史根目录。
- `python -m clot_sim.main inspect` 正常：1927 材料点、1561 锚点、153986 参考键、体积 0.240875 mm³、horizon 150.5 µm。
- `git diff -- clot_sim native configs inputs tests` 无改动；没有修改物理参数、损伤逻辑或现有测试。输入及原生库哈希核对通过。
- `.gitignore` 的 13 项路径检查通过；关键输入、C++、科研 restart 和日志保持可见，环境和明确缓存被忽略。`ldd` 未发现当前库缺失依赖或链接到旧项目目录。
- 删除后对其余 **119,930 个保留旧文件重新计算 SHA-256，全部与清理前一致**；没有修改旧科研结果。测试只更新新工程的可再生验证报告。
- 目前没有新工程运行所需的旧目录 Python/C++/数据依赖。外部系统库、独立虚拟环境和新工程本地原生库仍是显式依赖。
- 删除后的回归及文档中的 inspect 工作流已经通过。未执行全量历史模型复现、MPI 多进程或正式 3.6 ms 轨迹，不宣称这些验证已通过。

## 本轮临时审计索引

索引仅含路径、大小、哈希和只读 Git 状态，放在 `/tmp`，没有复制科研数组或项目源码。下面是用户批准附带清理的 **C 项**；只删除所列精确路径，释放量不计旧项目成果。

精确路径（23 个普通文件，当前分配 155521024 bytes）：

```text
/tmp/clot_cleanup_0930_du.tsv
/tmp/clot_cleanup_1011_du.tsv
/tmp/clot_cleanup_1015_du.tsv
/tmp/clot_cleanup_compare.py
/tmp/clot_cleanup_comparison.json
/tmp/clot_cleanup_extra_du.tsv
/tmp/clot_cleanup_git.json
/tmp/clot_cleanup_git_audit.py
/tmp/clot_cleanup_hash.py
/tmp/clot_cleanup_hashes.json
/tmp/clot_cleanup_inspect.json
/tmp/clot_cleanup_inventory.json
/tmp/clot_cleanup_legacy_du.tsv
/tmp/clot_cleanup_projects_du.tsv
/tmp/clot_cleanup_pytest.log
/tmp/clot_cleanup_recovery_git.json
/tmp/clot_cleanup_recovery_git.py
/tmp/clot_cleanup_recovery_symlinks.json
/tmp/clot_cleanup_remote_audit.py
/tmp/clot_cleanup_remote_blobs.json
/tmp/clot_cleanup_remote_heads.tsv
/tmp/clot_cleanup_staging_compare.json
/tmp/clot_cleanup_temp_du.tsv
```

## 最终执行记录

- A/B 两个目录均不存在；旧 `project-sync-20261020` 分支与共同 Git 对象库保留。
- C 的 23 个获批临时文件已删除，实际分配共 155,521,024 bytes（148.32 MiB）。没有删除清单之外的文件。
- 最终磁盘可用 **487,801,950,208 bytes（454.30 GiB）**；`temp_storage` 现占 **71.16 GiB**。
- 没有 push，没有 Git 历史重写，没有删除远端分支，没有清理其他项目或服务器。
- 日常维护只使用 `new_clot_simulation`；保留旧目录是为保存尚不可删除的科研数据与历史证据，不能误称历史磁盘内容已经全部清空。

## 第二批：旧 Clot 子目录及专用环境（2026-10-21，已批准并执行）

本节是在第一批完成后继续核对得到的精确清单。用户明确回复 **“批准第二批，保留必要原件后清理”**，下面七个路径已完成处理。上文的“本轮未访问服务器”“暂留”等是第一批执行记录；本节补充了只读服务器核验和逐文件保留结果。用户进一步要求“只处理与 clot 相关项目的文件”，本节据此缩小清理边界。

### 范围与必须保留的内容

- `clot_cleanup_recovery_snapshot_20261004T203008Z` 已在第一批批准后删除，不再计算回收空间。
- 10 月 11 日及 10 月 10 日同步目录是 **Git worktree**，但根部的 `model_generate.py`、`cfd_lumen_config.yaml`、`requirements.txt` 和原始 README 包含独立的血管生成项目。本批仅处理其中的 **`clot_clearance/` 子目录**；根部 24 个普通文件、`.git` 连接及共同 Git 对象库保留原位。不会把这两个父目录整棵删除。
- 保留 `temp_storage/ulm-3d-vascular-model-generation`、`ulm_particle_3d_particle0`、其他项目的同步目录及符号链接、BRAVA/声场/打印工程、共享 Python 解释器和环境。不遍历符号链接实施删除。
- 9 月 30 日共同 Git 工程约 **28.377 GiB**、10 月 15 日科研工作树约 **21.013 GiB**，本批仍保留。它们含大量尚未逐项分离的正式数据和原生科研状态，不能只因版本旧而删除。新工程运行不使用这些旧目录。
- `.venvs/clot_analysis` 约 **0.312 GiB** 暂留：包文件可以重建，但本次未找到直接使用它的 Clot 执行记录，不仅凭名字判定其用途。`.task_work`、历史迁移/审计证据、论文和未明确用途的文件仍保留。
- 服务器 `/workspace/y_clot/` 只读核对，全部保留。核对发现 `workspace_is_volume=false`，部分完整原生数据只有旧本地文件和服务器副本，没有在核对的 GitHub 分支 tip 找到原字节文件。因此不建议把服务器容器当作它们唯一的保留位置。

### 精确处理清单

以下为批准时的精确路径、实际分配字节及处理依据。第 1—3 项已经先完成下节的逐文件保留，再清除剩余旧内容；四个环境也已删除。重跑相应历史版本需重新安装依赖，这不影响新工程已经验证的独立环境。表内“建议处理”保留批准时的表述，执行结果见本节末尾。

| 原始路径 | 占用空间 | 文件用途 | 判断依据 | 建议处理 |
|---|---:|---|---|---|
| `/home/lzy/projects/temp_storage/clot_github_sync_20261011/clot_clearance` | 5,915,328,512 bytes / 5.509 GiB | 历史 Clot 源码、完整原生状态和结果 | 普通文件逐一 SHA-256/Git blob 比对；4 个未提交源码修改的实际内容在更新的远端 Clot 快照中存在；未直接在 GitHub 保存的必要原件列入保留表 | 保留清单中的原件后移除本子树，保留父工作树的非 Clot 文件 |
| `/home/lzy/projects/temp_storage/y_clot_github_sync_20261010/clot_clearance` | 1,831,718,912 bytes / 1.706 GiB | 迁移阶段 Clot 源码和验证数据 | 同样逐文件比对；独有修改 `AGENTS.md` 和本地文献搜索记录单独保留；原生状态纳入保留表 | 保留原件后移除本子树，保留父工作树的非 Clot 文件 |
| `/home/lzy/projects/clot_clearance` | 1,159,360,512 bytes / 1.080 GiB | 无 Git 的早期 Clot 工程 | 有旧 VTK 轨迹、图像和验证输出，不能直接视为整份重复；逐文件区分远端/本地已有内容和待保留内容 | 只转存必要历史原件，随后删除这个旧 Clot 工程 |
| `/home/lzy/projects/temp_storage/curved_clot_discrepancy_gpu_env_20261006` | 3,020,464,128 bytes / 2.813 GiB | 旧 Clot 曲壁 GPU 审计环境 | 历史 `INSTALLATION.json`、`RUN_COMMANDS.json` 明确引用；4,184 个包 RECORD 哈希匹配，无独有源码或 editable 安装 | 保存包版本记录后删除环境本身，共享解释器与 CUDA 系统环境不动 |
| `/home/lzy/projects/temp_storage/clot_literature_fem_env_20261007` | 156,057,600 bytes / 0.145 GiB | 旧 Clot 文献/FEM/渲染附加环境 | Clot 历史任务和复现说明明确引用；570 个包 RECORD 哈希匹配 | 保存包版本记录后删除环境本身，基础 dolfinx 环境不动 |
| `/home/lzy/projects/.task_envs/traction_failure_v1` | 333,651,968 bytes / 0.311 GiB | 旧 Clot 牵引映射分析环境 | Clot `DEVELOPMENT_NOTES.json` 明确记载；4,133 个包 RECORD 哈希匹配 | 保存版本后只删除该子目录，保留 `.task_envs` 父目录 |
| `/home/lzy/projects/.task_envs/material_audit` | 16,998,400 bytes / 0.016 GiB | 旧 Clot 材料论文核查环境 | Clot `EXECUTION_SUMMARY.json` 明确记载；538 个包 RECORD 哈希匹配 | 保存版本后只删除该子目录，保留 `.task_envs` 父目录 |

四个环境未发现未归属包文件的自定义二进制或可编辑安装。对其他项目脚本和文档的定向检索未发现对这些环境的引用；当前用户可读的进程中未发现相关打开文件、映射库或工作目录。26 个系统进程不可完整读取，不将此描述为全系统排他保证。

### 最小原件保留方案

[逐文件保留清单](archive/legacy_research/PRESERVATION_PLAN.csv) 有 **3,153 条源路径**，按 SHA-256 去重后为 **2,867 份原件、3,954,196,480 bytes（3.683 GiB）**。其中既有本地独有早期结果，也有已与服务器原件核对一致、但需要保留本地一份的完整状态。这里的“独有”是指在本次核对的保留目录和 7 个实时远端 Clot 分支 tip 中未找到同字节文件，不代表已穷尽所有 Git 历史或外部存储。

已经在同一文件系统中将这些原件移动到 `archive/legacy_research/data/`，逐一核对 inode、大小和 SHA，没有压缩、降采样或改写内容；286 条重复源路径引用同一份保留文件。原件按早期研究、10 月 11 日原生研究和 10 月 10 日研究分组，没有复制旧项目全目录。少量无法取代的历史构建/渲染脚本只作为实验来源记录保存，不被当前代码 import。旧 `AGENTS.md` 改为 `AGENTS.original.txt` 文件名，内容字节不变，不继承它的指令。含新归档目录元数据的实际占用为 **3,955,253,248 bytes（3.684 GiB）**。

包版本和基础解释器记录见 [ENVIRONMENTS.json](archive/legacy_research/ENVIRONMENTS.json)。两个随旧研究附带的 pypdf 安装也已按各自 RECORD 核对（各 65 项通过），可从包重新安装，无须保留一套重复的第三方 Python 源码。

保留数据已经单独列入 `.gitignore` 的精确路径 `/archive/legacy_research/data/`，清单和 [RETENTION.json](archive/legacy_research/RETENTION.json) 执行记录保留在 Git 中；没有忽略整个 `results/`。原件仅保留在本地，**未随源码推送**，不能因有清单就视为可删除。新工程不从 `archive/` 加载运行依赖。

### 空间、执行与验证

- 七个获批目录处理前合计 **11.580 GiB**。保留原件后，磁盘可用从 **487,820,865,536 bytes** 增至 **496,298,725,376 bytes**，观测净回收 **8,477,859,840 bytes（7.896 GiB）**；后台活动可能带来微小变化。下列临时审计索引的回收另计。
- 删除完成时磁盘可用 **462.214 GiB**；最终读数见末尾记录。
- 两个旧工作树已用 Git 稀疏检出，仅保留原有根文件（模式 `/*`、`!/*/`），各占 **225,280 bytes**；原有工作树注册、24 个根部普通文件、`.git` 指针、所有原分支和共同对象库保留，状态均干净。操作中 Git 曾物化其他已跟踪子目录，随后通过稀疏检出恢复到原有根文件集合；没有删除或修改原先存在的非 Clot 文件，也没有保留新增副本。
- 删除前重新核对 **31,198 个源文件 SHA、30 个符号链接**、Git 状态和四个环境的包 RECORD；已核对未提交内容的远端/归档保留位置。删除后再次核对全部归档映射、24 个根文件、Git refs/指针和当前两份原生库，全部通过。
- 当前 `clot_sim/`、`native/`、`configs/`、`inputs/`、`tests/` 相对重构基线 `29202fe` 无改动；运行源码和配置未引用七个待清理目录。通过临时 Python 文件审计钩子禁止打开这七个目录中的文件后，现有回归 **9 passed，24.52 s**，114 条既有 VTK/NumPy 废弃提示。未改测试、物理参数或损伤逻辑；没有重跑正式 3.6 ms 轨迹。保留清单的 3,153 条来源、2,867 个目标、去重内容一致性及目标路径边界检查通过。
- 删除后现有回归 **9 passed，18.89 s**，`inspect` 正常。当前物理源码、C++、配置、输入及测试仍与 `29202fe` 相同。没有重跑正式轨迹，没有修改其他项目或服务器。

用户批准附带清理本批专用的五个临时审计文件，以下精确路径已全部删除，合计 **19,804,160 bytes（18.89 MiB）**，不计入旧工程回收量：

```text
/tmp/clot_legacy_audit.py
/tmp/clot_legacy_audit.json
/tmp/clot_legacy_remote_check.py
/tmp/clot_legacy_env_check.py
/tmp/clot_legacy_envs.json
```

第二批最终磁盘可用 **496,318,365,696 bytes（462.232 GiB）**。七个获批路径和五个临时文件均已不存在；两个父 worktree 仍在，只保留原有根文件。归档数据的内容身份和执行结果已写入 `RETENTION.json`，原始逐文件清单保留不变。其余旧科研数据、共享依赖和其他项目没有附带清理。
