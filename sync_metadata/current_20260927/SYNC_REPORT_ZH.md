# 当前血管工作流 GitHub 同步报告

日期：2026-09-27。目标仓库：<https://github.com/OZAKI39/ulm-3d-vascular-model-generation>。

新分支：`sync/current-h0-dt1ms-wss-audit-20260927`，以 `37c40d7c08ca6527c8379c0ffa81436ae708f59f` 为父提交建立当前工作区快照。原分支与原始源码、配置、结果、服务器均未改写。本次提交的目录层级改为清理后的多个工作区，旧分支的根目录布局不再作为新分支入口。

## 范围和来源

完整导航见仓库根目录 README。当前 H0 场、dt=1 ms 微泡、独立 RBC 及最新 WSS 审计均保留。仍被当前依赖契约读取的旧名称文件也保留，不等同于重新启用旧阶段。

本地纳入 **13,238 个常规源文件，1,630,812,825 字节**，另保留 90 个转换后的相对链接及 950 条第三方上游 LFS 指针元数据。来源包括当前五个主工作区、共享几何代码及 A 数据、残差表格、当前目录和服务器路径说明、本轮清理/迁移回执。逐文件源路径、仓库路径、尺寸和 SHA256 见 [local_file_inventory.csv](local_file_inventory.csv)。本次新生成的说明、验证日志和服务器补充不计入该本地源文件数。

服务器通过已有 SSH 别名 `vast4090` **只读**核对。对当前 H0、微泡正式任务、出生采样来源、RBC 工程/结果、求解器构建和源码，共记录 **12,512 个文件，1,000,486,368 字节**。按内容去重后新增 **1,102 个唯一文件，56,918,505 字节（约 54.3 MiB）**；其余均映射到快照中已有字节。源文件计数按所选路径计，因此相同内容在不同运行目录重复出现，不能把它当成唯一数据量。

服务器补充内容位于 `server_evidence/current_20260927/`。`server_file_map.json` 是完整映射表，包含远端路径、标签、尺寸、SHA256、快照路径和收录方式；所有映射对象均实际做 SHA256 比较。其中 950 个上游 LFS 占位文件以 `upstream_lfs_pointer_metadata_only` 标明，校验的是清单中保存的完整指针文本，不能解释为上游真实大文件已上传。`server_inventory.json` 还记录被排除对象、链接、运行二进制哈希。`server_solver_provenance.json` 记录服务端 svMultiPhysics 源码提交 `c3f0bb892b765b718f61069ecd9726dbc6d177fd` 和三个已修改源码文件；其字节均已纳入比对，不以 Git 提交号代替实际源码校验。

没有保存 SSH 私钥、账户配置、已安装运行环境。大小与具体排除理由见 [EXCLUSIONS.md](EXCLUSIONS.md)。原工作区 Git 规则另存，根 `.gitattributes` 禁用换行和 LFS 转换以保留科学文件字节；这不代表补齐了第三方库原本未下载的 LFS 测试夹具。

首次推送因这些第三方 LFS 指针缺少目标仓库中的对象而被 GitHub 以 GH008 拒绝。随后仅在同步副本中将 950 个占位文件转为 `upstream_lfs_pointers.json` 可追溯清单，重做完整性检查并修订尚未发布的提交。没有从当前血管算例或轨迹中删减科学数据，也没有改写已经发布的分支。

## 实际执行的检查

| 检查 | 结果 | 证据 |
|---|---|---|
| 当前粒子/RBC 回归与短积分 | 78 项通过，0 失败/错误/跳过 | `validation/particle_current/summary.json`、`tests.log`、`tests.xml` |
| Network-H0 / 1D0D 当前测试 | 49 项通过 | `validation/network.log`、`network.xml` |
| 原工作区轨迹完整性 | 7,760 个受保护文件哈希全部一致；500 正式轨迹、500 配对点；dt=0.001 s | `validation/trajectory_integrity.json` |
| 同步副本轨迹完整性 | 使用副本自身代码和数据运行同一入口，全部通过 | `validation/snapshot_trajectory_integrity.json` |
| 当前管网输入保护 | 40 个输入映射至仓库副本后全部一致 | `validation/network_protected_inputs.json` |
| WSS 审计的独立附件重算 | 原始面片 WSS 与显示节点 WSS 的最大差均为 0 Pa | `validation/wss_portable.log` |
| 本地源文件与副本逐字节核对、服务器映射、内部链接 | 详见最后的机器可读检查结果；必须无失败后才提交 | `validation/snapshot_validation.json` |

同步时没有重新进行完整 3D FEM 求解，也没有重新积分全部 500 条轨迹。`check_current.py` 默认排除了 `worker_independent_fixed_subset` 与 `same_host_adapter_preserves_frozen_p9a1_samples` 两项完整重积分测试；不能把 78 项通过解释为本次重新完成生产求解。WSS 附件重算只证明保存和恢复一致，不证明 WSS 已达网格无关性；科学诊断和局限以 WSS 中文主报告为准。

同步副本第一次运行轨迹完整性检查时，发现受保护清单包含 RBC 来源中的 `USER_REQUEST.txt`。打包规则已经修正，保留该契约要求的原始字节，随后副本检查通过；未删改或重写科学保护清单。

## 复核命令

以下 `SNAPSHOT` 指克隆后的仓库根路径。`PYTHON` 需指向具有项目依赖的 Python 环境；本次实际使用 `/home/lzy/projects/temp_storage/ulm_particle_3d_particle0/.venv/bin/python`，WSS 使用 FEM 的 `.venv/bin/python`。依赖版本记录见 `validation/python_environment.json` 和 WSS `REPRODUCE.md`。

```bash
cd "$SNAPSHOT/ulm_particle_formal_p9a5"
PYTHONDONTWRITEBYTECODE=1 "$PYTHON" -B scripts/check_current.py --output "$SNAPSHOT/sync_metadata/current_20260927/validation/particle_current"
PYTHONDONTWRITEBYTECODE=1 "$PYTHON" -B scripts/verify_current_data.py
cd "$SNAPSHOT/ulm_3D_vascular"
PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 "$PYTHON" -B -m pytest -q -p no:cacheprovider tests/network_h0 tests/network_1d0d
cd "$SNAPSHOT"
PYTHONDONTWRITEBYTECODE=1 "$PYTHON" -B formal_3D_flow_solver/FEM_SimVascular/wss_audit/scripts/verify_bundle.py
```

此次 78/49 项测试在原工作区执行；随后又在同步副本中执行轨迹完整性与独立 WSS 检查。部分测试/全流程入口仍依赖原绝对目录，不能将上面的命令理解为所有测试都已在全新机器运行过。

`tools/package_local.py` 和 `tools/collect_server.py` 记录本次收集算法；前者依赖原 WSL 工作区，后者需要当前 SSH 权限。它们会重新生成快照文件，仅应在独立同步副本中运行。`tools/validate_snapshot.py` 会比对原工作区未变、快照完整性和服务器哈希映射，依赖本机原路径；`snapshot_manifest.csv` 则可供另一台机器仅对下载副本做逐文件 SHA256 复核（不包含自身和验证输出，避免自引用）。

科学源码/数据来自所列工作区的实际文件，而不是各自旧 HEAD 的检出。原 HEAD 和跟踪文件改动数量记录于 `source_git_provenance.json`；最终同步提交保存当前实际状态，不据旧 HEAD 宣称本地没有未提交修改。
