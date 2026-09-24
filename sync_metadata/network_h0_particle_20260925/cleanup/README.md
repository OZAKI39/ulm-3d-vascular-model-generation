# 血管流场 / Particle 代码清理记录

完成于 2026-09-24T22:12:45.461839+00:00。本次按用户要求清理本地相关工作树，同时梳理并删除服务器上的确认退役脚本。

共删除 **1,100 个代码文件、90,614 行**，另删除 26 个本地旧字节码缓存。删除文件共约 8.17 MiB；本次主要减少源文件与上下文噪声，未清理占空间较大的科学结果或第三方安装。

## 具体范围

| 本地范围 | 删除代码文件数 |
| --- | ---: |
| 原 FEM_SimVascular | 226 |
| Flow 工作树 FEM | 235 |
| Particle 工作树内 FEM 副本 | 235 |
| Particle 一次性脚本 | 33 |
| 本地合计 | 729 |

服务器：旧 SV 试装/探测/运行编排 63 个文件；P8.2A 历史部署副本中的过期部署/交付脚本 274 个；P9A2 副本 33 个；旧 P9A diagnosis 收尾脚本 1 个，共 371 个。

- 退役：早期 SV 阶段试装、失败重试、阶段准备和收尾脚本；旧 Particle 部署、迁移和一次性交付脚本。
- 一并归档：Flow 与 Particle 两份 FEM 副本中，各自的 `test_main_history_preserved.py` 和 `test_fem_sync_manifest.py`。它们断言整次历史复制交付的脚本/README 字节不变，与用户本次明确授权的代码清理不再相容。历史复制清单仍保留。这不是放宽当前物理测试。
- 保留：`sv_validation` 全部共享实现、Particle 全部科学包、当前 adapter 与 source snapshot、当前/历史数值回归、仍使用的 GPU 工具链构建与启动脚本、当前可视化接口。没有改动微泡公式、参数、入口采样、积分逻辑或可视化风格。
- 保留独立能力：Network v1 图/端口生成（H0 仍依赖）、旧 DOLFINx FEM、LBM、human/MeVO/TopBrain 预处理及尚不能确认无调用的代码。没有仅凭阶段编号或未被新主线调用就删除共享类和函数。
- 没有删除流场、轨迹、正式报告、网格、checkpoint、第三方源码/二进制或 Python/CUDA/MPI 环境。

逐文件原路径、SHA-256、字节数和删除理由见 [local_manifest.json](audit/local_manifest.json) 与 [server_manifest.json](audit/server_manifest.json)。实际 unlink 日志保存在 `audit/*deleted.jsonl`。静态调用/路径审计见 [local_plan.json](audit/local_plan.json)。

## 当前入口与上下文

[ACTIVE_VASCULAR_WORKFLOW.md](/home/lzy/projects/ACTIVE_VASCULAR_WORKFLOW.md) 是跨仓库入口；服务器对应 `/workspace/ACTIVE_VASCULAR_WORKFLOW.md`。
8 份原说明文件已备份，现有 README 已指向 Network H0 与当前 30 泡配对验证。Particle 过期的失败状态说明已更新为 `NETWORK_FLOW_MB_VALIDATION_PASS_WITH_STATIONARY`；500 条新流场生产尚未启动。
4 个仓库根目录新增 `.ignore`，只将生成数据、历史报告和第三方目录移出默认 `rg` 搜索；当前报告里的 adapter 目录明确保留。Python、Git 和文件路径读取不受影响。历史搜索用明确路径或 `rg --no-ignore`。

## 验证结果

- Particle：清理前后均 24 passed，包含 12 项当前配对验证和 12 项 P9A1 回归。
- Network A/H0：清理前后均 49 passed。
- 当前流场 `tests/flow_2mmps`：42 passed；当前主线共 **115 项通过，0 项失败**。
- FEM 历史全量：清理前 516 passed / 50 failed / 4 errors / 68 skipped；清理后 513 passed / 48 failed / 4 errors / 68 skipped。删除的 5 个测试属于上述旧复制交付合同，原有 2 个失败随其退役；其余失败/错误节点与清理前相同，**无新增失败**。历史失败包括旧服务器/源码绝对路径、未镜像产物与旧阶段冻结假设，不应表述为全量测试通过。
- 本地 929 个关键文件 SHA-256 不变，包括科学包、两份流场、当前配对结果/快照、共享 FEM 代码和旋转可视化目录。
- 服务器 122 个关键文件 SHA-256 不变，包括当前 Particle 科学包、入口、固定 GPU 求解器、PETSc、MPI 与必要 wrapper。
- 服务器当前 Particle 入口导入成功，H0 场实际载入 70,363 点 / 371,402 四面体；GPU 求解器 `ldd` 无缺失库。未重新积分或运行 CFD。
- 服务器 27,374 个保留文件的大小、mtime、mode 不变；本地对应审计见 [local_preservation_after.json](audit/local_preservation_after.json)。
- 归档逐文件哈希已校验，服务器归档另下载一份到本地并校验一致；本地恢复演练解包 763 个原文件且 SHA 全部匹配。
- 没有修改 Git HEAD 或暂存区字节，没有提交或重置。已有 staged/dirty 状态保留。

检测到本地 TopBrain/BRAVA 工作正在并发修改 `vascular_processing` 中的源码和 `outputs/topbrain_brava_transfer/environment_summary.json`。它们不属于删除/写入清单，均保持现场，未回滚或覆盖。元数据审计把这些外部并发变更单独列出，未将它们误报成“全部保留文件都完全未变”。

完整测试日志在 [logs](logs/)，机器摘要为 [SUMMARY.json](SUMMARY.json)。

## 归档与恢复

- 本地代码归档：`backups/local_retired_code.tar.gz`，SHA-256 `1603a51ae1adae746dfedea1b669316510e360b5c73916c82f4a0d822a208e13`。
- 服务器代码归档（本地副本）：`backups/server_retired_code.tar.gz`，SHA-256 `c1fddfd824880a9f320647af67c126f7a71e21ffafe014ecfdef6580e87c0941`。
- 服务器原归档目录：`/workspace/archives/vascular_workflow_cleanup_20260924T214356Z`。

只恢复缺失的本地退役文件，任何同名现存文件都会让工具拒绝操作：

```bash
python3 /home/lzy/archives/vascular_workflow_cleanup_20260924T214356Z/scripts/restore_retired.py --manifest /home/lzy/archives/vascular_workflow_cleanup_20260924T214356Z/audit/local_manifest.json --archive /home/lzy/archives/vascular_workflow_cleanup_20260924T214356Z/backups/local_retired_code.tar.gz --in-place
```

仅解包到新的审阅目录，并取回原 README：

```bash
python3 /home/lzy/archives/vascular_workflow_cleanup_20260924T214356Z/scripts/restore_retired.py --manifest /home/lzy/archives/vascular_workflow_cleanup_20260924T214356Z/audit/local_manifest.json --archive /home/lzy/archives/vascular_workflow_cleanup_20260924T214356Z/backups/local_retired_code.tar.gz --destination /tmp/vascular-cleanup-review --include-original-docs
```

服务器同一工具已放在原归档目录，具体命令见 `/workspace/ACTIVE_VASCULAR_WORKFLOW.md`。服务器源归档也可用本地副本先解包审阅。
历史报告保持原字节，因此其中有些退役脚本链接需要恢复文件后再访问；完整历史部署 source-manifest 校验也需要先恢复对应脚本。
