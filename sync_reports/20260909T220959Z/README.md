# 本次归档校验记录

分支：`sync/hemocell-mirheo-benchmark-20260909T220959Z`；UTC 名称：`20260909T220959Z`。基础提交：`48d53e9c86963a0d3520b8d262b79ec21fe76980`。

- [远程核查总入口](../../BENCHMARK_REVIEW.md)
- [SYNC_MANIFEST.json](SYNC_MANIFEST.json)：源路径、仓库路径、字节数、SHA256、执行位、复制/继承身份及 base commit。
- [EXCLUDED_FILES.json](EXCLUDED_FILES.json)：排除的二进制、完整环境、重复 ZIP、临时 profile、独立聊天答复文本及其复核影响。
- [PATH_MAPPING.json](PATH_MAPPING.json)：原 WSL 路径到实际仓库位置，不改写配置和 provenance。
- [SYNC_VALIDATION.json](SYNC_VALIDATION.json)：本次 CPU 只读语法、数据回读、引用、哈希、秘密扫描、Git 索引和来源保护检查。
- [SOURCE_GIT_STATE.json](SOURCE_GIT_STATE.json)：来源工程的原有 HEAD、分支、未提交状态和索引 SHA256。
- [SENSITIVE_SCAN_REVIEW.json](SENSITIVE_SCAN_REVIEW.json)：四处 Plotly URL 拼接表达式的非凭据判定，按精确 HTML 和表达式哈希限定，无原件脱敏。
- [verify_archive.py](verify_archive.py)：只读校验入口，无求解器导入、无编译、无 GPU/浏览器/物理测试调用。

在本分支的完整工作副本中可运行：

```bash
python3 sync_reports/20260909T220959Z/verify_archive.py
```

仅在原 WSL 机器上可加 `--check-sources` 比对原始文件。此脚本只核查归档完整性，输出的 PASS 不代表新的科学验证。
已有的物理回归与 Chrome 记录保留其原始日期和结果。科学状态仍为 PARTIAL，人工验收仍为 PENDING。

本清单不对 `SYNC_MANIFEST.json` 和 `SYNC_VALIDATION.json` 自身或互相递归计算哈希；这两份记录的身份由最终 Git 提交确定。
其他派生说明和校验脚本纳入清单。推送后的提交 SHA 与远程取回核查保存在独立工作区 `.git/sync_work/REMOTE_VERIFICATION.json`，无需反复 amend。

独立克隆采用 sparse checkout，仅减少本机物化的历史文件；Git 分支和索引完整继承基础提交，没有删除历史文件。
Mirheo 继承的 `.gitattributes` 已禁用文本转换和过滤器，新增 HemoCell 与同步目录另设局部规则。暂存采用逐文件白名单。
当前共享预算池另存为本轮 `archive_evidence/shared_budget_pool_at_sync.json`；基础分支原历史池文件和旧账本不被覆盖。

旧血管工作区的完整 git status 因已配置的 git-lfs 可执行文件不可用而标为 UNAVAILABLE。
该目录保持只读，未安装 LFS 或修改配置；其 HEAD、分支和索引另行比对。本轮源文件以复制前后 SHA256 与执行位核查。
普通 WSL 路径与随机种子保留。真实秘密扫描不打印匹配内容；如出现需处理的匹配，将先停止暂存并记录处理，不在原工程脱敏。
