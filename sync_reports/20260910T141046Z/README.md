# 本次 GitHub 增量归档核查

本次仅归档已完成内容，选定 `rbc_repair_20260910T131105Z`，从 `sync/hemocell-mirheo-single-rbc-benchmark-20260910T120318Z` / `7a0be010fa4cfeee4073bba5bab929216cf4df7b` 创建 `sync/hemocell-mirheo-single-rbc-repair-20260910T141046Z`。没有新仿真、编译、标定、物理测试或浏览器运行。根 README 与所有无关历史内容保留。

当前科学状态来自原件：CPU_REPAIR_PREPARED_AWAITING_NEW_AUTHORIZATION，修复运行 NOT_TESTED，材料 NOT_MATCHED，qualified_speedup=null，人工 PENDING。16 项 CPU 与 21 项浏览器 PASS 均为归档前已有记录。归档的完整性 PASS 只检验字节和引用。

- [远程核查总入口](../../BENCHMARK_REVIEW.md) · [中文原报告](../../mirheo_starter/data/single_rbc_repair/rbc_repair_20260910T131105Z/report_zh.md) · [原 HTML](../../mirheo_starter/test_code/outputs/single_rbc_repair/rbc_repair_20260910T131105Z/comparison_review.html)
- [清单](SYNC_MANIFEST.json) · [本地同步校验](SYNC_VALIDATION.json) · [路径映射](PATH_MAPPING.json) · [结论字段索引](RESULT_INDEX.json)
- [排除文件](EXCLUDED_FILES.json) · [大文件审查](LARGE_FILE_REVIEW.json) · [敏感信息扫描规则与复核](SENSITIVE_SCAN_REVIEW.json)
- [源 Git 状态](SOURCE_GIT_STATE.json) · [无写入检查](SOURCE_QUIESCENCE.json) · [原生补丁重建说明](NATIVE_SOURCE_ARCHIVE.md) · [补丁字节核查](PATCH_INTEGRITY.json)

源文件新增复制 87 个 / 38927152 字节；另有 185 个直接依赖核对为基础分支原有字节。最终新增/更新数含本目录派生说明与根入口，见提交和本地索引/远程回执。本轮新增单文件均小于 50 MiB，旧 57,568,143 字节 results.json 原样继承。核心 JSON/CSV 可直接读取，未转 LFS、未降采样或重新格式化。

排除清单逐文件记录 2518 个本轮相关 checkout/cache/compiled-library 文件的大小、SHA-256 与执行位；其他未选中的完整环境、vendor、浏览器临时配置、凭据和无关项目按目录策略排除。没有省略本轮已存在的科学观测或进行脱敏。原本未保存的准确崩溃状态、新候选尚未产生的结果继续标为缺失，不伪造占位。

原 JSON、报告和 HTML 内的 WSL 路径及同步前 `git_commit_or_push=false` 不改写。`comparison_results.json` 是浏览器检查前的结果快照；`delivery_receipt.json` 和 `review_latest.json` 记录稍后浏览器核查，HTML SHA 一致。旧 `*_pre_review_fix.json` 和 CPU/导出失败日志同样保留。同步链接核查发现原修复 HTML 的历史区域引用 `browser_delivery_verified/browser_check.json`，但原修复输出目录缺少该伴随文件；归档仅在新 HTML 相对路径放入[历史记录的逐字节副本](../../mirheo_starter/test_code/outputs/single_rbc_repair/rbc_repair_20260910T131105Z/browser_delivery_verified/browser_check.json)，精确来源见 PATH_MAPPING，未修改 HTML 或伪造新的浏览器核查。此记录的 HTML 哈希仍指向旧页面；当前页面的记录在 `browser_cpu_delivery`。

完整检出本分支后可运行：

```bash
python3 -B sync_reports/20260910T141046Z/verify_archive.py
```

该脚本只做 JSON/CSV/NPZ 读取、Python AST 语法、链接和哈希检查；不导入求解器，也不重跑原 CPU 科学测试或浏览器。`--check-sources` 仅适用于原 WSL 来源仍存在的机器。打开 HTML 请保留对应 data/output 目录。

清单不哈希自身或 validation，二者身份由最终 Git 提交确定。推送后独立从 GitHub 获取新分支，核对 HEAD、所有变更 blob/mode、继承树，以及实际读回的入口/结果/CSV/HTML。该最终回执保存在独立同步目录 `.git/sync_work/REMOTE_VERIFICATION.json`，避免为自引用反复提交；不能把推送前的 validation 当成远端已通过的证明。
