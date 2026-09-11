# 云端 RBC 与最新本地修复结果归档

本次只复制、校验、提交与推送，不运行仿真、编译、科学测试或浏览器。基础分支 `sync/hemocell-mirheo-single-rbc-repair-20260910T141046Z`，base commit `26c22ed53df1d867aa39443fd65173b0770409c0`；新分支 `sync/hemocell-mirheo-cloud-rbc-20260911T073240Z`。

[核查总入口](../../BENCHMARK_REVIEW.md) · [最终中文报告](../../cloud_results/cloud-rbc-full-20260911T000225Z-review/report_zh.md) · [原 HTML](../../cloud_results/cloud-rbc-full-20260911T000225Z-review/rbc_full_review.html) · [原始数据](../../cloud_results/cloud-rbc-full-20260911T000225Z/simulation)

[SYNC_MANIFEST](SYNC_MANIFEST.json) 记录精确来源、大小、SHA-256、执行位与新增/继承状态；[PATH_MAPPING](PATH_MAPPING.json) 保留原路径映射；[排除项](EXCLUDED_FILES.json) 保留省略原因及影响；[校验结果](SYNC_VALIDATION.json)、[大文件审查](LARGE_FILE_REVIEW.json)、[敏感信息复核](SENSITIVE_SCAN_REVIEW.json)、[原生源码说明](NATIVE_SOURCE_ARCHIVE.md)。

来源文件 4402 个，其中 110 个已由基础分支逐字节继承；4292 个新增或更新来源文件合计 851384410 字节。最大单文件 34542407 字节，未新增大于 50 MiB 的文件。最终提交计数还包括本次派生说明、清单与相对入口。所有 A0–A6 科学输出、云端原始封包和本地派生报告均保留，不按结果筛选。

当前云端结论是未完成 Γ=4、准备失稳，运行 NOT_COMPLETE / 数值 FAILED / 回传 PASS。旧构建失败、smoke 成功、CPU 准备态和最终正式失败属于不同阶段，保持原来的时间及状态。原始报告、CSV、JSON 和 HTML 没有改字节；本次仅新增归档入口与映射。缺测不补造，原始表头-only 文件保留。

完整环境、第三方 checkout、可执行库/对象、缓存、锁和 .env 排除。原生源码采用提交、原审计子集、完整修改文件/新头文件和补丁；构建不在归档中重跑。没有上传聊天全文或认证凭据。

```bash
python3 -B sync_reports/20260911T073240Z/verify_archive.py
```

该检查仅哈希、AST、JSON/CSV/JSONL/NPZ 读取、HDF5 签名与 XMF 文件引用，核对冻结依赖和已有浏览器记录，不进行几何重算或仿真。`--check-sources` 只用于仍拥有 WSL 原件的机器。

SYNC_MANIFEST 不包含自身或 SYNC_VALIDATION，避免递归哈希；这两个文件身份由 Git 提交确定。推送后的 HEAD、重新取回 Git 对象及字节核查记录保存在本工作树 Git 私有目录的 `sync_work/REMOTE_VERIFICATION.json`，不反复 amend 以记录自身提交。
