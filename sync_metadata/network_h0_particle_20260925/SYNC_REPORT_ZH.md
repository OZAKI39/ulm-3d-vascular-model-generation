# 当前血管流场与微泡工作流同步

分支：`sync/network-h0-particle-cleanup-20260925`。基于已发布工作流提交 `fab4cf0eb3c8f8ad1ce40dfe995181c147cb9a9d`，日期为 2026-09-25（柏林时间）。主分支不变。

本次在独立 Git 工作树整理后发布，保留各源工作树的未提交工作；包含清理后的当前 vascular A → ROI → H0 → 3D FEM → Particle 主线，以及依赖的历史科学模块、测试、必要输入、冻结流场、配对轨迹、旋转可视化、正式报告和日志。服务器保留源码、配置、原始执行记录与固定求解器补丁一并同步。

## 范围与来源

- 本地选择 29,159 条文件记录，合计 2.79 GiB，逐文件来源、SHA 与排除原因见 [local_file_inventory.csv](local_file_inventory.csv)。多个工作树合并时相同内容可能共享 Git 对象，条目数不等同于新增对象数。
- 服务器导出 17,213 个文件，发布 17,207 个原始文件和 6 份缺失 LFS 指针的 JSON 证据，原始大小约 1.04 GiB；逐项 SHA 已核对，见 [server_file_inventory.json](server_file_inventory.json)。源路径映射到 `server_evidence/`，当前工作流入口仍位于仓库主目录。
- 源 Git 提交、分支和未提交状态见 [source_git_provenance.json](source_git_provenance.json)。迁移只重写当前导航文档、Git 字节保存属性和默认检索排除规则，修改清单见 [packaging_edits.json](packaging_edits.json)。原导航文档另存于 `original_docs/`。
- 前一轮代码瘦身已体现在此分支，父分支中仍有的 265 个废弃路径已删除；清理清单和备份位置见 [cleanup/README.md](cleanup/README.md)。已归档的废弃源码没有重新纳入本分支。

不纳入虚拟环境、编译构建树、完整上游原始脑数据、可重建缓存、逐帧渲染、中间 CFD 状态、重复部署包和非必要大文件。一般新增文件阈值为 25 MiB；现有回归必需的流线候选池作为例外保留，约 29.5 MiB。最终没有超过 50 MiB 的文件。独立进行的 TopBrain/BRAVA 输出不属于本快照。

## 完整性与验证

当前科学源码和结果数组没有修改。服务器实际修改的 `main.cpp`、`petsc_impl.cpp`、`sv13q_reuse.h` 与仓库 vendor 文件逐字节一致。源文件对照结果见 [source_byte_verification.json](source_byte_verification.json)。

在独立临时目录仅使用已导出的文件验证：Particle 24 项、Network 49 项、Flow 42 项，共 **115 项通过，0 失败，0 错误，0 跳过**。原始 pytest 日志、JUnit XML 与 `summary.json` 一并保留。5 份 JSON 的临时迁移和 92 份解释器缓存条目的处理详见 [PORTABILITY_ZH.md](PORTABILITY_ZH.md)；全部 114 个科学源文件及 114 个服务器源码快照保持原期望 SHA。

发布前文本凭据扫描没有发现匹配项。历史 Taylor-Hood 官方示例中的 6 份缺失 LFS 占位文件已转换为 [JSON 证据](historical_lfs_pointer_evidence.json)，保留原始路径、哈希及指针全文，以通过 GitHub GH008 检查；当前 H0/Particle 依赖无此缺项。最大文件与扫描范围见 [preflight.json](preflight.json)。这不是完整安全审计。

`SNAPSHOT_SHA256.txt` 覆盖除清单自身外的所有发布文件。克隆后执行 `verify_snapshot.py --scientific-inputs` 可复核。仓库使用普通 Git 对象，不要求为当前工作流安装 Git LFS。

## 当前科学状态

配对验证仍为 `NETWORK_FLOW_MB_VALIDATION_PASS_WITH_STATIONARY`：30 个同初始条件，OLD/NEW 各 29 条完成、1 条 stationary、0 条数值失败；尚未开始 500 泡新流场生产运行。Taylor-Hood 验证状态为 `TAYLOR_HOOD_VALIDATION_STOPPED_RESOURCE_COST`，其代码和日志是停止阶段证据，不表示已有完整 P2/P1 结果。本次同步没有重新求解流场或积分轨迹。
