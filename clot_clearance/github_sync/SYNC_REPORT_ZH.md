# GitHub 同步说明

目标仓库：OZAKI39/ulm-3d-vascular-model-generation（同步时为公开仓库）。

新分支：`clot-clearance/energy-regularization-20260930`。基线：`main`，提交 `8a264014a6181d23d30f3ad4447931ce182e1ac4`。

此提交只新增 `clot_clearance/`。原有仓库文件和原本地工程均未修改；本地完整清单的 3812 条记录在复制后重新校验一致。科学文件不改路径字段、不转换换行、不改数值。

## 保留

- 1689 个文件/相对链接按原始身份复制。
- 1 个超限的核心 NPZ 文件按原始字节无损分为三片，实测还原 SHA-256 与原文件一致。
- 全部 18 组 regularization 算例的原始 states、配置、标定、源码快照、日志及诊断。
- 旧直管与强制碎裂的源码、结果、日志和正式可视化。
- 正式新动画、13 组 PNG/PDF、原报告和额外相对链接报告副本。

原始保留载荷 476,659,643 字节；Git 压缩后的传输大小不同。

## 排除

共 2122 个可再生成、重复或非源码依赖文件，672,036,429 字节。逐文件原因及原 SHA-256 在 SNAPSHOT_MANIFEST.json。

- Build products, installed tools, caches or temporary test files: 416 files, 59,481,765 bytes.
- Third-party paper full text; retain DOI and implementation mapping: 2 files, 4,876,892 bytes.
- Index of omitted regenerable VTK exports: 93 files, 96,163 bytes.
- Regenerable per-state VTK export; canonical NPZ states retained: 1545 files, 541,939,913 bytes.
- Superseded visualization layout draft: 61 files, 19,890,587 bytes.
- Derived display interpolation; raw solver states and frame maps retained: 3 files, 13,694,401 bytes.
- Large duplicate legacy GIF; corresponding MP4 and renderer retained: 2 files, 32,056,708 bytes.

## 验证

- `COPY_TESTS.log`：在隔离副本中运行当前 tests，48 passed。
- `SNAPSHOT_CHECK.json`：1690 条被保留源记录全部校验通过。
- `RESTORE_CHECK.json`：三片拼接还原的细网格 NPZ 与原件 SHA-256 一致。
- `PRE_PUSH_CHECK.json`：检查暂存范围、源文件未变、遗漏、文件大小及凭据特征。
- 增加的 gitignore 只作用于本子目录，允许已审查的科学数据，排除还原后的超限原文件与缓存。

本次只同步已经完成的开发内容，不运行新的科学算例。科学状态仍是 NOT CONVERGED、无主算例断键或脱落、碎片分辨率不可评估。原始校验记录中的本地绝对路径属于历史证据；GitHub 包采用本次 snapshot manifest 作为完整性标准。
