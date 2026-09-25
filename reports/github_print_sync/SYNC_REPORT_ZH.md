# 打印模型开发快照同步记录

本次仅同步当前开发状态。新分支 `sync/vascular-print-models-20260925` 以此前 BG001 审阅分支提交 `d8df343e129eb958e1860f5ece6f6b1da6e9f48d` 为基础，保留审阅历史，补齐当前模型生成代码和结果。未纳入其他 CFD 分支的求解输出。

## 范围与完整性

复制 4265 个源文件／链接条目，总文件字节数 1,114,247,613（约 1.038 GiB），最大单文件 32.48 MiB。文件与源工作区逐项 SHA-256 核对通过。普通 Git 对象直接保存文件，关闭文本换行转换，保留已有 CRLF 文件与受保护 UI 的原始哈希。

保留 s1-1 至 s1-4、s2、全部 vascular_processing、配套工具／配置／测试、共享 utils、VascularMD 与 Ultraliser 源码；保留语义迁移、精修、三个紧凑候选、历代端口盒体、STL/STEP/3MF、原始与拟合 SWC、节点／分支映射、日志、CSV、JSON、QC 图片和中文报告。小鼠单样本原始 SWC、图像及 mask 随包保存；派生目录中的数据链接改为包内相对链接。

不包含完整 TopBrain ZIP／全部原始影像、其他巨型数据集、嵌套 Git 历史、编译缓存、虚拟环境和无关 CFD 求解输出。文件级／目录级排除理由见 `excluded_files.csv`。TopBrain 的许可、实际 ITK-SNAP 标签表、25 份 MRA 标签与下载来源保存在 `data/topbrain_provenance/`。

## 真实验证结果

| 检查 | 结果 |
| --- | --- |
| 原项目目录下选定的 110 项回归 | 109 通过，1 失败 |
| 异地快照目录下同一组回归 | 104 通过，2 失败，4 个 setup 错误 |
| 4265 个源文件／链接条目完整性 | 通过 |
| 原工作区文件状态、分支及 Git 暂存区 | 保持不变 |
| BALANCED 原模型 STL | 闭合，1 个连通实体，体积 208.648126 mm³ |
| 最新血管与端口 STL | 闭合，1 个连通实体，体积 360.333048 mm³ |
| 盒体 STL／STEP 重新读取 | 通过 |
| 最新 assembly.3mf | 可读取，2 个物理网格实体 |
| 所列敏感凭据规则扫描 | 在已扫描文本中未发现匹配 |

原目录的唯一失败来自历史 `refined_roi/protection_before.json` 与后来更新过的 `s1-3_swc_roi_generate_MeVO.py` 哈希不一致。这是原项目已有的来源保护冲突；本次没有覆盖旧清单、伪造通过状态或回退最新可视化。详细回溯见 `source_workspace_pytest.log`。

异地目录另外暴露了历史绝对路径绑定、VascularMD 缺少自身 Git 元数据的迁移条件。README 提供恢复固定上游提交 Git 元数据的命令；外层项目 Git HEAD 不能代替上游 HEAD。旧清单不可简单批量替换路径后仍沿用其历史验收。初次异地验证日志完整保存于 `pytest.log`，没有标作全通过。

## 尚未完成的几何工作

最新版本是 `print_fixture_design_aligned`，已验证端面轮廓和轴向对齐。用户进一步指出 O3 接缝仍不光顺；周向侧壁切向／法向连续性的修复尚未完成。闭合、体积为正、截面对齐均不能替代光顺性验收。盒体和血管芯仍为两个实体，一体式 union 尚未实现。`prior_union_printer_probe.json` 仅是先前打印机发现记录。

## 复核入口

运行 `python3 tools/verify_print_snapshot.py` 对实际文件和相对链接进行只读核验。`SHA256SUMS` 覆盖快照文件但不包含校验表自身；`source_file_manifest.json` 另列源路径及原哈希。同步分支位于独立工作目录，原工作区无需切换分支。
