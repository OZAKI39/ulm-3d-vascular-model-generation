# 小鼠血管流场与配套微泡项目（2026-10-22 同步）

本分支保存小鼠血管项目；BraVa 项目见 `sync/brava-flow-project-20261022`。共同粒子/求解器代码是共享依赖，不代表两个几何或结果相同。祖先提交仍保留原仓库历史。

- [ROI-only 算例](computation_examples/mean-2p0-mmps-A-ROI-only-balanced-pressure-v1/CASE_README.md)：1500 条、dt 0.5 ms，1315 出流（787/61/467）、185 接触支持静止、0 执行失败。
- [best-feasible 算例](computation_examples/mean-2p0-mmps-A-best-feasible-balance-v1/CASE_README.md)：1500 条、dt 0.5 ms，1315 出流（0/133/1182）、185 接触支持静止、0 执行失败。
- [Network-H0 基线](computation_examples/mean-2p0-mmps-A-H0-pressure-v1/CASE_README.md)：500 条、dt 1 ms，430 出流（0/139/291）、70 接触支持静止。
- [原 2 mm/s 基线](computation_examples/mean-2p0-mmps/CASE_README.md)与[早期 FEM Stage-Q 参考](computation_examples/fem_reference_stage_q/CASE_README.md)。
- [WSS 验证报告](computation_examples/_shared/wss_validation_v2/WSS_VALIDATION_REPORT_V2.md)及 `validation_cases/`；细档血管和 O2 压力敏感性继续为 SKIPPED_BY_USER。
- [ROI-only 流场图表](computation_examples/mean-2p0-mmps-A-ROI-only-balanced-pressure-v1/visualization/OPEN_RESULTS.html)和[微泡图表](computation_examples/mean-2p0-mmps-A-ROI-only-balanced-pressure-v1/microbubble/OPEN_RESULTS.html)。

四套冻结场、求解配置与日志、WSS 表面、原始轨迹及验收记录均保留。`geometry_inputs/` 只补入此流场来源使用的源样本、ROI/表面和共享几何代码，不复制整个上游打印模型库。小鼠项目未迁入当前 BraVa 服务器；本次同步依据本地保留的权威结果，没有声称重新核验旧服务器。

本次是现有文件同步，不启动或恢复 CFD、轨迹、网格加密或敏感性计算，不修改科学数据。源目录、原有兼容链接和完整本地备份保留。

纳入 34,645 个源文件，约 4.97 GB。保留代码、关键日志、冻结网格/流场、原始正式状态、诊断表格和成品图表；排除可安装环境/二进制、缓存、重复中间帧、废弃批次大体积状态、可由已保留输入重建的全域审计数组。排除是发布筛选，没有删除原始数据；普通 Git 不上传单文件 ≥100 MiB 的派生大文件。没有 LFS，也没有将文件切片规避限制。

逐文件来源/大小/SHA256 见 [INCLUDED.json](sync_20261022/INCLUDED.json)；[EXCLUDED.json](sync_20261022/EXCLUDED.json) 记录排除原因、原位置及大文件哈希；[LINKS.json](sync_20261022/LINKS.json) 记录相对兼容链接。少数上游 `.gitattributes` 的 LFS 规则只在发布副本中停用，原文与前后哈希已保留。未在本地取得实际数据的上游测试资源 LFS 指针未发布为文件，路径和原始引用保存在 sync_20261022/UPSTREAM_TEST_RESOURCE_REFERENCES.json；本项目冻结 CFD 数据和原始接受状态均为实际文件。

历史配置和报告中的绝对路径保持原始字节；这不是空白机器一键运行包。HTML/视频可下载本分支后打开；GitHub 页面直接预览部分格式受限。克隆时可加 `--single-branch --depth 1` 仅获取此分支当前快照。

同步分支：`sync/mouse-vascular-flow-project-20261022`。最终提交和远端核验回执在本地 `github_sync_audits/flow_projects_20261022/`。
