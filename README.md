# BraVa CFD 与微泡项目（2026-10-22 同步）

本分支仅保存 BraVa 项目与必要共享依赖。小鼠血管项目见 `sync/mouse-vascular-flow-project-20261022`。两个项目不合并为同一个计算模型；祖先提交仍保留原仓库历史。

当前科学状态以 [computation_examples/brava_flow_roi_18mlmin/ACTIVE_DEVELOPMENT_SCOPE.md](computation_examples/brava_flow_roi_18mlmin/ACTIVE_DEVELOPMENT_SCOPE.md) 为准：P2/P1 路线在资源门槛处停止；Gmsh V2 场局部精度仍未通过，未运行新正式群体。

- [正式流场身份](computation_examples/brava_flow_roi_18mlmin/reports/ACTIVE_FLOW.json)：18 mL/min 历史 balanced_pressure_final；候选 0/15 为同一解的打印姿态。
- [旧正式群体账本](computation_examples/brava_flow_roi_18mlmin/population_v1/runs/realization_0000_contactfix/SUMMARY.json)：1565 条接受轨迹，458/527/437 出流，141 观察截止、2 数值失败。原始 states/chunk_*.npz 全部保留。
- [最新四视角](computation_examples/brava_flow_roi_18mlmin/population_v1/candidate_0/OPEN_RESULTS.html)与[流场可视化](computation_examples/brava_flow_roi_18mlmin/visualization/candidate_0/OPEN_RESULTS.html)。
- [P2 线性代数报告](computation_examples/brava_flow_roi_18mlmin/population_v1/analysis/p2_linear_algebra_rerun_20261006/P2_LINEAR_ALGEBRA_AND_RESOURCE_REVIEW_ZH.md)。
- [新服务器迁移与独立环境验证](deployment/newserver_20261009/MIGRATION_AND_VALIDATION_ZH.md)；当前入口和激活脚本在 `deployment/current_server/`，环境验证不等于局部物理精度认证。

当前固定灌注协议尚无新正式大批量结果；未来新边界按中心线阻力估算出口压力，历史校准场身份保持不变。服务器当前项目与本地已有对应文件 SHA 全部一致；服务器独有、与当前诊断相关的源码与执行证据从已核验备份补齐。

本次是现有文件同步，不启动或恢复 CFD、轨迹、网格加密或敏感性计算，不修改科学数据。源目录、原有兼容链接和完整本地备份保留。

纳入 8,749 个源文件，约 5.89 GB。保留代码、关键日志、冻结网格/流场、原始正式状态、诊断表格和成品图表；排除可安装环境/二进制、缓存、重复中间帧、废弃批次大体积状态、可由已保留输入重建的全域审计数组。排除是发布筛选，没有删除原始数据；普通 Git 不上传单文件 ≥100 MiB 的派生大文件。没有 LFS，也没有将文件切片规避限制。

逐文件来源/大小/SHA256 见 [INCLUDED.json](sync_20261022/INCLUDED.json)；[EXCLUDED.json](sync_20261022/EXCLUDED.json) 记录排除原因、原位置及大文件哈希；[LINKS.json](sync_20261022/LINKS.json) 记录相对兼容链接。少数上游 `.gitattributes` 的 LFS 规则只在发布副本中停用，原文与前后哈希已保留。未在本地取得实际数据的上游测试资源 LFS 指针未发布为文件，路径和原始引用保存在 sync_20261022/UPSTREAM_TEST_RESOURCE_REFERENCES.json；本项目冻结 CFD 数据和原始接受状态均为实际文件。

历史配置和报告中的绝对路径保持原始字节；这不是空白机器一键运行包。HTML/视频可下载本分支后打开；GitHub 页面直接预览部分格式受限。克隆时可加 `--single-branch --depth 1` 仅获取此分支当前快照。

同步分支：`sync/brava-flow-project-20261022`。最终提交和远端核验回执在本地 `github_sync_audits/flow_projects_20261022/`。
