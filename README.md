# 血管几何、3D流场、微泡与RBC工作流

2026-10-04：BraVa当前源码、冻结Gmsh V2、GPU steady field及误差诊断同步至 `sync/roi-flow-mb1500-wss-stop-analysis-20260928`。默认入口为 [BraVa active scope](brava_flow_roi_18mlmin/ACTIVE_DEVELOPMENT_SCOPE.md)，结果为 [报告与三张最终PNG](brava_flow_roi_18mlmin/production_workflow/OPEN_RESULTS.html)。

已完成 continuous population simulator、GPU particle backend、Gmsh V2 manufactured reference与GPU svMultiPhysics稳态场。科学状态仍为 **REFERENCE_PIPE_FIELD_FAIL_SCIENTIFIC + MIXED_FIELD_ERROR + P1_FIELD_TRAJECTORY_IMPACT_LARGE（DIAGNOSTIC_ONLY）**。V3、正式CFD-vs-oracle trajectory、dt convergence、新BraVa validation均NOT_RUN；不能写成全部轨迹已验证。

旧BraVa正式群体为1565条（458/527/437出流、141观察截止、2数值失败）。固定灌注1.240815660s+观察10s的新协议尚无正式大批量结果。旧`microbubble/`1500条独立任务保持停止；已有群体动画见 [最新四视角](brava_flow_roi_18mlmin/population_v1/candidate_0/views_0to1p2s_large_axes/OPEN_RESULTS.html)。本轮没有运行科学计算，没有更改production physics。

[同步与排除说明](brava_flow_roi_18mlmin/GITHUB_SYNC_MANIFEST_ZH.md) · [清理计划](brava_flow_roi_18mlmin/SYNC_CLEANUP_PLAN_ZH.md) · [大文件身份](brava_flow_roi_18mlmin/EXCLUDED_LARGE_ARTIFACTS.json)。其它小鼠、RBC、几何工作流保留。

| 内容 | 入口 |
|---|---|
| BraVa当前科学状态与结果身份 | [当前结果说明](brava_flow_roi_18mlmin/CURRENT_RESULTS_ZH.md) |
| BraVa真实网格、压力边界及求解日志 | [算例](brava_flow_roi_18mlmin/cases/balanced_pressure_final/) |
| candidate 0 / 15 的流场动画与4K图 | [candidate 0](brava_flow_roi_18mlmin/visualization/candidate_0/OPEN_RESULTS.html) · [candidate 15](brava_flow_roi_18mlmin/visualization/candidate_15/OPEN_RESULTS.html) |
| BraVa几何、ROI与四端口芯的生成来源 | [开发与使用指南](vascular_printing/开发与使用指南.md) |
| 小鼠血管A、ROI、1D/0D及边界设计 | [ulm_3D_vascular](ulm_3D_vascular/) |
| 已完成的小鼠ROI-only 3D流场 | [算例](ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-ROI-only-balanced-pressure-v1/) |
| 已完成的小鼠1500条微泡、四视角及停止原因分析 | [动画](ulm_particle_formal_p9a5/particle_3d/reports/microbubble_roi_only_dt0p5ms_n1500/OPEN_RESULTS.html) · [停止分析](ulm_particle_formal_p9a5/particle_3d/reports/microbubble_roi_only_dt0p5ms_n1500/stopping_analysis/MICROBUBBLE_STOP_ANALYSIS_ZH.md) |
| 网格、WSS审计、会议疑问图与旋转绘图代码 | [formal_3D_flow_solver](formal_3D_flow_solver/README.md) |
| RBC代码、轨迹及独立模型展示 | [保留结果索引](ulm_particle_formal_p9a5/CURRENT_RESULTS.md) |
| 当前服务器路径 | [路径表](CURRENT_SERVER_PATHS.md) |
| 此次同步范围、排除项、哈希核验 | [同步报告](sync_metadata/current_20260929/SYNC_REPORT_ZH.md) |

HTML需下载后在浏览器中打开。人脑与小鼠输入、压力、粒径样本和完成状态分别记录，不能混用；RBC结果不表示已完成人脑新场的RBC耦合计算。

本轮只收集和核验，不运行CFD或轨迹。未上传运行环境、可重建二进制和不重要的大型缓存；路径、大小、哈希与恢复方式见[排除说明](sync_metadata/current_20260929/EXCLUSIONS.md)。上游缺失的LFS载荷仅记录指针，不能当作实际数据。部分入口仍绑定原WSL/服务器路径和科学身份契约，快照不是可在任意机器直接启动的安装包。

BraVa流场的边界分流接近各三分之一；已确认内部截面存在约百分之几的局部速度通量缺陷。单网格、边界守恒与残差通过均不证明局部WSS或微泡输运精度。旧细档真实血管CFD、出口压力敏感性及本次BraVa微泡均保持各自的用户取消状态。

上次快照详细入口保存在[2026-09-28说明](sync_metadata/current_20260929/README_before_sync.md)。本次以`ca0ae424d01717f6a231bbea91405cbea8f095a8`为父提交，正常追加到指定分支，不改写历史。
