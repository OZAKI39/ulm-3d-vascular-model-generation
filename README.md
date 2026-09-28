# 血管几何、3D流场、微泡与RBC工作流

2026-09-29增量同步至 `sync/roi-flow-mb1500-wss-stop-analysis-20260928`。原有小鼠流程完整保留；最新新增BraVa人脑打印通道的18 mL/min流场及其动画。**BraVa微泡批次已按用户要求停止，未完成1500条，未生成微泡动画。**

| 内容 | 入口 |
|---|---|
| BraVa最新流场、数值限制与微泡停止状态 | [当前结果说明](brava_flow_roi_18mlmin/CURRENT_RESULTS_ZH.md) |
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
