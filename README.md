# 小鼠微血管：最新流场、微泡轨迹及数值证据

2026-09-28 同步快照，分支 `sync/roi-flow-mb1500-wss-stop-analysis-20260928`。当前主线为 **ROI-only-balanced-pressure-v1 流场 → 1500 条有限尺寸微泡轨迹，dt=0.5 ms → 动画及停止原因分析**。

这是多个实际开发目录及服务器证据的整合快照。最新管网和 CFD 开发来自上一同步工作树的实际文件，最新微泡、WSS 和动画来自各自当前工作区；未用旧 Git HEAD 覆盖未提交开发。

## 当前结果入口

| 内容 | 仓库路径 |
|---|---|
| ROI-only 边界设计、0D 与真实 3D 结果对照 | [中文报告](ulm_3D_vascular/reports/roi_only_boundary_balance_v1/ROI_ONLY_BOUNDARY_BALANCE_FORWARD_REPORT_ZH.md) |
| 血管 A、ROI、1D/0D 代码 | [ulm_3D_vascular](ulm_3D_vascular/)，重点 `network_1d0d/roi_only_hydraulics.py`、`roi_boundary_design.py` |
| 最新 3D FEM 网格、配置、日志及冻结场 | [ROI-only 算例](ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-ROI-only-balanced-pressure-v1/) |
| 表面导入、网格生成和求解器支持 | [formal_3D_flow_solver](formal_3D_flow_solver/README.md) |
| 流线、速度、压力和 WSS 动画 | [最新流场可视化](formal_3D_flow_solver/FEM_SimVascular/rotate_visualization/roi_only_balance_v1/OPEN_RESULTS.html) |
| 原始面片 WSS 旋转展示 | [raw_wss](formal_3D_flow_solver/FEM_SimVascular/rotate_visualization/roi_only_balance_v1/raw_wss/) |
| 应变率与剪切率展示 | [strain_shear_rate](formal_3D_flow_solver/FEM_SimVascular/rotate_visualization/roi_only_balance_v1/strain_shear_rate/) |
| J1 的法向、速度梯度、黏性牵引和切向 WSS 分步展示 | [j1_wss_pipeline](formal_3D_flow_solver/FEM_SimVascular/rotate_visualization/roi_only_balance_v1/j1_wss_pipeline/) |
| 1500 条微泡、英文动画、主体放大30% | [当前微泡入口](ulm_particle_formal_p9a5/particle_3d/reports/microbubble_roi_only_dt0p5ms_n1500/OPEN_RESULTS.html) |
| 微泡另外三个视角 | [additional_views](ulm_particle_formal_p9a5/particle_3d/reports/microbubble_roi_only_dt0p5ms_n1500/additional_views/OPEN_RESULTS.html) |
| 微泡在血管内停止的原因、逐条证据与图件 | [中文分析报告](ulm_particle_formal_p9a5/particle_3d/reports/microbubble_roi_only_dt0p5ms_n1500/stopping_analysis/MICROBUBBLE_STOP_ANALYSIS_ZH.md)、[图文入口](ulm_particle_formal_p9a5/particle_3d/reports/microbubble_roi_only_dt0p5ms_n1500/stopping_analysis/OPEN_RESULTS.html) |
| WSS 审计、两档网格敏感性及会议疑问位置图 | [wss_validation_v2](formal_3D_flow_solver/FEM_SimVascular/wss_validation_v2/)、[第一轮审计](formal_3D_flow_solver/FEM_SimVascular/wss_audit/WSS_AUDIT_REPORT.md) |
| RBC 代码、轨迹与独立模型展示 | [结果索引](ulm_particle_formal_p9a5/CURRENT_RESULTS.md)、[HemoCell 副本](ulm_particle_formal_p9a5/server_evidence/hemocell_restore/) |
| 当前服务器路径 | [CURRENT_SERVER_PATHS.md](CURRENT_SERVER_PATHS.md) |
| 同步范围、排除项、实际检查及复核方法 | [同步报告](sync_metadata/current_20260928/SYNC_REPORT_ZH.md) |

下载后用浏览器打开 HTML；GitHub 的文件页不会直接执行本地预览。

当前流场 SHA256：`fb3c6ad0815d156b09b289fc347766486f24ddb521b8e51e8bbc4cf1fd2ace12`。3D 出口流量比例为 O1=34.5037%、O2=27.2729%、O3=38.2234%，与 0D 设计的等分目标需分别解释。

微泡出口数量为 O1=787、O2=61、O3=467；1315 条完成出流，185 条接触支持静止，执行失败0。微泡数量比例不能替代体积流量比例。停止分析复核了全部185个末次位置，表明刚性球体近壁阻力与多面约束在当前模型中导致静止；尚非真实生理永久滞留证明。动画在终止后继续保持这些微泡的末坐标。

## 范围与使用限制

保留最新1500条的原始轨迹、完整逐轨迹审计、汇总导出、源码身份契约、主动画及附加视角。重要的 `gpu_mesh_input.npz` 虽超过50MiB也予以保留。旧 H0、上一版 best-feasible-balance 及独立 RBC 是有明确版本的参考，不是当前默认输入。

未上传运行环境、Git 元数据库、已安装 SimVascular 分发包、可重建流线候选缓存和已取消细档的三个大网格文件；排除对象的路径、原因及可得的哈希见 [排除说明](sync_metadata/current_20260928/EXCLUSIONS.md)。细档真实血管 CFD 和出口压力敏感性仍为用户取消，不会因本次同步恢复。

本次同步未运行新的 CFD 或微泡积分。服务端代码/日志只读收集；CPU/GPU 历史运行证据保留于原结果。RBC 为各自保留的独立模型，不宣称新流场上已有 RBC 耦合生产计算。

本快照不是自动安装包。部分科学入口仍包含原 WSL/服务器绝对路径和受保护身份契约；克隆到其他目录需按各模块复现说明配置路径、Python、MPI/PETSc/CUDA 等依赖，不能静默批量替换受保护源码。目录导航中的当前算例链接已调整为仓库内相对路径。

父提交：`69c6af94eabb87ba6db0823134c4b92b30a069db`。原分支和原工作区不因本次同步改变。
