# 当前三维血流与粒子研究目录

本索引于 2026-09-23 本地清理时建立。现有主工程路径、科学代码、数据、日志、环境及 Git 工作树保持原位置。

| 环节 | 文件夹 | 用途 |
|---|---|---|
| 原始血管数据与来源 | [ulm_3D_vascular](ulm_3D_vascular/) | 原始数据、血管重建、预处理、FEM 来源几何与参考条件 |
| 血管来源审计 | [vascular_migration_audit](vascular_migration_audit/) | 上游血管数据迁移与来源记录 |
| 几何与冻结交接仓库 | [ulm-3d-vascular-model-generation](ulm-3d-vascular-model-generation/) | 几何生成与冻结 FEM；流场、粒子工作树共用其 Git 存储，不可单独删除 |
| 原始 FEM 与 SimVascular | [formal_3D_flow_solver](formal_3D_flow_solver/) | 历史 FEM 来源、原生 SimVascular 开发、构建、环境、求解、日志、结果和服务器记录 |
| 新 2.0 mm/s FEM 流场 | [ulm_flow_mean_2p0_mmps](ulm_flow_mean_2p0_mmps/) | 新算例、原生输出、检查点、稳态场、压力/WSS/残差、Excel 数据 |
| 微泡与 RBC 全流程 | [ulm_particle_3d_particle0](ulm_particle_3d_particle0/) | Particle 源码、永久测试、环境、轨迹、阶段报告、动画与服务器作业记录 |
| SonoVue 粒径采样 | [sonovue_size_distribution_v0](sonovue_size_distribution_v0/) | 当前微泡粒径分布与采样依赖；根目录 FROZEN_SONOVUE_HISTOGRAM.csv 为来源数据 |
| SonoVue Git 同步仓库 | [github_sync/ulm_sonovue_20260920_111258](github_sync/ulm_sonovue_20260920_111258/) | 保留粒径采样模块对应的版本历史 |
| FEM 历史二维参考 | [ulm_microbubble_traj_gen_2D](ulm_microbubble_traj_gen_2D/) | 仅保留 FEM 来源清单引用的 7 个文件；不是完整二维运行工程 |

## 常用入口

- 粒子计算代码：[particle_3d/src/particle_3d](ulm_particle_3d_particle0/particle_3d/src/particle_3d/)
- 运行及渲染入口：[particle_3d/scripts](ulm_particle_3d_particle0/particle_3d/scripts/)
- 500 条微泡轨迹与 PPT：[particle8_2a_ppt](ulm_particle_3d_particle0/particle_3d/outputs/particle8_2a_ppt/OPEN_RESULTS.html)
- 最新正常 RBC–MB 共流旋转动画：[rbc_mb_flow_rotation](ulm_particle_3d_particle0/particle_3d/reports/rbc_mb_flow_rotation/OPEN_RESULTS.html)
- 最新流线、局部矢量、压力与 WSS：[rotate_visualization](formal_3D_flow_solver/FEM_SimVascular/rotate_visualization/OPEN_RESULTS.html)
- 新流场完整算例：[mean-2p0-mmps](ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps/)
- 残差 Excel：[field_diagnostics/excel](ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps/field_diagnostics/excel/)
- 学位论文研究整理：[THESIS_SUMMARY](ulm_particle_3d_particle0/THESIS_SUMMARY/)

## 流场版本

最新流场图采用 2.0 mm/s FEM 算例；500 条微泡轨迹仍采用约 0.352841 mm/s 的旧冻结 FEM 场；最新正常 RBC–MB 共流旋转示例使用理想化 Poiseuille 场。三者没有在本次清理中互相替换。

## 服务器记录

当前 FEM/Particle 的远程构建、求解、部署、运行、下载校验和服务器来源记录保留在上述主工程的 scripts、reports、logs、outputs 及 flow_cases 中。本次只删除本地目录；不连接或修改服务器端文件，也不改动 /home/lzy/.ssh。

保留根目录两份球体近壁流动参考 PDF，以及当前工作区 .vscode 设置。

## 清理证据

本次保留/删除范围、逐文件清单、SHA-256 复核与运行验证位于 [CLEANUP_AUDIT_20260923_174727](CLEANUP_AUDIT_20260923_174727/)。
