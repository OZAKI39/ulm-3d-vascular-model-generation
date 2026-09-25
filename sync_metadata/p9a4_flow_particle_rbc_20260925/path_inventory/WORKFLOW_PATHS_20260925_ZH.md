# 血管流场、微泡与 RBC 全流程路径索引（2026-09-25）

核查时间 UTC：本地 2026-09-25T00:11:38.182949+00:00；服务器 2026-09-25T00:11:37.211556+00:00。本地路径逐项存在性检查，服务器通过 SSH 只读扫描，主机 `f7c62a262077`，连接 `ssh -p 4159 root@50.115.148.16`（本地别名 `vast4090`）。本次仅新增路径索引和扫描记录，没有运行计算、修改科学代码或移动结果。

**先区分版本**：当前正式场是 NEW Network-H0 outlet-pressure 场，SHA `064cbd28f3efa72f426fc946b2f29da21f056c596609095e7283d39070aa55f4`。最新入口开发在独立 P9-A.4 worktree，100k source audit与30泡smoke已完成；本轮未推送的新增P9-A.4代码/结果只在新worktree及对应服务器目录，不能假定已包含于GitHub同步分支。P9-A.4状态为READY_FOR_LARGE_SAMPLE_REVIEW，正式500未运行。

原Particle工作树的默认frozen_reference仍是OLD equal-pressure 2 mm/s场（SHA `129ebb77396550b300b20ce29cde7b03919e6037b7a81de60fe75d6c6efb616d`）；Network配对runner及新P9-A.4入口显式选择NEW。更早P8.2A/PPT500使用更早Stage Q冻结场，不能与新H0结果混用。当前rotate_visualization的输入manifest已经绑定NEW H0，尽管其中兼容文件名仍写mean_2p0_mmps。

RBC与微泡共用particle_3d，不存在单独的当前RBC顶层仓库。已有RBC–MB flow_rotation演示采用理想化Poiseuille场；不是NEW Network-H0 RBC生产轨迹。rbc_mb_hydrodynamic有用户停止标记。服务器另有历史HemoCell RBC和LAMMPS/Palabos工作目录，下文单列；本次只核查其路径，不替这些历史模型作新的有效性验收。

**清单范围**：335条关键路径；本地目录扫描6471条、服务器目录扫描3289条。目录扫描最多向各工程根下展开5层，环境、Git对象、第三方源码/构建和大原始数据内部不作递归展开；不是全磁盘逐文件清单。RBC相关命名额外检查至8层，共712条，包含多个部署副本。

完整表格：[关键路径CSV](workflow_path_inventory_20260925/key_paths.csv)、[本地目录CSV](workflow_path_inventory_20260925/local_directories.csv)、[服务器目录CSV](workflow_path_inventory_20260925/remote_directories.csv)、[远端RBC相关路径](workflow_path_inventory_20260925/remote_rbc_paths.txt)。需要逐个文件时可查最新同步快照的local_file_inventory.csv / server_file_inventory.json，再加P9-A.4的delivery_manifest.json；这些是各自快照时间的清单，不冒充当前全盘实时清单。

目录中的 `src/`、`scripts/` 是源码/入口；`contracts/` 是模型与输入合同；`data/`、`outputs/`、`run/` 是数据/轨迹/原生结果；`reports/`、`logs/` 保存报告及执行记录；`figures/`、`animations/` 保存展示。具体以各阶段目录为准。


## 当前 P9-A.4

| 用途 | 已核对路径 | 状态 |
| --- | --- | --- |
| 全新开发worktree | [/home/lzy/projects/ulm_particle_population_inlet_p9a4](</home/lzy/projects/ulm_particle_population_inlet_p9a4>) | 存在 |
| 新入口源码 | [/home/lzy/projects/ulm_particle_population_inlet_p9a4/particle_3d/src/particle_3d/continuous_infusion.py](</home/lzy/projects/ulm_particle_population_inlet_p9a4/particle_3d/src/particle_3d/continuous_infusion.py>) | 存在 |
| NEW流场适配与worker | [/home/lzy/projects/ulm_particle_population_inlet_p9a4/particle_3d/src/particle_3d/population_inlet_p9a4.py](</home/lzy/projects/ulm_particle_population_inlet_p9a4/particle_3d/src/particle_3d/population_inlet_p9a4.py>) | 存在 |
| 新浓度与源合同 | [/home/lzy/projects/ulm_particle_population_inlet_p9a4/particle_3d/contracts/P9A4_CONTINUOUS_INFUSION_V1.json](</home/lzy/projects/ulm_particle_population_inlet_p9a4/particle_3d/contracts/P9A4_CONTINUOUS_INFUSION_V1.json>) | 存在 |
| 新永久测试 | [/home/lzy/projects/ulm_particle_population_inlet_p9a4/particle_3d/tests/particle9a4_population_inlet](</home/lzy/projects/ulm_particle_population_inlet_p9a4/particle_3d/tests/particle9a4_population_inlet>) | 存在 |
| P9-A.4全部报告及结果 | [/home/lzy/projects/ulm_particle_population_inlet_p9a4/particle_3d/reports/particle9a4_population_inlet](</home/lzy/projects/ulm_particle_population_inlet_p9a4/particle_3d/reports/particle9a4_population_inlet>) | 存在 |
| 总预览 | [/home/lzy/projects/ulm_particle_population_inlet_p9a4/particle_3d/reports/particle9a4_population_inlet/OPEN_RESULTS.html](</home/lzy/projects/ulm_particle_population_inlet_p9a4/particle_3d/reports/particle9a4_population_inlet/OPEN_RESULTS.html>) | 存在 |
| 100k全部source ledger | [/home/lzy/projects/ulm_particle_population_inlet_p9a4/particle_3d/reports/particle9a4_population_inlet/data/inlet100k/proposal_ledger.jsonl](</home/lzy/projects/ulm_particle_population_inlet_p9a4/particle_3d/reports/particle9a4_population_inlet/data/inlet100k/proposal_ledger.jsonl>) | 存在 |
| accepted births | [/home/lzy/projects/ulm_particle_population_inlet_p9a4/particle_3d/reports/particle9a4_population_inlet/data/inlet100k/accepted_births.json](</home/lzy/projects/ulm_particle_population_inlet_p9a4/particle_3d/reports/particle9a4_population_inlet/data/inlet100k/accepted_births.json>) | 存在 |
| 审计数据与机器摘要 | [/home/lzy/projects/ulm_particle_population_inlet_p9a4/particle_3d/reports/particle9a4_population_inlet/data](</home/lzy/projects/ulm_particle_population_inlet_p9a4/particle_3d/reports/particle9a4_population_inlet/data>) | 存在 |
| 生成、部署、审核与绘图脚本 | [/home/lzy/projects/ulm_particle_population_inlet_p9a4/particle_3d/reports/particle9a4_population_inlet/scripts](</home/lzy/projects/ulm_particle_population_inlet_p9a4/particle_3d/reports/particle9a4_population_inlet/scripts>) | 存在 |
| 本地/回传日志 | [/home/lzy/projects/ulm_particle_population_inlet_p9a4/particle_3d/reports/particle9a4_population_inlet/logs](</home/lzy/projects/ulm_particle_population_inlet_p9a4/particle_3d/reports/particle9a4_population_inlet/logs>) | 存在 |
| 30泡与point完整输出 | [/home/lzy/projects/ulm_particle_population_inlet_p9a4/particle_3d/reports/particle9a4_population_inlet/outputs/smoke30](</home/lzy/projects/ulm_particle_population_inlet_p9a4/particle_3d/reports/particle9a4_population_inlet/outputs/smoke30>) | 存在 |
| 30泡轨迹 | [/home/lzy/projects/ulm_particle_population_inlet_p9a4/particle_3d/reports/particle9a4_population_inlet/outputs/smoke30/trajectories](</home/lzy/projects/ulm_particle_population_inlet_p9a4/particle_3d/reports/particle9a4_population_inlet/outputs/smoke30/trajectories>) | 存在 |
| 对应point轨迹 | [/home/lzy/projects/ulm_particle_population_inlet_p9a4/particle_3d/reports/particle9a4_population_inlet/outputs/smoke30/point](</home/lzy/projects/ulm_particle_population_inlet_p9a4/particle_3d/reports/particle9a4_population_inlet/outputs/smoke30/point>) | 存在 |
| 8组PNG/PDF图件 | [/home/lzy/projects/ulm_particle_population_inlet_p9a4/particle_3d/reports/particle9a4_population_inlet/figures](</home/lzy/projects/ulm_particle_population_inlet_p9a4/particle_3d/reports/particle9a4_population_inlet/figures>) | 存在 |
| P9-A.4使用的NEW流场副本 | [/home/lzy/projects/ulm_particle_population_inlet_p9a4/particle_3d/reports/network_derived_flow_mb_validation_v1/server_bundle/inputs/NEW.vtu](</home/lzy/projects/ulm_particle_population_inlet_p9a4/particle_3d/reports/network_derived_flow_mb_validation_v1/server_bundle/inputs/NEW.vtu>) | 存在 |

## 当前 P9-A.4 集成工作树内容

| 用途 | 已核对路径 | 状态 |
| --- | --- | --- |
| vascular_network | [/home/lzy/projects/ulm_particle_population_inlet_p9a4/vascular_network](</home/lzy/projects/ulm_particle_population_inlet_p9a4/vascular_network>) | 存在 |
| particle_3d | [/home/lzy/projects/ulm_particle_population_inlet_p9a4/particle_3d](</home/lzy/projects/ulm_particle_population_inlet_p9a4/particle_3d>) | 存在 |
| formal_3D_flow_solver | [/home/lzy/projects/ulm_particle_population_inlet_p9a4/formal_3D_flow_solver](</home/lzy/projects/ulm_particle_population_inlet_p9a4/formal_3D_flow_solver>) | 存在 |
| sonovue_size_distribution_v0 | [/home/lzy/projects/ulm_particle_population_inlet_p9a4/sonovue_size_distribution_v0](</home/lzy/projects/ulm_particle_population_inlet_p9a4/sonovue_size_distribution_v0>) | 存在 |
| server_evidence | [/home/lzy/projects/ulm_particle_population_inlet_p9a4/server_evidence](</home/lzy/projects/ulm_particle_population_inlet_p9a4/server_evidence>) | 存在 |
| source_dependencies | [/home/lzy/projects/ulm_particle_population_inlet_p9a4/source_dependencies](</home/lzy/projects/ulm_particle_population_inlet_p9a4/source_dependencies>) | 存在 |
| THESIS_SUMMARY | [/home/lzy/projects/ulm_particle_population_inlet_p9a4/THESIS_SUMMARY](</home/lzy/projects/ulm_particle_population_inlet_p9a4/THESIS_SUMMARY>) | 存在 |

## 当前 vascular A → Network H0 → FEM

| 用途 | 已核对路径 | 状态 |
| --- | --- | --- |
| vascular A/ROI原始主工程 | [/home/lzy/projects/ulm_3D_vascular](</home/lzy/projects/ulm_3D_vascular>) | 存在 |
| 血管原数据 | [/home/lzy/projects/ulm_3D_vascular/vessel_model](</home/lzy/projects/ulm_3D_vascular/vessel_model>) | 存在 |
| SWC/ROI/表面数据 | [/home/lzy/projects/ulm_3D_vascular/outputs](</home/lzy/projects/ulm_3D_vascular/outputs>) | 存在 |
| 预处理共享代码 | [/home/lzy/projects/ulm_3D_vascular/vascular_processing](</home/lzy/projects/ulm_3D_vascular/vascular_processing>) | 存在 |
| Network 1D/0D模块 | [/home/lzy/projects/ulm_3D_vascular/network_1d0d](</home/lzy/projects/ulm_3D_vascular/network_1d0d>) | 存在 |
| H0求解入口 | [/home/lzy/projects/ulm_3D_vascular/scripts/run_a_network_h0_v2.py](</home/lzy/projects/ulm_3D_vascular/scripts/run_a_network_h0_v2.py>) | 存在 |
| H0 FEM构建 | [/home/lzy/projects/ulm_3D_vascular/network_1d0d/fem_h0_case.py](</home/lzy/projects/ulm_3D_vascular/network_1d0d/fem_h0_case.py>) | 存在 |
| H0远程求解 | [/home/lzy/projects/ulm_3D_vascular/scripts/solve_a_h0_fem_remote.py](</home/lzy/projects/ulm_3D_vascular/scripts/solve_a_h0_fem_remote.py>) | 存在 |
| Network v1报告 | [/home/lzy/projects/ulm_3D_vascular/reports/a_network_1d0d_boundary_v1](</home/lzy/projects/ulm_3D_vascular/reports/a_network_1d0d_boundary_v1>) | 存在 |
| H0 v2报告/数据/日志 | [/home/lzy/projects/ulm_3D_vascular/reports/a_network_1d0d_boundary_v2_idealized](</home/lzy/projects/ulm_3D_vascular/reports/a_network_1d0d_boundary_v2_idealized>) | 存在 |
| NEW H0完整FEM算例 | [/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-H0-pressure-v1](</home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-H0-pressure-v1>) | 存在 |
| NEW H0 frozen_flow | [/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-H0-pressure-v1/frozen_flow](</home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-H0-pressure-v1/frozen_flow>) | 存在 |
| NEW H0最终场 | [/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-H0-pressure-v1/frozen_flow/steady_flow_mean_2p0_mmps_A_H0.vtu](</home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-H0-pressure-v1/frozen_flow/steady_flow_mean_2p0_mmps_A_H0.vtu>) | 存在 |
| NEW H0网格 | [/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-H0-pressure-v1/SV_MESH](</home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-H0-pressure-v1/SV_MESH>) | 存在 |
| NEW H0原生运行/日志/检查点 | [/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-H0-pressure-v1/run](</home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-H0-pressure-v1/run>) | 存在 |
| NEW H0求解配置 | [/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-H0-pressure-v1/run/solver.xml](</home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-H0-pressure-v1/run/solver.xml>) | 存在 |
| NEW H0验收与通量报告 | [/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-H0-pressure-v1/reports](</home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-H0-pressure-v1/reports>) | 存在 |
| 共享FEM测量代码 | [/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/src/sv_validation](</home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/src/sv_validation>) | 存在 |
| 已停止Taylor-Hood审核 | [/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/reports/taylor_hood_p2p1_validation_v1](</home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/reports/taylor_hood_p2p1_validation_v1>) | 存在 |

## Network-H0 配对验证与冻结动力学

| 用途 | 已核对路径 | 状态 |
| --- | --- | --- |
| Network OLD/NEW配对验证 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/network_derived_flow_mb_validation_v1](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/network_derived_flow_mb_validation_v1>) | 存在 |
| 配对runner/analyze | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/network_derived_flow_mb_validation_v1/scripts](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/network_derived_flow_mb_validation_v1/scripts>) | 存在 |
| OLD/NEW输入与源码快照 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/network_derived_flow_mb_validation_v1/server_bundle](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/network_derived_flow_mb_validation_v1/server_bundle>) | 存在 |
| 流场身份与配对事件 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/network_derived_flow_mb_validation_v1/data](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/network_derived_flow_mb_validation_v1/data>) | 存在 |
| 配对轨迹输出 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/network_derived_flow_mb_validation_v1/outputs](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/network_derived_flow_mb_validation_v1/outputs>) | 存在 |
| 配对日志 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/network_derived_flow_mb_validation_v1/logs](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/network_derived_flow_mb_validation_v1/logs>) | 存在 |
| P9-A.1动力学 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/src/particle_3d/particle9a_motion.py](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/src/particle_3d/particle9a_motion.py>) | 存在 |

## 微泡和 RBC 共用计算工程

| 用途 | 已核对路径 | 状态 |
| --- | --- | --- |
| 整个工作树 | [/home/lzy/projects/ulm_particle_3d_particle0](</home/lzy/projects/ulm_particle_3d_particle0>) | 存在 |
| Particle 主目录 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d>) | 存在 |
| 全部运动与插值源码 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/src/particle_3d](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/src/particle_3d>) | 存在 |
| 入口、部署、计算、导出和渲染脚本 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/scripts](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/scripts>) | 存在 |
| 永久测试 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/tests](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/tests>) | 存在 |
| 粒径、RBC 几何、近场及人口模型契约 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/contracts](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/contracts>) | 存在 |
| 正式输出与轨迹 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/outputs](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/outputs>) | 存在 |
| 各阶段数据、图、审核报告及日志 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports>) | 存在 |
| 实际科学 Python 环境 | [/home/lzy/projects/ulm_particle_3d_particle0/.venv](</home/lzy/projects/ulm_particle_3d_particle0/.venv>) | 存在 |
| 解释器 | [/home/lzy/projects/ulm_particle_3d_particle0/.venv/bin/python](</home/lzy/projects/ulm_particle_3d_particle0/.venv/bin/python>) | 存在 |
| 学位论文研究整理 | [/home/lzy/projects/ulm_particle_3d_particle0/THESIS_SUMMARY](</home/lzy/projects/ulm_particle_3d_particle0/THESIS_SUMMARY>) | 存在 |

## 微泡计算与轨迹关键位置

| 用途 | 已核对路径 | 状态 |
| --- | --- | --- |
| 冻结 FEM 插值与梯度 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/src/particle_3d/field.py](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/src/particle_3d/field.py>) | 存在 |
| 微泡状态 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/src/particle_3d/microbubble.py](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/src/particle_3d/microbubble.py>) | 存在 |
| 微泡基础运动更新 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/src/particle_3d/integrator.py](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/src/particle_3d/integrator.py>) | 存在 |
| 水动力阻力 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/src/particle_3d/hydrodynamic_resistance.py](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/src/particle_3d/hydrodynamic_resistance.py>) | 存在 |
| 近壁正则化 V1 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/src/particle_3d/nearfield_regularization.py](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/src/particle_3d/nearfield_regularization.py>) | 存在 |
| Particle-6.5 运动组装 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/src/particle_3d/particle65_motion.py](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/src/particle_3d/particle65_motion.py>) | 存在 |
| 多出口轨迹积分 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/src/particle_3d/particle82_integration.py](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/src/particle_3d/particle82_integration.py>) | 存在 |
| 新准入阶段轨迹积分 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/src/particle_3d/particle82a_integration.py](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/src/particle_3d/particle82a_integration.py>) | 存在 |
| Particle-8.1 原始正式轨迹 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/outputs/particle8_1/trajectories](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/outputs/particle8_1/trajectories>) | 存在 |
| Particle-8.1 运行日志 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/outputs/particle8_1/logs](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/outputs/particle8_1/logs>) | 存在 |
| Particle-8.2 结果、服务器命令、日志和来源 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/particle8_2](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/particle8_2>) | 存在 |
| Particle-8.2A 正式轨迹 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/outputs/particle8_2a/trajectories](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/outputs/particle8_2a/trajectories>) | 存在 |
| Particle-8.2A 准入缓存与中间计算 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/outputs/particle8_2a/admission](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/outputs/particle8_2a/admission>) | 存在 |
| Particle-8.2A 审核数据 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/outputs/particle8_2a/review](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/outputs/particle8_2a/review>) | 存在 |
| Particle-8.2A 部署与服务器记录 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/particle8_2a](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/particle8_2a>) | 存在 |
| 500 条轨迹 PPT 交付目录 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/outputs/particle8_2a_ppt](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/outputs/particle8_2a_ppt>) | 存在 |
| 500 条轨迹的数据 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/outputs/particle8_2a_ppt/data](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/outputs/particle8_2a_ppt/data>) | 存在 |
| 500 条轨迹的动画 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/outputs/particle8_2a_ppt/animations](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/outputs/particle8_2a_ppt/animations>) | 存在 |
| 500 条轨迹的浏览入口 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/outputs/particle8_2a_ppt/OPEN_RESULTS.html](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/outputs/particle8_2a_ppt/OPEN_RESULTS.html>) | 存在 |

## RBC 源码与验证入口

| 用途 | 已核对路径 | 状态 |
| --- | --- | --- |
| normal_rbc_encounter.py | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/src/particle_3d/normal_rbc_encounter.py](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/src/particle_3d/normal_rbc_encounter.py>) | 存在 |
| rbc.py | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/src/particle_3d/rbc.py](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/src/particle_3d/rbc.py>) | 存在 |
| rbc_capillary_surrogate.py | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/src/particle_3d/rbc_capillary_surrogate.py](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/src/particle_3d/rbc_capillary_surrogate.py>) | 存在 |
| rbc_distribution.py | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/src/particle_3d/rbc_distribution.py](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/src/particle_3d/rbc_distribution.py>) | 存在 |
| rbc_integrator.py | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/src/particle_3d/rbc_integrator.py](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/src/particle_3d/rbc_integrator.py>) | 存在 |
| rbc_mb_encounter.py | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/src/particle_3d/rbc_mb_encounter.py](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/src/particle_3d/rbc_mb_encounter.py>) | 存在 |
| rbc_orientation.py | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/src/particle_3d/rbc_orientation.py](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/src/particle_3d/rbc_orientation.py>) | 存在 |
| generate_rbc_population.py | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/scripts/generate_rbc_population.py](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/scripts/generate_rbc_population.py>) | 存在 |
| run_particle2_validation.py | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/scripts/run_particle2_validation.py](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/scripts/run_particle2_validation.py>) | 存在 |
| run_particle3_real_validation.py | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/scripts/run_particle3_real_validation.py](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/scripts/run_particle3_real_validation.py>) | 存在 |
| prepare_rbc_mb_flow_rotation.py | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/scripts/prepare_rbc_mb_flow_rotation.py](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/scripts/prepare_rbc_mb_flow_rotation.py>) | 存在 |
| render_rbc_mb_flow_rotation.py | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/scripts/render_rbc_mb_flow_rotation.py](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/scripts/render_rbc_mb_flow_rotation.py>) | 存在 |
| verify_rbc_reproducibility.py | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/scripts/verify_rbc_reproducibility.py](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/scripts/verify_rbc_reproducibility.py>) | 存在 |
| 测试 rbc_mb_coflow | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/tests/rbc_mb_coflow](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/tests/rbc_mb_coflow>) | 存在 |
| 测试 rbc_mb_flow_rotation | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/tests/rbc_mb_flow_rotation](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/tests/rbc_mb_flow_rotation>) | 存在 |
| 测试 rbc_mb_hydrodynamic | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/tests/rbc_mb_hydrodynamic](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/tests/rbc_mb_hydrodynamic>) | 存在 |
| 测试 rbc_mb_interaction | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/tests/rbc_mb_interaction](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/tests/rbc_mb_interaction>) | 存在 |
| 测试 rbc_mb_normal | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/tests/rbc_mb_normal](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/tests/rbc_mb_normal>) | 存在 |

## RBC 运动、形状与展示

| 用途 | 已核对路径 | 状态 |
| --- | --- | --- |
| RBC 状态与几何 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/src/particle_3d/rbc.py](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/src/particle_3d/rbc.py>) | 存在 |
| RBC 几何分布 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/src/particle_3d/rbc_distribution.py](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/src/particle_3d/rbc_distribution.py>) | 存在 |
| RBC 姿态与四元数 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/src/particle_3d/rbc_orientation.py](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/src/particle_3d/rbc_orientation.py>) | 存在 |
| RBC 运动更新 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/src/particle_3d/rbc_integrator.py](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/src/particle_3d/rbc_integrator.py>) | 存在 |
| RBC 毛细血管变形代理 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/src/particle_3d/rbc_capillary_surrogate.py](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/src/particle_3d/rbc_capillary_surrogate.py>) | 存在 |
| 粒子形状 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/src/particle_3d/particle_shapes.py](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/src/particle_3d/particle_shapes.py>) | 存在 |
| 近壁运动 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/src/particle_3d/particle3_motion.py](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/src/particle_3d/particle3_motion.py>) | 存在 |
| 正常 RBC 和 MB 共流旋转模型 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/src/particle_3d/coflow_rotation.py](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/src/particle_3d/coflow_rotation.py>) | 存在 |
| RBC 几何分布契约 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/contracts/C57BL6_RBC_GEOMETRY_DISTRIBUTION_V0.json](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/contracts/C57BL6_RBC_GEOMETRY_DISTRIBUTION_V0.json>) | 存在 |
| 原血管内 RBC 运动验证和数据 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/particle2](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/particle2>) | 存在 |
| 近壁和变形代理验证、数据、图和日志 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/particle3](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/particle3>) | 存在 |
| 早期相遇/交互示例 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/rbc_mb_interaction](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/rbc_mb_interaction>) | 存在 |
| 正常形态 RBC 示例 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/rbc_mb_normal](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/rbc_mb_normal>) | 存在 |
| RBC 与 MB 共流示例 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/rbc_mb_coflow](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/rbc_mb_coflow>) | 存在 |
| 已有正常 RBC 与 MB 共流旋转完整结果 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/rbc_mb_flow_rotation](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/rbc_mb_flow_rotation>) | 存在 |
| 已有 RBC/MB 位置和姿态数组 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/rbc_mb_flow_rotation/data/motion.npz](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/rbc_mb_flow_rotation/data/motion.npz>) | 存在 |
| 已有 RBC/MB 轨迹及旋转 CSV | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/rbc_mb_flow_rotation/data/trajectories_and_rotation.csv](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/rbc_mb_flow_rotation/data/trajectories_and_rotation.csv>) | 存在 |
| 已有 RBC/MB 场景与流场参数 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/rbc_mb_flow_rotation/data/SCENE.json](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/rbc_mb_flow_rotation/data/SCENE.json>) | 存在 |
| 已有 RBC/MB 动画 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/rbc_mb_flow_rotation/animations](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/rbc_mb_flow_rotation/animations>) | 存在 |
| 已有 RBC/MB 测试日志 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/rbc_mb_flow_rotation/tests.log](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/rbc_mb_flow_rotation/tests.log>) | 存在 |
| 已有 RBC/MB 渲染日志 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/rbc_mb_flow_rotation/render.log](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/rbc_mb_flow_rotation/render.log>) | 存在 |
| 已停止的 hydrodynamic 尝试标记 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/rbc_mb_hydrodynamic/STOPPED_BY_USER.md](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/rbc_mb_hydrodynamic/STOPPED_BY_USER.md>) | 存在 |
| particle2 data | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/particle2/data](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/particle2/data>) | 存在 |
| particle2 figures | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/particle2/figures](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/particle2/figures>) | 存在 |
| particle2 logs | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/particle2/logs](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/particle2/logs>) | 存在 |
| particle3 data | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/particle3/data](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/particle3/data>) | 存在 |
| particle3 figures | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/particle3/figures](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/particle3/figures>) | 存在 |
| particle3 logs | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/particle3/logs](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/particle3/logs>) | 存在 |

## 服务器当前主线

| 用途 | 已核对路径 | 状态 |
| --- | --- | --- |
| flow_mean_2p0_mmps_A_H0_20260924T181140Z | `/workspace/flow_mean_2p0_mmps_A_H0_20260924T181140Z` | 存在 |
| particle_network_flow_mb_validation_v1_20260924T210047Z | `/workspace/particle_network_flow_mb_validation_v1_20260924T210047Z` | 存在 |
| particle9a4_population_inlet_20260924T233708Z | `/workspace/particle9a4_population_inlet_20260924T233708Z` | 存在 |
| P9-A.4源代码 | `/workspace/particle9a4_population_inlet_20260924T233708Z/particle_3d/src/particle_3d` | 存在 |
| P9-A.4部署manifest | `/workspace/particle9a4_population_inlet_20260924T233708Z/deployment_manifest.json` | 存在 |
| P9-A.4终端日志 | `/workspace/particle9a4_population_inlet_20260924T233708Z/smoke30.log` | 存在 |
| P9-A.4审计入口、数据及结果 | `/workspace/particle9a4_population_inlet_20260924T233708Z/particle_3d/reports/particle9a4_population_inlet` | 存在 |
| 30泡与point输出 | `/workspace/particle9a4_population_inlet_20260924T233708Z/particle_3d/reports/particle9a4_population_inlet/outputs/smoke30` | 存在 |
| H0算例 SV_MESH | `/workspace/flow_mean_2p0_mmps_A_H0_20260924T181140Z/mean-2p0-mmps-A-H0-pressure-v1/SV_MESH` | 存在 |
| H0算例 reports | `/workspace/flow_mean_2p0_mmps_A_H0_20260924T181140Z/mean-2p0-mmps-A-H0-pressure-v1/reports` | 存在 |
| H0算例 run | `/workspace/flow_mean_2p0_mmps_A_H0_20260924T181140Z/mean-2p0-mmps-A-H0-pressure-v1/run` | 存在 |
| H0算例 run/1-procs | `/workspace/flow_mean_2p0_mmps_A_H0_20260924T181140Z/mean-2p0-mmps-A-H0-pressure-v1/run/1-procs` | 存在 |
| Network验证 particle_3d/src | `/workspace/particle_network_flow_mb_validation_v1_20260924T210047Z/particle_3d/src` | 存在 |
| Network验证 inputs | `/workspace/particle_network_flow_mb_validation_v1_20260924T210047Z/inputs` | 存在 |
| Network验证 scripts | `/workspace/particle_network_flow_mb_validation_v1_20260924T210047Z/scripts` | 存在 |
| Network验证 outputs/NEW | `/workspace/particle_network_flow_mb_validation_v1_20260924T210047Z/outputs/NEW` | 存在 |
| Network验证 outputs/NEW_POINT | `/workspace/particle_network_flow_mb_validation_v1_20260924T210047Z/outputs/NEW_POINT` | 存在 |
| Network验证 logs | `/workspace/particle_network_flow_mb_validation_v1_20260924T210047Z/logs` | 存在 |
| Network验证 final_review | `/workspace/particle_network_flow_mb_validation_v1_20260924T210047Z/final_review` | 存在 |

## 服务器RBC代码副本

| 用途 | 已核对路径 | 状态 |
| --- | --- | --- |
| 当前Particle部署携带的RBC模块；不代表本次执行过RBC轨迹 | `/workspace/particle9a4_population_inlet_20260924T233708Z/particle_3d/src/particle_3d/rbc.py` | 存在 |

## 服务器历史运行/独立模拟工程

| 用途 | 已核对路径 | 状态 |
| --- | --- | --- |
| archives | `/workspace/archives` | 存在 |
| flow_mean_2p0_mmps_A_H0_TaylorHood_20260924T194503Z | `/workspace/flow_mean_2p0_mmps_A_H0_TaylorHood_20260924T194503Z` | 存在 |
| hemocell_restore | `/workspace/hemocell_restore` | 存在 |
| lammps_active | `/workspace/lammps_active` | 存在 |
| lammps_migration | `/workspace/lammps_migration` | 存在 |
| microbubble_lammps | `/workspace/microbubble_lammps` | 存在 |
| particle9a1_2mmps_20260924T001605Z | `/workspace/particle9a1_2mmps_20260924T001605Z` | 存在 |
| particle9a1_routing_stationary_audit_20260924T064607Z | `/workspace/particle9a1_routing_stationary_audit_20260924T064607Z` | 存在 |
| particle9a2_inlet_20260924T081657Z | `/workspace/particle9a2_inlet_20260924T081657Z` | 存在 |
| particle9a3_interior_inlet_20260924T100618Z | `/workspace/particle9a3_interior_inlet_20260924T100618Z` | 存在 |
| particle9a3b_flowfield_conservation_20260924T121004Z | `/workspace/particle9a3b_flowfield_conservation_20260924T121004Z` | 存在 |
| particle9a_2mmps_20260923T222942Z | `/workspace/particle9a_2mmps_20260923T222942Z` | 存在 |
| particle9a_2mmps_diagnosis_20260923T232212Z | `/workspace/particle9a_2mmps_diagnosis_20260923T232212Z` | 存在 |
| reference | `/workspace/reference` | 存在 |

## 服务器历史 HemoCell / LAMMPS 独立工程

| 用途 | 已核对路径 | 状态 |
| --- | --- | --- |
| hemocell_restore archive | `/workspace/hemocell_restore/archive` | 存在 |
| hemocell_restore logs | `/workspace/hemocell_restore/logs` | 存在 |
| hemocell_restore results | `/workspace/hemocell_restore/results` | 存在 |
| hemocell_restore toolchains | `/workspace/hemocell_restore/toolchains` | 存在 |
| hemocell_restore work | `/workspace/hemocell_restore/work` | 存在 |
| microbubble_lammps results | `/workspace/microbubble_lammps/results` | 存在 |
| microbubble_lammps work | `/workspace/microbubble_lammps/work` | 存在 |
| lammps_active build_cpu | `/workspace/lammps_active/build_cpu` | 存在 |
| lammps_active build_gpu | `/workspace/lammps_active/build_gpu` | 存在 |
| lammps_active build_provenance | `/workspace/lammps_active/build_provenance` | 存在 |
| lammps_active cases | `/workspace/lammps_active/cases` | 存在 |
| lammps_active contracts | `/workspace/lammps_active/contracts` | 存在 |
| lammps_active downloads | `/workspace/lammps_active/downloads` | 存在 |
| lammps_active logs | `/workspace/lammps_active/logs` | 存在 |
| lammps_active provenance | `/workspace/lammps_active/provenance` | 存在 |
| lammps_active scripts | `/workspace/lammps_active/scripts` | 存在 |
| lammps_active src | `/workspace/lammps_active/src` | 存在 |
| lammps_active upstream | `/workspace/lammps_active/upstream` | 存在 |
| lammps_active validation | `/workspace/lammps_active/validation` | 存在 |
| HemoCell RBC work/rbc_stage1_20260915_165632 | `/workspace/hemocell_restore/work/rbc_stage1_20260915_165632` | 存在 |
| HemoCell RBC results/rbc_stage1_20260915_165632 | `/workspace/hemocell_restore/results/rbc_stage1_20260915_165632` | 存在 |
| HemoCell RBC logs/minimal_migration_20260915_143631 | `/workspace/hemocell_restore/logs/minimal_migration_20260915_143631` | 存在 |

## 服务器归档

| 用途 | 已核对路径 | 状态 |
| --- | --- | --- |
| 同步/退役代码归档 | `/workspace/archives/github_sync_network_h0_20260924T221559Z` | 存在 |
| 同步/退役代码归档 | `/workspace/archives/vascular_workflow_cleanup_20260924T214356Z` | 存在 |

## 本地原生 FEM / SimVascular 工程

| 用途 | 已核对路径 | 状态 |
| --- | --- | --- |
| 原始完整工程 | [/home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular](</home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular>) | 存在 |
| 验证与后处理代码 | [/home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/src/sv_validation](</home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/src/sv_validation>) | 存在 |
| 求解、远程构建、部署及验证脚本 | [/home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/scripts](</home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/scripts>) | 存在 |
| 配置与运行策略 | [/home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/configs](</home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/configs>) | 存在 |
| 上游输入 | [/home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/inputs](</home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/inputs>) | 存在 |
| SimVascular / PETSc / MPI 等源码、依赖与构建 | [/home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/external](</home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/external>) | 存在 |
| 原生求解器补丁 | [/home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/patches](</home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/patches>) | 存在 |
| 求解结果及检查点 | [/home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/outputs](</home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/outputs>) | 存在 |
| 本地与服务器回传日志 | [/home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/logs](</home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/logs>) | 存在 |
| 阶段审核、性能和来源记录 | [/home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/reports](</home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/reports>) | 存在 |
| 永久测试 | [/home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/tests](</home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/tests>) | 存在 |
| Python 环境 | [/home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/.venv](</home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/.venv>) | 存在 |
| 早期 FEM 完整工程 | [/home/lzy/projects/formal_3D_flow_solver/FEM](</home/lzy/projects/formal_3D_flow_solver/FEM>) | 存在 |
| 早期 FEM 服务器环境和连接记录 | [/home/lzy/projects/formal_3D_flow_solver/FEM/remote](</home/lzy/projects/formal_3D_flow_solver/FEM/remote>) | 存在 |

## 历史 equal-pressure 2.0 mm/s 流场

| 用途 | 已核对路径 | 状态 |
| --- | --- | --- |
| 历史 equal-pressure 流场工作树 | [/home/lzy/projects/ulm_flow_mean_2p0_mmps](</home/lzy/projects/ulm_flow_mean_2p0_mmps>) | 存在 |
| 历史 equal-pressure 流场代码 | [/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/scripts/flow_2mmps](</home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/scripts/flow_2mmps>) | 存在 |
| 历史 equal-pressure 流场测试 | [/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/tests/flow_2mmps](</home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/tests/flow_2mmps>) | 存在 |
| 完整算例 | [/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps](</home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps>) | 存在 |
| 网格、壁面和出入口 | [/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps/SV_MESH](</home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps/SV_MESH>) | 存在 |
| 求解参数和边界条件 | [/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps/run/solver.xml](</home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps/run/solver.xml>) | 存在 |
| PETSc 参数 | [/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps/run/PETSC_OPTIONS.txt](</home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps/run/PETSC_OPTIONS.txt>) | 存在 |
| 原生运行目录 | [/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps/run](</home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps/run>) | 存在 |
| 原生求解日志 | [/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps/run/solver.log](</home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps/run/solver.log>) | 存在 |
| 原始场与检查点 | [/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps/run/1-procs](</home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps/run/1-procs>) | 存在 |
| 导出的最终稳态场 | [/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps/frozen_flow](</home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps/frozen_flow>) | 存在 |
| 最终流场 VTU | [/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps/frozen_flow/steady_flow_mean_2p0_mmps.vtu](</home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps/frozen_flow/steady_flow_mean_2p0_mmps.vtu>) | 存在 |
| 机器记录、守恒、收敛与服务器来源 | [/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps/reports](</home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps/reports>) | 存在 |
| 残差、压力与 WSS 诊断 | [/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps/field_diagnostics](</home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps/field_diagnostics>) | 存在 |
| 残差与 WSS 的 CSV/JSON 数据 | [/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps/field_diagnostics/data](</home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps/field_diagnostics/data>) | 存在 |
| 每幅残差图的 Excel 文件 | [/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps/field_diagnostics/excel](</home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps/field_diagnostics/excel>) | 存在 |
| 诊断静态图 | [/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps/field_diagnostics/figures](</home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps/field_diagnostics/figures>) | 存在 |
| 流线数据与原可视化 | [/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps/streamlines](</home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps/streamlines>) | 存在 |
| 流线数据目录 | [/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps/streamlines/data](</home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps/streamlines/data>) | 存在 |

## 历史粒子使用的冻结 FEM

| 用途 | 已核对路径 | 状态 |
| --- | --- | --- |
| Particle 工作树内的冻结输入 | [/home/lzy/projects/ulm_particle_3d_particle0/formal_3D_flow_solver/FEM_SimVascular/frozen_reference](</home/lzy/projects/ulm_particle_3d_particle0/formal_3D_flow_solver/FEM_SimVascular/frozen_reference>) | 存在 |
| 冻结速度、压力场文件 | `/home/lzy/projects/ulm_particle_3d_particle0/formal_3D_flow_solver/FEM_SimVascular/frozen_reference/flow/steady_flow_stage_sv1_3q.vtu` | **旧索引路径已失效** |
| 主 Git 工作树 | [/home/lzy/projects/ulm-3d-vascular-model-generation](</home/lzy/projects/ulm-3d-vascular-model-generation>) | 存在 |
| 主工作树内的冻结输入副本 | [/home/lzy/projects/ulm-3d-vascular-model-generation/formal_3D_flow_solver/FEM_SimVascular/frozen_reference](</home/lzy/projects/ulm-3d-vascular-model-generation/formal_3D_flow_solver/FEM_SimVascular/frozen_reference>) | 存在 |
| 当前默认frozen_reference内的历史equal-pressure场（OLD SHA） | [/home/lzy/projects/ulm_particle_3d_particle0/formal_3D_flow_solver/FEM_SimVascular/frozen_reference/flow/steady_flow_mean_2p0_mmps.vtu](</home/lzy/projects/ulm_particle_3d_particle0/formal_3D_flow_solver/FEM_SimVascular/frozen_reference/flow/steady_flow_mean_2p0_mmps.vtu>) | 存在 |

## 最新流场旋转可视化

| 用途 | 已核对路径 | 状态 |
| --- | --- | --- |
| 完整独立可视化工程 | [/home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/rotate_visualization](</home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/rotate_visualization>) | 存在 |
| 流线和局部速度矢量渲染 | [/home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/rotate_visualization/render_visualization.py](</home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/rotate_visualization/render_visualization.py>) | 存在 |
| 压力与 WSS 渲染 | [/home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/rotate_visualization/render_surface_fields.py](</home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/rotate_visualization/render_surface_fields.py>) | 存在 |
| 渲染输入数据 | [/home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/rotate_visualization/input_data](</home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/rotate_visualization/input_data>) | 存在 |
| 流线与局部矢量结果 | [/home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/rotate_visualization/results](</home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/rotate_visualization/results>) | 存在 |
| 压力与 WSS 结果 | [/home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/rotate_visualization/results/surface_fields](</home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/rotate_visualization/results/surface_fields>) | 存在 |
| 四种视图统一浏览入口 | [/home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/rotate_visualization/OPEN_RESULTS.html](</home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/rotate_visualization/OPEN_RESULTS.html>) | 存在 |
| 此前可视化版本记录 | [/home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/rotate_visualization_history](</home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/rotate_visualization_history>) | 存在 |

## 独立剪切升力审核

| 用途 | 已核对路径 | 状态 |
| --- | --- | --- |
| 独立审核根目录 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/shear_lift_audit](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/shear_lift_audit>) | 存在 |
| 审核代码 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/shear_lift_audit/code](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/shear_lift_audit/code>) | 存在 |
| 审核测试 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/shear_lift_audit/tests](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/shear_lift_audit/tests>) | 存在 |
| 输入清单及数据 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/shear_lift_audit/data](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/shear_lift_audit/data>) | 存在 |
| 服务器回传结果（本地完整版） | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/shear_lift_audit/data/remote_results](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/shear_lift_audit/data/remote_results>) | 存在 |
| P82A 正式轨迹导入副本 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/shear_lift_audit/data/imported_P82A](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/shear_lift_audit/data/imported_P82A>) | 存在 |
| 本地和远程日志 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/shear_lift_audit/logs](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/shear_lift_audit/logs>) | 存在 |
| 审核图 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/shear_lift_audit/figures](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/shear_lift_audit/figures>) | 存在 |
| 已有轨迹的展示动画 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/shear_lift_audit/animations](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/shear_lift_audit/animations>) | 存在 |
| 中文审核报告 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/shear_lift_audit/SHEAR_LIFT_MAGNITUDE_REVIEW.md](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/shear_lift_audit/SHEAR_LIFT_MAGNITUDE_REVIEW.md>) | 存在 |

## 上游输入、来源与历史记录

| 用途 | 已核对路径 | 状态 |
| --- | --- | --- |
| 血管原始数据、上游几何与来源 | [/home/lzy/projects/ulm_3D_vascular](</home/lzy/projects/ulm_3D_vascular>) | 存在 |
| 血管来源迁移记录 | [/home/lzy/projects/vascular_migration_audit](</home/lzy/projects/vascular_migration_audit>) | 存在 |
| 冻结 SonoVue 分布 | [/home/lzy/projects/sonovue_size_distribution_v0](</home/lzy/projects/sonovue_size_distribution_v0>) | 存在 |
| 原始冻结粒径直方图 | [/home/lzy/projects/FROZEN_SONOVUE_HISTOGRAM.csv](</home/lzy/projects/FROZEN_SONOVUE_HISTOGRAM.csv>) | 存在 |
| 保留的二维 FEM 参考依赖（7 个文件） | [/home/lzy/projects/ulm_microbubble_traj_gen_2D](</home/lzy/projects/ulm_microbubble_traj_gen_2D>) | 存在 |
| 本地清理证据 | [/home/lzy/projects/CLEANUP_AUDIT_20260923_174727](</home/lzy/projects/CLEANUP_AUDIT_20260923_174727>) | 存在 |

## GitHub 同步工作树与证据

| 用途 | 已核对路径 | 状态 |
| --- | --- | --- |
| 完整流程的精简同步工作树 | [/home/lzy/projects/github_sync/flow_microbubble_rbc_20260923](</home/lzy/projects/github_sync/flow_microbubble_rbc_20260923>) | 存在 |
| 同步说明、逐文件选择清单和 SHA256 | [/home/lzy/projects/github_sync/flow_microbubble_rbc_20260923/sync_metadata/flow_mb_rbc_20260923](</home/lzy/projects/github_sync/flow_microbubble_rbc_20260923/sync_metadata/flow_mb_rbc_20260923>) | 存在 |
| 发布日志、远端核对及未同步文件清单原始记录 | [/home/lzy/projects/github_sync/flow_mb_rbc_sync_20260923_audit](</home/lzy/projects/github_sync/flow_mb_rbc_sync_20260923_audit>) | 存在 |
| 远端发布核对回执 | [/home/lzy/projects/github_sync/flow_mb_rbc_sync_20260923_audit/publication_receipt.json](</home/lzy/projects/github_sync/flow_mb_rbc_sync_20260923_audit/publication_receipt.json>) | 存在 |
| SonoVue 历史同步工作树 | [/home/lzy/projects/github_sync/ulm_sonovue_20260920_111258](</home/lzy/projects/github_sync/ulm_sonovue_20260920_111258>) | 存在 |

## 服务器

| 用途 | 已核对路径 | 状态 |
| --- | --- | --- |
| 早期 FEM 工程 | `/workspace/formal_3D_flow_solver_FEM_lzy` | 存在 |
| 早期 FEM 的 Conda 环境 | `/workspace/formal_3D_flow_solver_FEM_lzy/remote/.env` | 存在 |
| SimVascular 原生完整工程及各阶段 | `/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy` | 存在 |
| 最终 SV1.3Q 阶段 | `/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3q` | 存在 |
| SV1.3Q 原生求解结果 | `/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3q/outputs/REAL_VASCULAR_GPU_ILU_REUSE_WINNER` | 存在 |
| SV1.3Q 服务器日志 | `/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3q/logs` | 存在 |
| SV1.3Q 构建与验证记录 | `/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3q/reports` | 存在 |
| 实际求解器源码 | `/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3q/external/svMultiPhysics-reuse` | 存在 |
| 实际求解器二进制 | `/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3q/external/build/svmp_gpu_reuse/svMultiPhysics-build/bin/svmultiphysics` | 存在 |
| PETSc 安装目录 | `/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3n/external/petsc325/install_gpu13` | 存在 |
| MPI 安装目录 | `/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3l/external/gpu_mpi_fortran` | 存在 |
| 2.0 mm/s 流场工程 | `/workspace/flow_mean_2p0_mmps_20260922` | 存在 |
| 2.0 mm/s 完整算例 | `/workspace/flow_mean_2p0_mmps_20260922/flow_cases/mean-2p0-mmps` | 存在 |
| 2.0 mm/s 服务器求解日志 | `/workspace/flow_mean_2p0_mmps_20260922/flow_cases/mean-2p0-mmps/run/solver.log` | 存在 |
| Particle-8.2 批处理、分片、检查点和验证 | `/root/particle8_2_runs` | 存在 |
| Particle 服务器 Python 环境 | `/root/particle8_2_runs/env` | 存在 |
| Particle 服务器 Python 解释器 | `/root/particle8_2_runs/env/bin/python` | 存在 |
| Particle 原冻结数据只读副本 | `/root/particle8_2_runs/readonly_inputs` | 存在 |
| P81 轨迹只读副本 | `/root/particle8_2_runs/readonly_inputs/p81/particle_3d/outputs/particle8_1` | 存在 |
| Particle-8.2 最终执行代码 | `/root/particle8_2_runs/2ec91fe4c6d142400a352bb9d949800e0cb3b489/final_execution/repo` | 存在 |
| Particle-8.2 完整审核结果 | `/root/particle8_2_runs/cbb45d8222b3569bce4bb4e2034f1c444b18584e/full_audit` | 存在 |
| 10 万点示踪参考结果 | `/root/particle8_2_runs/a4820609c9ed71f09bcea968b5819701632a2a1e/diagnostic_pipeline/point_basin_100000` | 存在 |
| Particle-8.2A 各提交代码与执行目录 | `/workspace/particle8_2a_20260922` | 存在 |
| Particle-8.2A 结果总目录 | `/workspace/particle8_2a_20260922/results` | 存在 |
| Particle-8.2A 准入计算数据 | `/workspace/particle8_2a_20260922/results/admission` | 存在 |
| Particle-8.2A 正式轨迹 | `/workspace/particle8_2a_20260922/results/trajectories/formal` | 存在 |
| 500 条轨迹 PPT 结果 | `/workspace/particle8_2a_20260922/results/ppt` | 存在 |
| 剪切升力审核完整目录 | `/workspace/shear_lift_audit_20260923` | 存在 |
| 升力审核最终代码 | `/workspace/shear_lift_audit_20260923/final_source/code` | 存在 |
| 升力审核最终测试 | `/workspace/shear_lift_audit_20260923/final_source/tests` | 存在 |
| 升力审核输入 | `/workspace/shear_lift_audit_20260923/input` | 存在 |
| 升力审核全量结果 | `/workspace/shear_lift_audit_20260923/results` | 存在 |
| 升力审核逐状态数据 | `/workspace/shear_lift_audit_20260923/results/states` | 存在 |
| 服务器上的冻结 SonoVue | `/home/lzy/projects/sonovue_size_distribution_v0` | 存在 |
| Particle 部署冻结 SonoVue 副本 | `/root/particle8_2_runs/readonly_inputs/sonovue_size_distribution_v0` | 存在 |

## Taylor-Hood 已停止的历史尝试

| 用途 | 已核对路径 | 状态 |
| --- | --- | --- |
| mean-2p0-mmps-A-H0-pressure-v1-taylor-hood-p2p1 | [/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-H0-pressure-v1-taylor-hood-p2p1](</home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-H0-pressure-v1-taylor-hood-p2p1>) | 存在 |
| mean-2p0-mmps-A-H0-pressure-v1-taylor-hood-p2p1-20260924T195557Z | [/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-H0-pressure-v1-taylor-hood-p2p1-20260924T195557Z](</home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-H0-pressure-v1-taylor-hood-p2p1-20260924T195557Z>) | 存在 |
| mean-2p0-mmps-A-H0-pressure-v1-taylor-hood-p2p1-20260924T201307Z | [/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-H0-pressure-v1-taylor-hood-p2p1-20260924T201307Z](</home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-H0-pressure-v1-taylor-hood-p2p1-20260924T201307Z>) | 存在 |

## P9历史验证报告

| 用途 | 已核对路径 | 状态 |
| --- | --- | --- |
| particle9a_2mmps | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/particle9a_2mmps](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/particle9a_2mmps>) | 存在 |
| particle9a_2mmps_diagnosis | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/particle9a_2mmps_diagnosis](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/particle9a_2mmps_diagnosis>) | 存在 |
| particle9a1_2mmps | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/particle9a1_2mmps](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/particle9a1_2mmps>) | 存在 |
| particle9a1_routing_stationary_audit | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/particle9a1_routing_stationary_audit](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/particle9a1_routing_stationary_audit>) | 存在 |
| particle9a2_inlet_sampling | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/particle9a2_inlet_sampling](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/particle9a2_inlet_sampling>) | 存在 |
| particle9a3_interior_inlet | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/particle9a3_interior_inlet](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/particle9a3_interior_inlet>) | 存在 |
| particle9a3b_flowfield_conservation | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/particle9a3b_flowfield_conservation](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/particle9a3b_flowfield_conservation>) | 存在 |

## Particle 各阶段原始输出

| 用途 | 已核对路径 | 状态 |
| --- | --- | --- |
| particle8_1 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/outputs/particle8_1](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/outputs/particle8_1>) | 存在 |
| particle8_2a | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/outputs/particle8_2a](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/outputs/particle8_2a>) | 存在 |
| particle8_2a_ppt | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/outputs/particle8_2a_ppt](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/outputs/particle8_2a_ppt>) | 存在 |
| particle8_full3d | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/outputs/particle8_full3d](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/outputs/particle8_full3d>) | 存在 |
| particle9a1_2mmps | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/outputs/particle9a1_2mmps](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/outputs/particle9a1_2mmps>) | 存在 |
| particle9a2_2mmps | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/outputs/particle9a2_2mmps](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/outputs/particle9a2_2mmps>) | 存在 |
| particle9a3_inlet_audit | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/outputs/particle9a3_inlet_audit](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/outputs/particle9a3_inlet_audit>) | 存在 |
| particle9a_2mmps | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/outputs/particle9a_2mmps](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/outputs/particle9a_2mmps>) | 存在 |

## 清理、同步与部署归档

| 用途 | 已核对路径 | 状态 |
| --- | --- | --- |
| github_sync_network_h0_20260924T221559Z | [/home/lzy/archives/github_sync_network_h0_20260924T221559Z](</home/lzy/archives/github_sync_network_h0_20260924T221559Z>) | 存在 |
| github_sync_network_h0_current.json | [/home/lzy/archives/github_sync_network_h0_current.json](</home/lzy/archives/github_sync_network_h0_current.json>) | 存在 |
| p9a4_smoke_20260924T233708Z | [/home/lzy/archives/p9a4_smoke_20260924T233708Z](</home/lzy/archives/p9a4_smoke_20260924T233708Z>) | 存在 |
| vascular_workflow_cleanup_20260924T214356Z | [/home/lzy/archives/vascular_workflow_cleanup_20260924T214356Z](</home/lzy/archives/vascular_workflow_cleanup_20260924T214356Z>) | 存在 |

## Git 同步快照

| 用途 | 已核对路径 | 状态 |
| --- | --- | --- |
| network_h0_particle_20260925 | [/home/lzy/projects/github_sync/network_h0_particle_20260925](</home/lzy/projects/github_sync/network_h0_particle_20260925>) | 存在 |
| flow_microbubble_rbc_20260923 | [/home/lzy/projects/github_sync/flow_microbubble_rbc_20260923](</home/lzy/projects/github_sync/flow_microbubble_rbc_20260923>) | 存在 |
| flow_mb_rbc_sync_20260923_audit | [/home/lzy/projects/github_sync/flow_mb_rbc_sync_20260923_audit](</home/lzy/projects/github_sync/flow_mb_rbc_sync_20260923_audit>) | 存在 |
| ulm_sonovue_20260920_111258 | [/home/lzy/projects/github_sync/ulm_sonovue_20260920_111258](</home/lzy/projects/github_sync/ulm_sonovue_20260920_111258>) | 存在 |
| 共享Git对象库（各worktree依赖） | [/home/lzy/projects/ulm-3d-vascular-model-generation/.git](</home/lzy/projects/ulm-3d-vascular-model-generation/.git>) | 存在 |

## 最新已同步快照内容

| 用途 | 已核对路径 | 状态 |
| --- | --- | --- |
| vascular_network | [/home/lzy/projects/github_sync/network_h0_particle_20260925/vascular_network](</home/lzy/projects/github_sync/network_h0_particle_20260925/vascular_network>) | 存在 |
| particle_3d | [/home/lzy/projects/github_sync/network_h0_particle_20260925/particle_3d](</home/lzy/projects/github_sync/network_h0_particle_20260925/particle_3d>) | 存在 |
| formal_3D_flow_solver | [/home/lzy/projects/github_sync/network_h0_particle_20260925/formal_3D_flow_solver](</home/lzy/projects/github_sync/network_h0_particle_20260925/formal_3D_flow_solver>) | 存在 |
| sonovue_size_distribution_v0 | [/home/lzy/projects/github_sync/network_h0_particle_20260925/sonovue_size_distribution_v0](</home/lzy/projects/github_sync/network_h0_particle_20260925/sonovue_size_distribution_v0>) | 存在 |
| source_dependencies | [/home/lzy/projects/github_sync/network_h0_particle_20260925/source_dependencies](</home/lzy/projects/github_sync/network_h0_particle_20260925/source_dependencies>) | 存在 |
| server_evidence | [/home/lzy/projects/github_sync/network_h0_particle_20260925/server_evidence](</home/lzy/projects/github_sync/network_h0_particle_20260925/server_evidence>) | 存在 |
| sync_metadata/network_h0_particle_20260925 | [/home/lzy/projects/github_sync/network_h0_particle_20260925/sync_metadata/network_h0_particle_20260925](</home/lzy/projects/github_sync/network_h0_particle_20260925/sync_metadata/network_h0_particle_20260925>) | 存在 |
| THESIS_SUMMARY | [/home/lzy/projects/github_sync/network_h0_particle_20260925/THESIS_SUMMARY](</home/lzy/projects/github_sync/network_h0_particle_20260925/THESIS_SUMMARY>) | 存在 |

## 最新已同步逐文件清单

| 用途 | 已核对路径 | 状态 |
| --- | --- | --- |
| local_file_inventory.csv | [/home/lzy/projects/github_sync/network_h0_particle_20260925/sync_metadata/network_h0_particle_20260925/local_file_inventory.csv](</home/lzy/projects/github_sync/network_h0_particle_20260925/sync_metadata/network_h0_particle_20260925/local_file_inventory.csv>) | 存在 |
| server_file_inventory.json | [/home/lzy/projects/github_sync/network_h0_particle_20260925/sync_metadata/network_h0_particle_20260925/server_file_inventory.json](</home/lzy/projects/github_sync/network_h0_particle_20260925/sync_metadata/network_h0_particle_20260925/server_file_inventory.json>) | 存在 |
| original_path_map.json | [/home/lzy/projects/github_sync/network_h0_particle_20260925/sync_metadata/network_h0_particle_20260925/original_path_map.json](</home/lzy/projects/github_sync/network_h0_particle_20260925/sync_metadata/network_h0_particle_20260925/original_path_map.json>) | 存在 |
| SYNC_REPORT_ZH.md | [/home/lzy/projects/github_sync/network_h0_particle_20260925/sync_metadata/network_h0_particle_20260925/SYNC_REPORT_ZH.md](</home/lzy/projects/github_sync/network_h0_particle_20260925/sync_metadata/network_h0_particle_20260925/SYNC_REPORT_ZH.md>) | 存在 |

## Particle 完整阶段目录

| 用途 | 已核对路径 | 状态 |
| --- | --- | --- |
| reports/particle0 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/particle0](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/particle0>) | 存在 |
| reports/particle1 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/particle1](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/particle1>) | 存在 |
| reports/particle4 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/particle4](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/particle4>) | 存在 |
| reports/particle5 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/particle5](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/particle5>) | 存在 |
| reports/particle6 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/particle6](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/particle6>) | 存在 |
| reports/particle6_5 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/particle6_5](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/particle6_5>) | 存在 |
| reports/particle7 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/particle7](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/particle7>) | 存在 |
| reports/particle8 | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/particle8](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/particle8>) | 存在 |
| reports/particle8_2a_ppt | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/particle8_2a_ppt](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/particle8_2a_ppt>) | 存在 |
| reports/rbc_mb_hydrodynamic | [/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/rbc_mb_hydrodynamic](</home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/rbc_mb_hydrodynamic>) | 存在 |

## 已知旧索引差异

旧索引中的 `/home/lzy/projects/ulm_particle_3d_particle0/formal_3D_flow_solver/FEM_SimVascular/frozen_reference/flow/steady_flow_stage_sv1_3q.vtu` 当前不存在；该默认冻结目录的现存文件为 `steady_flow_mean_2p0_mmps.vtu`，对应OLD equal-pressure场。原Stage Q文件副本与历史结果应按各自manifest查找，不能通过改文件名当作NEW场。
旧WORKFLOW_INDEX_ZH.md / WORKFLOW_PATHS_ZH.md保留原文，其中“最新”只代表其写作日期。旧ACTIVE_VASCULAR_WORKFLOW.md尚未涵盖刚完成的P9-A.4。本次新索引优先说明版本映射，未覆盖旧索引或受保护报告。
退役代码原路径与归档位置见 `/home/lzy/archives/vascular_workflow_cleanup_20260924T214356Z/audit/local_manifest.json` 和 `server_manifest.json`；代码已退役不表示其历史报告/轨迹已删除。

## Git工作树记录

```text
worktree /home/lzy/projects/ulm-3d-vascular-model-generation
HEAD c83dccea0f1fe5cd9fd3f5939e2eba883aeba9c2
branch refs/heads/sync/fem-simvascular-stage-q-particle-handoff-20260920

worktree /home/lzy/projects/github_sync/flow_microbubble_rbc_20260923
HEAD fab4cf0eb3c8f8ad1ce40dfe995181c147cb9a9d
branch refs/heads/dev/particle9a-2mmps-sync-snapshot

worktree /home/lzy/projects/github_sync/network_h0_particle_20260925
HEAD abb3ab05fb8dcddcf120765702b682f23a63a71d
branch refs/heads/sync/network-h0-particle-cleanup-20260925

worktree /home/lzy/projects/ulm_flow_mean_2p0_mmps
HEAD c83dccea0f1fe5cd9fd3f5939e2eba883aeba9c2
branch refs/heads/dev/flow-mean-2p0-mmps-20260922

worktree /home/lzy/projects/ulm_particle_3d_particle0
HEAD 775adc019536585a4ee2f4dbf2c083c959e0a339
branch refs/heads/dev/particle9a3b-flowfield-conservation-diagnosis

worktree /home/lzy/projects/ulm_particle_population_inlet_p9a4
HEAD abb3ab05fb8dcddcf120765702b682f23a63a71d
branch refs/heads/dev/p9a4-poisson-finite-size-flux-inlet-20260925
```

同一Git对象库被多个worktree引用；这些路径是不同版本/开发现场，不是可以直接互换的重复目录。
