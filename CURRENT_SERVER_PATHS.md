# 当前服务器路径

连接别名：`vast4090`（已有 SSH 配置）。本轮只清理 WSL，未修改服务器文件、SSH 配置或运行环境。以下保留的是当前计算所需路径和独立 RBC 入口。

| 用途 | 服务器路径 |
|---|---|
| 新边界条件 FEM 算例 | `/workspace/flow_mean_2p0_mmps_A_H0_20260924T181140Z/mean-2p0-mmps-A-H0-pressure-v1` |
| 当前微泡 dt=1 ms | `/workspace/particle9a5_formal_trajectories_20260925T080115Z_dt1ms` |
| 当前出生样本生成来源 | `/workspace/particle9a4_population_inlet_20260924T233708Z` |
| 粒子 Python 环境 | `/root/particle8_2_runs/env/bin/python` |
| FEM 运行库、MPI/PETSc 构建信息 | `/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3q/reports/svmp_reuse_build.json` |
| 当前 FEM 二进制 | `/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3q/external/build/svmp_gpu_reuse/svMultiPhysics-build/bin/svmultiphysics` |
| FEM 求解器源码 | `/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3q/external/svMultiPhysics-reuse` |
| FEM 环境脚本 | `/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3n/use_gpu13_env.sh` |
| HemoCell RBC 工程 | `/workspace/hemocell_restore/work/rbc_stage1_20260915_165632` |
| HemoCell RBC 结果 | `/workspace/hemocell_restore/results/rbc_stage1_20260915_165632` |

微泡目录下的 `particle_3d/reports/particle9a5_formal_trajectories/outputs/formal/` 是当前正式轨迹，`outputs/points/` 是配对点示踪，`logs/` 是本次运行日志。具体部署上下文保留在本地 P9A5 的 `data/run_context.json`。

RBC 路径表示当前保留的独立 RBC 模型/结果，不能据此认定已经完成 Network-H0 新场中的 RBC 耦合生产计算。
