# 2026-09-29 新增 BraVa 服务器路径

| 内容 | 路径 |
|---|---|
| BraVa当前独立算例 | `vast4090:/workspace/brava_flow_roi_18mlmin_20260928/` |
| 正式压力出口CFD | 该目录下 `cases/balanced_pressure_final/` |
| candidate 0 / 15流场动画 | 该目录下 `visualization/candidate_0/`、`candidate_15/` |
| 已停止的微泡批次 | 该目录下 `microbubble/`；`data/USER_STOP.json`为取消标记 |
| 隔离的GPU同步修复求解器 | 该目录下 `gpu_solver_fix/svmultiphysics` |

微泡服务为STOPPED、autostart=false、autorestart=false；本次同步不恢复计算。下列小鼠/RBC路径为独立保留流程：

# 当前服务器路径

连接别名：`vast4090`（已有 SSH 配置）。2026-09-27 的清理仅针对 WSL；2026-09-28 新增下列独立流场及微泡结果目录。旧依赖和独立 RBC 入口保留。

| 用途 | 服务器路径 |
|---|---|
| 当前 FEM：ROI-only-balanced-pressure-v1 | `/workspace/flow_roi_only_balance_20260928/mean-2p0-mmps-A-ROI-only-balanced-pressure-v1` |
| 当前微泡：1500 条，dt=0.5 ms | `/workspace/microbubble_roi_only_dt0p5ms_n1500_20260928` |
| 上一版 FEM：best-feasible-balance-v1 | `/workspace/flow_best_feasible_balance_20260928/mean-2p0-mmps-A-best-feasible-balance-v1` |
| 上一批微泡：1500 条，dt=0.5 ms | `/workspace/microbubble_best_balance_dt0p5ms_n1500_20260928` |
| CUDA/Torch 环境 | `/venv/main/bin/python` |
| 保留 H0 FEM 基线 | `/workspace/flow_mean_2p0_mmps_A_H0_20260924T181140Z/mean-2p0-mmps-A-H0-pressure-v1` |
| H0/CORE500 基线及本批只读源码依赖 | `/workspace/particle9a5_formal_trajectories_20260925T080115Z_dt1ms` |
| H0 入口样本来源 | `/workspace/particle9a4_population_inlet_20260924T233708Z` |
| 粒子 Python 环境 | `/root/particle8_2_runs/env/bin/python` |
| FEM 运行库、MPI/PETSc 构建信息 | `/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3q/reports/svmp_reuse_build.json` |
| 当前 FEM 二进制 | `/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3q/external/build/svmp_gpu_reuse/svMultiPhysics-build/bin/svmultiphysics` |
| FEM 求解器源码 | `/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3q/external/svMultiPhysics-reuse` |
| FEM 环境脚本 | `/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3n/use_gpu13_env.sh` |
| HemoCell RBC 工程 | `/workspace/hemocell_restore/work/rbc_stage1_20260915_165632` |
| HemoCell RBC 结果 | `/workspace/hemocell_restore/results/rbc_stage1_20260915_165632` |

最新批次的 `tracks/` 保存 1500 条原始轨迹和审计，`data/` 保存 0.5 ms 位置导出及统计，`OPEN_RESULTS.html` 为动画入口。运动积分为 8 CPU 进程，RTX 4090 用于 FP64 数据复核及 NVIDIA EGL 渲染。

H0 基线目录下的 `particle_3d/reports/particle9a5_formal_trajectories/outputs/formal/` 是该基线的正式轨迹，`outputs/points/` 是配对点示踪，`logs/` 是本次运行日志。具体部署上下文保留在本地 P9A5 的 `data/run_context.json`。

RBC 路径表示当前保留的独立 RBC 模型/结果，不能据此认定已经完成 Network-H0 新场中的 RBC 耦合生产计算。
