# 最终有效执行范围与复核命令

本轮已完成，后续命令仅用于复核已有数据。用户明确取消细档真实血管CFD和出口压力敏感性，均为SKIPPED_BY_USER；没有新的CFD启动权限。本轮发生过的原始命令、失败及执行方式变化按时间完整保留于`evidence/COMMANDS_CHRONOLOGICAL_HISTORY.md`，其中旧队列/细档/压力计划已失效，不应重新执行。

## 已实际完成的数值工作

- 生产入口修复和H0逐值回归、仿射/旋转/Couette/法向及配置拒绝测试：`stage1/`。
- 三档圆管实际CFD及最细圆管dt/2：`stage2/pipe_nr4`、`pipe_nr8`、`pipe_nr16_cpu_mpi8`、`pipe_nr16_cpu_mpi8_halfdt`。
- 血管原H0身份/时间序列复核、原网格CPU8复算、原网格GPU1半步长控制、中档CPU8/LU独立CFD。中档4079.429秒、第51步正常结束，118次线性校正无失败。
- 已对有效例运行实际生产WSS入口、日志/区域/截面分析、几何质量、壁面模板及控制体积分检查。四份J1/J2实际场子集恢复差均0Pa。
- 原—中两档敏感性汇总及五图；没有将细档或压力取消状态变成PASS。

每个实际启动例的`reports/launch.json`包含原始完整MPI命令、二进制/PETSc/代码/配置哈希；`run/solver.xml`与`run/PETSC_OPTIONS.txt`是生效输入；`reports/execution.json`和`run/solver.log`记录实际求解，而非任务队列状态。失败、未完成性能试验不属于科学结果。

## 在当前项目复核已有结果

工作目录 `/home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular`，不修改正式冻结H0：

```bash
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B rotate_visualization/prepare_surface_data.py --case /home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-H0-pressure-v1 --output wss_validation_v2/stage1/recomputed
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B wss_validation_v2/scripts/analyze_vessel.py --case wss_validation_v2/stage3/vessel_medium
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B wss_validation_v2/scripts/control_volume_checks.py --case wss_validation_v2/stage3/vessel_medium
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B wss_validation_v2/scripts/collect_results.py
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B wss_validation_v2/scripts/collect_geometry.py
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python -B wss_validation_v2/scripts/summarize_two_mesh.py
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python -B wss_validation_v2/scripts/make_figures.py --part vessel
```

以上命令已实际执行过；重做会更新本次验证目录的诊断文件，不会启动CFD。圆管已有验证不需为交付重复运行。

## 精简包独立复核

在解压后的`wss_validation_v2/`，只需NumPy：

```bash
python -B scripts/reproduce_local_wall_evidence.py
```

核对生产核心及数据哈希，在实际J1/J2子集上重算壁面牵引向量/模长；不能将其当作闭合域CFD或物理精度证明。完整原始数据及服务器路径见`ARTIFACT_INDEX.csv`。

## 自动流程的最终状态及取消保护

`watch_vessel_results.py`只处理已授权的原网格控制和中档，已完成退出。`finish_numerical_sequence.py`只检查/汇总两档数据及绘图，已完成退出；其中不存在生成、注册、启动或等待压力CFD的代码。细档等待程序全部停止并禁用autostart；取消标记由求解/初始化/等待/注册入口识别。不要启动新的细档/压力/性能/时间步任务。

只读服务器核查：

```bash
ssh vast4090 'supervisorctl -c /workspace/wss_validation_v2_20260927T1230Z/supervisord.conf status'
```

中档`vessel_medium_LU`为EXITED，细档为STOPPED，未注册O2压力程序。专用服务器根`/workspace/wss_validation_v2_20260927T1230Z`，不会操作系统supervisor或其他项目。若以后另获授权需要CFD复现，应使用新case目录和相同显式输入/初值/求解器，不覆盖本轮失败或成功日志；本次不执行该后续工作。
