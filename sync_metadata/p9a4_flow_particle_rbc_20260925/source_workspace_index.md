# 当前血管流场与微泡代码入口

更新：2026-09-25（柏林时间）。当前主线为 vascular A → ROI → Network H0 边界 → 3D FEM → 配对微泡轨迹验证。
先按下表读取对应入口，再沿实际 import 跟进。完整历史报告和原始数据按需读取。

| 功能 | 当前入口 |
| --- | --- |
| 小鼠血管 ROI、表面与端口预处理 | [ulm_3D_vascular/README.md](/home/lzy/projects/ulm_3D_vascular/README.md) |
| Network A H0 | [run_a_network_h0_v2.py](/home/lzy/projects/ulm_3D_vascular/scripts/run_a_network_h0_v2.py)，[idealized_h0.py](/home/lzy/projects/ulm_3D_vascular/network_1d0d/idealized_h0.py) |
| 新 H0 FEM 算例生成 | [fem_h0_case.py](/home/lzy/projects/ulm_3D_vascular/network_1d0d/fem_h0_case.py) |
| 3D 求解与验收 | [solve_a_h0_fem_remote.py](/home/lzy/projects/ulm_3D_vascular/scripts/solve_a_h0_fem_remote.py)，[validate_a_h0_fem.py](/home/lzy/projects/ulm_3D_vascular/scripts/validate_a_h0_fem.py) |
| FEM 共享测量、收敛与解析 | [src/sv_validation](/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/src/sv_validation)，[Stage Q flow_parser.py](/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/scripts/sv13q/flow_parser.py) |
| 当前流场、派生量与历史 2 mm/s 复现 | [FEM README](/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/README.md)，[flow_2mmps](/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/scripts/flow_2mmps) |
| 当前 OLD/NEW 配对轨迹入口 | [runner.py](/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/network_derived_flow_mb_validation_v1/scripts/runner.py) |
| 微泡积分与阻力模型 | [Particle README](/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/README.md) |
| 旋转流场图与动画 | [OPEN_RESULTS.html](/home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/rotate_visualization/OPEN_RESULTS.html)；原有两份 render 脚本保留 |

最新轨迹状态是 `NETWORK_FLOW_MB_VALIDATION_PASS_WITH_STATIONARY`：相同 30 个初始状态，OLD/NEW 各 29 条完成、1 条 stationary、0 条数值失败。详见[正式报告](/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/reports/network_derived_flow_mb_validation_v1/NETWORK_DERIVED_FLOW_MB_VALIDATION_ZH.md)。500 条新流场生产计算尚未启动。

NEW 冻结流场：`/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-H0-pressure-v1/frozen_flow/steady_flow_mean_2p0_mmps_A_H0.vtu`，SHA-256 `064cbd28f3efa72f426fc946b2f29da21f056c596609095e7283d39070aa55f4`。
Particle 默认 `frozen_reference` 仍绑定历史 2 mm/s 场；配对入口通过 `--label NEW` 显式加载 H0 场。

## 快速回归

在各仓库根目录运行，Python 使用 `/home/lzy/projects/ulm_particle_3d_particle0/.venv/bin/python`；设置 `PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VTK_SMP_MAX_THREADS=1`。

```bash
# ulm_3D_vascular
python -m pytest -q -p no:cacheprovider tests/network_h0 tests/network_1d0d
# ulm_particle_3d_particle0
python -m pytest -q -p no:cacheprovider particle_3d/tests/network_flow_mb_validation_v1 particle_3d/tests/particle9a1_2mmps
# ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular
PYTHONPATH=src python -m pytest -q -p no:cacheprovider tests/flow_2mmps
```

## 服务器与历史代码

SSH：`vast4090`。当前 Particle 目录为 `/workspace/particle_network_flow_mb_validation_v1_20260924T210047Z`；H0 求解目录为 `/workspace/flow_mean_2p0_mmps_A_H0_20260924T181140Z`。
当前 GPU 求解器位于 `/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3q/external/build/svmp_gpu_reuse/svMultiPhysics-build/bin/svmultiphysics`，仍依赖 `sv1_3n/use_gpu13_env.sh` 和 Stage L 的 MPI/Fortran 运行库。

`sv_validation/sv13*.py` 和 `particle_3d/src/particle_3d/particle*.py` 中的旧阶段名称仍有实际调用；核心科学代码及当前报告内的 source snapshot 保留。
Taylor–Hood 任务已因资源成本停止，代码与结果保留供追溯；旧 DOLFINx FEM、LBM 和其他血管预处理分支保留为独立能力，不能仅因不在本次主线上判定无效。

已退役脚本、原 README、删除理由、依赖审计及回归日志位于 [vascular_workflow_cleanup_20260924T214356Z](/home/lzy/archives/vascular_workflow_cleanup_20260924T214356Z)。历史报告中的旧脚本路径需要时可按清单恢复。
各根目录 `.ignore` 只缩小默认 `rg` 文件搜索范围；不会影响 Python、Git 或数据读取。历史查找可显式指定路径，或使用 `rg --no-ignore`。
