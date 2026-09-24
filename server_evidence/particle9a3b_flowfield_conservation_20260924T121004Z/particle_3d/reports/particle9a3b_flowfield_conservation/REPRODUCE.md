# 复现 P9-A.3B 诊断

这些脚本只做 inventory、静态速度场诊断和画图，不生成微泡。原场、policy、旧 P9-A.3 证据均只读。输出会写入当前 checkout 的新 `particle9a3b_flowfield_conservation` 报告目录；若保留本次交付，先复制到新的 checkout/新运行目录，不覆盖本报告。

本地仓库根目录：`/home/lzy/projects/ulm_particle_3d_particle0`。

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=particle_3d/src
P9A3B_CASE=/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps
P9A3B_SOURCE=/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/vendor/svMultiPhysics_stage_q/Code/Source/solver
.venv/bin/python particle_3d/scripts/inventory_particle9a3b.py --repo "$PWD" --case "$P9A3B_CASE" --source "$P9A3B_SOURCE"
.venv/bin/python particle_3d/scripts/run_particle9a3b_diagnosis.py --repo "$PWD" --case "$P9A3B_CASE" --p9a3 particle_3d/reports/particle9a3_interior_inlet --workers 6
.venv/bin/python -m pytest particle_3d/tests/particle9a3b_flowfield_conservation -q -p no:cacheprovider
.venv/bin/python particle_3d/scripts/render_particle9a3b_diagnosis.py
.venv/bin/python particle_3d/scripts/finalize_particle9a3b_diagnosis.py
```

本次实际计算在 `/workspace/particle9a3b_flowfield_conservation_20260924T121004Z`，Python 为 `/root/particle8_2_runs/env/bin/python`。服务器命令把 `--case` 换成 `/workspace/flow_mean_2p0_mmps_20260922/flow_cases/mean-2p0-mmps`，`--source` 换成 `/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3q/external/svMultiPhysics-reuse/Code/Source/solver`，`--p9a3 inputs/p9a3`。执行时两脚本位于运行根目录；交付后亦保留规范 `particle_3d/scripts` 副本。测试需 export 同一 `P9A3B_CASE`。

六个位置在计算前按旧 candidate CSV 选定：最近合法几何位置到 0.07/1/10/30/50 µm，再加最后一个合法位置。不是按新流量寻找有利截面。裁切只做这6个真实 slab，解析控制复用其中一个体积。全域散度逐单元解析求解，无 CFD 时间推进。

`inventory_particle9a3b.py` 读取原生单进程 step71 checkpoint 的源代码规定格式，若节点、DOF、step、文件长度不符则拒绝。`velocity_representation_inventory.json` 必须先存在才允许运行主诊断。库依赖沿用既有 NumPy/SciPy/PyVista/VTK/Matplotlib/Pytest，没有安装 GPU 依赖。

报告、容差 provenance 和交付保护清单依赖本次源审计与人工阅读，不由数值 runner 自动编造。`protected_before.json` 来自计算前；复跑不能重新生成一份“之前清单”来假装已保护本次旧状态。
