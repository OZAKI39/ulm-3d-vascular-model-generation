# Stable rigid-sphere near-field V0

科研结论、近似范围和未完成的 twist 请先阅读 `STABLE_LAMMPS_RIGID_SPHERE_NEAR_FIELD_V0_REPORT.md` 与 `THEORY_DERIVATION.md`。

## 已冻结环境

唯一上游版本为 LAMMPS 22Jul2025 Update 6，commit `9c5ab448c78a14fd534619622162ba418d6a1fb1`。CPU/GPU 核心二进制身份见 `LAMMPS_ACTIVE_INSTALL_AUDIT.json`。新外部 embedding 的实际二进制 SHA256、动态链接和编译配置见 `provenance/PROJECT_BUILD_PROVENANCE.json`。构建依赖 C++17、MPI、HDF5、PNG/JPEG/Zlib；Kokkos 分支使用原冻结 GPU 库。

```
cmake -S RESULT_ROOT -B NEW_BUILD_DIR -G Ninja \
  -DENGINE_ROOT=/workspace/lammps_active -DCMAKE_BUILD_TYPE=Release
cmake --build NEW_BUILD_DIR -j4
```

`cases/*/case.cfg` 使用相对输入路径。每个案例的初始粒子、固定数值设置和驱动定义保存在 `CASE_CONTRACT.json`。如需将来复现实验，应先复制为新的结果目录，以保留本次原始证据；从该复制目录的单个 case 目录启动相应 embedding 即可：

```
HWLOC_COMPONENTS=-gl OMP_NUM_THREADS=1 \
mpirun --allow-run-as-root --bind-to none -np 1 \
  /ABSOLUTE/NEW_BUILD_DIR/rigid_lmp_cpu case.cfg
```

I_mpi4 对应 MPI4。Kokkos 兼容测试使用 `rigid_lmp_gpu case.cfg kokkos`。本轮实际执行命令、返回码和耗时在各 case 的 `RUN_RECEIPT.json` 中，不应由示例命令替代。

## 独立复核

`src/reference_rigid_sphere_lubrication.py` 不依赖生产 C++。NumPy/SciPy/h5py 的 Python 环境可以执行单案例或全局 finalizer：

```
OPENBLAS_NUM_THREADS=1 python -B scripts/validate_cases.py K_real_8
OPENBLAS_NUM_THREADS=1 python -B scripts/finalize_rigid_v0.py
```

全局 finalizer 的旧 passive 回归还读取报告所列本地/远端原始 baseline 路径。WSL 中安装 VTK 的环境使用 `--full-local`，重新检查真实 STL、内部点与连续线段。`scripts/visualize_rigid_v0.py` 生成 QA 图与精确回读的 ASCII VTP。

`provenance/control_scripts/` 是本轮操作归档，包括已经完成的候选退役；不属于新求解器启动流程。

## 核查文件

- `SHA256SUMS`：全量结果校验。
- `FINAL_TERMINAL_SUMMARY.txt`、`FINAL_SUMMARY.json`：最终状态。
- `LOCAL_INDEPENDENT_FINALIZER.json`：WSL 独立数值复核。
- `validation/`：100000 组系数、10000 组组合运动、矩阵、梯度、C_gap、MPI、排列、回归及真实重放证据。
- `PAIR_HYDRODYNAMIC_EVENTS.csv`：带 case 索引的全部球对事件；各 case 原始 CSV 保留。
- `visualization/CASE_K_TRAJECTORIES.vtp` 与 `FROZEN_LUMEN.vtp`：可一起在 ParaView 中打开。

Twist、壁面水动力与黏附均未完成；RBC、形变关闭。自动核验和生成图形没有代替人工审阅。
