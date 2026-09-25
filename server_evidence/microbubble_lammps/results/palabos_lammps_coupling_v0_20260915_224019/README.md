本目录是 PALABOS_LAMMPS_COUPLING_V0 的独立工程验证归档。

阅读 PALABOS_LAMMPS_COUPLING_V0_REPORT.md 和 FINAL_TERMINAL_SUMMARY.txt。
冻结物理与阈值：contracts/PALABOS_LAMMPS_COUPLING_V0_CONTRACT.json。
真实场只有一份：fields/FROZEN_FLOW_FIELD_V0.h5；所有算例相对路径引用同一文件。

不运行任何求解器的独立复核：

```bash
sha256sum -c SHA256SUMS
/usr/bin/python3 -B src/finalize_palabos_lammps_coupling.py --root "$PWD" --output /tmp/coupling_reaudit.json
/usr/bin/python3 -B scripts/audit_geometry_local.py "$PWD"
```

Python reference 使用 NumPy、h5py、SciPy；精确 STL 后处理另需 VTK。
源快照身份复核需要原始 WSL PBS/BSA smoke 或其原 Vast 路径；归档没有重复复制完整原快照。
几何后处理命令会重新写本目录的两个几何审计文件；需要保持整包字节不变时请在独立副本运行。

重新编译只读 LAMMPS baseline 的外部适配器：

```bash
cmake -S . -B /path/to/new/build -G Ninja -DENGINE_ROOT=/path/to/validated/lammps/work -DCMAKE_BUILD_TYPE=Release
cmake --build /path/to/new/build -j 4
```

CMake 使用当前已验证平台的 CUDA 13.2 与 Kokkos 静态库布局；迁往其他环境需配置依赖路径。
不修改 LAMMPS core。CPU/GPU embedding 在启动时向现有 public fix_map 注册 frozen/flow/drag。
`binaries/` 是实际执行二进制的证据，不能假定跨机器直接可用。
运行使用 `HWLOC_COMPONENTS=-gl OMP_NUM_THREADS=1 CUDA_VISIBLE_DEVICES=0`；前者只排除 GL 显示拓扑探测。
`scripts/run_validation.py` 是本轮固定远端一次性执行记录：检查启动标记，拒绝覆盖或二次求解，不是可直接重复提交的模板。
任何新科研 case 应另建目录并明确授权；本阶段结束后停止。
