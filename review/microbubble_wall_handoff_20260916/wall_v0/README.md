本目录保存 MICROBUBBLE_WALL_HYDRODYNAMICS_V0 的代码、冻结输入、实际结果和独立核查。

总状态为 BLOCKED_LOCAL_PLANE_VALIDITY；技术实现通过，Case I/J 未获得完整真实时域结果。先读 [中文报告](MICROBUBBLE_WALL_HYDRODYNAMICS_V0_REPORT.md)。

复核只读取结果，不执行求解器：

```bash
cd /home/lzy/projects/compre_output/wall_hydrodynamics_v0/20260916_130547
sha256sum -c SHA256SUMS
OPENBLAS_NUM_THREADS=1 /home/lzy/projects/microbubble_wall_hydrodynamics_v0_20260916_130547/env/bin/python -B src/finalize_wall_hydrodynamics_v0.py --phase offline
OPENBLAS_NUM_THREADS=1 /home/lzy/projects/microbubble_wall_hydrodynamics_v0_20260916_130547/env/bin/python -B src/finalize_wall_hydrodynamics_v0.py --phase runtime
OPENBLAS_NUM_THREADS=1 /home/lzy/projects/microbubble_wall_hydrodynamics_v0_20260916_130547/env/bin/python -B src/finalize_wall_hydrodynamics_v0.py --phase supplemental
```

独立 finalizer 会重写 validation 结果文件；若要保留原始归档，请先复制目录再复核。环境版本见 provenance/LOCAL_VALIDATION_REQUIREMENTS.txt。离线 RMBW 生成使用之前冻结的 Python3.8 参考环境；不是这个新 Python3.12 核查环境。生成脚本有 FIRST 输出保护，不覆盖原始 RMBW 输出。

生产驱动位于 src/rigid_lammps_main.cpp；仅链接既有 LAMMPS 静态库，不重新编译其核心。远端 cmake 配置/编译命令见 provenance/control_scripts/build_remote.py，源与二进制 SHA 见运行 receipt 和 provenance/RUNTIME_SOURCE_SHA256.json。全部现有 case 已到终态；run_wall_phase.py 拒绝隐式重跑。不要直接重启现有 case 来冒充新证据。

静态表与 contract 在 tables/、contracts/；几何/PCA10k×2 在 raw/、validation/；全部实际 timestep 和两阶段记录在 cases/；12 张图和 ParaView 数据在 visualization/。Case I VTP 为空；Case J VTP 只有初始点，NaN 表示未求解量。

远端原始数值结果：/workspace/microbubble_lammps/results/wall_hydrodynamics_v0_20260916_130547
远端独立构建：/workspace/microbubble_lammps/work/wall_hydrodynamics_v0_20260916_130547
本地独立项目：/home/lzy/projects/microbubble_wall_hydrodynamics_v0_20260916_130547

RMBW 源码/GPL 来源和 Eigen 头文件许可分别在 provenance/reference_audit/、provenance/EIGEN_COPYRIGHT.txt 中。生产未复制 RMBW 实现源码。
