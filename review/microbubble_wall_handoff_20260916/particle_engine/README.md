# Minimal LAMMPS particle engine V0

本目录是一套独立的有限尺寸刚性球软件工程验证。直径来自冻结的 SonoVue 连续经验 CDF sampler；1000 和 50000 颗只用于 engine/GPU smoke，不代表浓度。

- 源码：官方 LAMMPS 22 Jul 2025 Update 6。两套 CMake build 共用同一份已核验源码；CPU GRANULAR，GPU GRANULAR + KOKKOS CUDA ADA89；LATBOLTZ 关闭。
- 单位：SI。转换器只执行 `diameter_um * 1e-6`，ID 从 sampler 的零基编号映射为 LAMMPS 的一基编号。密度 1000 kg/m³ 是技术测试值。
- `contracts/`、原始 population、data/input 及验证源码在第一次 `run 0` 前封存于 `FROZEN_INPUT_SHA256SUMS`。
- 执行顺序：直径 round-trip → A → B → C → D → E MPI1/MPI4 → F GPU → CPU/GPU 无接触比较。逐关独立 finalizer 检查；失败即停止依赖算例，无自动重试。
- F 仅一次执行，使用已有 Nsight Systems 记录 CUDA kernel；输出包括原始 `.nsys-rep` 和只读 SQLite 导出。带 profiler 的 wall time 只作该次运行记录，不作优化或硬件速度结论。
- 所有算例只输出 initial/final dump，浮点使用 17 位有效数字，坐标比较使用 unwrapped coordinates。
- `build_provenance/` 保存下载、构建、编译器、哈希、帮助输出和对应版本的官方文档摘录；源码包、build tree 和 binary 留在远端 work 目录，不纳入下载结果。

冻结 sampler 本体和 population 在原 WSL NumPy 1.26.4 环境中生成。远端转换/验证只读取这些 population，不用另一个 NumPy 版本重新抽样。`scripts/prepare_stage.py` 是该次生成过程的完整记录，不能在已封存目录上覆盖运行。

运行控制器：`python3 -B scripts/run_stage.py --root RESULT_DIR --work WORK_DIR`。控制器拒绝已有 stage/运行标记，不能当作自动续跑或重复测试入口。

独立复核（不启动 LAMMPS）：

```bash
sha256sum -c SHA256SUMS
python3 -B src/finalize_lammps_particle_engine.py --root . --output /tmp/lammps_independent_audit.json
```

数值门槛全部读取冻结合同。技术接触通过不等于 SonoVue 微泡相互作用被验证：bubble interaction、lubrication、transport、drag、wall、adhesion、Palabos–LAMMPS coupling 均待另立科学合同；RBC 关闭。
