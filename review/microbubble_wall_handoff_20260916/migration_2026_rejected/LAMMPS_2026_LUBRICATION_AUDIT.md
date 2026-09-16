# LAMMPS 2026 不等半径修正审计：FAIL

本报告是 Phase A 的上游源码身份/数学审计，不是 Phase C 的 production interaction 实现。目标始终冻结为 `patch_2Sep2026` / `d71abe6102c44577442ba7f03b7378a83166b9fd`。

## 文档与实际源码不一致

目标随附的 `doc/src/pair_lubricate.rst`（117 行起）及 `pair_lubricateU.rst`（126 行起）将不等半径的对称化修正标为 4Jul2026。实际下载源码经过与同一完整 commit 的官方 raw 文件逐字节比较。

然而 `pair_lubricate_poly.cpp` 使用 `h_sep = gap/r_i` 进入对数项，squeeze 对数项也缺少旧版实际源码中的 beta 前因子；`pair_lubricateU_poly.cpp` 的对应三个系数路径具有同样差异。旧 stable_22Jul2025_update6 源码则具有 symmetric gap `xi = 2*gap/(r_i+r_j)` 与相应 beta 因子。原始 diff 保存在 `provenance/upstream_lubrication/OLD_TO_TARGET.diff`。这里不依据版本日期推断修正，直接审计文件与计算结果。

## 实测

`source_correction_audit/source_coefficient_probe.cpp` 将上游系数片段原样提取编译；提取文件 SHA、行号和字节内容已记录。10,000 组不等半径分别交换 i/j：

| 上游路径 | flaglog | squeeze 最大相对交换误差 | shear 最大相对交换误差 |
|---|---:|---:|---:|
| 目标 lubricate/poly | 0 | 8.9915e-16 | 0（关闭） |
| 目标 lubricateU/poly | 0 | 8.6834e-16 | 0（关闭） |
| 目标两条路径 | 1 | 0.7935591195 | 0.8558137020 |
| 旧版 official poly 片段 | 1 | 1.1729e-15 | 1.0536e-15 |

所有 10,000 组目标 flaglog=1 系数均超出 1e-12 对称性门槛。独立的真实 LAMMPS 检查使用未修改的新 CPU binary、`pair_style lubricate/poly`、run 0，只读实际 atom force 数组：半径 0.65/1.5 µm，gap 0.05 µm，轴向速度 ±0.001 m/s、mu=0.001 Pa·s、flagfld=0、flagHI=1、flagVF=0。

- flaglog=1：Fx1 = -1.6313089638965207e-10 N；Fx2 = +2.1054918721048569e-10 N；合力 +4.7418290820833615e-11 N；归一化合力缺陷 **0.22521241449120213**。
- flaglog=0：归一化合力缺陷 **2.0004983091242208e-15**。
- 交换粒子 ID 后结论相同。四次检查均正常退出、0 timestep，没有惯性/transport 运行。

## 判定与适用范围

纯 1/gap normal-only 分支的对称性通过。本结果**不等于** normal-only V0 方案不可行；也不能把这一较窄的通过结果替代用户明确要求的“4Jul2026 完整不等半径修正已存在”硬门。

`POLYDISPERSE_LUBRICATION_FIX_PRESENT = NO_AS_REQUIRED_BY_DOCUMENTED_SYMMETRY`。因此 `LAMMPS_2026_MIGRATION = FAIL`，不提升目标版，不删除旧版，不开始 Phase C。没有补丁修改上游数学，也没有换用其它版本绕过冻结目标。既有 pipeline 回归单独记录，不掩盖源码门失败。

## 参考和复现

- [冻结目标官方发布](https://github.com/lammps/lammps/releases/tag/patch_2Sep2026)
- [冻结 commit 的 lubricate/poly 源码](https://github.com/lammps/lammps/blob/d71abe6102c44577442ba7f03b7378a83166b9fd/src/COLLOID/pair_lubricate_poly.cpp)
- [冻结 commit 的 lubricateU/poly 源码](https://github.com/lammps/lammps/blob/d71abe6102c44577442ba7f03b7378a83166b9fd/src/COLLOID/pair_lubricateU_poly.cpp)
- [官方说明](https://docs.lammps.org/pair_lubricateU.html)
- 复算：`python scripts/finalize_source_correction_audit.py --root .`（NumPy；不调用 solver）。
- `source_correction_audit/actual_lammps_flaglog*_swap*/` 包含完整输入、原始 forces.dump 和日志。
- 源码片段保留上游 GPL 身份；许可证随 provenance 归档。
