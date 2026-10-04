# PFPD_AT1_PETSC_VI_SOLVER_V1

## 1. 为什么换 solver 而不换 physics

上一轮已证明 C5 正相场 stationary roots 存在。本轮唯一问题是成熟 VI 求解器能否从合法历史态自动识别约束并找回同一分支。SVK spectral split、AT1、几何、horizon、Yc、Gc、sc=0.95 和全部科学阈值未改变。结论仅针对 stationary solver candidate，不是稳定性、裂纹或真实血栓预测。

## 2. 什么是 variational inequality

相场下界为上一接受载荷的 s，上界为 sc。下界要求 q/D≥0，内部要求 q/D=0，上界要求 q/D≤0，其中 D=max(Ψ+,Yc)>0；力自由度无界且要求归一化残差为零。项目仍独立检查 projected KKT≤1e−8、force≤1e−9。实际 PETSc 标量下界、上界、内部根测试和向量往返全部通过，见 mapping/SCALAR_AND_VECTOR.json。

## 3. RSLS 做了什么

[官方 RSLS](https://petsc.org/main/manualpages/SNES/SNESVINEWTONRSLS/) 在边界符号约束下建立自由子系统并执行 Newton 和线搜索。本轮使用 vinewtonrsls、默认 backtracking、PREONLY + MUMPS，无参数扫描；SNES atol/rtol=1e−12、stol=1e−14、max_it=30。三个 VI 监控及终止原因全部启用，原始文本/二进制保存在各 load 文件夹及运行日志。

初次接线的只读 Vec 访问错误和一次监控集合生命周期崩溃均已修复并保留开发日志。匹配版本 virs.c:509 在 SNESMonitor(:584) 前释放 IS_inact，因此不在 monitor 中读取 getVIInactiveSet；改用公开 C 函数 SNESVIGetActiveSetIS 在当前 (z,F) 上生成集合。lambda 通过公开 SNESLineSearchGetLambda 读取。此改动不影响 Newton 算法。

## 4. 与上一轮 explicit active-set 为什么本质相似

两者都求解约束下的 stationarity，并对未固定的相场分量施加 q=0。差别是本轮分类与迭代由 PETSc 实现；旧诊断仅作为事后 oracle。正式求解代码没有 active bond IDs，也不从已知正相场根初始化。Jacobian gate C 使用正根仅为导数检查，和求解入口隔离。

## 5. PETSc/petsc4py 安装与版本

PETSc 3.25.3，petsc4py 3.25.3，MPI MPICH Version:      5.0.1。Python `/home/lzy/.cache/codex-dolfinx-env/bin/python`；scalar=float64、PetscInt=int32，MUMPS 可用。复用现有环境，未修改其包。缺失的项目导入依赖以 binary wheels 安装到本轮 external/python_deps；完整安装日志和固定实际版本清单见 provenance。初次 llvmlite 版本不匹配已依据 numba metadata 修正为 0.49.0。没有源码构建或 CUDA 安装。

## 6. residual / bounds 映射

z 前 225 项沿用旧 scaled free displacement 坐标，后 5392 项为全部 s，共 5617 项。几何固定自由度不入向量。机械 residual 为原始 positive force/force_floor；phase residual 直接为 q/D，没有叠加 clip。稀疏 AIJ 包含 Jxx、Jxs、Jsx 和局部对角 Jss；没有全尺寸 dense 方阵。

保存的结构 nnz=2,482,417，CSR 占用 28.43 MiB；PETSc getInfo 的 memory 字段此构建为 0（不是实际无内存），因此另报 CSR 字节和进程 peak RSS。

## 7. Jacobian 验证

采用 inspect/AST 复用上一轮 ResidualBlocks 的原始导数语句：仅把独立位移列循环交给 ProcessPoolExecutor，去掉 projected C 行置零，以及把 projected diagonal 换回原有 Qs/D。生成的可读构造函数已归档；单测检查只有这三处变化。没有第二套导数推导。

A 最后接受 s=0、B 首个 weak activation、C 已知正相场根，分别 10 个固定随机方向 × 3 档 ε=[1e−5,3e−6,1e−6]，共 90 次。PETSc AIJ Jv 对中心差分最大相对误差 4.186862e-06，通过上一轮同量级门槛 5e−5。完整方向误差见 jacobian/GATE.json。

在 s=0 做双侧差分需要负 s 的代数探针。仅诊断函数在原始退化函数的私有 AST 副本中去掉 s-domain guard，保留完全相同 (1−s)² 公式和 split。正式 callback 仍调用原函数、拒绝负 s；从未把差分探针当作物理解或 initial guess。单测验证延拓在合法域与原函数一致。

## 8. first activation result

ε=0.004498390736212，初值来自 C5_W_004 加边界增量，相场全零。11 次 Newton 后科学接受，max s=0.000585374862705，force=3.705e-15，KKT=3.708e-14。事后内部键 [4488, 5059]，precision=1.0、recall=1.0。

首个载荷先只调整机械位移；接近激活时 BT 缩短步长，随后 PETSc 释放两个相场约束并进入正相场根。这是局部、夹持关联的 activation，不是裂纹预测。

## 9. original failure load result

ε=0.004514117555921，从前一个 PETSc 接受态继续。force=2.435e-15，KKT=6.722e-14，max s=0.00970408183757；完整科学门槛通过。旧 projected solver 源码和原失败结果未修改。

## 10. 七个 oracle roots regression

| ε | force | KKT | max s | Δs∞ | scaled Δx₂ | relative Δenergy | relative Δreaction | accepted |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| 0.004498390736 | 3.71e-15 | 3.71e-14 | 0.000585374863 | 1.12e-13 | 4.87e-13 | 4.44e-16 | 2.84e-14 | True |
| 0.004499390736 | 4.05e-15 | 1.41e-13 | 0.00117007874 | 1.96e-13 | 4.52e-13 | 1.62e-14 | 2.00e-15 | True |
| 0.004502390736 | 3.03e-15 | 3.41e-13 | 0.00292017816 | 4.87e-13 | 4.52e-13 | 1.15e-14 | 4.44e-16 | True |
| 0.004507390736 | 1.47e-15 | 1.07e-14 | 0.00582371806 | 9.59e-14 | 2.51e-13 | 1.11e-15 | 1.55e-15 | True |
| 0.004514117556 | 2.43e-15 | 6.72e-14 | 0.00970408184 | 1.83e-13 | 2.96e-13 | 1.12e-14 | 1.11e-15 | True |
| 0.004517390736 | 4.73e-15 | 3.76e-13 | 0.0115815539 | 6.42e-13 | 5.57e-13 | 1.35e-14 | 2.28e-14 | True |
| 0.004547390736 | 2.92e-15 | 8.01e-14 | 0.0284727189 | 1.24e-13 | 3.71e-13 | 1.12e-14 | 1.34e-14 | True |

预声明 agreement：phase Linf≤1e−8；自由位移差除 length_scale 后的 L2≤1e−8；energy relative≤1e−9；reaction relative≤1e−8。能量与反力分母分别使用对应 oracle 的绝对基准值。七点均满足，当前无多根证据；不据此主张根唯一。

## 11. continuation 是否自动维护 active set

七点按真实保存的 ε 排序，将原失败载荷置于 offset_1e−05 与 offset_2e−05 之间。每一点保存 z_initial、lower、accepted state；测试逐元素验证初值等于上一 PETSc 解的位移预测及上一 s，未逐点喂 oracle。所有点 inactive inequalities、严格下界、s<sc、不可逆性、固定边界、finite positive J 和旧安全检查通过。

## 12. RSLS vs SSLS

[官方 SSLS](https://petsc.org/main/manualpages/SNES/SNESVINEWTONSSLS/) 为另一成熟 VI formulation。本轮仅两点，读取对应 RSLS 的 z_initial 与 lower；同一 residual、Jacobian、tol 和线性配置，无调参。

| case | method | force | KKT | SNES | linear | seconds | scientifically accepted |
|---|---|---:|---:|---:|---:|---:|---|
| first | RSLS | 3.705e-15 | 3.708e-14 | 11 | 11 | 19.769 | True |
| first | SSLS | 3.864e-05 | 0.000e+00 | 30 | 30 | 73.621 | False |
| original failure | RSLS | 2.435e-15 | 6.722e-14 | 3 | 3 | 5.020 | True |
| original failure | SSLS | 3.171e-15 | 7.973e-14 | 5 | 5 | 13.788 | True |

SSLS 状态：DIFFICULT；科学接受与 PETSc reason 独立。详细 reason/error/门槛在 ssls/*/SUMMARY.json。

## 13. pfm-cracks / 2015 PDAS 对照

固定 `tjhei/cracks`、`paper-2015`、`342678db32e1736b3f1eae601f55322907f584c9`，GPL-2.0-or-later。源码行号、不可逆相场方向转换、FE 专属部分和停止/集合变化逻辑见 [PFM_CRACKS_PDAS_NOTES.md](references/PFM_CRACKS_PDAS_NOTES.md)。旧源码要求残差阈值和集合稳定共同满足，并有最大迭代限制；没有把其 FE 框架移植进来。

[Kolditz–Mang–Wick 2023](https://arxiv.org/html/2306.03501v2) 仅作为效率参考：若未来出现反复集合切换或过高组合 Newton 成本，再单独评估其改进。本轮不实现其算法。

## 14. CPU / MPI / GPU 实际效率

硬件复核为 12 physical / 24 logical；BLAS/OMP/MKL 均单线程；5 个 spawn 进程并行位移导数列。按可用内存、父进程预算和至少 20% RAM 预留选择，不使用 24 MPI ranks。continuation 串行，两个 SSLS 也串行执行以避免各自列工作池叠加。

| MPI ranks | same-root seconds | accepted |
|---:|---:|---|
| 1 | 19.930 | True |
| 2 | 23.284 | True |
| 4 | 24.268 | True |

最快 1 rank，故正式 continuation 使用该模式。每种模式只测一次，主张限于本机本轮观测；rank-zero BAQP 并未变成分布式 physics。平均 SNES/linear 均 3.714；PREONLY 的 linear iterations 表示直接线性求解次数，不能解释为 Krylov 收敛难度。continuation 路径 wall time=44.802 s，各 solve 时间合计=44.639 s。

RTX 4060 存在，但此 PETSc 未配置 CUDA/Kokkos/HIP，VecCUDA/AIJCUSPARSE 不可用；GPU_BACKEND=UNAVAILABLE，没有 GPU 性能结论。C7 optional smoke 仅完成 46794 DOF 的 residual 映射与 AIJ 初始化：低载荷 ε=0.00025，附近预测载荷 ε≈0.00314907。后三点 cheap scan 使用低载荷松弛位移的线性外推，未机械平衡，也不构成新 activation bracket。该状态 max Ψ+/Yc≈1.0093。完整 Jacobian、CSR/AIJ 临时数组与列工作池的保守估算为 6.49 GB，超过保留 20% RAM 后的预算，故未装配 C7 Jacobian、未执行 Newton，不能称 C7 coupled solver 已验证。执行证据见 c7/SMOKE.json；未运行任何 C7 Hessian、critical cone、mode 或空间稳定性工作。

## 15. tests

初始一次 full suite：350 passed、3 known failed、0 skipped。结束一次 full suite：368 passed、3 failed、0 skipped。原失败 ID/message/trace 是否逐字不变：True。开发期间仅运行本轮 targeted tests。

历史保护：23650 files、417 symlinks；modified=[]，missing=[]，changed_links=[]。生产文件修改数=0。

## 16. 当前 scientific status

```
PETSC_RUNTIME: AVAILABLE
PETSC_RSLS_ADAPTER: VERIFIED
PETSC_RSLS_C5: VERIFIED
FIRST_ACTIVATION_ROOT: RECOVERED_AUTOMATICALLY
ORIGINAL_FAILURE_LOAD: RECOVERED
C5_ORACLE_ROOTS: 7 / 7 PASSED
ACTIVE_SET_DISCOVERY: AUTOMATIC
RSLS_VS_DIAGNOSTIC_ORACLE: AGREES
PETSC_SSLS: DIFFICULT
PFM_CRACKS_PDAS_REVIEW: COMPLETE
MPI_BENCHMARK: 1 rank fastest
GPU_BACKEND: UNAVAILABLE
FORMAL_C5_STATIONARY_SOLVER: VERIFIED_CANDIDATE
C5_STABILITY: NOT ASSESSED THIS ROUND
C7_STABILITY: NOT RUN
PRODUCTION: UNCHANGED
```

## 17. 唯一下一步

等待人工审核本轮 stationary solver candidate。只有审核接受 VERIFIED_CANDIDATE 后，下一轮才能另行启动 C5 stability；本轮在此停止。

图件沿用 dark background、Arial regular、English labels、无箭头和无双 y 轴。图件见 `results/pfpd_at1_petsc_vi_solver_v1/figures/`：6 组、7 张图，每张有 PNG/PDF。所有分支图仅使用本轮 accepted roots。
