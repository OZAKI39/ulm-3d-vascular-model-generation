# PFPD_AT1_C5_STABILITY_V1

## 1. 为什么现在终于可以研究 stability

PETSc RSLS 已在上一轮从零相场自动恢复并验证七个 C5 stationary roots。本轮直接读入这些根，不重跑七点回归；另读入紧邻 activation bracket 下侧的 intact 状态和 diagnostic weak point。九个已有状态的 checksum、force、KKT、bounds、irreversibility、energy 与 max phase 均通过复核。七个 PETSc 根的 reaction 也与原摘要一致；两个旧 intact 摘要未存 reaction，此项明确记录为无历史标量参考，未静默补造。

## 2. stationary 与 stable 的区别

stationary 表示自由力平衡并满足相场 KKT；stable 还要求允许的小扰动不能降低总能量。求解器收敛、相场 activation 和稳定性变化是三件不同的事。本轮结果：激活弱点和首个正相场根都稳定，稳定正相场分支存在；在规定窗口的全部审计状态中未找到失稳。

## 3. 什么是 Hessian / curvature

固定上下端位移后，考察原 BAQP elastic energy 加原 AT1 crack energy。Hessian 表示总能量在微小扰动方向上的弯曲程度。正曲率像谷底，负曲率意味着存在向下方向。固定边界不参加 variation，也不添加外力势。这里没有使用 PETSc residual Jacobian 代替物理 Hessian。

统一使用历史尺度：Ls=2.5e-06 m，Ss=1，Es=7.930274414062e-13 J。所有表中曲率无量纲。三档差分 [1e−5,3e−6,1e−6]；复用原 full_hessian 中的 physical-gradient 差分循环，仅限制变量集合并分派独立列。

Hxx 最大原始非对称误差 5.505e-09，cross reciprocity 最大误差 2.228e-08；均先记录，再取对称部分。Hss 用原 pair weights 的 2WΨ+/Es 对角，AT1 线性项二阶导数为零；analytic-vs-gradient FD 最大误差 9.394e-11。Hxx 三尺度符号和幅值门槛通过。

直接 scalar-energy 相减的最小幅度起初出现约 5%–20% 误差，未放宽 2% 阈值。保留 FLOAT64_PRECHECK 和 FLOAT64_MECHANICAL_ORACLE。新增诊断接口只重绑定原 BAQP、spectral split、State.energy、AuditProblem.parts 的同一 code objects，使用 long-double 浮点构造和求和；不复制能量公式。NumPy linalg 仍在 float64 做小矩阵运算，但在高精度 Green strain 形成后才转换。该接口从不用于 stationary solver。

扩展精度与正式总能量逐状态相对一致性门槛为 1e−12，全部通过；10 个确定随机机械方向 × 三幅度的物理能量曲率复核，最大相对误差 2.413e-04。因此改善来自减小浮点抵消误差，不是调物理、步长或验收阈值。

每个状态有 68 个近零 principal Green strain 的键。它们全部属于 intact、strong-lower、ds=0 的键；在 critical subspace 上 g=1，Ψ+ + Ψ− 恢复原始光滑 SVK 能量。两个可动相场键不在谱零面，因此限制后的能量具有本轮需要的二阶导数。仍逐探针记录所有 raw sign changes；任何发生切换的幅度不作普通 centered 验证，而使用可行的单侧物理能量变化。没有引入谱平滑。

## 4. 为什么只分析真正可动的 phase directions

沿用 s_tol=1e−8、q_tol=1e−7，strong bound 的 ds 固定为零，weak lower 只允许 ds≥0，weak upper 只允许 ds≤0，interior 不限符号。初始 intact 状态 m=0，weak point 为两个弱下界键，全部 13 个正相场状态均只有两个 interior 键；没有构造 5392×5392 dense phase Hessian。

## 5. 为什么 mechanical condensation 能省很多计算

先检查 225×225 Hxx 正定。给定相场方向后，位移可以选择最有利的松弛；用一次 Cholesky 分解和 cho_solve 计算这次松弛，再从固定机械状态的相场刚度中减掉松弛收益，得到只有 m×m 的 Schur 矩阵 S。S 正定意味着即使机械位移尽量帮助相场降能量，所有相场方向仍向上弯曲。

计算采用 [SciPy Cholesky 接口](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.cho_factor.html)，不显式求逆。正定块与 Schur 的关系可参见 [Schur complement 推导](https://www.cis.upenn.edu/~jean/schur-comp.pdf)。本轮各 m=2 的 unconstrained S 均严格正定，所以整个允许 cone 稳定，无需运行受约束优化来证明稳定。SLSQP 仅保留为遇到不定 S 时的负 witness 搜索工具；没找到负 witness 不能证明稳定。

## 6. activation point 是否稳定

S_WEAK，ε=0.004497390736212：Hxx 最小值 0.0338400871，S 最小值 0.00212248883。两个 weak-lower 键满足单侧约束；正定 S 证书和单侧物理能量验证通过，分类 STATIONARY_STABLE。activation bracket 本身沿用旧结果，未重算。

## 7. first positive root 是否稳定

S1，ε=0.004498390736212，max s=0.000585374862705，两个 interior 键。S 最小值 0.00212468092，分类 STATIONARY_STABLE。故本轮证据排除了“刚激活就已失稳”这一解释。

## 8. positive-phase branch 能稳定多久

先完成已有 9 个状态，确认七个正根均稳定，才按固定 Δε=2.5e−5 新增 E01–E06。每点只继承前一 PETSc 接受态，1 MPI rank，原 RSLS 和全部阈值不变。共 15 个状态均科学 stationary 且通过稳定性证书。

最后 ε=0.004697390736212，max s=0.105282593434。E05 的 max s=0.0932857291；下一个规定粗步 E06 首次跨过 0.10，随即停止，没有为了精确碰到阈值修改 load step，也没有继续到 0.006。结论针对这 15 个审计点与规定步长，不是整个连续载荷区间的解析证明。

完整状态表：

| state | strain | force | KKT | strong lower / weak lower / interior | λmin Hxx | m | λmin S | best feasible phase Rayleigh | classification |
|---|---:|---:|---:|---|---:|---:|---:|---:|---|
| S0_LAST | 0.004497343750 | 9.93e-16 | 0.00e+00 | 5392/0/0 | 0.03384008 | 0 | N/A | 0.03384008 | STATIONARY_STABLE |
| S_WEAK | 0.004497390736 | 1.22e-15 | 1.34e-13 | 5390/2/0 | 0.03384009 | 2 | 0.002122489 | 0.002204956 | STATIONARY_STABLE |
| S1 | 0.004498390736 | 3.71e-15 | 3.71e-14 | 5390/0/2 | 0.03384014 | 2 | 0.002124681 | 0.002124681 | STATIONARY_STABLE |
| S2 | 0.004499390736 | 4.05e-15 | 1.41e-13 | 5390/0/2 | 0.0338402 | 2 | 0.002126873 | 0.002126873 | STATIONARY_STABLE |
| S3 | 0.004502390736 | 3.03e-15 | 3.41e-13 | 5390/0/2 | 0.03384037 | 2 | 0.002133446 | 0.002133446 | STATIONARY_STABLE |
| S4 | 0.004507390736 | 1.47e-15 | 1.07e-14 | 5390/0/2 | 0.03384066 | 2 | 0.002144395 | 0.002144395 | STATIONARY_STABLE |
| S5 | 0.004514117556 | 2.43e-15 | 6.72e-14 | 5390/0/2 | 0.03384105 | 2 | 0.002159112 | 0.002159112 | STATIONARY_STABLE |
| S6 | 0.004517390736 | 4.73e-15 | 3.76e-13 | 5390/0/2 | 0.03384124 | 2 | 0.002166268 | 0.002166268 | STATIONARY_STABLE |
| S7 | 0.004547390736 | 2.92e-15 | 8.01e-14 | 5390/0/2 | 0.03384296 | 2 | 0.002231705 | 0.002231705 | STATIONARY_STABLE |
| E01 | 0.004572390736 | 2.85e-15 | 6.48e-14 | 5390/0/2 | 0.0338444 | 2 | 0.002286043 | 0.002286043 | STATIONARY_STABLE |
| E02 | 0.004597390736 | 2.08e-15 | 9.03e-14 | 5390/0/2 | 0.03384584 | 2 | 0.002340222 | 0.002340222 | STATIONARY_STABLE |
| E03 | 0.004622390736 | 4.89e-15 | 2.10e-14 | 5390/0/2 | 0.03384728 | 2 | 0.002394258 | 0.002394258 | STATIONARY_STABLE |
| E04 | 0.004647390736 | 5.06e-15 | 7.93e-14 | 5390/0/2 | 0.03384871 | 2 | 0.002448165 | 0.002448165 | STATIONARY_STABLE |
| E05 | 0.004672390736 | 5.17e-15 | 8.30e-14 | 5390/0/2 | 0.03385015 | 2 | 0.002501956 | 0.002501956 | STATIONARY_STABLE |
| E06 | 0.004697390736 | 3.00e-15 | 4.94e-14 | 5390/0/2 | 0.03385159 | 2 | 0.002555644 | 0.002555644 | STATIONARY_STABLE |

所有 weak/strong upper 计数均为零。S0 的 Rayleigh 为 mechanical mode；其余为 ||p||=1 的 pᵀSp。恢复的 full mode 按统一 [u,p] norm 重新归一化，真实 full-mode 曲率另存在 MODE_ENERGY.json，不能把 Schur eigenvalue 当作 full-mode eigenvalue。

## 9. 第一次 instability interval

FIRST_STABILITY_LOSS=NOT_FOUND_IN_WINDOW。没有 verified unstable upper endpoint，因此没有合法失稳 bracket、宽度、失稳时 max s 或 active count。没有运行 midpoint refinement，也没有用 brentq。STABILITY_LOSS_RELATIVE_TO_ACTIVATION=INCONCLUSIVE 仅表示首次失稳位置未知；activation point 和已审计正分支稳定性已经解决。

## 10. unstable mode 在哪里

没有 verified unstable mode，故不提供 critical mode PR、top10 fraction、模式中心或位置图。保存的正曲率验证方向不是失稳模式。

## 11. 是否靠近 clamp

历史首次 activation 的两个键位于夹持附近；该事实不能替代本轮尚不存在的 unstable-mode localization。没有把 activation 键或正曲率验证方向称为物理裂纹。

## 12. C3 vs C5 有什么本质变化

只引用历史 C3：约 1% strain 时 185/185 键同时开启，历史失稳区间 [0.0099921875,0.01] 接近 activation。C5 约 0.449739% 时仅 2/5392 键开启，随后存在稳定正相场状态。本轮没有 C7，也不构成空间收敛结论。

## 13. tests

初始 full suite：368 passed / 3 known failed / 0 skipped。结束 full suite：390 passed / 3 failed / 0 skipped。原 3 项失败的 IDs、messages、tracebacks 逐字一致：True。未改旧测试、未加 xfail/skip。开发中仅 targeted tests；一次在分支尚未写完时启动的 targeted run 因最终汇总文件尚不存在出现两项失败，完成后复跑已保留日志。

PRE/POST manifests 与 PRESERVATION_CHECK 保护 27961 个历史文件及 449 个链接；历史修改=[]，生产修改=[]，旧测试修改=[]。

## 14. CPU/GPU/计算时间

12 physical / 24 logical，9 个 ProcessPoolExecutor worker，各 BLAS/OMP/MKL=1，按初始可用内存扣除至少 20% 总 RAM 和父进程预算。Hessian 与 PETSc 使用同一列工作池，PETSc 始终单 rank，没有 ranks×workers 嵌套扩张。唯一 Hxx/cross 矩阵构造总 wall time=54.854 s；修复能量复核精度时全部复用 hash 匹配缓存，未重新构造这些列。

父进程 peak RSS=509.7 MiB；最大 worker peak RSS=345.9 MiB。各进程峰值保守相加=3.538 GiB，不冒充同步测得的系统峰值。GPU=NOT USED，没有安装 GPU stack。

## 15. 当前 scientific status

```
C5_STABILITY_AUDIT: COMPLETE
PHYSICAL_HESSIAN: VERIFIED
MECHANICAL_SUBSYSTEM: STABLE
ACTIVATION_POINT_STABILITY: STABLE
FIRST_POSITIVE_ROOT_STABILITY: STABLE
STABLE_POSITIVE_PHASE_BRANCH: EXISTS
FIRST_STABILITY_LOSS: NOT_FOUND_IN_WINDOW
STABILITY_LOSS_INTERVAL: NOT_AVAILABLE
STABILITY_LOSS_RELATIVE_TO_ACTIVATION: INCONCLUSIVE
CRITICAL_MODE: NOT_AVAILABLE
MODE_LOCALIZATION: INCONCLUSIVE
PETSC_RSLS: UNCHANGED_VERIFIED_CANDIDATE
C7: NOT RUN
SPATIAL_CONVERGENCE: NOT YET ASSESSED
PRODUCTION: UNCHANGED
```

## 16. 唯一下一步

等待人工审核本轮能量曲率、精度接口和审计窗口结论。本轮在 max(s) 首次达到 0.10 后停止；未启动 C7、后失稳分支或任何后续模型工作。

本轮只生成 4 组有数据支撑的图：stability vs strain、activation vs stability、stationarity quality、历史 C3 vs C5。由于没有 bracket 或 unstable mode，不生成 energy_along_mode 的 stable/unstable 对比或 critical_mode_map。
