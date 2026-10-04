# PFPD_AT1_C7_EARLY_BRANCH_V1

## 1. 为什么现在做 C7

C5 已有可复核的 stationary roots 和稳定性证据，现在检查同一个物理立方体在 C7 采样下何时开启相场，以及求解器能否连续找到开启后的平衡态。C7 为 343 粒子、46059 键、735 个自由机械 DOF；读取历史 corrected geometry，未生成另一套几何。

L=0.25 mm，h=0.0416666666667 mm，δ/h=6.03。体积和上下表面权重分别通过 L³、L² 检查，几何 SHA256=25bc8e2580066b7baaac5292015e0ea3d8c659d90970be51ebee99e702fbbc05。

## 2. C5 已经告诉我们什么

历史 C5 约 0.449739% strain 时两个键首次开启；激活点和已有正相场状态稳定，已审计至 ε=0.004697390736、max s=0.105282593434，窗口内未找到失稳。本轮只读取这些结果。

## 3. 为什么本轮不直接做 C7 stability

相场下界 KKT 失效表示 s=0 不再满足相场平衡，stationary 表示力与相场方程同时满足；二者都不能代替局部变分稳定性。本轮只构造求解 raw residual 所需的 reduced Jacobian，未构造稳定性 Hessian、Schur 补、critical cone 或负曲率模式。

## 4. C7 activation bracket

真正机械平衡的 s=0 状态给出 [0.003137573242188, 0.003137634277344]，宽度 6.103515625e-08。宽扫在首次 min q 变号后停止，再按 midpoint 二分；未使用 brentq。随后仅为获取 weak set 做诊断二分，ε=0.003137607455254，原激活区间保留。

所有 26 个 intact 求解点均检查自由力、固定边界、positive J、安全范围，并对完整 46059 根键计算 q。使用原 SVK spectral 物理力和原 PETSc RSLS 设置，phase pool 为空使 s 精确固定为零；负主 Green 应变只作诊断。扫描从历史已接受 ε=0.00025 状态接续；宽扫按递增加载，二分从已接受 lower 状态接续。旧 cheap predictor ε≈0.00314907 未被当作激活结果。

## 5. 第一次多少 bonds 激活

首个诊断 weak set 为 2 / 46059，比例 4.34225667079e-05；IDs=[40174, 44110]。最终 bracket 上端的 first violated set IDs=[40174, 44110]。所有 q 分位数、分类 IDs 与完整数组随状态保存。

## 6. activation 在哪里

按预声明位置规则分类为 BOUNDARY_ASSOCIATED。质心 [0.0, 0.0, -4.811544270433496e-20] m；RMS 半径 104.166667 μm；键中点包围盒尺寸 [  0.           0.         208.33333333] μm；最近 clamp 距离的 min/Q1/median/Q3/max 为 [20.83333333 20.83333333 20.83333333 20.83333333 20.83333333] μm。

方向 |nx|,|ny|,|nz| 加权均值 [0.0, 0.0, 1.0]；0–90° 对加载轴夹角的 15° 分箱计数 [2, 0, 0, 0, 0, 0]。逐键端点、中点、长度、q、Ψ+ 与全部边界距离见 activation/SPATIAL.json。等效粒子参与体积 1.832619811e-13 m³，只是离散空间比较指标，不能称为裂纹宽度。对称的上下夹持区会使全局 RMS 半径较大；它不等同于单个局部区域的宽度。

## 7. 与 C5 是否相似

C5 激活 ε=0.004497390736212，C7 ε=0.003137607455254，差值 -1.359783281e-03，相对变化 -30.235%。激活载荷标签为 SHIFTING，不能由两级数据判断收敛。C5/C7 物理区域标签 CONSISTENT；C5 clamp 距离中位数 31.25 μm，C7 为 20.8333333 μm。

## 8. reduced phase pool 是什么

只将 735 个自由位移和少量候选相场变量送入 PETSc。初始池包含正相场历史、前一步 weak 集合、predictor 上 q_scaled≤1e−3 的键及共享任一端点的一环邻键。1e−3 直接沿用历史 near-active 阈值，在计算前声明；一环不作递归扩张。

## 9. 为什么这不改变 physics

池外相场精确固定在前一个接受态的下界。每次残差仍从完整物理状态计算；池内相场残差是未经投影截断的 q/max(Ψ+,Yc)，PETSc 自己判断其下界、内部或上界状态。原 BAQP 权重、谱分裂、相场律、几何、物理尺度、阈值和 PETSc 设置均不变。

导数通过执行原已验证 ResidualBlocks 的相同源代码构造，仅将 W0、相场权重和未知向量限制到当前 pool；Jxx、Jxs、Jsx、Jss 以 sparse AIJ 装配。没有 46059² 相场矩阵，也没有 46794² 总矩阵。参见 [PETSc RSLS 算法说明](https://petsc.org/release/manualpages/SNES/SNESVINEWTONRSLS/)。本机实际版本沿用 3.25.3。

## 10. full-KKT 如何保证被排除 bonds 真正可以固定

每次 solve 后重算全部键的 KKT；池外 q_scaled<−1e−8 的键及其一环邻键全部加入后重解，最多 8 次。池超过 25% 即停止。两个 C5 回归均从空池发现池外违例并自动扩至 426 键，随后恢复完整 RSLS 参考根；不是只在人工选择的活动键上检验。C7 初始池 1218，最大池 1218，实际扩池 0 次。

C5 两点回归均满足 phase Linf≤1e−8、scaled x L2≤1e−8、energy relative≤1e−9、reaction relative≤1e−8。C7 first weak state 的 10 个确定方向 × 三档步长 Jv 复核最大误差 1.936279511e-05，门槛 5e−5；使用可行的二阶单侧 raw residual 差分。

## 11. PETSc RSLS 是否跨过 C7 activation

FOUND。首个正根 ε=0.003137634277344，force=1.431281502e-15，全键 projected KKT=9.672792013e-14，max s=2.05673374e-05。从 bracket 下端真实 s=0 状态出发，不输入活动 IDs 或未来正根；初始下界全零。

## 12. early positive branch 是否存在

EXISTS，共 20 个顺序接受的正相场状态。后续每步只使用前一接受态及其 s 下界；从 Δε=1e−6 起步，连续五点满足预声明容易收敛条件后才升至 2e−6、5e−6。最终 ε=0.003197634277344，max s=0.0441886508。停止原因 20_ACCEPTED_POSITIVE_STATES。没有稳定性分类。

| index | strain | force | full KKT | max s | weighted mean s | interior | weak | pool | expansions | wall s | PETSc iterations |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 0.003137634277 | 1.431e-15 | 9.673e-14 | 2.0567337e-05 | 3.2930177e-09 | 2 | 0 | 1218 | 0 | 472.717 | 10 |
| 1 | 0.003138634277 | 3.002e-15 | 2.087e-14 | 0.00078634657 | 1.2590123e-07 | 2 | 0 | 1218 | 0 | 97.427 | 2 |
| 2 | 0.003139634277 | 3.092e-15 | 1.298e-13 | 0.0015510772 | 2.4834156e-07 | 2 | 0 | 1218 | 0 | 88.478 | 2 |
| 3 | 0.003140634277 | 2.665e-15 | 1.313e-13 | 0.0023147617 | 3.7061439e-07 | 2 | 0 | 1218 | 0 | 88.707 | 2 |
| 4 | 0.003141634277 | 2.742e-15 | 1.430e-13 | 0.0030774023 | 4.9272008e-07 | 2 | 0 | 1218 | 0 | 88.885 | 2 |
| 5 | 0.003142634277 | 4.565e-15 | 4.188e-14 | 0.0038390014 | 6.1465902e-07 | 2 | 0 | 1218 | 0 | 92.893 | 2 |
| 6 | 0.003144634277 | 3.328e-15 | 1.930e-13 | 0.0053590845 | 8.5803815e-07 | 2 | 0 | 1218 | 0 | 90.442 | 2 |
| 7 | 0.003146634277 | 2.271e-15 | 1.081e-13 | 0.0068750297 | 1.1007548e-06 | 2 | 0 | 1218 | 0 | 92.787 | 2 |
| 8 | 0.003148634277 | 2.988e-15 | 2.248e-13 | 0.0083868555 | 1.3428118e-06 | 2 | 0 | 1218 | 0 | 100.000 | 2 |
| 9 | 0.003150634277 | 2.600e-15 | 1.902e-13 | 0.0098945802 | 1.5842122e-06 | 2 | 0 | 1218 | 0 | 99.969 | 2 |
| 10 | 0.003152634277 | 4.082e-15 | 6.547e-14 | 0.011398222 | 1.824959e-06 | 2 | 0 | 1218 | 0 | 96.942 | 2 |
| 11 | 0.003157634277 | 1.071e-15 | 1.115e-13 | 0.015139585 | 2.4239851e-06 | 2 | 0 | 1218 | 0 | 141.227 | 3 |
| 12 | 0.003162634277 | 1.120e-15 | 6.956e-14 | 0.018855823 | 3.0189887e-06 | 2 | 0 | 1218 | 0 | 148.623 | 3 |
| 13 | 0.003167634277 | 1.204e-15 | 8.649e-14 | 0.022547215 | 3.6100141e-06 | 2 | 0 | 1218 | 0 | 139.186 | 3 |
| 14 | 0.003172634277 | 1.121e-15 | 5.646e-14 | 0.02621403 | 4.1971045e-06 | 2 | 0 | 1218 | 0 | 132.060 | 3 |
| 15 | 0.003177634277 | 1.377e-15 | 1.264e-13 | 0.029856536 | 4.7803028e-06 | 2 | 0 | 1218 | 0 | 134.047 | 3 |
| 16 | 0.003182634277 | 9.968e-16 | 6.667e-14 | 0.033474995 | 5.3596509e-06 | 2 | 0 | 1218 | 0 | 134.288 | 3 |
| 17 | 0.003187634277 | 1.593e-15 | 1.257e-13 | 0.037069665 | 5.9351903e-06 | 2 | 0 | 1218 | 0 | 133.427 | 3 |
| 18 | 0.003192634277 | 1.214e-15 | 3.005e-14 | 0.0406408 | 6.5069615e-06 | 2 | 0 | 1218 | 0 | 135.282 | 3 |
| 19 | 0.003197634277 | 1.068e-15 | 9.721e-15 | 0.044188651 | 7.0750047e-06 | 2 | 0 | 1218 | 0 | 132.417 | 3 |

## 13. damage zone 是否增长

使用原 phase_detection_tol=1e−8，首点 active count=2，末点=2。逐点保存 phase 加权质心、RMS 半径、包围盒、夹持距离和等效参与体积；幅值增长与空间支持区域扩张分开判断。真实结果见 activation_region_growth 图及各状态 SPATIAL.json。

## 14. C5/C7 spatial comparison

激活比例 C5=0.000370919881，C7=4.34225667e-05；RMS 半径 C5=93.75 μm，C7=104.166667 μm。统一 physical coordinates 的图用于比较位置。另从历史 C5 正分支读取与 C7 相对激活载荷最接近的状态，仅保留 offset 差≤1e−6 的配对，共 8 对，详见 C5_C7_EARLY_SPATIAL_COMPARISON.json；未重解、未插值，也未比较 C7 稳定性。

激活诊断状态的反力名义应力：C5=5.63843821 Pa，C7=2.94749468 Pa。最初两个接受正根之间的 Δmax(s)/Δε 分别为 584.703876、765.779229；加权平均相场斜率分别为 0.881731359、0.122608216。这是有限步长割线，不是激活点精确导数。在 C7 相对激活载荷窗口内的已存正根上，两级 active count 均保持 2，激活跳变之后的活动键增长率均为 0 bonds/strain。逐级应变窗口和计算定义见 C5_C7_RESPONSE_COMPARISON.json。C5 反力仅在归档 S_WEAK 状态上用原物理力重新求值，未求解新的平衡态。

## 15. CPU/GPU 实际使用和 speedup

重新检测到 12 physical / 24 logical cores；根据可用 RAM 选择 6 个 ProcessPoolExecutor worker，BLAS/OMP/MKL 每进程单线程，PETSc 1 rank。科学求解与 Jacobian gate 的累计 wall time=4546.593 s；这是各实际阶段 wall time 相加，包含并行等待，不是各 CPU 核时间总和。

外部同步采样观察到进程树 RSS 峰值 2.873 GiB，最低 MemAvailable 2.307 GiB；任务进程 VmSwap 最大 0 bytes。采样从 C7 宽扫中段开始，故报告为观察峰值而非全程精确峰值；逐阶段预算及父/worker 峰值另存。

GPU_BACKEND=RUNTIME_UNAVAILABLE。RTX 4060 与 Numba CUDA 驱动访问可用，但缺少已验证的 GPU 数组/批量特征分解后端；PETSc 无 CUDA/Kokkos/HIP；隔离目录的 CuPy/CUDA binary wheel 下载在 6 次不完整下载后失败，未源码编译或修改 system Python。CPU 全 46059 键 constitutive+q batch 十次热运行中位数 62.853 ms。GPU correctness/speed gates 未执行，speedup 和峰值 VRAM 均为 N/A；所有科学结果为 CPU。二进制安装方式参考 [CuPy 官方说明](https://docs.cupy.dev/en/stable/install.html)。

## 16. tests

开始 full suite：390 passed / 3 known failed / 0 skipped。结束：412 passed / 3 failed / 0 skipped。原三项失败的 IDs、messages、tracebacks 完全一致：True。开发中仅 targeted tests；保留空池索引和元数据字段冲突的早期调试日志，修复只涉及新增代码，未调整物理或验收门槛。

历史保护检查：{'files': 29403, 'symlinks': 481, 'modified': [], 'missing': [], 'symlink_changes': [], 'production_modified': [], 'existing_tests_modified': [], 'passed': True}。原生产、旧测试和历史结果不改写。PRE/POST manifests、geometry/state/physics/previous-state hashes 与自动 resume 测试提供可追溯记录。

## 17. 当前 scientific status

```
C7_ACTIVATION_AUDIT: COMPLETE
C7_S0_ACTIVATION: BRACKETED
C7_ACTIVATION_INTERVAL: [0.0031375732421875004, 0.0031376342773437502]
C7_FIRST_ACTIVATION_SET: 2 / 46059
C7_ACTIVATION_PATTERN: BOUNDARY_ASSOCIATED
C5_TO_C7_ACTIVATION_STRAIN: SHIFTING
C5_TO_C7_ACTIVATION_REGION: CONSISTENT
REDUCED_PHASE_POOL: VERIFIED
C7_PETSC_RSLS: VERIFIED_FOR_EARLY_BRANCH
C7_FIRST_POSITIVE_ROOT: FOUND
C7_EARLY_POSITIVE_BRANCH: EXISTS
C7_MAX_PHASE_REACHED: 0.04418865075454984
C7_STABILITY: NOT_ASSESSED
SPATIAL_STABILITY_CONVERGENCE: NOT_YET_ESTABLISHED
GPU_BACKEND: RUNTIME_UNAVAILABLE
PRODUCTION: UNCHANGED
```

## 18. 唯一下一步

等待人工审核 C7 激活区间、空间结构、reduced pool 与全键 KKT 证据。本轮完成后停止，不自动运行 C7 stability、C9 或后续模型工作。
