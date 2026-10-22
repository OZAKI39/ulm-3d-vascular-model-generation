# Frozen Gmsh V2 GPU manufactured field 验证

2026-10-04。本轮结论：**REFERENCE_PIPE_FIELD_FAIL_SCIENTIFIC**。正式 GPU 执行路线已得到合格稳态，但冻结 particle-band 与 bulk 速度门槛未全部通过，按用户第 21 条停止。此失败发生在科学 field gate，运行没有超时、发散或回退到 CPU。

## 执行与科学身份

执行 authority 为 `cases/balanced_pressure_final/reports/execution.json`；直接使用其 solver binary、环境 wrapper、MPI launcher、PETSc options 与 RA adaptive preconditioner reuse。8 MPI ranks 共享 RTX 4090 的 GPU 0，OMP=1，历史配置 `use_gpu_aware_mpi=0` 保持。PETSc 原生 Mat 为 `mpiaijcusparse` / `seqaijcusparse`，Vec 为 `mpicuda`；进程实际加载 CUDA runtime、cuSPARSE 与指定 PETSc，原生日志观察到 RA reuse。ASM 子域 CPU LU 属于同一正式执行路线。

科学 authority 为原 `manufactured_pipe/R1_GMSH_V2/cfd`：mesh、caps、25/50/75 µm shells、精确入口速度、精确 Cauchy 出口 traction、wall no-slip、rho=1056 kg/m³、mu=0.00345312 Pa·s、Umean=0.1 m/s、CFD dt=0.01 s、nonlinear tolerance=1e-10 与 KSP rtol=1e-10 / atol=1e-24 全部保持。solver.xml 和 mesh/BC 输入逐字节一致。没有将 BraVa 历史校准边界用于直管。

旧 CPU step10 native checkpoint 本身完整有限、binary/mesh/physics/rank 相同，但缺少独立确认 rank-local node-order 对应的记录，因此按授权从 **GPU step0 fresh run** 开始。旧 CPU step10 没有用于科学认证或伪造 restart。

- mesh SHA256：`3a16731d70547ab4b657000c1c99d9af6aa4f0581f205cd045c37e410c1bbbd7`。
- solver SHA256：`43cd7d7b9ae9bde91275982225c6bf084859141cdaf2ec1030e1b83c3b94691d`。
- PETSc options SHA256：`00fc38aa688187f52f07cc866d597c4feedb34ddc11a1384ede9e64e8d8ab994`。
- 完整运行身份：[GPU_MANUFACTURED_RUN_IDENTITY.json](GPU_MANUFACTURED_RUN_IDENTITY.json)。

## 稳态与原生保存证据

连续五个合格保存区间如下；阈值 Eu≤1e-5、EQ≤1e-6，质量与入口流量误差≤1e-6，finite/no-slip 等原 gate 不变。

| 保存步区间 | Eu | EQ | 判定 |
|---|---:|---:|---|
| 20→30 | 9.946883163e-07 | 9.527429935e-07 | PASS |
| 30→40 | 3.070658744e-09 | 9.304132311e-10 | PASS |
| 40→50 | 1.698538689e-11 | 9.084229269e-13 | PASS |
| 50→60 | 1.058131794e-13 | 8.293935684e-16 | PASS |
| 60→70 | 7.368270390e-16 | 0.000000000e+00 | PASS |

第 70 步触发正常 STOP_SIM，原生求解器完成当前步骤后在 **第 71 步**正常退出并保存最终 VTU/native checkpoint；分析使用该最终保存态。最终质量不平衡 `2.369695910e-16`，壁面节点速度为 0，全部场值有限。原生及监督退出码均为 0；318 次线性求解收敛，没有 divergence 或 nonlinear-limit warning。运行身份保存最终场重新测量结果及五个区间，保存节律仍为每 10 步。

- [最终 VTU](manufactured_pipe/R1_GMSH_V2/cfd_gpu/run/8-procs/result_071.vtu)：SHA256 `095753aeb22f6c5684bcd14cb8c14e6f50c2334a57ded9805b848611ed0f16e2`。
- [最终 native checkpoint](manufactured_pipe/R1_GMSH_V2/cfd_gpu/run/8-procs/stFile_071.bin)：SHA256 `fda5282e4a7f3365f8b28f70c36417775de260dcc8b6644850783e7d86c615f5`。
- [求解器 execution](manufactured_pipe/R1_GMSH_V2/cfd_gpu/reports/execution.json)、[监督 receipt](monitoring/manufactured_gmsh_v2_gpu/EXECUTION.json)。

GPU CFD wall time **5334.385 s（88.91 min）**，低于 4 h 硬上限。GPU RTX 4090，driver 595.84，实际 CUDA runtime 13.2；5 s 采样的整卡峰值显存 4264 MiB、平均 GPU 利用率 2.45%，进程树平均 CPU 占用 5.88 核、峰值 RSS 10.525 GiB。整卡显存包含已有约 827 MiB 分配，不能当作本任务独占显存。这些仅为运行记录，没有 CPU/GPU 性能比较。

## 冻结探针科学门槛

960/960 原始探针严格包含于原 mesh，全部保留。五个 particle bands 各 96 点，使用原始 P1 barycentric u、逐单元完整 G；`G_ij=du_i/dx_j`，`s=(I-nnᵀ)Gn`，n 为解析圆柱向内法向。没有 gradient recovery、平滑、移动或外推。向量采用相对 L2，G 采用完整 Frobenius 相对 L2。

| 距解析壁面 µm | u（≤5%） | G（≤8%） | s（≤8%） | curl（≤8%） | 判定 |
|---:|---:|---:|---:|---:|---|
| 0.377 | 8.519477% | 4.424559% | 4.366459% | 4.411926% | FAIL（u） |
| 1 | 4.965430% | 4.440457% | 4.382446% | 4.427685% | PASS |
| 2 | 4.428310% | 4.471119% | 4.413316% | 4.458140% | PASS |
| 5 | 4.346174% | 4.661532% | 4.541800% | 4.640587% | PASS |
| 10 | 4.456588% | 4.995713% | 4.845931% | 4.971915% | PASS |

| Bulk 距壁 µm | u 相对 L2（≤5%） | 判定 |
|---:|---:|---|
| 29.838961 | 5.000639% | FAIL |
| 59.677923 | 5.591184% | FAIL |
| 149.194807 | 5.748240% | FAIL |
| 298.389615 | 4.784138% | PASS |

29.838961 µm 组的真实值为 5.000639%，超过 5%，即使保留两位小数显示 5.00% 也必须 FAIL。四组 bulk 均保留各 96 点，按原逐组规则验收。

三个内部截面 Q 相对误差分别为 2.653820%, 2.850242%, 2.582880%，均通过 5%；dp/dz 误差 **2.586847%**，通过 8%。WSS 误差 3.191699% 和 0.1 µm 带均仅作诊断，没有作为 hard gate。

## 结果含义与适用范围

本轮已经排除“尚未达到冻结稳态”和“GPU 后端没有真正启用”作为本次 field 未通过的解释。G/s/curl 的五个粒子带通过，但速度在最靠壁带及部分 bulk 带超限；全局守恒不能替代这些局部精度门槛。

既有诊断将误差写为 `u_h-u_exact=(u_h-Iu_exact)+(Iu_exact-u_exact)`。在 0.377 µm 带，解析节点值的同网格 P1 插值误差为 **9.758956%**，CFD 与该插值之差为 **4.688421%**，实际 CFD 总误差为 **8.519477%**。这说明有限网格/P1 表示本身有显著误差，同时 CFD 离散解与解析插值也有差异。各项范数不能相加，存在方向抵消，不能据此将总误差唯一归因于 mesher 或某个求解器环节。诊断插值没有替换 CFD 场。

此结论只认证本次冻结 V2 manufactured case 的失败，不等于已解释或修复全部 BraVa 异常轨迹。用户规定失败即停，因此没有改 mesh、BC、solver、threshold 或继续数值研究。尚无合格 particle field；oracle、CFD-vs-oracle 与 dt 收敛均未运行。

机器结果：[FIELD_VALIDATION.json](FIELD_VALIDATION.json)。原始采样：[FIELD_SAMPLES.npz](manufactured_pipe/R1_GMSH_V2/cfd_gpu/analysis/FIELD_SAMPLES.npz)。

可视化：[300 dpi PNG](gmsh_reference/figures/gpu_v2_field_validation.png) / [PNG](gmsh_reference/figures/gpu_v2_field_validation.png) / [HTML 入口](OPEN_RESULTS.html)。图同时显示稳态、全部 particle bands、bulk 与积分门槛，EQ=0 以文字明确标记，未伪造对数轴上的正数。

<!-- V2_ERROR_DIAGNOSIS -->
## 后续冻结场诊断（2026-10-04，fd6c068c）

本轮没有重跑CFD或改原gate。误差分解结论为 **MIXED_FIELD_ERROR**：0.377/149.195 µm支持P1表示限制，29.839/59.678 µm支持CFD−P1离散解差异；各项存在负交叉项，不能相减百分数。0.377 µm的P1误差以切向项为主，不能仅凭缩小wall altitude保证通过。

新授权的exact-vs-P1四场景诊断已完成，**P1_FIELD_TRAJECTORY_IMPACT_LARGE（DIAGNOSTIC_ONLY）**，由A的10 µm局部通过时间误差5.148892%触发；四组无径向分离、gap差或contact类别改变。它未使用CFD场，不是正式oracle/CFD trajectory certification；正式dt收敛仍未运行。

详见[误差分解与敏感性报告](V2_P1_REPRESENTATION_DIAGNOSIS_ZH.md)。[V3条件性策略](PROPOSED_V3_MESH_POLICY.json)仅为PROPOSAL_ONLY，未生成网格。原 `REFERENCE_PIPE_FIELD_FAIL_SCIENTIFIC` 继续有效。
