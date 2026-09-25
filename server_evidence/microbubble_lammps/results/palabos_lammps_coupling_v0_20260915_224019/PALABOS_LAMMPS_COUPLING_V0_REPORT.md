# PALABOS_LAMMPS_COUPLING_V0 核查报告

本阶段 **PASS**：在冻结速度场与技术性 Stokes 阻力条件下，已验证单位/坐标映射、三线性插值、力注入、完整粒子推进、MPI 一致性和 Kokkos host-fix 兼容性。独立 Python 复核在 Vast 和下载后的 WSL 均通过。

这些证据只支持工程耦合链正确。真实场是 PBS/BSA 5000 步瞬态 smoke，尚未证明长期稳态；入口倍率在新介质中未重新校准，原 smoke 观察到瞬态回流。技术密度 1000 kg/m³ 和 Stokes 力仅供数值验证，不能据此声称 SonoVue 水动力学、真实输运、margination、壁面停留或分支选择已验证。

本轮没有执行 Palabos 求解，没有改变流体参数、STL、入口倍率、出口条件或原有工程源码。bubble-bubble、lubrication、wall hydrodynamics、adhesion 与最终 transport model 均为 PENDING；two-way coupling 与 RBC 为 OFF。

## 冻结输入与实现

- Palabos GPU commit：`4127697e90169bbef982295f1d1c933cf6e90caa`；沿用已核验 smoke 的构建 provenance。
- LAMMPS：22Jul2025 Update 6，commit `9c5ab448c78a14fd534619622162ba418d6a1fb1`。复用原 CPU/GPU 静态库，外部 embedding 程序注册 custom Fix，未修改 core 或 baseline binary。动态 plugin 未使用；共享库、Fix 源文件和可执行文件均有 SHA256。
- SonoVue histogram SHA256：`2c9f783c6169421b06c057e16653878e7385139a179682c0648519e72a9f5198`。原 sampler 和 contract 逐文件核验，所有 population 保留原 seed、CSV、metadata 与 source IDs。
- 唯一真实场：`fields/FROZEN_FLOW_FIELD_V0.h5`，SHA256 `7828ed32b0cdd1e0e85836a564b6737128e7dbe30181ddd43110233b37b64a7f`，3,957,891 bytes，182,694 个流体节点。稀疏索引排序且唯一；完整 fluid lookup 区分 SOLID 与 MISSING。
- 原始 `fields_5000.bin` SHA256：`abf6d649abad347a5ae43ec10978c8db6f3330dfc723b47ad4dfe078b773a87f`。原快照保存 PackedNode 的三分量 LU 速度；按冻结 dx/dt 一次转换为 `Velocity_m_s`。独立 finalizer 对全部 ID 和向量重新比较。没有读取或替代使用无效 `VelocityMagnitude_m_s`。
- `dx=1.9989918081065344e-07 m`；origin `[8.00903149195763e-05, 5.5701359566515553e-05, 0.00011143341591952393] m`；dims 493×497×280。公式 `x_m=origin+(ix,iy,iz)*dx`；展平 `id=ix+nx*(iy+ny*iz)`。原 VTI CellData 几何中心对应 Palabos node，不另加半格偏移。
- synthetic uniform/linear/missing-corner 与真实场使用同一 HDF5 schema 和同一 C++ sampler 路径。无 nearest、clamp、silent extrapolation；八角必须全部合法流体节点。
- `F=3*pi*mu*d*(u_f-v_p)`，mu=0.001 Pa·s；SI 单位。每步在 POST_FORCE 内存注入 atom->f。独立核验 dump 中实际力，以及采样流速、直径与 force-evaluation velocity。
- 保留原 `fix nve/sphere`。速度依赖力在其标准半步速度阶段计算；轨迹同时保存 force_v、applied_F 和完整步状态的诊断 F，避免把两种速度时相混为一谈。
- 冻结 particle dt=6.2499999999999991e-11 s。支持下限给出 tau_min=3.125e-8 s、tau_max=1.53125e-6 s；采用 tau_min/500，满足 <=tau_min/50。该选择在任何正式数值输出产生前完成，未按结果改动。
- 冻结清单包含 86 个输入/代码文件，所有正式运行与后处理均重新验证。初次 supervisor 启动遇到 Unix socket 路径过长，求解器启动次数为 0；只缩短控制 socket 路径后启动。本轮 9 个算例各执行一次，没有数值重试。

## 插值、坐标和力硬门

| 检查 | 最大误差 | 阈值 | 结果 |
|---|---:|---:|---|
| 256 个原 VTI 节点坐标 | 2.71050543e-20 m | 1e-12 m | PASS |
| 10,000 个 uniform 查询 | 4.33680869e-19 m/s | 1e-13 | PASS |
| 10,000 个 linear 查询 vs 解析 | 1.77635684e-15 m/s | 1e-12 | PASS |
| linear C++ vs Python | 8.8817842e-16 m/s | 1e-12 | PASS |
| 1,000 个真实节点 identity | 9.48676901e-18 m/s | 1e-12 | PASS |
| 10,000 个真实 interior C++ vs Python | 3.25260652e-19 m/s | 1e-12 | PASS |
| 10,000 组力函数相对误差 | 0 | 1e-12 | PASS |
| 其中 100 组零力绝对误差 | 0 N | 1e-25 N | PASS |

outside x/y/z、NaN、缺失角点、真实 solid 六类查询均返回明确 invalid 状态；正式 Fix 遇到 invalid 会中止，不使用 CLI invalid 行中的数值占位作为流速。

## 完整轨迹结果

| 算例 | N / MPI | steps | 物理时长 (s) | 归一速度误差 | 归一位置误差 | 结果 |
|---|---:|---:|---:|---:|---:|---|
| case_a_uniform | 1 / 1 | 55035 | 3.4396875e-06 | 1.67112139e-05 | 4.35898938e-05 | PASS |
| case_b_linear | 1 / 1 | 76653 | 4.7908125e-06 | 1.20814192e-05 | 3.35601116e-05 | PASS |
| case_c_real | 1 / 1 | 17331 | 1.0831875e-06 | — | — | PASS |
| case_d_d10 | 1 / 1 | 7590 | 4.74375e-07 | 0.000121174537 | 0.000316099501 | PASS |
| case_d_d50 | 1 / 1 | 16670 | 1.041875e-06 | 5.51717499e-05 | 0.000143915598 | PASS |
| case_d_d90 | 1 / 1 | 41130 | 2.570625e-06 | 2.23607733e-05 | 5.83265371e-05 | PASS |
| case_e_mpi1 | 1000 / 1 | 1000 | 6.25e-08 | — | — | PASS |
| case_e_mpi4 | 1000 / 4 | 1000 | 6.25e-08 | — | — | PASS |
| kokkos_compatibility | 1000 / 1 | 1000 | 6.25e-08 | — | — | PASS |

A/D 使用整个 trajectory 与解析解比较：速度误差除以 |u_f-v0|，位置误差除以 |u_f-v0|·tau_p；所有初速为零。位置误差没有除以全局坐标或大 box 长度。A/D 门为 1e-3。B 使用 DOP853（rtol=1e-12，位置 atol=1e-20 m、速度 atol=1e-15 m/s），相同归一方式，门为 2e-3。

Case A 第一颗 seed=20260918 的直径为 3.51891177935619 µm，tau=6.87930006160651e-07 s；末态 slip/initial slip=0.00673604978，覆盖 >=5 tau。

Case D 根据 seed=20260920 的 10,000 个原样本选择最近 d10/d50/d90：

| 分位目标 | 实际直径 (µm) | 公式 tau (s) | 数值 1/e crossing (s) |
|---|---:|---:|---:|
| d10 | 1.30676175165544 | 9.48681264216e-08 | 9.48368881867e-08 |
| d50 | 1.93664587354855 | 2.08366513307e-07 | 2.08335268603e-07 |
| d90 | 3.04206446360301 | 5.14119788929e-07 | 5.14088541118e-07 |

三者解析轨迹均通过；公式 tau 与 d² 一致，测得的 1/e crossing 严格递增。crossing 使用相邻轨迹记录插值，未拟合或更改 drag 参数。

Case E 使用 seed=20260921 的 1000 粒子，同一初态分别执行 MPI1/MPI4；最终 ID 和直径逐项完全相同，最大位置差 0 m，最大速度差 0 m/s。初始周期中心间距和保守运动上界证明整个短程无接触；表面间隙下界 2.83553315 µm。没有粒子间物理作用。

## Case C 的几何与科学边界

seed=20260919 的 100 个样本中，先按最接近冻结 d50 的规则选定 1.97471353346382 µm 粒子，之后才做几何放置。初始 center-wall 距离 2.38282920197874 µm；要求 radius+3dx=1.58705430916387 µm。粒径、seed、STL、安全 margin 均未为几何适配而改变。

冻结短程 17331 步，1.0831875e-06 s。总位移 5.74441972854e-11 m，说明真实场链路实际驱动了粒子，但路径非常短，仅用于工程查询和推进检查。

每一步均检查有限值、完整八角查询和 Lipschitz 保守壁距下界。下载后又对全部 17332 个记录用原闭合 STL 计算精确三角面距离及 enclosure：全部在腔内；最小 center-wall 距离 2.38282920198 µm；最小球表面到壁间隙 1.39547243525 µm，大于 3dx=0.599697542432 µm。终端 WALL_CLEARANCE 指球表面间隙，不是中心到壁距离。

该安全规则是防止当前无壁面模型的工程算例进入未定义区域；它不是 wall force，不提供真实 bubble wall-interaction 结论。

## Kokkos 兼容性与资源证据

同一 1000 粒子短例由 GPU/Kokkos binary 执行，非有限值、丢粒子、invalid query 均为 0。NVML 100 ms 采样记录的 utilization 峰值 33%，VRAM 峰值 492 MiB。采样可能遗漏短峰值，不据此评价性能。

保留原始 `kokkos_host_fix_trace.nsys-rep`、SQLite 和 CUDA API/memcpy 明细。trace 涵盖整次 invocation（含初始化和 run 0），不能把这些总量直接标成纯每步性能：

| 观测 | 次数 | 数据/累计时长 |
|---|---:|---:|
| H2D | 27399 | 4,831,070,612 bytes；0.40145 s |
| D2H | 22191 | 4,490,300,876 bytes；0.350662 s |
| cudaDeviceSynchronize_v3020 | 94094 | 0.563311 s API 时间 |
| cudaStreamSynchronize_v3020 | 1180 | 0.0393838 s API 时间 |
| cudaEventSynchronize_v3020 | 2012 | 0.00160747 s API 时间 |

这些计数与保守 host Fix 的 ALL_MASK 同步路径一致。CUDA kernel 累计时间约 0.0579732 s；API、copy、kernel 时间可能重叠，不能相加当总 wall time。没有改写 Kokkos kernel 或进行 GPU 优化；GPU_PERFORMANCE_READY=NO，COUPLING_GPU_ACCELERATION=PENDING。

## 完整性、大小与复核入口

远端与 WSL 独立 finalizer 均为 PASS；全部数值记录中 NONFINITE_COUNT=0、LOST_ATOMS=0。原 LAMMPS 源码/二进制在前后 inventory 中一致，原 LAMMPS 与 PBS/BSA 结果清单远端再次通过；WSL 原 sampler、LAMMPS 结果、PBS/BSA 结果也逐文件通过。详情见 PALABOS_LAMMPS_COUPLING_PROVENANCE.json。

远端数值证据共 233 个校验项，下载后 `sha256sum -c SHA256SUMS` 全部通过；NUMERICAL_EVIDENCE_SHA256SUMS 保留第一次数值封存，最终 SHA256SUMS 另覆盖本地审计与报告。结果总量约 0.281 GiB（逻辑文件字节数，四舍五入），小于 2 GiB；真实 HDF5 仅一份。

- `contracts/`：物理、选择、时步、geometry safety 和误差门。
- `src/`：C++ library/Fix、Python reference 与独立 finalizer。
- `cases/`：原输入、初末 dump、全/稀疏轨迹、每步计数器、执行日志和资源 trace。
- `validation/`：原始查询、C++ 输出、原 VTI 对照、单元审计。
- `LOCAL_INDEPENDENT_FINALIZER.json`、`LOCAL_EXACT_GEOMETRY_AUDIT.json`：本地独立最终核查。
- `FINAL_TERMINAL_SUMMARY.txt`：完整用户要求的终端字段。

远端：`/workspace/microbubble_lammps/results/palabos_lammps_coupling_v0_20260915_224019`。
本地：`/home/lzy/projects/compre_output/palabos_lammps_coupling_v0/20260915_224019`。

本阶段到此停止。NEXT_STAGE=USER_REVIEW_BEFORE_BUBBLE_BUBBLE_INTERACTION_V0。
