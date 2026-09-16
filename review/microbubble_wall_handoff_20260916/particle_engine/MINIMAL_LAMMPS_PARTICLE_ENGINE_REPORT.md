# MINIMAL_LAMMPS_PARTICLE_ENGINE_V0 核查报告

阶段状态：**PASS**。本阶段只验证有限尺寸刚性球的软件工程能力。真实 SonoVue 微泡间相互作用、lubrication、流体阻力、壁面、黏附和流固耦合均未验证；RBC 关闭。

运行采用 Vast RTX4090（CC 8.9），独立 work/result 目录。未运行 Palabos，未修改正式纯流体基线、Stage4、PBS/BSA、Step3C 或 RBC 资产。

LAMMPS 版本：Update 6 for Stable release 22 July 2025；tag `stable_22Jul2025_update6`；commit `9c5ab448c78a14fd534619622162ba418d6a1fb1`。官方来源：[release](https://github.com/lammps/lammps/releases/tag/stable_22Jul2025_update6)。源码压缩包 SHA256：`34a2526440d52f220f86d9c6537c471184edd25e1fbaa655ffa0c9d41c0ab78a`。两套 binary 共用同一份源码。官方压缩包中 15013 个普通文件已逐一比对 SHA256，零差异。

CPU binary SHA256：`188ebc166c967cb004af3df6e75b1903e4f5df7123eca46264899831a845bee4`。GPU binary SHA256：`82dd028a0fd648b5e2355c6b955a367d7bfac4847b0ed69640e5d101fed04f65`。CMake 构建均启用 GRANULAR；GPU 另外启用 KOKKOS CUDA，`Kokkos_ARCH_ADA89=ON` 来自该版本源码的架构声明。LATBOLTZ 显式关闭。编译使用 gcc/g++ 13.3.0、CUDA nvcc 13.2.86、OpenMPI 4.1.6、CMake 3.28.3；并行编译限制 8。完整 flag、cache、帮助和日志在 `build_provenance/`。

单位为 SI；原始直径 µm 仅乘 `1e-6` 一次。`atom_style sphere` 保存 diameter、radius、mass、位置、速度、omega。密度 1000 kg/m³ 是技术测试值，**不是 SonoVue 微泡真实密度**；质量按球体积计算，半径固定。

粒径由原冻结 sampler 在原 WSL NumPy 1.26.4 环境生成；远端只读取结果。每组输入绑定 sampler source、contract、histogram、N、seed 和 population CSV 哈希。没有新采样逻辑、重新加权、裁尾或粒径调整。Case D 的 1000 颗为 ENGINE_TEST_ONLY；Case F 的 50000 颗为 GPU_ENGINE_SMOKE_ONLY，均不是实验浓度或正式科研 population。

首个 run 0 前完成单位、阶段和技术接触合同封存。技术参数：Kn=0.001 N/m，Kt=gamma_n=gamma_t=friction=0，dampflag=0；不启用 limit_damping。初始 overlap=2.7898202830022029e-08 m，dt=1e-08 s，200 步。预检 dt·ω=0.010508498 < 0.05。CPU `gran/hooke/history` 与 GPU counterpart 的语法和 half-neighbor 要求均据冻结源码核查。该接触参数不表示 SonoVue 材料特性。

| 算例 | 状态 | N | steps | MPI / engine | wall s | loop steps/s | peak RSS KiB |
|---|---|---:|---:|---|---:|---:|---:|
| diameter_roundtrip | PASS | 1000 | 0 | 1 / cpu | 1.00796 | 0 | 104292 |
| case_a_single_ballistic | PASS | 1 | 200 | 1 / cpu | 0.805831 | 47621.4 | 99424 |
| case_b_two_separated | PASS | 2 | 200 | 1 / cpu | 0.906378 | 7192.54 | 102932 |
| case_c_technical_contact | PASS | 2 | 200 | 1 / cpu | 0.704363 | 7419.11 | 107704 |
| case_d_polydisperse_1000 | PASS | 1000 | 200 | 1 / cpu | 1.00692 | 7704.78 | 104000 |
| case_e_mpi_migration/cpu_mpi1 | PASS | 1000 | 1000 | 1 / cpu | 0.804404 | 7376.05 | 104984 |
| case_e_mpi_migration/cpu_mpi4 | PASS | 1000 | 1000 | 4 / cpu | 1.10868 | 11832.1 | 162592 |
| case_f_kokkos_gpu | PASS | 50000 | 1000 | 1 / gpu | 4.24988 | 1575.82 | 1.26064e+06 |
| cpu_gpu_comparison/cpu | PASS | 1000 | 200 | 1 / cpu | 1.00834 | 4933.34 | 103692 |
| cpu_gpu_comparison/gpu | PASS | 1000 | 200 | 1 / gpu | 1.21171 | 897.759 | 264680 |

wall time 包括 MPI 启动、初始化、稀疏输出；F 还包含 profiler 开销，不能用于 CPU/GPU 性能比较。RSS 为每 0.1 s 采样的进程树合计峰值，短时峰值可能漏采。loop 时间来自对应实际 steps/count 的 LAMMPS 日志，逐项记录在 runtime CSV。

直径 round-trip 最大误差：4.440892098500626e-16 µm（门槛 1e-10）。A 解析位置最大误差：4.811147140404426e-19 m。B 检查静止位置、速度和零接触力。C 初始作用力对称相对误差：0.0，初始 overlap=2.789820283002137e-08 m，最终 overlap=0 m；另外检查初始 Hooke 力、排斥方向和能量尺度速度上界。

E 使用相同输入、固定 4×1×1 分解，初始 ownership={'0': 200, '1': 300, '2': 200, '3': 300}，最终 ownership={'0': 200, '1': 300, '2': 300, '3': 200}，实际 owner 改变粒子数=1000。证据来自 dump 的 proc 字段，而非只凭 mpirun rank 数。MPI1/MPI4 最大位置差=0.0 m；使用 unwrapped 坐标。

CPU/GPU 无接触比较：位置最大差=4.1338612898456646e-23 m，速度最大差=0.0 m/s，ID 与 diameter 按 ID 排序精确匹配。没有使用接触算例强行要求 CPU/GPU 逐点一致。

F 的 CUDA kernel 证据：{'status': 'PASS', 'kernel_count': 28205, 'nve_sphere_kernel_count': 1010, 'gran_hooke_history_kernel_count': 1002, 'evidence_scope': 'Case F complete launch including initialization and sparse outputs; no production performance claim'}。原始 trace 在 `cases/case_f_kokkos_gpu/gpu_trace.nsys-rep`，SQLite 导出同目录，按 kernel 汇总在 `validation/GPU_KERNEL_SUMMARY.csv`。GPU VRAM 峰值=576.0 MiB；GPU-resident 样本利用率中位数=3.0%，全部采样峰值=61.0%。排除 context 建立/退出端点后，驻留样本中央一半共 7 个，VRAM 范围=2.0 MiB（冻结门槛 256 MiB）。该稳定性只覆盖本次短 smoke，不代表长程无内存增长。

日志中的 Dangerous builds 为 not checked：因为冻结配置采用 every 1、delay 0、check no，每一步无条件重建 neighbor list；本轮不把缺失计数记作零。独立补充核查会验证该日志声明和实际 200/1000 次重建次数。冻结 finalizer 中对应布尔值只表示没有发现非零计数，不是已测得危险重建次数为零。

独立 finalizer 不导入 converter，也不以 LAMMPS 的成功退出代替数值判断。它重新读取 population、实际 data、initial/final dumps、log、MPI ownership 和 CUDA activity，逐项核验。采样输出与 thermo 日志的 nonfinite 检查、particle count/ID 和 lost error 联合使用；本阶段没有高频 full dump。

独立数值核查=PASS，记录到 nonfinite=0，lost atoms=0。远端到 WSL 哈希核验=PASS。下载只包含紧凑证据与代码，不含源码包、build tree、checkpoint 或 restart。

构建后的帮助启动曾在 hwloc GL 显示扫描中等待：strace 显示尝试 X0–X6 并在 localhost:6006 的 X11 握手阻塞。移除 DISPLAY 的排查无效（该环境本身未设 DISPLAY）；仅在任务子进程设置 HWLOC_COMPONENTS=-gl 后，MPI1 GPU 帮助在 0.95 s 完成。此设置只关闭显示拓扑探测，不关闭 CUDA。原始失败/超时排查证据保留，所有排查均为零 timestep，未重新编译 binary。依据：[hwloc 官方组件选择说明](https://www.open-mpi.org/projects/hwloc/doc/v2.0.4/a00324.php)。

预检期间修正过与冻结版本打印格式不匹配的 KOKKOS 日志识别器，发生于首次正式 run 前，保留修改前后哈希；未改数值门槛。不存在为通过结果而调整接触参数、粒径或容差的操作。

TECHNICAL_CONTACT_ENGINE = PASS；BUBBLE_BUBBLE_PHYSICS_VALIDATED = NO；BUBBLE_BUBBLE_INTERACTION_MODEL = PENDING；BUBBLE_BUBBLE_LUBRICATION = PENDING；SONOVUE_CONTACT_MATERIAL_PROPERTIES = NOT_DEFINED。

MICROBUBBLE_TRANSPORT_MODEL、PALABOS_LAMMPS_COUPLING、HYDRODYNAMIC_DRAG、MICROBUBBLE_WALL_MODEL、MICROBUBBLE_ADHESION 均为 PENDING。本阶段完成后停止，等待用户审核；没有自动进入耦合开发。

远端结果：`/workspace/microbubble_lammps/results/lammps_particle_engine_20260915_214759`。本地结果：`/home/lzy/projects/compre_output/lammps_particle_engine/20260915_214759`。
