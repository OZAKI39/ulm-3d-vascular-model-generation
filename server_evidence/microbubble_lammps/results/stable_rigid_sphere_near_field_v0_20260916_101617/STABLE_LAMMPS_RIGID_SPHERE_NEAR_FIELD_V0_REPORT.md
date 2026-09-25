# STABLE_LAMMPS_RIGID_SPHERE_NEAR_FIELD_V0

两个刚性微泡靠近时，不仅会因为正面挤压液体而减速，还会因为擦肩滑动和旋转受到液体阻力。本模型考虑这些刚性球运动，但不允许微泡改变形状。

本轮结果为 **PASS_WITH_TWIST_PENDING**。normal、tangential、横向旋转、平移—旋转耦合、硬防重叠及 A–K 中除 C3 外的验证已完成。**轴向 twist 尚未通过匹配与实现验证，不能称完整旋转模型 PASS。**

## 实际完成的结果

- 真实旧 8 微泡使用相同直径、初始位置和冻结 Palabos 场，从零运行至 **0.02 s**。旧 passive 于 0.005675914432167097 s 停止，新结果最小球间隙 **0.0227822247886 µm**，最大角速度 **18.9253889548 rad/s**，未触发壁面停止。
- 强收敛 F 算例接触后继续完成，C_gap=0.4 下有 381 个实际约束步；按 gap>=-1e-12 m 门槛，全部算例接受态和连续位移段的重叠数为 0。约 -1e-21 m 的接触舍入误差保留在原始结果中，没有改写为零。
- MPI1/4 的位置、速度、角速度差异均为 0；MPI4 中记录到 100 次跨 rank 活跃球对观察。G 初始三个球对同时作用，H 交换 ID 后恢复物理对应仍一致。
- 另外复用旧均匀流案例的两个原始粒子和记录轨迹，完成非空球对均处于 cutoff 外的 J2 回归；共同步位置差异 2.11758e-22 m、速度差异 2.1684e-19 m/s，终点时间相同。独立重算最小间隙 5.28302e-06 m，大于 cutoff 1.35839e-07 m；结果见 `validation/NONEMPTY_PAIR_PASSIVE_REGRESSION.json`。
- 旧 passive 单球案例在共同的 754 个输出步上中心位置差异为 0，速度差异约 2.03e-20 m/s。新积分器累加时间产生一个终点舍入补齐小步，记录为 3011 个接受步；共同步仍对应旧 3010 步，最终时间相同，最终位置差异 2.17e-19 m。
- Kokkos 兼容短测完成 50 步。新增阻力解仍在主机执行；本轮没有进行 GPU 性能优化，GPU_PERFORMANCE_READY=NO。

## 稳定版、退役及完整性

唯一 active LAMMPS 为 **22Jul2025 Update 6**，commit `9c5ab448c78a14fd534619622162ba418d6a1fb1`。

```
CPU SHA256 = 188ebc166c967cb004af3df6e75b1903e4f5df7123eca46264899831a845bee4
GPU SHA256 = 82dd028a0fd648b5e2355c6b955a367d7bfac4847b0ed69640e5d101fed04f65
active link = /workspace/lammps_active
```

先保存候选版 16200 个文件、1727916384 字节的退役清单及 SHA256，再删除候选源码、编译目录、CPU/GPU 和 custom coupling 二进制、临时运行树。原来的 998 项失败证据完整保留。退役清单记录删除前事实；删除后结果见 `LAMMPS_ACTIVE_INSTALL_AUDIT.json` 及 `LAMMPS_2SEP2026_CANDIDATE_RETIREMENT_RECEIPT.json`。

稳定版 15195 项源码/构建资产重新校验，差异为 0；旧 engine、coupling、passive、PBS/BSA 和上轮迁移失败结果清单全部通过。WSL SonoVue sampler 原清单也通过。正式核心、Palabos 场、几何、旧 passive 和 sampler 未修改。

## 为什么修复原生数学

详细证据见 [源码审计](STABLE_LAMMPS_LUBRICATION_SOURCE_AUDIT.md)、[独立推导和近似范围](THEORY_DERIVATION.md)、[冻结 reference contract](contracts/RIGID_SPHERE_LUBRICATION_REFERENCE_CONTRACT.json)。

独立提取运行的原生 sq/sh 标量与修正的 XA/YA 一致。但异径球的 YB 不能用球半径乘 YA 替代；原生 pump 也不能直接作为包含 shear 转动力矩后的剩余阻力。`lubricateU/poly` 的完整旋转段还有自身系数差异。这些原始失败保留，未更换版本、放宽门槛或关闭必需模式。

新增 project-owned 模块从独立横向矩阵的 Schur 补构造 pump，并保留核对过的 sq/sh 标量。有限间隙使用明确冻结的工作共轭、客观性补全，使整体刚体运动不产生虚假内部阻力。**这是一套声明了闭合与截断的 V0 近场近似，并非精确双球全距离 Stokes 解。** 系数门槛验证实现与冻结公式一致，不是微泡实验精度声明。

Twist 的公开 XC 理论项包含非奇异集体/远场阻力，与保留的孤立 Stokes 转动阻力如何匹配仍未验证。Ball–Melrose 作者仓储无全文，出版社访问返回 403；未声称读到完整原文。可读取的独立推导、Jeffrey–Onishi 和 Townsend 修正资料已保留。因此 C3 明确为 BLOCKED_TWIST_NOT_IMPLEMENTED。

## 数学与数值证据

| 项目 | 实测 | 门槛/结果 |
|---|---:|---|
| 独立系数比较 | 100000 组；最大相对误差 4.53765e-14 | <=1e-10，PASS |
| 力反对称 | 10000 组合状态；4.31837e-15 | <=1e-12，PASS |
| 力矩和力臂总平衡 | 4.43716e-15 | <=1e-12，PASS |
| 负耗散违规 | 0 | PASS |
| 仿射场 C++ 梯度误差 | 9.82254e-11 s^-1 | <=1e-10，PASS |
| C++ / Python 矩阵差异 | 3.93158e-14 | <=1e-10，PASS |
| 小系统 N=2/3/5/10 | 80 组；条件数最大 196.068 | 正定；性质审计后选 PCG |
| PCG / dense 解差异 | 4.91845e-13 | <=1e-9，PASS |
| 求解器记录的最大残差 | 4.49768e-14 | <=1e-10，PASS |
| C_gap=0.4 vs 0.1，F 轨迹/间隙差 | 3.77193e-11 / 6.13418e-11 m | 两项 <=1e-8，选择 0.4 |
| 独立 STL 核查 | 17346 个状态/阶段点；无违规 | VTK 与 C++ 距离差 4.70103e-20 m |

C++ 使用模态 B^T D B；独立 Python 使用显式 3×3 矩阵块与 dense 解。约束由独立主动面枚举/KKT 核验。WSL finalizer 不加载 C++ 库，重新检查每个积分阶段、每个输出状态和 swept-sphere 距离，并用 VTK 检查内部点、真实 STL 距离与连续线段证书。

## 全部实际算例

| 算例 | 独立核查 | 接受步 | 终点 s | 最大阶段速度误差 m/s | 最大阶段角速度误差 rad/s |
|---|---|---:|---:|---:|---:|
| A_normal_cg04 | PASS | 400 | 0.02 | 1.57e-18 | 1.65e-12 |
| F_contact_cg04 | PASS | 400 | 0.0002 | 1.3e-15 | 1.62e-09 |
| A_normal_cg02 | PASS | 400 | 0.02 | 1.57e-18 | 1.65e-12 |
| F_contact_cg02 | PASS | 405 | 0.0002 | 1.03e-15 | 1.43e-09 |
| A_normal_cg01 | PASS | 400 | 0.02 | 1.57e-18 | 1.65e-12 |
| F_contact_cg01 | PASS | 430 | 0.0002 | 1.12e-15 | 1.63e-09 |
| B_tangential | PASS | 200 | 0.002 | 2.06e-18 | 1.25e-12 |
| C1_equal_transverse | PASS | 100 | 0.001 | 4.88e-19 | 9.09e-13 |
| C2_opposite_transverse | PASS | 100 | 0.001 | 3.56e-20 | 2.27e-13 |
| D_pass_by | PASS | 500 | 0.005 | 2.06e-18 | 1.25e-12 |
| E_arbitrary_pair | PASS | 100 | 0.001 | 5.08e-19 | 3.84e-13 |
| G_three_particle | PASS | 200 | 0.002 | 1.5e-17 | 1.95e-11 |
| H_permuted | PASS | 200 | 0.002 | 1.7e-17 | 1.57e-11 |
| I_mpi1 | PASS | 200 | 0.002 | 2.06e-18 | 1.25e-12 |
| I_mpi4 | PASS | 200 | 0.002 | 2.06e-18 | 1.25e-12 |
| J_passive | PASS | 3011 | 0.174331657559 | 5.42e-20 | 3.58e-15 |
| K_real_8 | PASS | 346 | 0.02 | 2.42e-18 | 6.88e-12 |
| kokkos_smoke | PASS | 50 | 0.0005 | 1.84e-18 | 1.25e-12 |
| J2_passive_separated_pair | PASS | 801 | 0.0376098661923 | 4.34e-19 | 0 |

C1/C2 的动态测试由外加横向力矩驱动，允许完整六自由度响应；另有纯指定运动的分模态矩阵核验。B 初始法向速度为 0；D 实际产生平移和转动。C3 未执行。总计 19 个求解案例，没有重复运行已经完成的 A 案例。

首个 A 求解完成后，报告 JSON 写入曾因 NumPy 布尔值类型失败。修复只转换报告类型，保留失败日志并对已有数据重做 finalizer；未重启该求解。源码提取器的一次开发期编译错误也已保留。

## 数据、图形与 ParaView

- [所有案例球对事件表](PAIR_HYDRODYNAMIC_EVENTS.csv)：normal/shear/pump 系数、力、力矩和约束标记；twist 明确写为 NOT_IMPLEMENTED。
- 每案例 `TRAJECTORIES.csv`、`INTEGRATION_STAGES.csv`、`SOLVER_HISTORY.csv`、运行日志、命令/二进制 SHA256 收据与独立核查均保留。
- 原始 nearest_gap 是 8 µm LAMMPS 搜索范围内的最近球对；`1e100` 表示该范围没有邻居。K 的最终最近间隙由离线坐标核查重新计算，存于 `validation/CASE_K_NEAREST_PAIR_GAPS.csv`，并写入 VTP，未把 sentinel 当作真实距离。
- [阻力模式](visualization/VIS_RESISTANCE_MODES.png)、[法向接近](visualization/VIS_NORMAL_APPROACH.png)、[切向滑动](visualization/VIS_TANGENTIAL_SLIDE.png)、[横向转动](visualization/VIS_ROTATIONAL_LUBRICATION.png)、[平移—旋转耦合](visualization/VIS_TRANSLATION_ROTATION_COUPLING.png)、[真实血管旧/新轨迹](visualization/VIS_REAL_8BUBBLE_RIGID_HYDRO.png)。法向图中 passive 和 lubrication-only 两条对照为明确标注的独立参考曲线；正式约束轨迹来自实际 LAMMPS 运行。
- [CASE_K_TRAJECTORIES.vtp](visualization/CASE_K_TRAJECTORIES.vtp) 含 ID、直径、时间、速度、角速度、最近间隙和活跃球对数，精确数组回读 PASS。可同时打开 `visualization/FROZEN_LUMEN.vtp`。轨迹没有平滑处理。

**HUMAN_VISUAL_REVIEW=PENDING**。自动生成与数组回读不代表人工科学审阅通过。

## 实现与限制

LAMMPS 负责 sphere storage、8 µm full neighbor list、ghost、MPI ownership/迁移和 ID；新增 embedding/custom fix 负责瞬时 6N 阻力解和 RK2 更新。正式粒子对来自 LAMMPS 邻居表，没有用 O(N²) 穷举替代生产邻居搜索。当前小规模验证采用复制的 dense 全局阻力存储和 PCG，不构成大粒子数的性能证明。

数值正则化 h_min=0.001*Reff、作用区 h<0.2*Reff、局部仿射背景闭合和有限间隙客观性补全均在测试前冻结。它们不是 shell、弹性或实测粗糙度。壁面仍仅作 3dx 安全停止，没有实现壁面水动力。adhesion、RBC 和形变均未进入本阶段。

## 归档

远端：`/workspace/microbubble_lammps/results/stable_rigid_sphere_near_field_v0_20260916_101617`

WSL：`/home/lzy/projects/compre_output/stable_rigid_sphere_near_field_v0/20260916_101617`

REMOTE_TO_WSL_INTEGRITY = **PASS**。最终以根目录 `SHA256SUMS` 和控制目录中的双端核验收据为准。没有 GitHub 推送，也没有自动进入下一阶段。

NEXT_STAGE = USER_REVIEW_BEFORE_MICROBUBBLE_WALL_HYDRODYNAMICS_V0
