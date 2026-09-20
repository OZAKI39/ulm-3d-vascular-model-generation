# Frozen FEM → Particle-0 Handoff

FEM DEVELOPMENT: **FROZEN**。Frozen FEM performance baseline: **Stage SV1.3Q**，状态 **ILU_REBUILD_POLICY_WINNER_FOUND**。

Next work: **Particle-0 → Particle-8**。从现在开始 FEM 是只读背景流；本分支完成冻结交接，尚未开始任何粒子代码开发。

## 项目已接受的验证边界

Mesh convergence: **NOT PERFORMED BY USER DECISION / NOT PLANNED BEFORE PARTICLE DEVELOPMENT**。

Time-step sensitivity: **NOT PERFORMED BY USER DECISION / NOT PLANNED BEFORE PARTICLE DEVELOPMENT**。

CPU/GPU field-level equivalence: **DEFERRED BY USER DECISION**。

这些是 USER-ACCEPTED PROJECT DECISIONS，不是通过的验证。FEM performance tuning: frozen after Stage Q。Stage Q is the frozen particle-development background-flow baseline；不能解释为网格收敛、时间步收敛、CPU/GPU 科学等价或新的 FEM production approval。

## 单向耦合

第一版仅允许 FEM → RBC / MB。RBC / MB → FEM 不允许，粒子运动不会重新求流场。不要重跑 FEM、调整网格、时间步、PETSc 或 ILU 来启动 Particle-0。

## Frozen Field Contract

以下路径全部相对于本文件所在目录，都是分支中的实际文件，无需原 WSL 或 GPU 服务器。

| 用途 | 相对此目录的路径 | SHA256 |
|---|---|---|
| 背景流 | [frozen_reference/flow/steady_flow_stage_sv1_3q.vtu](frozen_reference/flow/steady_flow_stage_sv1_3q.vtu) | `373ff549430cab710d57eb4406f2e7588f7a5bb68ebaa68cf32f8da5662e26f1` |
| 正式体网格 | [frozen_reference/SV_MESH/mesh-complete.mesh.vtu](frozen_reference/SV_MESH/mesh-complete.mesh.vtu) | `1a7a69dcba52f0475b2280775921e4dedcdfefa67b7965096d86cb9dba38bdf9` |
| WALL | [frozen_reference/SV_MESH/mesh-surfaces/WALL.vtp](frozen_reference/SV_MESH/mesh-surfaces/WALL.vtp) | `06c5d1d1975cf9470096d9b7e560e2173c9733c847d5bc364aab7355ef12a11b` |
| OUTLET_03 | [frozen_reference/SV_MESH/mesh-surfaces/OUTLET_03.vtp](frozen_reference/SV_MESH/mesh-surfaces/OUTLET_03.vtp) | `9e1d8d387eaff15ea3421c85051dd2f194d8adf57c8a839b4d894c5957e8546e` |
| OUTLET_01 | [frozen_reference/SV_MESH/mesh-surfaces/OUTLET_01.vtp](frozen_reference/SV_MESH/mesh-surfaces/OUTLET_01.vtp) | `1a051e6e2313c9be3985b16bd2a060fa111e8e8f8035682ab29eaa1961acdfed` |
| INLET | [frozen_reference/SV_MESH/mesh-surfaces/INLET.vtp](frozen_reference/SV_MESH/mesh-surfaces/INLET.vtp) | `fba5b444da15ebb7b881612acde9de1056395ff4c2903dd0b52e45e57b3fa9f5` |
| OUTLET_02 | [frozen_reference/SV_MESH/mesh-surfaces/OUTLET_02.vtp](frozen_reference/SV_MESH/mesh-surfaces/OUTLET_02.vtp) | `7dac9192aafbe935a595f4424c952e6fb6922e879b8e3821e3dc314f5b48cf33` |

正式体网格：70363 节点，371402 四面体，每单元 4 个节点（VTK type 10）。外表面见 [frozen_reference/SV_MESH/mesh-complete.exterior.vtp](frozen_reference/SV_MESH/mesh-complete.exterior.vtp)；其哈希及全部边界见 [boundary manifest](frozen_reference/boundary_manifest.json)。

坐标使用原始 Cartesian x/y/z（右手系），单位 m；未平移、旋转或重新缩放，未赋予解剖轴标签。边界框 `[xmin,xmax,ymin,ymax,zmin,zmax]` = `[7.894550799392164e-05, 0.00017674970149528235, 4.1103659896180034e-05, 0.00013389615924097598, 7.978476787684485e-05, 0.0001598619855940342]` m。流场和正式网格节点坐标及顺序完全一致；逐行比较确认 tetra cell rows 也一一对应，差异只在单元内部的局部顶点排列（从 volume 到 flow 的排列为 `[1, 0, 2, 3]`）。`tetra_id` 约定为**正式 volume mesh 的零起始 cell index**；其 `GlobalElementID` 是从 1 开始。本冻结文件对的 cell index 对应关系已验证，未来替换文件时必须重新核对，不能只凭相同单元数推断。

冻结时刻为 step **71**，`TimeValue = 1.0678826650986389e-05` s；step 70 首次通过连续 5 区间稳态，完成当前时间步后在 step 71 安全停止。背景流按稳态固定使用；该时间是 FEM 状态的生成时刻，不是粒子必须采用的初始时刻。

`Velocity`：POINT association，Float64，3 分量，m/s。`Pressure`：POINT association，Float64，1 分量，Pa。VTU 无 cell arrays。导出的是四节点四面体上的节点表示；这说明导出数据格式，不额外断言求解器所有内部空间的离散阶次。完整信息见 [flow manifest](frozen_reference/flow/flow_field_manifest.json)。

velocity gradient、vorticity、strain rate：**NOT STORED IN FROZEN VTU**。Particle-0 必须从 tetra velocity field 推导梯度，并明确单元内插值、跨面梯度和边界处的取值规则。本次未实现这些功能。

## FlowSampler contract（未来接口，尚未实现）

输入 `position[3]`：原坐标系中的 SI 位置，单位 m。

| 返回值 | 来源与约定 |
|---|---|
| `velocity[3]` | 由节点 `Velocity` 在所定位四面体内插值，m/s |
| `pressure` | 由节点 `Pressure` 插值，Pa |
| `velocity_gradient[3,3]` | Particle-0 推导，约定 G[i,j] = ∂u_i/∂x_j，s⁻¹ |
| `vorticity[3]` | Particle-0 推导 curl(u) = (G[2,1]−G[1,2], G[0,2]−G[2,0], G[1,0]−G[0,1])，s⁻¹ |
| `strain_rate[3,3]` | Particle-0 推导 (G+Gᵀ)/2，s⁻¹ |
| `inside_lumen` | Particle-0 的点定位结果；边界容差必须显式定义 |
| `tetra_id` | 正式体网格零起始 cell index；域外建议为 -1 |

接口建议：域外明确返回 `inside_lumen=false`，不得默默外推并伪装为域内结果；其余域外返回值及异常规则由 Particle-0 固定并测试。不得把节点数组读取成功等同于已完成 FEM sampling 验证。

## 单位与物理条件

geometry/length = m；time = s；velocity = m/s；pressure = Pa；density = kg/m³；dynamic viscosity = Pa s；Q = m³/s。

实际配置：ρ = 1056 kg/m³，μ = 0.00345312 Pa s，Q target = 2.7369132390905703e-15 m³/s，dt = 1.5040600916882215e-07 s。此 Q 是参考数值条件，不是实验泵流量。运行 XML 的入口有符号值为 -2.736913239876193e-15，保留了既有入口离散归一化校正，不能替换为手抄的 Q target。

建议 3D Particle solver 内部核心计算全部保持 SI：length m、time s、velocity m/s、force N、torque N m、pressure Pa、viscosity Pa s。若显示 µm，只在 I/O / visualization 层显式用 `1 m = 10^6 µm` 转换，禁止隐式混用。

## Boundary semantics 与法向

| 语义 | 实际 face ID | 面积（m²） |
|---|---:|---:|
| WALL | 1 | 2.024190790127075e-09 |
| OUTLET_03 | 2 | 5.049058023679602e-12 |
| OUTLET_01 | 3 | 4.355871419037442e-12 |
| INLET | 4 | 7.756795804427715e-12 |
| OUTLET_02 | 5 | 4.238794052117287e-12 |

所有 accepted VTP 的三角形绕序已逐面核对：`n_out = cross(x1-x0,x2-x0)/|cross|` 指向流体域外，验证依据是所属四面体中心与面中心。入口 inward normal = `-n_out`。`∫u·n_out dA` 在入口为负、出口为正，正入口流量 Qin 是入口有符号积分的负值。壁面使用各三角形局部法向；面积加权平均法向不能代替局部壁面法向。

VTP 未存储 normals 数组，Particle-3 / Particle-7 应按上述约定计算；不要猜测方向。边界和外表面的三角形集合构成完整、不重叠的外边界。编号来自 [face_map.json](configs/face_map.json)，与实际 XML 和 Stage Q 输入哈希交叉核对。

## Solver provenance

PETSc **3.25.5**（commit `a0b507cc1cc9eac739c34d6e71521f4854e84e36`），CUDA **13.2**，Open MPI **4.1.6**（Fortran bindings=all），1 RTX 4090、MPI ranks=1、OMP=1。svMultiPhysics upstream commit `c3f0bb892b765b718f61069ecd9726dbc6d177fd`。

正式配置为 GMRES/right、restart100、rtol1e-10、atol1e-24、max_it2000、ASM overlap2、ILU(2)。RA：首次 fresh ILU 健康求解迭代数作为 I_ref；严格超过 1.5× 标记下次重建，age≥5 强制重建。实际完整运行中增长触发与失败恢复均未触发，不能据此证明 adaptive trigger 本身比固定 R5 更快。

实际 XML 中有旧 `petsc-jacobi` 初始标签，运行时 PETSC_OPTIONS 将其覆盖为 ASM/ILU；以 [实际 options](frozen_reference/run/PETSC_OPTIONS.txt) 和已核对的 [runtime manifest](frozen_reference/run/runtime_manifest.json) 为准。配置没有被“清理”成另一份新配置。环境证据见 [GPU manifest](frozen_reference/run/gpu_environment_manifest.json)。

源码在 [vendor/svMultiPhysics_stage_q](vendor/svMultiPhysics_stage_q)，生命周期、safe-stop、Stage P/Q ILU 修改链见 [PATCH_PROVENANCE.md](patches/PATCH_PROVENANCE.md)。未提交 solver binary、PETSc/CUDA/MPI 安装树或 build cache。最终 checkpoint 已提交，仅用于 FEM 可复现证据，不是 Particle-0 背景流输入。

## 新聊天第一件事

读取本文件和 [FROZEN_FEM_BASELINE.md](FROZEN_FEM_BASELINE.md)，运行下述只读完整性检查，再为 Particle-0 的采样接口、点定位和梯度验证制定实现计划；粒子实现留到新聊天。

在此 FEM 子目录创建自己的 Python≥3.11 环境并 `python -m pip install -e .`。已验证的读取依赖版本见 [python_reader_versions.json](sync_metadata/python_reader_versions.json)。

```bash
python -B scripts/fem_freeze_sync/validate_frozen.py
python -B scripts/fem_freeze_sync/verify_manifest.py
python -B -m pytest -q tests/test_fem_sync_manifest.py tests/test_frozen_*.py tests/test_particle_handoff_contract.py tests/test_no_particle_implementation_yet.py tests/test_main_history_preserved.py
```

这些命令不运行 CFD，不连接远端；新校验脚本从自己的路径定位仓库，运行目录无需是原 WSL 路径。历史开发脚本的移植限制、Stage Q 原测试命令及历史失败见 [PORTABILITY_NOTES.md](PORTABILITY_NOTES.md)。

## 2D reference（只读链接）

参考分支：[sync/ulm-microbubble-traj-gen-2d-20260917](https://github.com/OZAKI39/ulm-3d-vascular-model-generation/tree/sync/ulm-microbubble-traj-gen-2d-20260917)。WSL provenance：`/home/lzy/projects/ulm_microbubble_traj_gen_2D`。

可参考 `particle_field_sampling.py`、`particle_mobility.py`、`particle_mobility_transport.py`、`particle_collisions.py`、`particle_constrained_step.py`、`particle_inlet_flux.py`、`particle_counter_rng.py`。`red_blood_cell_transport.py` 是旧 reduced-order RBC model，Particle-0→8 不继续把它作为正式 RBC 动力学。此次没有复制、修改或 merge 2D 分支。

后续路线见 [PARTICLE_RESEARCH_ROADMAP.md](PARTICLE_RESEARCH_ROADMAP.md)。
