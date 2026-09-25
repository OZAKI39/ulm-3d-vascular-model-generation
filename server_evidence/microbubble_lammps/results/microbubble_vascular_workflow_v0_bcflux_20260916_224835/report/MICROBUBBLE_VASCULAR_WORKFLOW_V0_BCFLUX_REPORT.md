# MICROBUBBLE_VASCULAR_WORKFLOW_V0 — BCFLUX

本轮已实现并实际运行自适应入口注入、stable rigid RK2/LAMMPS 和独立验证。**总体技术状态 PARTIAL**：LOW/MEDIUM 遇到有限半径球的几何侧壁停滞，没有自然出口事件；这些结果保持原样，没有加入壁面力或删除滞留粒子。原来的 16.4% 入口 sampler 缺口不再作为 total Q 的阻塞。

**NOT EXPERIMENTAL CONCENTRATION**。全部动力学用例采用 TEST_ONLY_DIRECT_NUMBER_FLUX；生产微泡浓度仍是 UNSPECIFIED。

1. **入口总流量现在从哪里获得？** 使用当前 CFD 保存并传入 PBS/BSA Palabos 的物理边界目标 **2.7369132390905703×10⁻¹⁵ m³/s**。实际 solver 参数、运行 readback、入口命令和冻结场来源哈希一致。它是权威的 case 边界目标，不是实验实测或瞬态全截面实测流量。
2. **为什么不再从 sparse Frozen Flow 积分整个入口？** 该采样器在约 16.4% 的真实近边界区域缺少完整八角点支持。全截面速度没有定义，不能假装补齐后积分。总 Q 已有独立的已核验边界来源。
3. **16.4% 缺口现在影响什么？** 影响可选注入中心的位置支持和空间分布可信度；不改变 total Q，也不缩减目标 number flux。有效面积约 83.6% 还要再满足 inward velocity、有限半径净空及无重叠要求，实际可选中心区域会更小。
4. **微泡怎样自适应进入？** rank0 累计 λ×dt，先重试 FIFO，再生成新的整数需求。准入发生在 timestep 起点并立即参与两次 RK2 solve。每步可为 0、1 或多颗；没有 fixed cadence、one-per-step 或 constant-N controller。
5. **怎样处理过大 bubble？** 使用冻结 SonoVue histogram/inverse CDF，直径转半径。按当前“一格入口净空”契约，半径超过约 1.370125 µm 的认证上界永久记为 SIZE_INADMISSIBLE_AT_INLET，不缩小、不改原 draw、不进入临时 backlog。
6. **怎样处理临时拥堵？** 保留 ID、source uniform draw、半径和出生时间，FIFO 重新找位置。容量档在 pending=21 超过门限20时显式停止；不存在 silent drop。若随机位置搜索没有测试到任何几何可行点，则显式报搜索预算失败，不能冒充拥堵或尺寸不合格。
7. **有没有真正运行 RK2/LAMMPS？** 有。LOW 2418步、MEDIUM 814步、容量档143步、加权比较1243步，均执行实际 native LAMMPS insertion/storage/neighbors 和 stable rigid 两阶段 midpoint RK2；每个用例都跑 MPI1/MPI4。MEDIUM 有122条真实 BB 水动力记录。
8. **微泡有没有从 outlet 出去？** 正常入口输运的出口和 inlet-backflow 计数都是 **0**。另有独立几何/LAMMPS fixture 验证四个 port 的 outward crossing、反向排除和实际删除；fixture 不是流场输运出口，也未计入 flux_timeseries。
9. **MPI1/MPI4 是否一致？** 四组配对运行的 source draws、ID、注入时间、尺寸拒绝、pending、RK2 位置、near-wall/port 日志及终态逐项一致，位置差为0。当前短轨迹都留在同一子域；独立 fixture 另验证粒子从 owner 2 迁到 3 后 ID、radius、omega、位置和删除计数不变。

## 流量来源链：必须保留的替代关系

[AUTHORITATIVE_INLET_FLOW_AUDIT.json](../validation/AUTHORITATIVE_INLET_FLOW_AUDIT.json) 给出原始绝对路径、保存副本、run identity、port、单位和 SHA256。链路如下：

- Global 1D 的 global edge814 → anchor003274 的 cut_000，ROI `flow_rate_m3_s=7.693508475538942e-16`，原 profile=PARABOLIC。
- CFD surface preparation 的 original/final BC 保留该旧 Q。axis-aligned geometry 的端口身份也保持一致；workflow inlet boundary_index=0、CellEntityId=4。
- **当前 `configs/cfd_flow.yaml` 显式将 target_volume_flow 改为 2.7369132390905703e-15**。保存的 accepted reference-scaled runtime contract 与 Stage4 `physical_bc_contract.json` 一致；后者明确写明 `current target supersedes earlier s3/s4 Q`。因此不能声称旧 1D Q 原封不动传到了 PBS/BSA。
- Stage4 lattice target → PBS/BSA NUMERICS_INPUT → generated solver_parameters → actual runtime readback 均保留当前 Q。PBS/BSA 换密度/黏度时保持体积流量，不把旧介质的质量流量直接当新介质目标。
- actual run_5000 的 numerics contract SHA 为 `05b1abaa757535c6d62006d09823343c97912c8f0777f10eb17b6f658b141936`；fields_5000.bin SHA 为 `abf6d649abad347a5ae43ec10978c8db6f3330dfc723b47ad4dfe078b773a87f`，与 Frozen Flow contract 对应。

Palabos 实际命令是 uniform normal nominal Q/native-cap-area，附带 ramp 与数值倍率 **1.1197286861799598**。数值 profile 积分 **3.0646002653954225e-15 m³/s** 与物理 target 分开存储；本轮 provider 不把倍率再乘到物理 Q 上。新介质倍率的实验/稳态校准仍 NOT_PERFORMED，不升级原 transient field 的科学状态。

**A：TOTAL INLET FLOW RATE = KNOWN authoritative case target。**

**B：LOCAL FROZEN-FIELD COVERAGE = PARTIAL，83.623545% supported / 16.376455% unsupported。**

旧积分完整保存在 provenance；[局部交叉检查](../validation/FROZEN_FIELD_LOCAL_CROSSCHECK.json) 保留 VALID-area Q_positive=2.471611096426182e-15 m³/s，角色严格是 FROZEN_FIELD_LOCAL_CROSSCHECK_ONLY。没有将 unknown ring 设为零、nearest、clamp、外推，或放大局部积分来凑 total Q。

## 实现与冻结契约

[WORKFLOW_CONTRACT.json](../contracts/WORKFLOW_CONTRACT.json) 冻结顺序、单位、FIFO、两种 population 语义、rollback、端口和终止规则。源码入口为 [workflow_main.cpp](../src/workflow_main.cpp)、[adaptive_injection.hpp](../src/adaptive_injection.hpp)、[workflow_geometry.hpp](../src/workflow_geometry.hpp)。

`AuthoritativeBoundaryFlowProvider.volume_flux(time)` 当前返回常数 Q；NumberFluxProvider 同时支持 C_b×Q 和 TEST_ONLY_DIRECT_NUMBER_FLUX。生产浓度未指定，因此没有启用真实 C_b。SOURCE_POPULATION 的整数目标控制源 draw，永久 size rejection 消耗源需求；ADMITTED_POPULATION 的目标控制 admitted+pending，拒绝后产生新的候选身份。二者都保留 `target-admitted` 原始欠额，不把永久拒绝偷偷变成 backlog。

在每次 trial step 中，rank0 复制 controller、位置 RNG 和 FIFO；更新累计后实际插入 LAMMPS，再组装 neighbor list 和执行 RK2。拒绝则删除本次新增 tags，核对原 active state 未变，丢弃 trial 的抽样与累计并用更小 dt 重试。只有 accepted step 提交 controller 和注入事件。试算日志留在 RETRY_HISTORY；重新尝试不会双计。注入时间是 t_n，对应的需求分配窗口为 [t_n,t_n+dt]，连续整数阈值最多提前一个 accepted timestep。

所有准入点都要求原生 sampler VALID、u·n_in>0、WALL_ONLY 距离≥a+dx、与 ACTIVE 和当步新增粒子无重叠。默认按实际切割截面的三角形面积均匀提议后拒绝；加权模式再按正向 u_n 做拒绝采样，包络来自 native nodes 投影速度上界，因此不重构速度。两者都是有效支持与有限半径约束下的条件分布，不是 FULL_PHYSICAL_INLET_SAMPLING。

尺寸认证使用 distance-to-WALL_ONLY 的 1-Lipschitz 性质，对完整截面三角形分支细分。最大距离下/上界为 1.5700142623 / 1.5700241672 µm；扣除 dx=0.1998991808 µm 得上述半径界，间隔约1e-11 m。落入未决窄带会显式停止，不能靠有限次随机失败认定永久过大。[认证记录](../validation/INLET_SIZE_CAPACITY.json) 保留 witness 与计算预算。原 histogram、sampler、rigid_math、frozen_flow 和几何 STL 均未改。

## 实际测试结果

典型半径0.96836 µm，Q/section-area=0.35138 mm/s，典型入口清空时间约6.65 ms。LOW=0.12/清空时间；MEDIUM=0.65/清空时间；压力试跑=8/清空时间，后续容量分支=16/清空时间。几何包围盒尺度给出的约0.426 s只是名义 transit estimate。原计划与压力档调整原因见 [test_flux_cases.yaml](../configs/test_flux_cases.yaml)。全部为软件测试，**NOT EXPERIMENTAL CONCENTRATION**。

| Case | TEST_ONLY λ (bubble/s) | 实际时间(s) | accepted steps | source / rejected / admitted / pending | BB记录 | 停止原因 |
|---|---:|---:|---:|---|---:|---|
| LOW | 18.0465 | 0.235279082 | 2418 | 4 / 1 / 3 / 0 | 0 | WALL_CONSTRAINT_STALL |
| MEDIUM | 97.7517 | 0.0709474101 | 814 | 8 / 2 / 6 / 0 | 122 | WALL_CONSTRAINT_STALL |
| CAPACITY_STRESS_LIMIT | 2406.2 | 0.01044351 | 143 | 30 / 5 / 4 / 21 | 322 | INJECTION_CAPACITY_EXCEEDED |
| FLUX_WEIGHTED | 18.0465 | 0.12 | 1243 | 2 / 0 / 2 / 0 | 0 | MAX_PHYSICAL_TIME |

这张表是 MPI1；MPI4 完全对应。最终验收使用 `*_release_mpi1/4`，见 [CASE_INDEX.json](../contracts/CASE_INDEX.json)。早期 `*_final_*`、不带 release 的 pilot 和 launcher 诊断保留作开发记录，不混入最终统计。

压力初跑1203.10 bubble/s在0.0196432 s达到17个 pending后先触发 WALL_CONSTRAINT_STALL；提高至2406.20 bubble/s的独立压力 case 保持20个 pending门限不变，在0.0104435 s达到21个并触发预期 INJECTION_CAPACITY_EXCEEDED。四颗已准入、21个pending、5个size rejected，正好等于30次 source draws。没有无限循环、重叠或漏计。

LOW/MEDIUM 未跑到配置的0.25/0.15 s。它们的最小 wall gap 接近冻结的1e-10 m safety margin，减小 dt 到下限后显式停止。**停滞粒子仍保留在最后 accepted ACTIVE state**；没有把全局 stop 冒充每颗粒子的 terminal removal。因此 terminal_failure_after_admission=0，守恒式仍成立。不得把这个结果描述为已完成长时间贯通输运。

## 独立验证与证据边界

[validate_workflow.py](../scripts/validate_workflow.py) 从CSV重新计算累计通量、draw/rejection/FIFO与两条数量恒等式，核对每步 native LAMMPS atom count，逐点检查 sampler support、入口三角形和净空；逐步重算 RK2 midpoint/end position，使用独立 VTK 距离查询和连续段距离界检查 wall，并穷举实际小群体 pair 检查 swept gap。它没有只读取 engine PASS。

所有真实 BB stage 及确定性抽选的无相互作用 stage 使用独立 Python resistance assembly/dense/KKT 解核对 PCG；最终用例的最大按半径缩放的速度差 **6.29e-17 m/s**。每一步 RK2 位置恒等式误差为0。BB twist 仍 PENDING；本轮没有改原物理矩阵。[总验证](../validation/INDEPENDENT_WORKFLOW_VALIDATION.json)、[MPI逐文件比较](../validation/MPI_CONSISTENCY.json)。

13项 injector unit tests 覆盖 λdt<1、=1附近与>1、多步累计、变 dt/变 rate、SOURCE/ADMITTED、永久拒绝、FIFO清空后恢复、ID唯一与 trial rollback。6000个 native position samples 分别验证均匀和 VALID_SUPPORT_FLUX_WEIGHTED 接口，独立截面求积的三个边际KS误差均<0.05；这些是固定中位半径的空间诊断，不是6000颗实际入流。见 [POSITION_DISTRIBUTION_COMPARISON.json](../validation/POSITION_DISTRIBUTION_COMPARISON.json)。

四个端口的 outward crossing、reverse/no crossing、WALL_ONLY开放 cap、侧壁段拒绝、tag保留/删除，以及一次显式 owner migration 均在 MPI1/MPI4 native LAMMPS fixture中验证；MPI4覆盖4个port owners。见 [PORT_LIFECYCLE_MPI4.json](../validation/PORT_LIFECYCLE_MPI4.json)。fixture位移是预设的几何测试，**不是实际流场出口**。正常输运目前没有验证到跨rank近接BB pair；这部分仍保留 stable baseline 的既有证据，不夸大本次短轨迹的MPI覆盖。

near-wall记录只涉及BIOLOGICAL_CORE_WALL最近子集、gap/a≤0.2的离散观测时间，记录ENTER/LEAVE与右删失residence；extension/collar不作为生物近壁事件。它们是几何诊断，不表示adhesion。共享region边缘的等距最近点由独立validator核验其距离，不强求两个BVH选中同一个triangle ID。

首次 launcher 没继承 stable baseline 的 HWLOC_COMPONENTS=-gl，worker尚未启动便停在硬件探测；该尝试已终止并保存 LAUNCHER_DIAGNOSTIC。随后所有实际运行继承该环境设置。没有改系统配置、安装依赖或运行/修改 Palabos。

## 科研与功能状态

| 字段 | 值 |
|---|---|
| AUTHORITATIVE_INLET_FLOW_STATUS | PASS |
| AUTHORITATIVE_INLET_FLOW_M3_S | 2.7369132390905703e-15 |
| FROZEN_FIELD_INJECTION_PLANE_COVERAGE | 83.62354517646728% |
| INJECTION_POSITION_POLICY | VALID_SAMPLER_SUPPORT_ONLY |
| INLET_POSITION_DISTRIBUTION_STATUS | WORKFLOW_V0_APPROXIMATION |
| ADAPTIVE_INJECTION_CONTROLLER_STATUS | PASS |
| SIZE_ADMISSIBILITY_STATUS | PASS |
| TRANSIENT_BACKLOG_STATUS | PASS |
| DYNAMIC_PARTICLE_LIFECYCLE_STATUS | PASS |
| FLOW_SAMPLER_RUNTIME_STATUS | PASS |
| RK2_STATUS | PASS |
| BUBBLE_BUBBLE_STATUS | PASS_WITH_TWIST_PENDING |
| HARD_WALL_STATUS | PASS |
| OUTLET_TRACKING_STATUS | PASS |
| NEAR_WALL_EVENT_STATUS | PASS |
| PARTICLE_ACCOUNTING_STATUS | PASS |
| MPI_STATUS | PASS |
| WORKFLOW_TECHNICAL_STATUS | PARTIAL |
| PRODUCTION_BUBBLE_CONCENTRATION_STATUS | UNSPECIFIED |
| FLOW_PHYSICS_STATUS | ENGINEERING_TRANSIENT_FIELD_ONLY |
| WALL_HYDRODYNAMICS_STATUS | NOT_IMPLEMENTED |
| RMBW_USED | NO |
| CF2003_USED | NO |
| FALADE_BRENNER_USED | NO |
| WSS_USED | NO |
| ADHESION_STATUS | OFF |
| REAL_SCIENTIFIC_PRODUCTION_READY | NO |

`DYNAMIC_PARTICLE_LIFECYCLE_STATUS=PASS` 包含实际注入/回滚以及隔离fixture删除/migration；`OUTLET_TRACKING_STATUS=PASS` 只证明端口实现与fixture，没有自然输运出口。整体保持PARTIAL以反映这个差别。浓度未指定、transient冻结场、部分位置支持、壁面物理关闭均未升级。Rallabandi、Higdon–Muldowney、wall lubrication、adhesion、RBC、buoyancy、lift和particle inertia也均未启用。

## 文件与复现

本地：`/home/lzy/projects/compre_output/microbubble_vascular_workflow_v0/bcflux_20260916_224835`

远端结果：`/workspace/microbubble_lammps/results/microbubble_vascular_workflow_v0_bcflux_20260916_224835`

隔离构建：`/workspace/microbubble_lammps/work/microbubble_vascular_workflow_v0_bcflux_20260916_224835`

原 stable rigid CPU binary和六个物理/采样/几何源文件的哈希重核见 [PROJECT_BUILD_PROVENANCE.json](../provenance/PROJECT_BUILD_PROVENANCE.json)。旧 BLOCKED stage、vascular源项目和其Git状态由本轮前后清单验证；最终镜像核验见 [REMOTE_LOCAL_SHA256_CHECK.json](../validation/REMOTE_LOCAL_SHA256_CHECK.json)。未git add/commit/push/checkout/reset/clean。

远端构建后使用 `python3 scripts/run_suite.py 1 release` 与 `python3 scripts/run_suite.py 4 release`；已完成同hash运行会复用，不静默覆盖。单case入口是 `scripts/run_case.py CASE RANKS NEW_LABEL`。本地复核：

```bash
OPENBLAS_NUM_THREADS=1 /usr/bin/python3 -B scripts/validate_workflow.py LOW_release_mpi1 LOW_release_mpi4 MEDIUM_release_mpi1 MEDIUM_release_mpi4 CAPACITY_STRESS_LIMIT_release_mpi1 CAPACITY_STRESS_LIMIT_release_mpi4 FLUX_WEIGHTED_release_mpi1 FLUX_WEIGHTED_release_mpi4
/usr/bin/python3 -B scripts/validate_mpi.py
OPENBLAS_NUM_THREADS=1 /usr/bin/python3 -B scripts/validate_position_sampling.py
OPENBLAS_NUM_THREADS=1 /usr/bin/python3 -B scripts/plot_workflow.py
```

每个run有 flux_timeseries、injection_events、TRAJECTORIES、INTEGRATION_STAGES、SOLVER_HISTORY、PAIR_HYDRODYNAMIC_EVENTS、PORT_EVENTS、NEAR_WALL_EVENTS、RETRY_HISTORY、PENDING_FINAL、case快照及execution hash。另提供0.05 s rolling-window flux统计；离散短窗误差不作单步通量FAIL。

![累计注入](../visualization/target_vs_actual_injection.png)

![Active count](../visualization/active_bubbles_vs_time.png)

![FIFO backlog](../visualization/pending_backlog_vs_time.png)

![真实输运出口计数](../visualization/outlet_counts_vs_time.png)

![源与准入尺寸](../visualization/source_vs_admitted_size_distribution.png)

![实际轨迹](../visualization/trajectories_3d.png)

![注入位置及独立位置诊断](../visualization/injection_positions.png)

![粒子守恒](../visualization/particle_accounting.png)
