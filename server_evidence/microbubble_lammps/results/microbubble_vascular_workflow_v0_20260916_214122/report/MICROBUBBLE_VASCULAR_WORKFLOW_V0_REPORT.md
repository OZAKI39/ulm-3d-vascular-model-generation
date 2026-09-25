# MICROBUBBLE_VASCULAR_WORKFLOW_V0：入口初始化阻塞报告

**这个 workflow 目前到底能做什么？** 当前新 stage 能恢复真实几何来源、分离壁面和四个端口、复现冻结 SonoVue 尺寸，并在入口内移平面成功初始化 W1 的一颗微泡。它尚不能运行完整微泡输运：W2/W3 的固定尺寸微泡无法同时放在这个入口平面内，因此已按用户第 57 节停止动力学开发。

`WORKFLOW_TECHNICAL_STATUS = BLOCKED`

`STOP_INLET_BATCH_INITIALIZATION_INFEASIBLE`

该结论不是求解器失败，也不是运行结果。没有启动本阶段 LAMMPS、RK2、MPI1/MPI4 轨迹实验。预检通过不等于 workflow 通过。

## 停止原因和可审核证据

用户第 18–21 节指定入口向内偏移的 `INJECTION_PLANE`、考虑有限半径、固定 seed 和有界位置重采样。本阶段把它明确实现为所有中心同处一个平面的 `BATCH_AT_T0_SINGLE_PLANE`，冻结配置 seed=42；W1/W2/W3 分别取同一尺寸序列的前 1/8/16 个。用户第 57 节要求：“当前 inlet injection plane 无法稳定采样”时，“使用 STOP_<REASON>”，记录 EXPECTED / ACTUAL / IMPACT / NEXT_OPTIONS，“然后停止等待用户”。

| 项目 | EXPECTED | ACTUAL |
|---|---|---|
| W1 | 1 颗通过 fluid、墙距、流场查询检查 | 206 次位置候选后成功；仅初始化预检 |
| W2 | 8 颗同时满足单平面、无穿墙、无互穿 | 必要几何条件失败，初始化失败 |
| W3 | 16 颗满足相同条件 | 必要几何条件失败，初始化失败 |

| Case | 内移距离 (µm) | 实际截面面积 (µm²) | Σπa² (µm²) | 需求/面积 | 最大半径 (µm) | 可容纳半径上界 (µm) |
|---|---:|---:|---:|---:|---:|---:|
| W1 | 1.458258 | 7.806134 | 4.974610 | 0.637269 | 1.258359 | 1.576666 |
| W2 | 2.059054 | 7.803644 | 40.264844 | 5.159749 | 1.859155 | 1.577746 |
| W3 | 2.059054 | 7.803644 | 68.880294 | 8.826683 | 1.859155 | 1.577746 |

计算用真实生产 STL 的内移切面，选取靠近 inlet 中心的闭合流体轮廓；没有用另一个支管的切面，也没有拿 cap 面积替代内移截面。球心共面时，不互穿球在该平面上产生半径 a 的不相交圆盘；不穿墙又要求圆盘包含于入口流体截面。因此 `Σπa² ≤ A_section` 是必要条件。W2/W3 大幅违反此条件，再多随机位置重试也不能成功。

独立 validator 用三维叉积重算截面面积，并使用更宽松的凸包面积 7.803868 µm² 复核；结论相同。截面两条正交切向宽度约为 3.155761 和 3.155491 µm，给出内切圆半径上界 `min(widths)/2 = 1.577746 µm`。最大球半径 1.859155 µm 超过此上界，**即使 clearance=0 也放不下**。这不是由额外安全余量造成的失败。该上界只是必要条件，不声称是精确最大内切半径。

0.5、1、2、3、4、5 µm 六个内移位置的敏感性检查给出截面面积 7.791714–7.809669 µm²、半径上界 1.574186–1.578272 µm；在这些已检查平面中，改变小幅内移距离仍不能消除两项障碍。不据此声称已验证任意深度的三维初始化方案。

- **IMPACT：** W2/W3 无有效初始坐标；完整 workflow 的 G1–G5 原生几何套件、RK2、安全步进、出口生命周期、事件统计及 MPI 对比均未开展。W1 的初始化结果不能替代其输运验证。
- **NEXT_OPTIONS A：** 用户授权三维入口缓冲区/薄层中的 t=0 批量放置，并明确超径微泡的准入规则。三维排布仅解决共面拥挤问题；它本身不能让超径微泡通过当前狭窄入口。若采用尺寸筛选，需披露筛选后分布及拒绝数。
- **NEXT_OPTIONS B：** 用户指定能容纳目标尺寸与批量规模的更宽真实入口/几何，并配套匹配的背景流；不缩放现有血管冒充原始几何。
- **NEXT_OPTIONS C：** 用户先将下一步范围限定为 W1，或明确批准有披露的几何条件化尺寸群体。仅改尺寸也仍须检查多泡单平面的总容量。

本轮没有缩小微泡、重抽尺寸、挑选有利 seed、改成连续注入或默默采用三维缓冲区。W2/W3 在数学必要条件失败后直接结束，未消耗无意义的 5000 次位置重试预算。结构化证据见 [STOP.json](../validation/STOP.json)、[初始化预检](../validation/INITIALIZATION_PREFLIGHT.json) 和 [独立验证](../validation/PREFLIGHT_VALIDATION.json)。

## 几何来源与入口/出口

来源已独立核对为 PASS，详见 [GEOMETRY_PROVENANCE.json](../GEOMETRY_PROVENANCE.json) 和 [输入清单](../provenance/INPUT_MANIFEST.tsv)。

| 环节 | 实际来源 |
|---|---|
| ROI | raw-analysis__fMOST_0_5_6_0_0_6_0001_02_01__anchor_003274 |
| 保存的 ROI run | sampling/20260825_133201_radius_plus_structure_k5 |
| Ultraliser | ultraliser_anchor003274_20260825_133350 |
| CFD preprocess | global_to_roi_anchor003274_20260825_183628 |
| CFD surface prepare | vmtk_tps_boundarynormal_crossseam_finalized_recovery_anchor003274_20260826_221611 |
| CFD 坐标系 | axis_aligned_inlet_geometry_anchor003274_20260829_111451 |

生产闭合 STL SHA256：`840da5e1c43ec31ac70ec781cbf75b32940bc73538b5dba83c58af1ddf5e36fb`。

带标签的 axis-aligned VTP SHA256：`83801f326d3266fee6ea1d6d03bcd85fe54191577b0a7421b0cac338cbfd59c8`。

刚体变换后 VTP 顶点误差最大 8.04e-14 µm；生产 float32 米制 STL 与高精度 VTP 的最大量化差为 1.20537e-5 µm。45,512 个 biological core 面与原 Ultraliser 面在统一坐标表示下精确对应；标签、拓扑、保存的 ROI 与 SWC/manifest 链条均核对。旧脚本名在当前工作区缺失的部分以只读 `git show HEAD` 恢复到新 stage 的 provenance，并未恢复或修改原工作区。

实际读取到的数组包括 `CellEntityIds`、`boundary_index`、`boundary_origin_code`、`boundary_origin`、`port_id`、`SurfaceRegionId`、`SurfaceRegion`、`RemeshEntityId`、`boundary_type_code`。区域选择依据实际标签及来源代码，未猜测数组名。

| 派生表面 | source boundary index | CellEntityIds | 三角形数 | cap 面积 (µm²) |
|---|---:|---:|---:|---:|
| INLET | 0 | 4 | 56 | 7.819753 |
| OUTLET_0 | 1 | 3 | 44 | 4.416141 |
| OUTLET_1 | 2 | 5 | 42 | 4.315071 |
| OUTLET_2 | 3 | 2 | 49 | 5.110621 |

OUTLET_0/1/2 是本阶段明确映射的 ID，不等同于源 boundary index。各端口中心、外法向、文件路径和面积见 [BOUNDARY_MANIFEST.json](../geometry/BOUNDARY_MANIFEST.json)。

`WALL_ONLY` 有 67,071 个侧壁三角形：45,512 个 `BIOLOGICAL_CORE_WALL`、21,349 个 `ARTIFICIAL_EXTENSION_WALL`、210 个 `ARTIFICIAL_BOUNDARY_SURGERY_WALL`。这三类均保留为未来几何壁约束；只有第一类可进入 biological near-wall 统计。人工 collar 单独标记，UNKNOWN=0。191 个 inlet/outlet cap 全部排除于 WALL_ONLY，壁与端口分区互斥且覆盖全部 67,262 个生产面。

所有派生 STL 都是原始 50 字节三角形记录的精确子集，单位米。VTP 使用同一生产坐标，仅合并完全相同坐标的索引以建立邻接；已独立读回验证。原始 vascular geometry 未修改。

## 微泡生成、入口放置与水流读取

沿用冻结 `sonovue_sampler.py` 与 `FROZEN_SONOVUE_HISTOGRAM.csv`，其输入量是 **diameter_um**，通过 `radius_m = diameter_um × 0.5e-6` 转换。采用数量加权直方图和 bin 内均匀逆 CDF，保留 uniform 随机数、seed、diameter 与 radius。原冻结尺寸支持对应半径 0.375–2.625 µm。每个 case 的 `POPULATION_DRAWN.csv` 是完整尺寸抽样记录；W2/W3 不伪造初始位置文件。

[workflow_v0.yaml](../configs/workflow_v0.yaml) 冻结单平面 t=0 batch、seed=42、5000 次/泡位置候选预算。偏移规则为 `max(4*dx, a_max+dx)`，初始化 clearance 为 `dx=1.9989918081065344e-7 m`，取自真实流场网格。位置候选按 inlet 三角形面积均匀采样后内移，再检查闭合参考几何内部、三维 wall-only 距离减半径、流场八角点有效性、泡间无重叠。只使用纯几何余量，不引入壁力。尺寸与位置使用分离的确定性随机序列。

W1 保存于 [BUBBLES_INITIAL.csv](../cases/W1/BUBBLES_INITIAL.csv)：球心约 `(103.729603, 57.542546, 157.348141) µm`，半径 1.258359 µm，最小壁面 gap=0.258527 µm，大于 0.199899 µm 的余量。预检背景速度约 `(-1.673376e-6, 1.712158e-7, -6.500126e-4) m/s`。重复生成的尺寸及已接受初始坐标文件逐字节一致。

流场文件 SHA256：`7828ed32b0cdd1e0e85836a564b6737128e7dbe30181ddd43110233b37b64a7f`；几何 hash 与 field contract 匹配，单位为 SI。它仍是 PBS/BSA 5000-step smoke field：

`FLOW_PHYSICS_STATUS = ENGINEERING_TRANSIENT_FIELD_ONLY`

W1 查询由独立 Python 预检按原八角点 trilinear 规则检查，八角点都必须存在且为 fluid；不是本阶段原生 C++ `FlowFieldSampler` 的运行验证，也未验证沿轨迹持续可查询。原生 `frozen_flow.cpp/.hpp` 从 stable baseline 原样复制到 src，尚未编译执行。冻结场保留历史 transient backflow 限制；不能升级为 validated steady vascular flow。

## 动力学与事件功能的实际完成度

| 用户要求 | 本轮实际状态及后续接入位置 |
|---|---|
| LAMMPS、MPI ownership、neighbors/ghosts | 已只读审计 stable rigid 主干并保留源码副本；未建立新运行入口 |
| bubble-bubble | 历史 baseline 状态保留 `PASS_WITH_TWIST_PENDING`；拟继承 bulk/excess resistance、PCG；本 workflow 未运行、未重验 twist |
| RK2 | 已审计原 midpoint、两次流场采样/solve 及现有自适应约束；未移植/执行新端口感知循环 |
| 墙面防穿透 | WALL_ONLY 已正确派生；本轮只验证初始化球体距离，没有实现/验证动态 swept-wall、dt halving、stall |
| 出口/入口反流 | caps 与 wall 已分离，ID/外法向已冻结；crossing 与粒子移除生命周期未实现/运行 |
| near-wall events | 区域标签与 1.0/0.5/0.2/0.1 阈值已入配置；enter/exit/residence 与独立重算未实现/运行 |
| 轨迹、粒子最终 summary | 无轨迹实验；未输出虚假的 events、终态、停留时间或出口计数 |
| MPI1 与 MPI4 | NOT_RUN，不能声称一致或无 silent particle loss |

拟实现的完整版本描述仍为：**Frozen-flow vascular microbubble transport with bubble-bubble hydrodynamics and geometric hard-wall exclusion.** 这描述目标功能，不表示当前 stopped stage 已具备完整功能。

stable rigid baseline 的原始源码、CMake 与 build receipt 保存在 `provenance/stable_rigid_source/`。没有以 wall_hydrodynamics 主程序作为主干。wall stage 的 `wall_distance.cpp/.hpp` 仅存于 `provenance/wall_geometry_donor/` 供审计；其中发现 PCA/local-plane validity 内容，故未整文件复制到 runtime src。停止前没有完成纯几何 C++ 提取。没有复制或加载 RMBW HDF5 table、wall_model resistance、wall shear RHS。

预期动力学仍为 `R_total=R_bulk+R_bubble_bubble_excess`，但本轮没有执行此系统。RMBW、CF2003、Falade–Brenner、Rallabandi、Higdon–Muldowney、WSS、BIE/QBX、wall lubrication、adhesion、RBC、buoyancy、lift、particle inertia 均未启用。没有壁面 penalty/spring/repulsion。`WALL_HYDRODYNAMICS_STATUS=NOT_IMPLEMENTED`；即使未来几何事件可统计，也不能据此解释最终 near-wall residence、sliding/rolling velocity 或 adhesion probability。

## 状态、验证与文件

| 状态字段 | 值 | 范围 |
|---|---|---|
| GEOMETRY_PROVENANCE_STATUS | PASS | 来源、坐标、hash 和 lineage |
| PORT_CLASSIFICATION_STATUS | PASS | 侧壁与 1 inlet/3 outlet 资产分离 |
| SONOVUE_INITIALIZATION_STATUS | FAIL | W1 预检 PASS；W2/W3 INITIALIZATION_FAILED |
| FLOW_SAMPLER_STATUS | NOT_RUN | 原生 workflow sampler；仅 W1 Python 查询预检通过 |
| BUBBLE_BUBBLE_STATUS | PASS_WITH_TWIST_PENDING | 仅保留已有 baseline 状态 |
| WORKFLOW_BUBBLE_BUBBLE_VALIDATION_STATUS | NOT_RUN | 本阶段未开展 |
| HARD_WALL_STATUS | NOT_RUN | 动态路径约束未验证 |
| OUTLET_TRACKING_STATUS | NOT_RUN | crossing/lifecycle 未验证 |
| NEAR_WALL_EVENT_STATUS | NOT_RUN | 时序事件未验证 |
| RK2_STATUS | NOT_RUN | 新 workflow 未执行 |
| MPI_STATUS | NOT_RUN | 无 MPI1/4 对比 |
| WORKFLOW_TECHNICAL_STATUS | BLOCKED | 初始化 STOP |
| FLOW_PHYSICS_STATUS | ENGINEERING_TRANSIENT_FIELD_ONLY | PBS/BSA smoke field |
| WALL_HYDRODYNAMICS_STATUS | NOT_IMPLEMENTED | 本轮无 wall physics |
| ADHESION_STATUS | OFF | 无黏附模型 |
| REAL_SCIENTIFIC_PRODUCTION_READY | NO | 不能作科研 production 结论 |
| LARGE_N_READY | NO | 无大规模验证 |
| GPU_PERFORMANCE_READY | NO | 无 GPU 开发/性能宣传 |

用户状态模板中的部分 PASS/FAIL 字段在强制 STOP 后并未执行，特明确使用 **NOT_RUN**，避免伪造通过或把未执行误写为数值失败。终端使用 `MICROBUBBLE_VASCULAR_WORKFLOW_V0_STOPPED`，不输出会误导完整性的 COMPLETE。机器可读状态见 [WORKFLOW_STATUS.json](../WORKFLOW_STATUS.json)。

[validate_workflow.py](../scripts/validate_workflow.py) 当前是**独立预检 validator**：核对完整互斥边界分区、STL 记录精确对应、切面面积及凸包必要条件、直径/半径与 seed、W1 球体 gap、重复初始化一致性、没有伪造轨迹、输入保护及 Git 状态。它尚未实现用户第 43 节要求的完整轨迹 validator；原生 G1–G5 与所有 runtime 验证在 JSON 中明确列为 NOT_RUN。

预检检查结果：PASS；40 个记录的本地原始输入 hash 未改变；15 个历史远端/本地输入 hash 匹配；两个 Git 工作区的 branch、commit 与完整 porcelain 状态前后相同。保留工作区原有 dirty/untracked 状态，没有 git add/commit/push 或任何历史 stage 修改。

Vast hostname 重新确认为 `f7c62a262077`，microbubble 与 HemoCell roots 均存在，无 STATE_DRIFT_DETECTED。当前 CPU 可执行文件及 rigid_math 库只读重新 hash，与历史 build receipt 一致；未执行。LAMMPS 记录为 `22Jul2025 Update 6`、commit `9c5ab448c78a14fd534619622162ba418d6a1fb1`，编译器/MPI/本地依赖版本见 [BUILD_CURRENT_IDENTITY.json](../provenance/BUILD_CURRENT_IDENTITY.json)。没有新 build objects、安装或 GPU 任务。

本地 stage：`/home/lzy/projects/compre_output/microbubble_vascular_workflow_v0/20260916_214122`

远端 stage：`/workspace/microbubble_lammps/results/microbubble_vascular_workflow_v0_20260916_214122`

远端 work：`/workspace/microbubble_lammps/work/microbubble_vascular_workflow_v0_20260916_214122`

同步及原始远端输入保护的最终证据见 [REMOTE_LOCAL_SHA256_CHECK.json](../validation/REMOTE_LOCAL_SHA256_CHECK.json)。该回执检查清单内文件逐个 SHA256、文件集合，并单独校验同步清单/回执，避免自引用 hash。未同步巨大 build objects 或 cache。

复现已完成预检（在 stage 目录内，使用已记录的本地依赖环境）：

```bash
OPENBLAS_NUM_THREADS=1 /usr/bin/python3 -B scripts/initialization_preflight.py
OPENBLAS_NUM_THREADS=1 /usr/bin/python3 -B scripts/validate_workflow.py
OPENBLAS_NUM_THREADS=1 /usr/bin/python3 -B scripts/plot_preflight.py
```

这些命令仅重建/检查预检，不运行输运。`prepare_inputs.py` 读取 `configs/source_paths.yaml` 的只读上游路径，需原始 WSL 上游项目；`probe_remote.py` 是首次创建 stage 的历史预检，不应在同一目录直接重跑。`workflow_v0.yaml` 的 simulation 数值是尚未执行和资格验证的计划参数，当前没有可消费它们的完整 C++ engine。

仅生成三张有真实数据支持的图。未生成要求中的后八张轨迹/事件/MPI 图及 WORKFLOW_TRAJECTORIES.vtp，因为不存在相应运行数据。

![真实壁面、入口、三个出口及 W1 内移平面](../visualization/workflow_geometry.png)

![W1 实际初始球截面；不是轨迹](../visualization/initial_bubbles.png)

![W2/W3 共面容量及超径障碍](../visualization/injection_feasibility.png)

当前等待用户决定上述 NEXT_OPTIONS，之后才能改变初始化契约或缩小任务范围。不会自动进入近壁水动力、黏附、RBC、GPU 或 large-N 开发。
