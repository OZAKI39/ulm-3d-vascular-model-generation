# ADAPTIVE_FLUX_INJECTION_AND_OPEN_BOUNDARY_V0：入口流量覆盖阻塞报告

本轮已按新要求废止 `BATCH_AT_T0_SINGLE_PLANE`，不再尝试把 8/16 颗微泡同时铺在入口。**当前仍未形成可运行的自动注入 workflow**：新的入口流量审计发现，现有冻结采样器在真实截面约 16.4% 的面积上没有定义可用速度，不能可靠给出整个入口的正向流量。因此按本次请求第 73 节停止。

`WORKFLOW_TECHNICAL_STATUS = BLOCKED`

`STOP_INLET_FLUX_INCOMPLETE_SAMPLER_COVERAGE`

这是本轮新发现的入口流量覆盖问题，和上一轮单平面 population 容量不足是不同的阻塞。上一 BLOCKED stage 保持原样。

## 六个问题的直接回答

1. **微泡现在怎样自动进入血管？** 当前尚未实现自动注入。新契约要求在 timestep boundary 按通量需求逐候选准入；真实入口流量的前置检查阻塞，尚未接入 LAMMPS。
2. **程序怎样决定此刻进入几颗？** 目标是累计 `∫λ(t)dt` 后处理整数需求，不采用固定间隔、one-per-step 或固定 N_active。累计控制器本轮尚未实现/测试。
3. **入口暂时堵住怎么办？** 目标是保持候选尺寸与 source draw 的 FIFO pending queue，并检查 backlog 数量、年龄和通量欠额。本轮尚未实现，不声称拥堵恢复已通过。
4. **尺寸本身进不去怎么办？** 新要求中的 SOURCE_POPULATION 与 ADMITTED_POPULATION 语义必须分开测试；前者记录源群体拒绝，后者拒绝后继续抽样以满足准入通量。原 SonoVue 分布不改，production control basis 未指定；本轮尚未开展尺寸准入/两种模式验证。
5. **怎样保证数量不凭空增加/消失？** 目标是 source、size rejection、pending、admitted、active、exit 和 admission 后终态的逐步恒等式审计。本轮没有动态粒子运行，未验证守恒。
6. **怎样查看实时入口/出口 flux 和粒子数？** 本轮只有入口积分审计与覆盖图，没有运行产生的 flux_timeseries、injection_events、轨迹或粒子计数。报告不输出虚构时间序列。

## 继承与来源保护

[BASELINE_CHECK.json](../provenance/BASELINE_CHECK.json) 已核验上一 stage 清单内 93 个文件，以及其 40 个记录的原始输入。另将上一 stage 全部 96 个文件冻结为本轮前后保护清单。远端对同样的 96 个旧 stage 文件和 15 个原始输入再次核对，结果 PASS；hostname 仍为 `f7c62a262077`。后续独立 validator 核实本地旧 stage、原始输入及两个 Git 工作区的 branch/commit/完整 status 均未变化。

直接继承 `GEOMETRY_PROVENANCE_STATUS=PASS`、`PORT_CLASSIFICATION_STATUS=PASS`，没有重做 geometry provenance。WALL_ONLY 的 67,071 个侧壁面与 191 个 cap 保持互斥，1 inlet、3 outlets 的既有 ID 映射不变；biological core、extension 和 boundary surgery 三类仍单独标记。没有修改 vascular source。

| 输入 | SHA256 |
|---|---|
| 原生产 closed STL | `840da5e1c43ec31ac70ec781cbf75b32940bc73538b5dba83c58af1ddf5e36fb` |
| 带标签生产 VTP | `83801f326d3266fee6ea1d6d03bcd85fe54191577b0a7421b0cac338cbfd59c8` |
| Frozen Flow | `7828ed32b0cdd1e0e85836a564b6737128e7dbe30181ddd43110233b37b64a7f` |
| SonoVue sampler | `c3f5ab75e20503a495cd107cef9b63240bda890e90f601d89b1f28fa8f292ece` |
| SonoVue histogram | `2c9f783c6169421b06c057e16653878e7385139a179682c0648519e72a9f5198` |

完整清单见 [INPUT_MANIFEST.tsv](../provenance/INPUT_MANIFEST.tsv)。stable rigid 仍是未来动力学主干：原始 source/build receipt 保存在 provenance，`frozen_flow` 和 `rigid_math` 的原样副本在 src。只有原生 flow query probe 被编译执行，没有运行旧 wall stage 主程序、旧 RK2 或 LAMMPS 动力学。

浓度检索覆盖 SonoVue 项目、handoff review、stable baseline 和 vascular configs/docs。[搜索证据](../provenance/CONCENTRATION_SEARCH.txt) 中，SonoVue README 明确写明正式数量与浓度尚未定义；其它命中为测试粒子群、PBS/BSA 介质浓度或上游软件示例，不能当作本项目实验 bubble number concentration。

`PRODUCTION_BUBBLE_CONCENTRATION_STATUS = UNSPECIFIED`

`CONTROL_BASIS_STATUS = UNSPECIFIED`

`TEST FLUX IS NOT EXPERIMENTAL CONCENTRATION.` 本轮因入口流量 STOP，尚未选取或运行 TEST_ONLY flux。

## 入口流量审计

[配置](../configs/inlet_flux_audit.yaml) 的所有长度/流量使用 SI。固定内移规则为 `max(4*dx, source_radius_support_upper+dx)`，其中 dx=0.199899 µm，冻结 SonoVue 支持上界 radius=2.625 µm，因此内移距离为 **2.824899 µm**。这是一次性选择固定开放入口平面，不随每次 draw 或 active 数变化；未来可准入的有限球完全位于 cap 内侧。该平面实际面积 **7.788988 µm²**，由原始闭合几何切割并三角化，未用 cap area × mean speed。

法向采用 `n_in=-n_out`。在每个查询有效的面积积分点计算：

- `Q_positive = ∫ max(u·n_in,0) dA`
- `Q_backflow = ∫ max(-u·n_in,0) dA`
- `Q_net = ∫ u·n_in dA`

四个三角形细化等级均使用三点面积求积。下表的流量**仅覆盖采样器 VALID 区域**，不是完整入口流量。

| 细化级别 | 三角形数 | 积分点数 | 有效面积占比 | 局部 Q_positive (m³/s) |
|---|---:|---:|---:|---:|
| 2 | 1,072 | 3,216 | 82.2423% | 2.4509541714e-15 |
| 3 | 4,288 | 12,864 | 83.7136% | 2.4740722597e-15 |
| 4 | 17,152 | 51,456 | 83.5125% | 2.4708675266e-15 |
| 5 | 68,608 | 205,824 | 83.5817% | 2.4716110964e-15 |

最高两级局部 Q_positive 相对差 **0.0301%**，满足 1% 的求积细化建议。但这个结果仅说明**有定义区域上的积分收敛**，不能弥补其余面积上的未知速度。最高级 VALID 面积约 6.510166 µm²，SOLID-query 面积约 1.278822 µm²；没有 MISSING、OUTSIDE 或 NONFINITE_POSITION 查询。

局部 Q_backflow=0，局部 Q_net=Q_positive，恒等式残差为 0。这只涉及当前有效区域，不证明未知环带无反流，也不否认该 transient field 的历史 backflow。完整的 Q_positive、Q_net、Q_backflow 和 backflow fraction 在 [INLET_FLUX_AUDIT.json](../validation/INLET_FLUX_AUDIT.json) 中均为 **null/未确定**，没有用局部值冒充。

## 为什么 SOLID 不能在这里解释成零速度

原始 C++ `FlowFieldSampler` 要求八个插值角点全部是 fluid 且均有 velocity；只要一个角点不是 fluid，就返回 SOLID 状态和数值占位 `{0,0,0}`。它不保证该查询点本身位于连续 STL 的固体区域。

[既有采样契约证据](../provenance/INHERITED_SAMPLER_POLICY.json) 明确规定 `HARD_FAIL; no clamping, nearest fallback or extrapolation`。此次没有把无效返回的占位零当作物理速度，没有删除未知面积后重归一化，也没有外推或自动附加壁面零速重构。这些都将改变已冻结的速度采样含义。

独立核验采用不同于细化面积求积的方法：把原截面三角形逐一裁剪到 Cartesian 插值 cell，按该 cell 的八角点支持分类，得到有效面积占比 **83.623545%**、未定义面积占比 **16.376455%**。裁剪面积和等于实际截面积；它与最细求积有效占比仅差 0.04188 个百分点。

另外检查了 128 个返回 SOLID 的积分位置，它们全部处于闭合 STL 的内部，到表面的距离约 0.000090–0.191670 µm。这证明缺口涉及真实流体内部的近边界带，而非可以直接丢弃的外部固体。

在 Vast 新 work 目录编译原样的 `frozen_flow.cpp`，由新增只读 CLI 查询 2,147 个实际积分位置。原生 C++ 与 Python 的状态逐点一致，速度绝对差最大 **0 m/s**。这一测试验证采样实现的一致性，没有声称它覆盖完整截面。见 [NATIVE_SAMPLER_EQUIVALENCE.json](../validation/NATIVE_SAMPLER_EQUIVALENCE.json)、[编译记录](../provenance/NATIVE_PROBE_BUILD.json)。

独立 [validate_workflow.py](../scripts/validate_workflow.py) 还以常数、x、x²、xy 的解析面积积分验证求积公式及细化，重算 raw 数据中的局部流量、输入 hash 和 Git 状态；[预检验证](../validation/PREFLIGHT_VALIDATION.json) 为 PASS。它当前只验证本次阻塞证据，没有升级为已完成的动态 workflow validator。

## STOP：EXPECTED / ACTUAL / IMPACT / NEXT_OPTIONS

- **EXPECTED：** 在完整真实 injection section 上可靠得到正向、净向和反向体积流量，并通过最高两级 <=1% 收敛门。
- **ACTUAL：** 局部积分满足收敛要求，但约 16.4% 的截面上速度查询无定义；原生与独立实现都证实这个缺口。无法给出完整入口正向通量。
- **IMPACT：** 触发本次请求第 73 节“无法可靠得到 inlet positive flux”时停止的要求。没有构建使用局部 Q 冒充全截面 Q 的 `FrozenFieldInletFluxProvider`，没有开发/验证 controller、source/admitted semantics、pending queue、动态 insertion/removal、RK2、出口或 MPI workflow。
- **NEXT_OPTIONS A：** 提供或授权从 Palabos 导出完整物理入口截面的边界感知速度/面积积分数据，明确包含 positive/net/backflow 分量。
- **NEXT_OPTIONS B：** 授权单独开发并验证用于 flux audit 的壁面边界速度重构，明确插值契约、覆盖范围与误差门；不能默默采用 nearest、clamp 或补零。
- **NEXT_OPTIONS C：** 明确把下一轮缩小为已披露有效内部支持区域上的直接 TEST_ONLY number flux 软件测试，同时保留完整入口 Q 未确定。此方案需要用户修改当前完整入口流量的前置条件，不属于本轮 flux audit PASS，也不定义真实实验浓度。

见机器可读 [STOP.json](../validation/STOP.json)。本轮没有通过改几何、改尺寸、假速度或忽略无效返回绕过条件，也不继续修复已废止的旧 batch 初始化方案。

## 功能与科研状态

| 字段 | 状态 / 范围 |
|---|---|
| GEOMETRY_PROVENANCE_STATUS | PASS，继承并重核 hash |
| PORT_CLASSIFICATION_STATUS | PASS，继承既有 mapping |
| INLET_FLUX_INTEGRATION_STATUS | BLOCKED；局部求积收敛，完整截面覆盖不成立 |
| ADAPTIVE_INJECTION_CONTROLLER_STATUS | BLOCKED；前置 STOP，尚未实现 |
| SIZE_ADMISSIBILITY_STATUS | NOT_RUN |
| TRANSIENT_BACKLOG_STATUS | NOT_RUN |
| DYNAMIC_PARTICLE_LIFECYCLE_STATUS | NOT_RUN |
| FLOW_SAMPLER_STATUS | PASS，仅原生/独立状态与速度一致性测试 |
| RK2_STATUS | NOT_RUN |
| BUBBLE_BUBBLE_STATUS | PASS_WITH_TWIST_PENDING，仅保留历史 baseline 状态 |
| WORKFLOW_BUBBLE_BUBBLE_VALIDATION_STATUS | NOT_RUN |
| HARD_WALL_STATUS | NOT_RUN |
| OUTLET_TRACKING_STATUS | NOT_RUN |
| NEAR_WALL_EVENT_STATUS | NOT_RUN |
| PARTICLE_ACCOUNTING_STATUS | NOT_RUN |
| MPI_STATUS | NOT_RUN |
| WORKFLOW_TECHNICAL_STATUS | BLOCKED |
| CONTROL_BASIS | UNSPECIFIED |
| PRODUCTION_BUBBLE_CONCENTRATION_STATUS | UNSPECIFIED |
| FLOW_PHYSICS_STATUS | ENGINEERING_TRANSIENT_FIELD_ONLY |
| WALL_HYDRODYNAMICS_STATUS | NOT_IMPLEMENTED |
| ADHESION_STATUS | OFF |
| REAL_SCIENTIFIC_PRODUCTION_READY | NO |
| LARGE_N_READY / GPU_PERFORMANCE_READY | NO / NO |

Validation C 已完成审计并判为 BLOCKED；A/B/D–L、flux-weighted histogram 和 near-wall 时序重算均 NOT_RUN。使用 NOT_RUN 如实区分“未执行”与数值 FAIL，不伪造模板中尚未验证的 PASS。`ADAPTIVE_FLUX_INJECTION=NOT_IMPLEMENTED`，目标模式虽已授权，但当前不能写 YES。终端打印 STOPPED，不宣称完整 workflow COMPLETE。

近壁水动力保持关闭：RMBW、CF2003、Falade–Brenner、Rallabandi、Higdon–Muldowney、WSS、wall lubrication/penalty、adhesion、RBC、buoyancy、lift、particle inertia 均未启用；twist 未修改。任何未来正常轨迹的目标仍为 bulk + bubble-bubble excess + geometric hard wall。当前没有此类运行结果。

## 产物、同步与复现

本地：`/home/lzy/projects/compre_output/microbubble_vascular_workflow_v0/adaptive_flux_20260916_222104`

远端 results：`/workspace/microbubble_lammps/results/microbubble_vascular_workflow_v0_adaptive_flux_20260916_222104`

远端 work：`/workspace/microbubble_lammps/work/microbubble_vascular_workflow_v0_adaptive_flux_20260916_222104`

[WORKFLOW_STATUS.json](../WORKFLOW_STATUS.json) 保存完整分项状态；[REMOTE_LOCAL_SHA256_CHECK.json](../validation/REMOTE_LOCAL_SHA256_CHECK.json) 保存最终同步核验。只读 native probe 的二进制保留在远端新 work 目录，没有同步大型 build cache、没有安装包、没有 GPU 或 large-N 任务。

已实现的本地审计入口：

```bash
OPENBLAS_NUM_THREADS=1 /usr/bin/python3 -B scripts/audit_inlet_flux.py
OPENBLAS_NUM_THREADS=1 /usr/bin/python3 -B scripts/validate_workflow.py
OPENBLAS_NUM_THREADS=1 /usr/bin/python3 -B scripts/plot_audit.py
```

validator 读取保存的原生查询回执；需要重建远端只读 probe 时可单独运行 `scripts/verify_native_sampler.py`。`probe_environment.py` 是首次创建新远端 stage 的历史操作，不应在同名目录重跑。没有调用上一 stage 的写入型脚本。

本轮图像是实际流速支持诊断、局部积分收敛，以及标明“未实现”的目标流程图；没有 injection position density、累计注入/出口、active count、backlog、尺寸准入分布或轨迹图，因为不存在相应数据。

![真实入口正向速度及未定义区域](../visualization/inlet_positive_flux_map.png)

![局部积分收敛与速度覆盖](../visualization/inlet_flux_convergence.png)

![目标流程与当前停止位置；后续模块尚未实现](../visualization/adaptive_injection_workflow.png)

上一 BLOCKED stage、所有历史物理 stage 和 vascular source 保持原样；未 git add、commit、push、checkout、reset 或 clean。当前按第 73 节等待用户决定后续数据/插值契约，不自动进入下一开发阶段。
