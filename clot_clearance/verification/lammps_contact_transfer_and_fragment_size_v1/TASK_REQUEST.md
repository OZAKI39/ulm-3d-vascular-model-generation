# CODEX TASK — LAMMPS_CONTACT_TRANSFER_AND_FRAGMENT_SIZE_V1

你已经知道如何进入服务器和运行项目。**本轮必须实际执行代码改动、少量真实 LAMMPS 求解、碎块统计和可视化，而不是仅写计划。** 结束后停止，等待人工审查。

## A. 任务所处阶段与最终目标

项目现在的最终链路是：

`microbubble cluster morphology → streaming/shear force field → LAMMPS clot-like fragmentation → fragment outcome → surrogate → inverse PZT control`。

高保真 B1/B2 phase-field 路线已由项目政策终止，完整 PD2 fracture 开发停止。现役碎裂器为 **LAMMPS BPM/PERI LPS**；只要求趋势、机制和可见碎块在 `EXPLORATORY_ONLY` 范围内可用，不做患者定量预测。

前一轮 `LAMMPS_CONTACT_CASCADE_DIAGNOSIS_AND_MINIMAL_FIX_V1` 已显示：在同一半球、Liu 材料、100 Pa 局部 Gaussian、2T_cross 条件下，原生 1.35dx 接触使 bulk 断键达约 99.95%、debris 达约 99.60%；关掉接触或改为 **1.0dx** 后，bulk 断键约 5.3%，debris 约 3.5%–3.8%，最大主体保留约 95.6%。该结论仅限已测工况。

**本轮不再问“LAMMPS 能否碎裂”，而是回答三个实际问题：**

1. `contact_distance_factor=1.0` 在**另一种 streaming pattern**上是否仍能避免明显的接触诱发粉碎？
2. 同一个 Gaussian pattern 的载荷从 100 Pa 加到 **200 Pa**，是否能产生有实际体积的游离碎块，或仍会全局粉碎？
3. 现有结果到底出现了多大的**真正游离碎块**、出现多久、是否持续存在？旧 force–duration map 哪些边界需要最少量更新？

`ROUND_NAME = LAMMPS_CONTACT_TRANSFER_AND_FRAGMENT_SIZE_V1`

`PRIMARY_USE = QUALITATIVE_FRAGMENTATION_SCREENING`

`BIOLOGICAL_QUANTITATIVE_VALIDATION = NOT_ESTABLISHED`

`PATIENT_SPECIFIC = NOT_ESTABLISHED`

## B. 权威输入、保护历史与源码

先读取以下旧任务产物，**以实际工作区已有文件为准**；若字段与此提示词略不同，在 `STARTING_STATE.md` 记录真实值，不悄悄重定材料：

- `REVIEW_LAMMPS_COHERENT_FAILURE_FORCE_TIME_MAP_V1.md`
- `REVIEW_LAMMPS_CONTACT_CASCADE_DIAGNOSIS_AND_MINIMAL_FIX_V1.md`
- `FINAL_STATUS.json`、contact patch/构建记录、材料卡、原生 broken bond events、保存帧与脚本
- CONTACT_ON 原生基线、CONTACT_OFF 诊断与 CONTACT_MODIFIED (1.0dx) 的原始轨迹

保护所有尚未提交的历史成果。开始记录分支、HEAD、`git status --short`、工作区 `diff`、可执行文件和 patch 的 SHA、关键目录哈希。**不得 `git reset --hard`、`git clean`、`git stash`、覆盖未提交文件、删除历史数据、改写上一轮 37-case map。不得不经授权 commit/push/PR。**

优先复用上一轮独立变体：

- 来源 `tools/lammps-contact-onset-v1/`
- 构建 `tools/lammps-contact-onset-v1-build/lmp`

这些是上轮报告中的相对路径；使用前确认实际存在、`contact_distance_factor` 可用、源码身份和版本匹配，不要盲用或覆盖官方原装安装。

LAMMPS 官方版本：`stable_30Sep2026`，原报告 SHA `8de817dd79bfe4525d5d39246a212d833e6dee07`。

官方参考：

- https://github.com/lammps/lammps/tree/stable_30Sep2026
- https://docs.lammps.org/bond_bpm_peri.html
- https://docs.lammps.org/pair_bpm_peri.html
- https://github.com/lammps/lammps/blob/stable_30Sep2026/src/BPM/bond_bpm_peri.cpp
- https://github.com/lammps/lammps/blob/stable_30Sep2026/src/BPM/pair_bpm_peri.cpp

本轮**不改**官方 LAMMPS 本构、不重写 contact、不追加 progressive softening、VES、PeriLab 或 PD2。只复用已验证的 `1.0dx` 工程候选。

## C. 全部力学配置冻结

除明确命名的 `pattern_id` / `tau_budget` 之外，所有科学输入保持与上一轮 `CONTACT_MODIFIED` Gaussian 100 Pa / 2T_cross **字节可比**：

- 3D hemisphere：`R = 0.5625 mm`
- `dx = 50 µm`；`horizon = 3.01 dx`
- 2971 clot nodes，三层 1203 fixed wall anchors（如实际文件相同）
- 现有 bulk/interface/wall bond topology，固定 3 层锚点，其余 clot 全部可动
- 来源：`Liu 2024 whole-blood bulk + porcine-skin interface proxy`；非真实人体血管壁
- `K = 33333.333... Pa`，`G = 940 Pa`，`rho = 1000 kg/m³`（密度为工程假定）
- `s00_bulk = 0.15`、`s00_interface = 0.140624859375`、`alpha = 0`
- 数值阻尼 `500 s^-1`
- `contact_distance_factor = 1.0`，其余 contact stiffness / law 不变
- `dt = 0.5 µs`，`T_cross = 0.580174451 ms`；脉宽 `2 T_cross`
- 同一 smooth pulse、同一后续观察窗口，物理终止时间约 `2.33 ms`，对应约 `4660` steps；如实际输入有精确取整，以同一冻结脚本为准
- 同一 CURRENT_POSITION_STREAMING_PROVIDER_V2，力随当前坐标离场衰减；沿用原始 exposed-node area weights
- 同样的时间采样/原生 bond events，优先保存足够的逐帧信息用于碎块生命周期
- MPI ranks 先沿用已有每 case 2 ranks，不做并行规模测试

所有 `tau_budget` 指**初始参考构型等总绝对切向力除以初始表面权重面积之和**。200 Pa 是**工程探索预算**，不是已测微泡载荷或临床破裂阈值。不同 pattern 的峰值局部 traction 未必相同，必须单独报告。

不得改 `s00`、`K/G`、密度、阻尼、界面校准、锚层、接触因子、几何、网格、horizon、仿真结束时刻以追求某种动画结果。

## D. Stage 0：只读重分析历史轨迹（不增加力学计算）

重点复用上一轮 `CONTACT_MODIFIED` Gaussian 100 Pa / 2T_cross 完整原生轨迹。

已有报告指出：0.500 ms 时曾出现 **49-node detached component**，随后继续破碎。把这条记录作为 **碎块检测器的交叉检查**，但不能只相信报告，需从 surviving-bond graph 与 interface-bond history 独立读出。

对于每一帧：

1. 用原生断键事件重建 surviving **bulk** bond graph；wall anchors 不计入 clot component。
2. 对每个 clot-only connected component，检查其成员是否仍通过 surviving interface bonds 连到任一 wall anchor。
3. `ATTACHED`、`DETACHED` 身份分开；仅拓扑断开不等于已经远离壁面，另保存 COM 位移和速度。
4. 以参考 nodal volume `V_node = dx³` 累加体积：
   `V_fragment = Σ V_i`；`mass = ρ V_fragment`。
5. 计算**体积等效球直径**：
   `d_eq = (6 V_fragment / π)^(1/3)`。
   它不是最大粒径、最长几何尺寸，也不是直接测量的真实碎片直径。

### 两档尺寸尺度与节点门槛（必须执行精确取整）

**文献关注尺度：`d_eq ≥ 200 µm`**。机械取栓文献使用 `>200 µm` 碎片分组，本轮的 ≥200 µm 统计是接近该实验分组的参考量，**不是普适临床危险阈值**。

**工程可辨认尺度：`d_eq ≥ 300 µm`**。这是考虑当前 `dx=50 µm` 的项目工程定义，**不是原论文测得的血栓破裂/栓塞阈值**。

本轮 `dx=50 µm` 时：

- 1 node `V=125000 µm³`
- `d_eq ≥ 200 µm` ⇔ **至少 34 nodes**
- `d_eq ≥ 300 µm` ⇔ **至少 114 nodes**（注意不是 113 nodes；门槛应向上取整）
- 49-node component `d_eq ≈ 227 µm`，应进入 ≥200 参考候选，但不达 ≥300 工程大块

禁止继续把 `≥3 nodes` 简称“大碎块”；可保留旧分类器定义以保证历史可比性，但另设新的尺寸统计字段。

### 空间形态的最小提示

大于 200/300 µm 的 component 仍可能是细长的一串节点，不能只依靠等效球直径就说是“立体大块”。对所有 ≥114-node 候选额外报告：

- 参考构型的 PCA/主方向包围盒三轴尺度；
- 最短轴是否至少约 `2dx`；不满足则标 `GEOMETRICALLY_THIN_CANDIDATE`；
- 不删除或隐藏这些 component，只如实标记。

**不在本轮追求更细网格或物理准确的 fragment-size distribution。** 小于当前粒子分辨率的碎屑必须标 `NUMERICALLY_UNRESOLVED_DEBRIS`，不能称“真实 50 µm 血栓微粒”。

### 生成、持久性和消失

记录每一帧：

- `N_detached_ge200`, `N_detached_ge300`
- `mass_detached_ge200`, `mass_detached_ge300`，及相对于初始 clot 质量的比例
- `largest_detached_nodes / d_eq / mass`
- `largest_attached_component_fraction`
- `debris_1_2node_fraction`
- fragment COM、velocity、interface attachment flag

汇总：

- `EVER_OBSERVED`：任一保存帧出现过
- `MAX_SIMULTANEOUS_COUNT`
- `PEAK_DETACHED_MASS_FRACTION`
- `AT_FINAL_FRAME`
- `FIRST_APPEARANCE_TIME`
- `LAST_SEEN_TIME`
- `FRAME_PERSISTENCE_COUNT`

如果做跨帧身份跟踪，用 particle-ID 集合最大交叠并保存分裂关系；明确标记 tracking ambiguity，不能因为最大 ID 相同就声称是同一完整碎块。可将 ≥3 连续保存帧作为“持续出现的工程标签”，但同时展示精确持续时间，并注明不是生物学阈值。

## E. Stage 1：另一个 streaming pattern 的接触修正迁移（必做一例）

新增正式算例：

`OPPOSING_SHEAR | tau_budget=100 Pa | duration=2T_cross | CONTACT_MODIFIED=1.0dx`

- 所有材料/几何/时间参数按 C 冻结；
- 与历史 `OPPOSING_SHEAR 100 Pa 2T_cross / native 1.35dx` 只读对照；
- 与历史 `GAUSSIAN 100 Pa 2T_cross / modified 1.0dx` 一起看接触修正是否跨 pattern 可用；
- 不为实现“类似 Gaussian”而调 Opposing 自己的强度。

记录真实首次 bulk/interface break、bulk/interface 断键率、debris、最大附着主体、≥200/≥300 游离块、动能、接触活跃对、最小非键合距离、深度重叠与全历程分类。

原 native 1.35dx case 不要重新跑。新 case 只能改变 contact configuration 一项相对于自己的旧同 pattern 输入。

**接触修正的跨 pattern 工程判据：**

- 模型包含真正的非键合接触（不等于 contact OFF）；
- 相比同 pattern 旧 native 结果，断后非物理级联明显下降，或至少没有无原因地增加；
- 对保存帧无明显严重互穿（沿用上轮 `distance <0.5dx` 的关键诊断线，并报告实际比例和持续时间）；
- 无 lost atoms、missing bonds、NaN、时间提前截断。

若仍然粉碎，记录 `CROSS_PATTERN_CONTACT_TRANSFER=NOT_ESTABLISHED`；**不得自动改接触距离、s00 或阻尼**。

## F. Stage 2：一个更强载荷（必做一例）

新增正式算例：

`LOCAL_GAUSSIAN_SHEAR | tau_budget=200 Pa | duration=2T_cross | CONTACT_MODIFIED=1.0dx`

只改变 `tau_budget`，其它与历史 `GAUSSIAN 100 Pa / modified` 相同。

这是**单个强载荷探针**，不是振幅扫描。即使结果仍是 PULVERIZATION，也按真实结果验收；**不允许为了画出 200/300 µm 大块而追加 300、400、500 Pa 或重新调 s00**。

特别记录：

- 有没有**大块曾经出现、随后再碎裂**；
- 有没有**持续存在到最后的游离大块**；
- ≥200 和 ≥300 µm 两个尺度下的峰值、终值、生命周期；
- 是否整块脱粘、局部碎裂、保持主体、还是再度进入颗粒化；
- 200 Pa 下若又粉碎，是在何时出现、与 contact event/动能是否同步，不把增强外载造成的粉碎直接归咎于接触修正失败。

## G. Stage 3：旧 force–duration map 只做**有条件**的局部更新

上一轮 37-case map 属于 `CONTACT_1.35dx` **历史物理模型**。本轮 `CONTACT_1.0dx` 是另一物理模型配置，必须版本分离：

`contact_model_id = NATIVE_1P35DX | MODIFIED_1P0DX`

原 37 个结果、分类和图片一个也不能改。不得拿旧 1.35dx 标签当作新 1.0dx 标签。

**只有在 Stage 1、2 的正式结果都正常且 Stage 1 支持新接触候选后**，允许最多 **2 个**新模型近边界补点，优先：

1. `LOCAL_GAUSSIAN_SHEAR / 31.6227766 Pa / 2T_cross / 1.0dx`
2. `OPPOSING_SHEAR / 31.6227766 Pa / 2T_cross / 1.0dx`

这两点在旧 1.35dx map 中靠近此前的 NO_FAILURE→PULVERIZATION 转换区域。先在 `CONDITIONAL_CASE_DECISION.json` 明确选择/拒绝原因和冻结的输入，然后再运行。

如果 Stage 1 接触修正未跨 pattern 起效、Stage 2 出现明显数值异常，或现有两例不足以支持区分边界，则**不运行**这两个补点，标 `MAP_REFRESH=DEFERRED`。

最多正式新增求解数：**2 必做 + 2 条件执行 = 4**。不跑 37 点全图、不重做旧 contact ON、不开展网格收敛、无其他强度/持续时间/载荷范式扫描。

给出 `LOCAL_CONTACT_MODEL_TRANSITION_COMPARISON`：仅以**真实执行的点**展示 1.35 与 1.0 模型的差异。不得凭这两三点平滑拟合完整新边界，未测区域标 UNKNOWN。

## H. 最小必要验证与统计边界

只检查和本轮直接相关的内容：

1. 新碎块统计：至少 1/2/3/5/34/49/113/114 nodes 的阈值判别与 d_eq 计算；附着/游离的 interface graph 判别；避免墙节点计数；
2. 正确读回历史保存帧与当前两个新 case 的真实断键事件；一条已有 49-node 游离块记录能被重现或若不一致给出原因；
3. 新 case 与冻结旧输入做 diff，仅允许 pattern、amplitude 与接触配置这些命名维度变化；
4. 原生非键合接触仍生效，不存在持续严重重叠；
5. 各新 case 实际 timestep、时间窗口、质量、运行状态完整，旧历史哈希不变。

**禁止多余验证：** 不做 synthetic stress benchmark、MPI/GPU 评测、PD2/FEM 比较、材料/强度重新校准、长时间疲劳、VES、软化 C++、contact factor sweep、horizon/mesh convergence、压力-only、重复 37-case map、微泡 CFD、PZT 或 surrogate。

参考文献与解释边界：

- Chueh, Puri & Gounis (2017), mechanical thrombectomy distal emboli grouping `>200 µm`：https://pubmed.ncbi.nlm.nih.gov/26919972/ ，DOI https://doi.org/10.1136/neurintsurg-2015-012256 。这是**取栓碎屑尺寸分组**，不是血栓碎裂普适阈值。
- Microtripsy clot-debris experiment：https://pmc.ncbi.nlm.nih.gov/articles/PMC4528908/ 。其细碎屑通常远小于当前 50 µm 模型能分辨的尺度，不能由 LAMMPS 单粒子推断真实碎屑分布。
- Liu et al. 2024 whole-clot fracture/interface：https://doi.org/10.1016/j.jmps.2024.105858 。仍是来源身份与现有校准约束；不借本轮改动重新拟合。

## I. 图像、短片和数据

输出 3 张英文科学图（PNG/PDF），优先复用现有绘图程序，不触发新 solver：

1. `contact_transfer_and_stronger_load.png/pdf`

   - 历史 Gaussian 100 Pa / 1.0dx
   - 新 Opposing 100 Pa / 1.0dx
   - 新 Gaussian 200 Pa / 1.0dx
   - 并列各自对应旧 native 1.35dx（存在的才展示）
   - bulk/interface broken bonds、debris fraction、attached bulk、真实 contact、分类时间线
2. `detached_fragment_size_and_lifetime.png/pdf`

   - 时间轴上 ≥200、≥300 µm **游离** fragment 数及质量占比
   - `ever / peak / persistent / final` 清楚分离
   - 标记历史 49-node ~227 µm 例子，若再次识别成功
   - 图下注明 `200 µm = literature reporting scale; 300 µm = dx-based engineering display threshold`。
3. `targeted_contact_map_refresh.png/pdf`

   - 老接触与新接触分图或不同图层，不能当同一个 map
   - 只放真实已运行点；未采样区域明确为空/UNKNOWN
   - 若 Stage 3 不执行，改成小范围旧/新类别比较图，不制造插值热图。

生成一个 3D 同步短片 `contact_transfer_fragment_size.mp4` 和 GIF preview：

- 至少展示两个本轮新实际案例和历史 Gaussian100 1.0dx；
- 同一相机、同一物理时刻；只用 solver 保存帧，绝不插值造断键；
- 粒子和 surviving bonds 全保留，anchored base 与 detached components 有明确区分；
- ≥200/≥300 µm 的游离 component 可在侧边统计/描边强调，但**不能把小碎屑隐藏、放大成大块、用颜色假切**；
- fps 不代表真实实验或超声时间。

不需要为每个补点生成独立 600dpi 四图或高帧数动画；可压缩渲染成本。

## J. 接口与结果结构

继续扩展当前 `pd_clot/lammps_fragmentation/`，不另起完整求解器。

建议最小模块：

- `detached_fragment_metrics.py`
- `fragment_lifetime.py`（仅在需要时）
- `targeted_contact_transfer.py`
- `targeted_map_refresh.py`

每个正式 case 保存：

`contact_model_id, material_id, interface_id, pattern_id, tau_budget_Pa, pulse_duration_s, contact_distance_factor, dx_m, horizon_m, time_end_s, n_steps, first_break_times, broken_bulk_fraction, broken_interface_fraction, debris_fraction, largest_attached_fraction, detached_ge200_count_peak, detached_ge200_count_final, detached_ge300_count_peak, detached_ge300_count_final, detached_ge200_mass_fraction_peak, detached_ge300_mass_fraction_peak, largest_ever_detached_d_eq_um, largest_final_detached_d_eq_um, component_lifetimes, dynamic_failure_classification, min_nonbonded_distance_over_dx, wall_seconds, provenance, EXPLORATORY_ONLY`。

注意：之前的 `RESOLVED_FRAGMENTATION` 工程状态（≥3 nodes）与本轮 `LITERATURE_SIZE_200UM_DETACHED` / `ENGINEERING_SIZE_300UM_DETACHED` **分别保存**，不可覆盖历史定义。

## K. 交付与验收状态

建议目录：

`verification/lammps_contact_transfer_and_fragment_size_v1/`

`results/lammps_contact_transfer_and_fragment_size_v1/`

至少交付：

- `TASK_REQUEST.md`, `STARTING_GIT_STATE.md`, `FROZEN_CONFIGURATION.json`
- `SOURCE_REGISTER.md`, `SIZE_THRESHOLD_POLICY.md`, `CONTACT_VARIANT_PROVENANCE.md`
- `HISTORICAL_GAUSSIAN_REANALYSIS.json`
- `CASES/OPPOSING_100_D2/` 原生输入、日志、broken-bond events、轨迹、后处理
- `CASES/GAUSSIAN_200_D2/` 同上
- 若执行 Stage 3：`CASES/GAUSSIAN_31P622_D2/` 与 `CASES/OPPOSING_31P622_D2/`
- `PER_FRAME_DETACHED_COMPONENTS.csv`, `FRAGMENT_LIFETIMES.csv`, `FRAGMENT_SIZE_SUMMARY.csv`
- `CROSS_PATTERN_CONTACT_DECISION.json`, `CONDITIONAL_CASE_DECISION.json`
- `TARGETED_OLD_NEW_MAP_COMPARISON.csv`, `PERFORMANCE_REAL_CASES.json`
- 三张科学图、一个 MP4/GIF、绘图数据与 QA
- `FINAL_STATUS.json`, `REVIEW_LAMMPS_CONTACT_TRANSFER_AND_FRAGMENT_SIZE_V1.md`, `evidence.zip`

`FINAL_STATUS.json` 至少区分：

- `CONTACT_1P0_TRANSFER = SUPPORTED_FOR_TESTED_PATTERNS | NOT_ESTABLISHED | FAILED`
- `CONTACT_1P0_PRODUCTION_DEFAULT = CONDITIONAL_EXPLORATORY_CANDIDATE | NOT_PROMOTED`
- `SECOND_PATTERN_TEST = PASS | FAIL`
- `STRONGER_LOAD_TEST = PASS | FAIL`
- `OLD_MAP_REFRESH = TARGETED_POINTS_ONLY | DEFERRED`
- `DETACHED_SIZE_200UM = OBSERVED | NOT_OBSERVED`
- `DETACHED_SIZE_300UM = OBSERVED | NOT_OBSERVED`
- `PERSISTENT_LARGE_DETACHED_FRAGMENT = OBSERVED | NOT_OBSERVED`
- `WHOLE_CLOT_DETACHMENT = OBSERVED | NOT_OBSERVED`
- `PULVERIZATION_AT_200PA = OBSERVED | NOT_OBSERVED`
- `LPS_PROGRESSIVE_SOFTENING = NOT_IMPLEMENTED`
- `VES = DEFERRED`
- `B1_B2 = TERMINATED_BY_PROJECT_POLICY_CHANGE`
- `PD2_FRACTURE_DEVELOPMENT = STOPPED`
- `REAL_BUBBLE_STREAMING = NOT_CONNECTED`
- `BIOLOGICAL_QUANTITATIVE_VALIDATION = NOT_ESTABLISHED`

### 何时允许“接触候选升级”？

只有 Stage 1 在另一个 pattern 显示稳定的真实断键/碎块运动、活跃接触、无持续严重重叠、避免不必要级联，而且参数冻结，才允许写：

`CONTACT_1P0_PRODUCTION_DEFAULT = CONDITIONAL_EXPLORATORY_CANDIDATE`

意思仅是：**可作为以后探索性 screening 的候选默认接触配置**，不覆盖旧 map、不声称真实 clot 接触半径已校准、也不表示已在所有载荷上保证稳定。若 200 Pa 再次粉碎，作为超载工况如实记录，不自动推翻 100 Pa 的独立证据。

### 停止规则

两例必做完成、尺寸分析完成、条件补点按已声明分支执行或跳过、图件/短片完成后立即停止。即使没有出现 ≥300 µm 大块，也算诚实结果，**不能继续加力或调材料强行得到预想动画**。

## L. REVIEW 必须用中文回答的关键问题

1. Opposing 100 Pa / 2T 的 `1.0dx` 接触是否保持抑制异常粉碎的效果？与旧 `1.35dx` 有什么差别？
2. Gaussian 200 Pa 下是更多可辨认游离碎块、整块脱落、仍保持主体，还是再次粉碎？何时发生？
3. 过去 49-node 游离块的等效尺寸与生命周期是多少？为什么不能把 `≥3 nodes` 叫大块？
4. 本轮真实出现过多少个 `d_eq≥200 µm` 和 `≥300 µm` 的游离块？峰值、最终、持续时间分别是多少？
5. 有没有细长结构被等效球直径误认为三维大块？
6. 旧 37-case map 哪些失效边界仍可解释、哪些已与 `1.0dx` 模型不兼容？新增了哪几个接触修正点？
7. `1.0dx` 是否有资格升为有限范围内的 exploratory screening 默认候选？为什么？
8. 单个真实 case 的 wall time、后处理和渲染成本是多少？
9. 当前是否已经有证据要求 progressive softening、VES、PeriLab？若没有，明确不要启动。
10. 将来 microbubble→streaming→fragmentation 时，怎样分别优化“碎裂程度”和“大游离块带来的栓塞风险”，而不是简单最大化碎块体积？

完成后停止，等待人工审查。