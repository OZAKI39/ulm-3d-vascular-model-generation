# PROJECT_ACCURACY_AND_MODEL_HIERARCHY_SPEC_V1

**项目名称**：ULM / 3D vascular model generation — clot mechanics branch  
**规范版本**：V1  
**日期**：2026-10-07  
**规范性质**：项目级开发与验收政策。用于后续 Codex / 人工开发、仿真数据生成、模型升级、代理模型训练与结果审查。  
**适用阶段**：从当前可调用完整态血栓 FEM 模块开始，覆盖材料、损伤/断裂、微泡/CFD 载荷映射、代理模型及 inverse PZT control。  

---

## 1. 项目目标重新定义

本项目的默认目标不是“获得单一数值方法的最高保真结果”，而是构建一个：

> **机制导向（mechanism-oriented）、控制导向（control-oriented）、计算成本可承受的近似仿真器。**

它必须优先支持：

1. 不同微泡形态 / 声学输入 / 表面载荷之间的**趋势判断**；
2. 不同候选治疗方案之间的**排序**；
3. 位移、应变、能量、破坏风险、碎裂程度等结果的**正确数量级**；
4. 以足够低的单样本成本生成 surrogate / inverse-control 所需数据；
5. 只有在简单模型出现明确风险时，才升级到高成本模型。

因此，项目以 **CONTROL-ORIENTED SUFFICIENCY V2** 为最高开发原则。

历史上已经完成的 PD2 / 边界 / 曲面验证继续作为科学证据保留，不因本规范失效；但它们不再决定所有后续工况必须使用 PD2。

---

## 2. 三条不可放宽的物理底线

任何模型，即使属于“低保真”或“探索级”，都必须满足以下三条底线。

### 2.1 方向 / 符号必须正确

对于有明确方向意义的响应量：

- 压力不能产生与约定相反的主要法向响应；
- 剪切主方向不能反转；
- 反力、功、能量符号必须符合定义；
- 损伤不能导致已知非物理的“越坏越储能”一类反向行为。

若简单方法在方向或符号上与可信参考、解析核算或实验事实矛盾，则该模型不得用于该用途，即使其标量相对误差较小。

向量方向可用

\[
\cos\theta=
\frac{\mathbf q\cdot \mathbf q_{\rm ref}}
{\|\mathbf q\|\,\|\mathbf q_{\rm ref}\|}
\]

衡量。探索性默认建议 \(\theta\le 30^\circ\)；更严格用途可另行规定。

### 2.2 数量级必须正确

如果已有可信参考或实验范围，则结果不得跨越明显错误的数量级。

对正标量可定义

\[
r_{\rm mag}
=
\frac{\max(|q|,q_{\rm floor})}
{\max(|q_{\rm ref}|,q_{\rm floor})}.
\]

当没有更强的用途标准时，探索性结果至少应满足

\[
\frac13 \le r_{\rm mag}\le 3.
\]

这只是“数量级正确”的最低底线，不代表已经达到定量验证。

### 2.3 排序必须稳健

如果任务的真正问题是“A 是否优于 B”，则排序稳定性优先于单点高精度。

若某正响应量的默认相对不确定度记为 \(\varepsilon\)，只有当

\[
Q_A(1-\varepsilon)>Q_B(1+\varepsilon)
\]

时，才把 \(A>B\) 视为稳健排序。

对本规范默认探索容差 \(\varepsilon=0.30\)：

\[
\frac{Q_A}{Q_B}>
\frac{1+0.30}{1-0.30}
\approx 1.857
\]

时，最坏情形下排序仍不反转。

不满足该条件时，输出应为：

`RANKING_UNCERTAIN`

而不是强行给出优劣结论。

---

## 3. 精度等级

不同输出和用途采用不同精度等级，不再使用“一项 FAIL 阻断全部模型”的策略。

### Tier A — Screening / Exploration

用途：

- 参数扫描；
- 微泡形态或载荷候选筛选；
- surrogate 初始数据；
- 机制趋势图；
- 快速 sensitivity。

默认连续量相对差：

\[
e_q=
\frac{|q-q_{\rm ref}|}
{\max(|q_{\rm ref}|,q_{\rm floor})}
\le 30\%.
\]

30% 是**项目默认探索容差**，不是材料定律、实验置信区间或论文定理。

同时必须满足：

- 方向/符号正确；
- 数量级正确；
- 排序满足 §2.3，或明确标 `RANKING_UNCERTAIN`。

### Tier B — Threshold / Damage Classification

用途：

- “会不会坏”；
- “是否接近脱落”；
- “是否超过损伤阈值”。

不能简单套用 30%。

若使用阈值 \(Q_c\)，定义

\[
m=\frac{Q-Q_c}{Q_c}.
\]

若模型相对不确定度为 \(\varepsilon\)，则：

**稳健超阈值**

\[
Q(1-\varepsilon)>Q_c
\]

**稳健低于阈值**

\[
Q(1+\varepsilon)<Q_c.
\]

对默认 \(\varepsilon=0.30\)：

\[
Q/Q_c > 1.429
\]

才稳健判为超阈值；

\[
Q/Q_c < 0.769
\]

才稳健判为低于阈值。

中间区间一律标：

`NEAR_THRESHOLD_UNCERTAIN`

并升级到更高模型层级或实验核对。

### Tier C — High-value Quantitative Reference

用途：

- 论文定量主结果；
- surrogate 临界区校准；
- fracture onset 附近；
- 实验对照；
- 控制器关键约束。

默认目标应显著严于 Tier A。若没有实验规定，数值方法间可以暂以 10% 左右作为高价值参考目标，但**不得把这个数自动推广到所有量**。

必须逐输出定义：

- 允许相对差；
- 允许绝对差；
- reference scope；
- 数据/实验来源。

---

## 4. 模型层级

以后所有仿真必须显式声明 `model_level`。

### Level 0 — Reduced / Analytical Load Model

职责：

- 将微泡 / CFD / 声学结果压缩为低维载荷描述；
- 快速解析核算；
- 载荷幅值、位置、宽度、方向和时间尺度描述。

典型变量：

\[
z_{\rm load}
=
[p,\tau_x,\tau_y,c_x,c_y,c_z,\ell,\ldots].
\]

不用于裂纹拓扑。

### Level 1 — Fast FEM Screening Model（默认）

**默认日常求解器。**

职责：完整态大变形、局部压力/剪切、wall kinematics、位移/应变/能量/反力、damage-risk proxy、大规模参数扫描和数据集生成。

当前已有 callable FEM 接口可直接作为 Level 1 基础。

原则：

- 优先 native FEM；
- 优先成熟 DOLFINx/PETSc；
- 允许 Tier A 级近似；
- 不要求与 PD2 全场一致；
- 现有历史 PD2–FEM 约 10–11% 位移差不再自动触发新审计。

### Level 2 — FEM Phase-field / Cohesive Fracture

仅用于：

- Level 1 判定为高风险；
- 接近破坏阈值；
- 需要裂纹路径、峰值载荷和承载下降；
- 需要 fracture surrogate 的少量标签。

不是默认 dataset generator。

优先采用成熟 phase-field / cohesive implementation，不为了项目重新发明断裂数值方法。

### Level 3 — PD2 High-fidelity Spot-check

PD2 新角色：

> **reference expert / spot-check model**

而不是默认数据生成器。

触发条件：

1. Level 1 给出方向异常；
2. Level 1 与 Level 2 对临界结论冲突；
3. wall / surface / nonlocal 区域明显支配结果；
4. surrogate 不确定度高且该样本价值高；
5. 需要检查非局部效应；
6. 论文中需要展示方法敏感性。

不得因为 10–20% 的正常连续量差异自动触发 PD2。

建议数据集中 PD2 只占少量高价值样本；具体比例由成本和相关性实测决定，不预设固定百分比。

### Level 4 — Explicit Dynamic / Fragmentation

仅在真正需要以下物理时启用：

- 微泡 collapse / jet impact；
- 高频冲击；
- 快速惯性响应；
- 裂后碎块拓扑；
- fragment size distribution；
- 动态抛射 / 碰撞。

此层可采用 simplified bond-based PD、显式 FEM、已发表 clot-PD fragmentation 路线；GPU 仅在规模足够且实测收益明确时使用。

不得为了平均 microstreaming / radiation-force 问题启用这一层。

### Level 5 — Fluid / Fragment Transport

职责：fragment advection、aspiration capture、downstream escape、血流/吸力竞争。

输入来自 Level 2 / 4 的 fragment geometry / release state。固体断裂求解器不负责把碎块一直运输到下游。

---

## 5. 准静态与动态的分流规则

定义特征时间

\[
T_{\rm clot}\sim \frac{L}{c_s}.
\]

若载荷变化时间尺度 \(T_{\rm load}\) 满足

\[
T_{\rm load}\gg T_{\rm clot},
\]

优先按准静态问题处理。

适合准静态/隐式：周期平均 microstreaming、acoustic radiation force、缓慢 wall stretch、长时间平均 pressure/shear。

需要动态模型：microbubble collapse、jet impact、单周期高速应力波、高频循环累积损伤（若不能用有效疲劳模型表示）、fragment release 后惯性运动。

**显式/隐式不按“哪个更快”统一选择，而按物理时间尺度选择。**

---

## 6. Damage-risk proxy 与真实 fracture 必须分开

Level 1 可定义快速 proxy，例如：

\[
R_{\Psi}=\frac{\Psi_{\rm driver}}{\Psi_c}
\]

或

\[
R_\varepsilon=\frac{\varepsilon_{\rm eq}}{\varepsilon_c}.
\]

它只能回答：哪里风险更高、哪个方案更危险、哪些样本值得升级。

输出名必须包含 `RISK_PROXY`，不得称 crack / fracture probability / detached volume / true damage。

只有 Level 2 / 4 才允许输出真实 fracture / fragmentation 类型结果。

---

## 7. 升级规则

所有工况从最低足够层级开始。

### Level 1 → Level 2

满足任一：risk proxy 落入 `NEAR_THRESHOLD_UNCERTAIN`、需要裂纹路径、需要峰后承载、需要实验 fracture 对照。

### Level 1/2 → Level 3

满足任一：方向/排序与预期或实验冲突、surface/wall/nonlocal effect 是结论核心、FEM 与已有 PD2 spot-check 差异超过该用途允许值、surrogate 在高价值区预测不确定。

### 任意准静态层 → Level 4

仅当

\[
T_{\rm load}\lesssim T_{\rm clot}
\]

或明确研究 impact / fragmentation dynamics。

### Level 2/4 → Level 5

只有当需要回答碎块捕获/逃逸和 downstream embolization。

---

## 8. 计算成本政策

1. 重用已有状态；
2. Level 0 / Level 1 优先；
3. 独立工况 case-level 并行；
4. Level 2；
5. Level 3；
6. Level 4 / 5。

每个昂贵求解器先 profiling：geometry/reference preprocessing、assembly、nonlinear iteration、linear solve、postprocessing、IO。

未经 profiling，不启动 GPU rewrite、FFT 改造或大规模算法重构。

数据集阶段优先 case-level CPU parallelism。

GPU 仅在 workload 足够大、数据搬运可接受、端到端实测收益 >10%、数值 parity 通过时采用。否则 CPU 是正式路径。

---

## 9. 数据集与 surrogate 策略

大部分标签来自 Level 1 FEM；少量临界标签来自 Level 2；少量 high-value spot checks 来自 Level 3；真正碎块拓扑来自 Level 4。

每个标签至少保存：

- `model_level`
- `method_id`
- `material_id`
- `geometry_id`
- `load_id`
- `quantity_definition`
- `reference_scope`
- `estimated_model_uncertainty`
- `validation_status`

不能把不同模型标签当作无差别 ground truth。

建议最终 surrogate 输出

\[
(\hat y,\hat\sigma_{\rm model})
\]

即预测值及模型不确定度。临界控制区不确定度过大时，回到高层级 physics oracle。

---

## 10. 项目停止规则

以下情况不再自动开展新数值方法研究：连续量跨方法差 <30% 且方向/量级/排序稳定；误差不改变治疗方案排序；误差不改变 damage/fracture 分类；surrogate 已满足用途；只有“数学上还能更精确”但无明确控制收益。

只有真实 application blocker 才重新打开 PD2 boundary、更多 RK、GPU solver rewrite、FFT、额外 mesh-objectivity/stability campaign。

---

## 11. 规范与历史结果关系

本规范不改写任何历史结论。PD2 kinematics/variational mechanics 的验证、历史曲面 V1 的 NOT_ESTABLISHED、历史 5/10/20% gates、当前 callable intact FEM module 的限定 PASS 都继续保留。

变化的是：**这些历史数值不再自动成为所有新用途的开发阻断条件。**

---

## 12. 规范性状态关键词

- `PASS_FOR_SCREENING`
- `PASS_FOR_RANKING`
- `PASS_FOR_QUANTITATIVE_REFERENCE`
- `NEAR_THRESHOLD_UNCERTAIN`
- `RANKING_UNCERTAIN`
- `EXPLORATORY_ONLY`
- `SPOT_CHECK_REQUIRED`
- `DYNAMIC_MODEL_REQUIRED`
- `NOT_SUPPORTED_FOR_FRACTURE`
- `NOT_SUPPORTED_FOR_CLINICAL_CLAIMS`

禁止用单一 `PASS` 覆盖所有用途。

---

## 13. 来源与政策依据

本规范的调整依据包括：当前项目已具备 callable intact FEM 接口、局部载荷、壁面运动学和 P2/P3 限定检查；真实下游 control tolerance 目前仍为空；前几轮已证明约 10–11% 的 FEM/PD2 差异不应自动触发更多方法审计；用户明确项目只需要“大概正确”的机制、趋势、排序和数量级。

加速综述同时指出：当前约 3000 粒子的隐式高阶 RK PD 与大量 GPU/FFT 文献适用场景明显不同，第一步应 profiling；数据生成最有效的并行优先是独立工况并行；动态和准静态问题应按时间尺度分流。

本规范中的 30% Tier-A 默认值、排序 margin、阈值不确定区和模型层级，是本项目新的工程政策，不是附件中论文或对话直接给出的定理。
