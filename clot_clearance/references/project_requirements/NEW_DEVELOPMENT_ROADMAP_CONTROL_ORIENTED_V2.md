# NEW_DEVELOPMENT_ROADMAP_CONTROL_ORIENTED_V2

**日期**：2026-10-07  
**路线版本**：V2  
**最高原则**：低成本优先、用途驱动、分层升级、只在关键样本使用高保真。

---

# 0. 当前起点

已经具备：曲面半球血栓几何、native FEM P2/P3、局部 pressure/shear/mixed load、prescribed wall kinematics、callable request→solve/load→evaluate→export、保存的 PD2 spot-check 能力、已知历史跨方法差异与适用边界。

当前尚缺：实验来源血栓材料、快速 damage-risk proxy、少量真实 fracture、bubble/CFD→reduced load、surrogate/inverse control、快速 impact/fragment transport 支路。

---

# Phase 1 — 项目政策落地与接口收尾

1. 更新 PROJECT_ACTIVE_CONTEXT：
   - `CONTROL_ORIENTED_SUFFICIENCY_V2=ACTIVE`
   - `DEFAULT_SOLVER=FEM`
   - `PD2_ROLE=SPOT_CHECK_REFERENCE`
   - `DEFAULT_SCREENING_TOLERANCE=30_PERCENT_WITH_DIRECTION_MAGNITUDE_RANKING_GATES`
2. 修复接口：REUSE_RESULT 用途上下文；质量失败 `response=null`，诊断放 `diagnostic_response`。
3. 结果增加 `model_level/intended_use/accuracy_tier/estimated_uncertainty/escalation_status`。

不新增科学平衡求解。

---

# Phase 2 — 实验来源材料 + Fast FEM Screening

默认 Level 1 FEM。

首轮只接一套来源明确的超弹性材料，不跨文献拼参数。

代表性工况：global pressure、global shear、localized pressure、localized shear、mixed、wall stretch。

输出：displacement、F/J、stress、energy、wall reaction、damage-risk proxy。

Tier A 验收：数值求解通过；方向/量级合理；不要求与 PD2 一致到 10%。

---

# Phase 3 — Fast Damage-Risk Proxy

建立低成本“哪里更可能先坏”的指标，例如

\[
R_\Psi=\Psi_{\rm driver}/\Psi_c
\]

或实验支持的 stretch/stress/energy threshold。

明确 `RISK_PROXY != FRACTURE`。

输出：risk map、high-risk volume fraction、hotspot regions、robust/uncertain threshold classification。

明显低风险不升级；明显高风险进入 fracture candidates；near-threshold 进入 Level 2。

---

# Phase 4 — Selected FEM Phase-field Fracture

只对 Phase 3 高价值/临界样本计算真实裂纹。

首轮建议 3–10 个代表性 cases，而不是整个参数空间。

输出：crack/damage field、crack extension、reaction-displacement、load drop、fracture energy、damaged volume。

目标是方向/区域/分类正确，不追求所有局部场极高精度。

---

# Phase 5 — Bubble / CFD → Reduced Load Model

将完整 CFD/bubble 字段压缩为

\[
z_{\rm load}=[p,\tau,c,\ell,\text{orientation},\ldots].
\]

若需要时间变化，则使用 \(z_{\rm load}(t)\)。

Mean/slow loading（streaming、radiation force）进入准静态 Level 1/2；collapse/jet/shock 标记 `DYNAMIC_MODEL_REQUIRED`，进入 Phase 8。

---

# Phase 6 — Dataset Generation

大部分：Level 1 FEM。  
临界：Level 2。  
spot-check：Level 3 PD2。

参数空间覆盖 geometry、material、load amplitude/location/width、wall condition、bubble-derived latent inputs。

并行优先 case-level CPU parallelism。

---

# Phase 7 — Surrogate + Active Sampling

输入：

\[
z_g+z_m+z_{\rm load}
\]

输出：

\[
z_{\rm response}
\]

包括 displacement metrics、energy/stress summaries、risk proxy、fracture classification、damaged volume、crack/fragment descriptors。

先比较低成本模型；不预先锁定 FNO。主动补样区域：高 surrogate uncertainty、ranking 可能反转、near-threshold、control high-value region。

---

# Phase 8 — Optional Dynamic / Fragmentation Branch

仅在 microbubble collapse/jet、高频循环损伤、碎块 topology/size 或 launch velocity 重要时启动。

优先 simplified explicit PD、published clot-PD failure model 或 explicit FEM；GPU 只在规模实测值得时使用。PD2 不自动成为此支路默认方法。

---

# Phase 9 — Fragment Transport / Aspiration

断裂模块只输出 fragment geometry/size/release time/release velocity，之后交给 CFD/particle tracking/aspiration model，回答 capture vs escape。

---

# Phase 10 — Inverse PZT Control

最终控制链：

\[
z_{\rm PZT}\rightarrow z_{\rm bubble}\rightarrow z_{\rm load}\rightarrow z_{\rm clot}.
\]

逆问题：

\[
z_{\rm clot}^{target}\rightarrow z_{\rm load}^{target}\rightarrow z_{\rm bubble}^{target}\rightarrow z_{\rm PZT}.
\]

高保真 solver 不进入在线优化循环。

---

# 新路线默认决策树

```text
Start
  |
  v
Reduced/CFD load descriptor available?
  |
  v
Level 1 FEM screening
  |
  +-- direction/sign wrong? --------> SPOT_CHECK / DEBUG
  |
  +-- order of magnitude wrong? ----> SPOT_CHECK / MATERIAL REVIEW
  |
  +-- ranking margin weak? ---------> MORE SAMPLES / HIGHER LEVEL
  |
  +-- near damage threshold? -------> Level 2 phase-field
  |
  +-- surface/nonlocal conclusion
  |   is scientifically central? ---> Level 3 PD2 spot-check
  |
  +-- fast impact/collapse? --------> Level 4 dynamic
  |
  +-- fragment transport needed? ---> Level 5 CFD/transport
  |
  `-- otherwise --------------------> accept Level 1 result
```

---

# 近期三轮开发建议

## Round A — POLICY_AND_FAST_MATERIAL_SCREENING_V1

- 规范落地；
- API 两项修复；
- 新材料卡；
- 6 个代表性 FEM screening；
- damage-risk proxy 初版；
- 四张验证图；
- 不做完整 fracture campaign。

## Round B — SELECTED_FRACTURE_DEMONSTRATION_V1

- 从 Round A 选 3–5 个样本；
- FEM phase-field；
- reaction–displacement；
- crack/damage evolution；
- 输出动画；
- 对 near-threshold 情况定义 uncertainty。

## Round C — BUBBLE_TO_LOAD_REDUCED_INTERFACE_V1

- 从 microbubble / CFD 结果提取 \(z_{\rm load}\)；
- 连接现有 FEM solver；
- 做第一组 bubble-shape→clot-response ranking；
- 为 surrogate 数据生成做准备。

---

# 不再自动进行的工作

除非真实应用阻断：更多 PD2/FEM 11% 差异审计、额外 cube polynomial tests、surface correction variants、大规模 Hessian/stability、PD2 GPU rewrite、FFT-PD 改造、为方法接近而调物理参数、所有样本都做 fracture。

---

# 新路线每轮只问四个问题

1. 这个开发会不会改变治疗方案的趋势或排序？
2. 会不会改变 damage/fracture 分类？
3. 会不会降低数据集生成成本？
4. 会不会让 bubble→load→clot→control 更接近闭环？

若四个答案都是“不会”，则优先级降低。


## QUALITATIVE_FRAGMENTATION_SCREENING_V1 — 2026-10-07 policy change

The new user-authorized fragmentation route supersedes the prior B1/B2 production route; historical intact FEM and PD2 results remain preserved.

- QUALITATIVE_FRAGMENTATION_SCREENING_V1 = ACTIVE
- PRIMARY_FRAGMENTATION_SOLVER = LAMMPS_BPM_PERI_LPS
- PRIMARY_FRAGMENTATION_USE = QUALITATIVE_PATTERN_SCREENING
- PHASE_FIELD_FRACTURE_ROLE = ARCHIVED_OPTIONAL_REFERENCE_ONLY
- B1_B2_ROUTE = TERMINATED_BY_PROJECT_POLICY_CHANGE
- PD2_FRACTURE_DEVELOPMENT = STOPPED
- PD2_ROLE = OPTIONAL_INTACT_SPOT_CHECK_REFERENCE
- BUBBLE_CFD_STREAMING = FUTURE_PROVIDER
- FRAGMENTATION_SCOPE = EXPLORATORY_ONLY
- PATIENT_SPECIFIC_FRAGMENTATION = NOT_ESTABLISHED

V1 uses the existing hemisphere geometry, equal-spaced cubic material points, LPS real bond topology and three equal-budget analytic tangential loading patterns. Fragment identity comes from surviving-bond connected components. No new phase-field, PD2, CFD, GPU kernel or convergence campaign is authorized in this task.


### COHERENT_FAILURE_AND_FORCE_TIME_MAP_V1 completed — 2026-10-07

COHERENT_FAILURE_AND_FORCE_TIME_MAP_V1 = COMPLETE_STOPPED_FOR_HUMAN_REVIEW. The user explicitly authorized changing one anchor layer to three and four additional interface calibrations; exactly four were run, with the final 1957.41 Pa interface peak passing the 1700–3500 Pa gate. Frozen s00_bulk=0.15 (6780.9 Pa), s00_interface=0.140624859375. Liu 2024 SI Table S2 established G=940 Pa, K=33333.333 Pa; the interface inferred ranges belong to porcine-skin peeling, not the separate inner-aorta adhesion observation.

Completed 37 formal cases (27 coarse + 10 log-midpoint refinements), campaign wall 788.09 s. Final categories: 16 NO_FAILURE, 21 PULVERIZATION; whole detachment/resolved fragmentation/hybrid not observed as final states. Transient resolved/hybrid phases exist before pulverization. All 10 Pa cases had zero broken bonds. Equal budgets refer to initial/reference mean traction; actual loads follow current-position spatial decay. Source strength mapping does not validate fracture energy or post-break native contact dynamics. All results EXPLORATORY_ONLY.

Four PNG/PDF figures, three native-state MP4/GIF pairs, native checks, Chinese 12-question review and evidence archive are under clot_clearance/verification/lammps_coherent_failure_force_time_map_v1/ and clot_clearance/results/lammps_coherent_failure_force_time_map_v1/. No automatic further calibration, amplitude grid, VES, fatigue, erosion or bubble/CFD work. PRIMARY_FRAGMENTATION_SOLVER=LAMMPS_BPM_PERI; B1/B2 remain terminated, PD2 fracture remains stopped and phase-field remains archived optional reference. Await human review.
