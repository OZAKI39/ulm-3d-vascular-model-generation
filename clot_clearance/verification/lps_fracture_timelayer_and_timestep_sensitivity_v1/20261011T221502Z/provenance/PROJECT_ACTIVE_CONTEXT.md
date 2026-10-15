## 2026-10-07 — CLOT_LITERATURE_MATERIAL_AND_CONTROLLED_FRACTURE_V1

本轮已完成接口修复、文献 Ogden 材料接入和 MAT_LM/MAT_W 原生 P2 完整态检查。独立 FEM AT2 候选在唯一三维厚度半试样上实际推进；最终 `FRACTURE_OUTCOME=NOT_OBSERVED_WITHIN_BUDGET`，已接受态物理分级 `DAMAGE_INITIATED`。接受 12 个状态，最大 Δ/H=0.197656，新高损伤扩展=0 mm；停止原因为 `WALL_BUDGET`。

先读 `clot_clearance/verification/clot_literature_material_and_controlled_fracture_v1/REVIEW_CLOT_LITERATURE_MATERIAL_AND_CONTROLLED_FRACTURE_V1.md`、`FINAL_STATUS.json`、`INTERFACE_AND_MODEL.md`。五图、MP4/GIF、原生体网格序列在同名 results 目录。资源审计记录了主机可用内存低于目标、首轮绘图预估不足及绘图修正，详见 `RESOURCE_POLICY_AUDIT.json`。所有失败/切步与源码版本保留；独立检查点复算结果 `True`。没有 commit、push、PR 或 cleanup。

参数来自牛全血体外文献；isochoric Ogden 与全能退化 AT2 是明确适配，Gc 未对新模型重标定，ell=0.5 mm 未标定。试样为资源受限的 native P1/P1（半球保持 P2）；没有网格客观性、体积锁死排除或定量实验强度认证。先于峰值的低载点不计为承载下降，初始相场尾部不计为新增损伤。

PD2 phase force 未实现，历史约11%差异与原FAIL不变；患者验证、adhesion、fragment transport、CFD/bubble/PZT、surrogate 未开展。未知用途容差仍未获批准。到本轮候选路径和交付为止，停止自动扩展物理或追加仿真。

---

## 2026-10-07 — CLOT_INTACT_PHYSICS_INTERFACE_AND_LOCALIZED_LOADING_V1

Completed callable candidate `clot_clearance/pd_clot/pfpd_candidate/intact_physics/` with Python + JSON CLI. Branch: clot-clearance/intact-physics-interface-localized-loading-v1; handoff HEAD remains c9815e8e2b9f2b66fb84719529a572c1d72cf233; no commit/push/PR. Prior uncommitted work retained.

Exactly five new native FEM equilibria: LP, LS, LM, W on P2 and LM_P3 on the same fine P1 mesh. All numerical quality gates passed. Saved GLOBAL M P2 and saved PD2 M read only; zero new PD2, mesh or random samples. LM_REUSE exercised Python no-solve reuse. Demo SVK mu=1000 Pa, lambda=0 remains explicitly uncalibrated; Psi+ diagnostic only. Prescribed wall uses all native boundary DOFs. Traction is reference dead GLOBAL/PATCH with c=(0,0,R), ell=R/2, not force-normalized.

LM same-geometry P2/P3 differences: apex vector 0.014433%, native u_RMS 0.008921%, original transverse field 0.012350%, common-volume ROI Psi+ P99 0.016410%. Numerical finite-check scope only; LP/LS/W remain single-discretization exploratory results. Reused 261405 fine-own-domain common points (ROI 254396); P99 finite estimate, no strict tail bound.

New tests: 26 passed. Five inspected English 600 dpi PNG/PDF figures and VTK/CSV display exports. Main review/config/results: `clot_clearance/verification/clot_intact_physics_interface_and_localized_loading_v1/`; figures: `clot_clearance/results/clot_intact_physics_interface_and_localized_loading_v1/`. Read REVIEW, FINAL_STATUS, INTERFACE_AND_MODEL and RESPONSE_SUMMARY first. `evidence.zip` is lightweight; native checkpoints/common arrays stay server-side with hashes. No prior state/manifests/tests altered; context history retained.

EXPLORATORY_MODULE_READINESS: ready for the declared interface and five-case scope. BIOLOGICAL_MATERIAL_CALIBRATION: NOT_ESTABLISHED. ADHESION/FRACTURE/DYNAMIC_RESPONSE/CFD_BUBBLE_PZT/SURROGATE: NOT_IMPLEMENTED_THIS_ROUND. HISTORICAL_V1_APPLICATION_ACCEPTANCE: UNCHANGED_NOT_ESTABLISHED. DOWNSTREAM_CONTROL_TOLERANCES: NEEDS_USE_REQUIREMENT. Historical 5/10/20% gates and PD2/FEM discrepancy untouched.

STOP: completed declared module, limited localized-loading cases and five figures. Do not automatically start discrepancy audits, material fitting or adhesion/phase-field development. Select next work from intended physical quantities and available data.

---

# 当前增补：固定态虚功定位与用途评估（2026-10-07）

`CURVED_CLOT_WALL_OPERATOR_LOCALIZATION_AND_USE_ASSESSMENT_V1` 已完成：**定位闭合 PASS；归属敏感性 CONVENTION_DEPENDENT；历史应用 NOT_ESTABLISHED；新用途 NEEDS_USE_REQUIREMENT。** 先读[本轮中文 REVIEW](clot_clearance/verification/curved_clot_wall_operator_localization_and_use_assessment_v1/REVIEW_CURVED_CLOT_WALL_OPERATOR_LOCALIZATION_AND_USE_ASSESSMENT_V1.md)、[用途表](clot_clearance/verification/curved_clot_wall_operator_localization_and_use_assessment_v1/USE_AND_ERROR_LEDGER.csv)与[四张图](clot_clearance/results/curved_clot_wall_operator_localization_and_use_assessment_v1/README.md)。下方旧上下文正文原样保留。

- 使用原 fine P2 顶点 ustar 与保存 v3；421739键/14422单元独立闭合，Delta_operator=8.902954998486e-13 J，最大闭合误差/W_ref=2.047e-14。零新平衡、零新样本。
- 两归属均为ring净差最大；A全局正负取消49.0224%，ring占2.5550%体积、23.4479%节点绝对份额；remainder净差A正/B负。不能转成11.2%位移误差因果分摊或唯一边界根因。
- 当前M的native FEM位移/共同体积能量有限参考保留；PD2探索须带差异，原10%替代FAIL与节点P99参考失败保留。未知实际用途容差留空，无新15/20%许可。
- 新6项测试通过，四图600dpi PNG/PDF且实际QA。短FP64 GPU contraction/scatter一致但较CPU慢，正式路径CPU；原源码/测试/manifest/结果不改，无commit/push/PR/cleanup。
- **本轮停止。下一步只建议按输出保留方法与差异标签，准备后续完整态模块接口需求；没有明确实现错误或用途精度要求时，不继续追求方法完全一致。** 本轮没有启动后续模块。

---

# 当前增补：曲面应用差异审计（2026-10-07）

`CURVED_CLOT_WALL_APPLICATION_DISCREPANCY_AUDIT_V1`：**AUDIT_COMPLETION=PASS；CURVED_APPLICATION_ACCEPTANCE=NOT_ESTABLISHED。** 先读[本轮中文 REVIEW](clot_clearance/verification/curved_clot_wall_application_discrepancy_audit_v1/REVIEW_CURVED_CLOT_WALL_APPLICATION_DISCREPANCY_AUDIT_V1.md)、[逐输出状态](clot_clearance/verification/curved_clot_wall_application_discrepancy_audit_v1/FINAL_STATUS.json)及[五张新图](clot_clearance/results/curved_clot_wall_application_discrepancy_audit_v1/README.md)。下面旧上下文原文保留，属于截至 2026-10-06 的历史记录。

- 原 127 项指标复算一致；完整 native P2 DOFs 已验证恢复。修复仅限新审计适配器的实际单元包含检查：旧 DOLFINx SI 点定位容差会允许约 1 um 外推；原平衡状态与节点恢复不改。
- 旧 16.0439% ROI 节点能量 P99 变化仍有效。共同体积 ROI main/fine P99 变化 0.302573%，同细几何 P2/P3 为 0.076628%；原节点 P2/P3 仍变化 5.546193% >5%，原节点能量参考未建立。
- 原横向位移差 11.209514%，纯点定位修正后 11.209516%，对 P3 为 11.197767%，均 >10%。同几何 P2/P3 横向位移变化 0.015488%。只给 M / 本轮输出的有限位移与共同体积参考支持，不能推广为全应用验证。
- 三方向虚功闭合；外功差在舍入范围。v3 的主要差异为内部离散功，伴随 P2→P1 transfer 项，不能等同纯边界误差。
- 一个同几何 P3 科学案例收敛；两次 SNES 调用（首次 monitor 第0步失败，从同一初态恢复），零新增 PD2 平衡。新小 pytest 7 pass，旧文档 hash test 仍原样 fail；被短路的实际力/能量断言已独立补测。
- 隔离 CuPy FP64 实测通过，但真实端到端收益约1.6%，正式后处理选CPU。五图均600dpi PNG/PDF并实际逐张检查。
- 当前本地分支 `clot-clearance/curved-clot-wall-discrepancy-audit-v1`，HEAD仍为handoff c9815e8e；保留上一轮和本轮未提交内容，没有commit/push/PR/cleanup。

**本轮已停止。** 唯一建议下一步是固定状态的 v3 内部虚功空间定位审计；不自动开展任何新计算或新物理。

---

# clot 项目当前上下文

更新：2026-10-06，`CURVED_CLOT_WALL_INTACT_MECHANICS_V1` 完成。新会话先读本文与[本轮中文 REVIEW](clot_clearance/verification/curved_clot_wall_intact_mechanics_v1/REVIEW_CURVED_CLOT_WALL_INTACT_MECHANICS_V1.md)，再按需读历史报告；不递归展开历史树。

**当前科学终点：曲面完整态 candidate 的预算内开发、计算与六图已完成，但 `CURVED_CLOT_WALL_INTACT_MECHANICS=NOT_ESTABLISHED`，不可宣布曲面应用已验证。** 旧平面基准 `HIGHER_ORDER_BOUNDARY_MECHANICS=VERIFIED_FOR_INTACT_BENCHMARKS` 保持历史含义。

## 本轮曲面血栓结果与停止点

分支 `clot-clearance/curved-clot-wall-intact-20261006`；从 handoff commit `c9815e8e2b9f2b66fb84719529a572c1d72cf233` 开始。新增 candidate 包 `clot_clearance/pd_clot/pfpd_candidate/curved_clot_wall/`；数据/配置/日志位于 `clot_clearance/verification/curved_clot_wall_intact_mechanics_v1/`；[六图、PDF 与 VTU/VTP](clot_clearance/results/curved_clot_wall_intact_mechanics_v1/README.md)。没有 commit/push/PR、cleanup 或生产替换。

R=0.5625 mm 的半球，底部圆盘三分量零位移；Gmsh P1 四面体顶点为粒子，正体积按单元四等分。主/细 h=62.5/50 µm，固定 δ=251.25 µm，SVK μ=1000 Pa、λ=0。粒子/键数为 1771/132594 和 3086/421739。全节点 rank 9，最大条件数 177.594/186.741，P2 梯度误差约 2e-14。几何、曲面名义死载功装配、独立理想半球积分、零载与非求解整体平移检查通过。无完整力依赖闭合节点。

| 已完成检查 | 结果 |
|---|---|
| 一个扰动起步二次剪切制造恢复 | PASS；残差 3.078e-13，Fbar/Ftilde 与 Ψ⁺ 恢复通过，仅 MANUFACTURED_ONLY |
| P/S/M_HALF/M 与 finer M | 全部真实平衡通过；finer M 残差 6.971e-11，外载前后不变，node/bond detF>0 |
| PD2 主/细 M | PASS；apex 变化 1.407%，u RMS 2.187%，ROI node/bond Ψ⁺ P99 0.492%/0.0228% |
| 原生 P2 DOLFINx 的两个 M 参考 | 均收敛；apex 变化 0.0200%、u RMS 0.321%，但 ROI node Ψ⁺ P99 变化 16.044% > 5% |
| 独立应用参考 | NOT_ESTABLISHED；没有追加第三套 FE 或调整门槛 |
| 暂以细 FE 比较细 PD2 | u RMS 差 10.054%、横向 u 剖面差 11.210%、ROI node Ψ⁺ P99 差 23.444%，数值门槛亦超限；参考未可信，不能解释成已知真实误差 |

本轮复用原 PD2 能量/力、投影、TRF/LSMR 与 callback。42 项直接接受输入最终哈希不变；旧源码、tests、manifests 与 scientific results 未修改。新 8 项实时测试通过；三套旧定向 tests + 新 tests 合计 43 passed / 1 failed。唯一失败是旧保护清单要求的文档哈希已在 handoff HEAD 不同，证据见 `PREEXISTING_HISTORICAL_HASH_DIFFERENCES.json`；没有改旧断言使其通过。未运行全量 pytest。

最多三路 PD2 CPU 求解，BLAS/OMP 单线程；FEM 与 finer M 按内存余量并行。主应用约 6–8 分钟/例，细 M 约 12.47 分钟；两套 FE 连输出约 18.6/36.6 秒。实际设备为 RTX 4060 Laptop 8 GiB；没有已验证 GPU 力学 backend，使用 CPU / Matplotlib Agg。完整成本与采样局限见本轮 `PERFORMANCE_SUMMARY.json`。

**到此停止，等待人工审阅。** 不自动继续 FE 加密、方法研究、adhesion、phase、fracture、CFD、microbubble、PZT 或 surrogate。材料未标定；本轮只形成未获完整应用验证的 intact physics-oracle 候选输入/输出接口。

## 项目方向与身份

最终目标是 PZT → bubble → morphology → streaming/shear → clot response，以可验证模块和紧凑响应关系支持未来控制/代理模型。当前只验证血栓固体力学基础；没有已训练代理或端到端 PZT—血栓结果。

当前主线是 **完整 P2/RK 梯度 + 完整态变分 BAQP + 已验证基准上的非局部边界转换**。历史 `PFPD_AT1_SPECTRAL_CANDIDATE_V1` 的 C1 BAQP、AT1、SVK 谱分裂和 PETSc RSLS 保留；PD2 尚无 phase force。两者不能合称已经验证的完整 PFPD 2026 reproduction。

Git 根：`/home/lzy/projects/temp_storage/clot_clearance_github_sync_20260930`；Python 子项目 `clot_clearance/`。交接来源分支 `clot-clearance/pfpd-progress-20261005`、HEAD `695e75166aa864809460ea9a1845fe72abce229c`；交接分支 `clot-clearance/pd2-boundary-handoff-20261006`。当前提交 SHA 以 Git 为准。

## 最近三个已完成里程碑

| 阶段 | 当前可接受结论 | 首读报告 |
|---|---|---|
| RK oracle | 全节点九项基满秩，最大条件数 218.069；九项再现、轴向/剪切/XZ、Julia 代表节点对照通过；省略交叉项的负对照失败 | [RK REVIEW](clot_clearance/verification/higher_order_baqp_rk_oracle_v1/REVIEW_HIGHER_ORDER_BAQP_RK_ORACLE_V1.md) |
| 变分能量—力 | 同一离散能量的负一阶变分；C1 回归、客观性、全局力/当前位置力矩、方向差分和选定 DOF 检验通过 | [能量—力 REVIEW](clot_clearance/verification/higher_order_baqp_energy_force_candidate_v1/REVIEW_HIGHER_ORDER_BAQP_ENERGY_FORCE_V1.md) |
| 边界回接 | 独立 PD2 bulk oracle、闭合图、合力/参考力矩/虚功及四个扰动起步制造恢复通过；唯一经典剪切控制与 bulk 模型差一致 | [边界 REVIEW](clot_clearance/verification/higher_order_baqp_boundary_reintegration_v1/REVIEW_HIGHER_ORDER_BAQP_BOUNDARY_REINTEGRATION_V1.md) |

共同主基准 N19=6859 节点、690020 条无向键，h=62.5 µm、δ=251.25 µm、L=1.125 mm，二次径向核与物理梯形体积，SVK μ=1000 Pa、λ=0。演示参数未生物标定。

当前边界定义 `f_BC=q_bulk*V-f_FD` 保留非仿射体域物理差异。D0/D1/D2/D3 累计为 1/257/1797/5821；最大坐标偏移 0/4/8/12h。N19 不闭合，N25 仍读取半体积边界，N27 的 13h 中心充分闭合。D1=257 包含观察点；与旧 256 直接邻点计数不矛盾。

四个恢复态机械残差分别为 5.191387e−10、8.658278e−11、1.150552e−12、5.778522e−11，均通过 1e−8；梯度和 Ψ⁺ 分别通过 1e−5、1e−4 门槛。经典剪切机械残差 3.780054e−12，目标位移偏移 1.693642e−14 m；相对模型载荷响应闭合误差约 1.40%，不能赋予微小偏移全部数字物理意义。

边界载荷 d50/d90/d95/d99：PD2 0/2/2/4h，C1 1/2/3/4h，是 **转换支持，不是物理 BC 厚度**。仿射实际外载与保存 fBC 的舍入量级 −qV 差异有独立核算。C1–PD2 场图比较历史 C1 目标态与 PD2 求解态。原求解器、生产 BAQP、已接受 PD2 力、旧结果与旧 manifest 均未被边界回接改写。

三个原定向套件各 12 项；历史 handoff 记录为 **36 passed in 4.42 s**。当前曲面轮次的实跑结果及 handoff 文档哈希差异见上文，不能将旧记录当作今天的测试结果。只用同步文件的 20 项可移植子集见 [环境与命令](HANDOFF_ENVIRONMENT.md)。

## 直接相关代码

以下相对 `clot_clearance/`：

- `pd_clot/pfpd_candidate/higher_order_rk/`：完整 P2 基/矩阵/梯度、冻结参考与变分能量—力伴随。
- `pd_clot/pfpd_candidate/nonlocal_boundary/`：几何、制造场、转换、PD2/C1 backend、bulk oracle、原求解器 adapter 与诊断。
- `pd_clot/pfpd_candidate/{baqp,kernel,energy,svk_spectral_split}.py`：保留的基础模块；`pd_clot/regularization/equilibrium.py`：复用的 TRF/LSMR。
- `verification/higher_order_baqp_*/run_audit.py`、`tests/`、`plot_figures.py`：完整复现入口；部分会写结果或运行长求解，不为阅读项目而自动执行。
- `scripts/` 下有限域、完整模板、非仿射与近边界审计入口；全部相关源和 tests 在分支中。

## 保留的历史结论

C1 截断支持引起非零三阶奇矩，进而导致二次 Fbar、投影 Ftilde 和 Ψ⁺ 偏差。既有轴向/剪切 max Fbar 偏差 /0.004≈0.0297584，剪切 Ψ⁺ 偏差约 6.43%。完整 P2 已解决所测场的重构问题；不再停留于“尚未选择修正”的旧状态。详见 [C1 因果 REVIEW](clot_clearance/verification/baqp_near_boundary_kinematic_energy_audit_v1/REVIEW_BAQP_NEAR_BOUNDARY_KINEMATIC_ENERGY_V1.md)。

旧有限域几何 core patch FAIL、完整模板原 raw Z 门槛超限保持原判；后续系数参考审计解释浮点传播，没有改写这些历史结论。[历史导航](HISTORICAL_INDEX.md) 按轮次保留证据；[交接前上下文原文](github_handoff_20261006/PRE_HANDOFF_ACTIVE_CONTEXT.md) 仅供溯源，其中“当前/下一任务”是旧时间点用语。

| 历史相场状态 | C5 | C7 |
|---|---:|---:|
| 粒子 / 自由机械 DOF / 键 | 125 / 225 / 5392 | 343 / 735 / 46059 |
| 激活应变 | 0.004497390736212 | 约 0.003137607455254 |
| 首激活 | 夹持关联、2 键、d/h≈0.5 | 夹持关联、2 键、d/h≈0.5 |
| 最大已存 s | 0.105282593434 | 0.044188650755 |
| 已有验证 | 15 态稳定窗口；未找到首次失稳 | 20 个驻定接受态；稳定性 NOT_ASSESSED |

激活应变偏移约 −30.2%；mesh objectivity NOT_ESTABLISHED。Phase A 确认 affine 无热点、真实平衡第一层应变集中；axial-only 未解决。不能把新的完整态 PD2 结果当成修复了旧相场激活网格依赖。

## 尚未验证及禁止自动重启的范围

本轮两套指定半球的局部满秩、条件数、体积与边界检查已通过；对其它曲面/点云的泛化仍未建立。整体曲面应用因独立参考与暂时跨模型超限而未获验证。PD2 二次再现和变分一致不证明积分精确、稳定性、mesh objectivity、phase、fracture、adhesion、Ogden、CFD 或微泡加载。后续阶段未启动。

历史 cleanup 最终仍有 72 项旧字节依赖未恢复，用户已说明无备份。清理前全量记录 412 passed / 3 failed；恢复审计 407 passed / 8 failed，不是今天的全量测试结果。原 3 项 staggered 失败见 `tests/test_pfpd_at1_staggered.py`。用户已取消历史 byte-hash/论文 exact-equivalence 的科学开发硬门，但未允许修改旧断言、旧 manifest 或继续清理。保留[恢复报告](clot_clearance/verification/cleanup_recovery_audit_v1/REVIEW_CLEANUP_RECOVERY.md)。

不依据 Git `??`、ignore 或 cache 名称删除文件。C5 stability cache、C7 branch 和全部 accepted states 仍是科学证据。旧 scientific manifests 描述完整服务器产物，不是本轮 GitHub 文件清单。

## 交接范围与停止点

这是选择性 ADD/UPDATE/SYNC。大数组、重复图/论文/环境保留本地；精确路径/大小/哈希与替代摘要见 [LOCAL_ARTIFACTS_NOT_SYNCED.md](LOCAL_ARTIFACTS_NOT_SYNCED.md)，提交范围见 [HANDOFF_SYNC_MANIFEST_20261006.md](HANDOFF_SYNC_MANIFEST_20261006.md)。没有删除历史 tracked 大文件、重写 Git 历史或创建 PR。

以上交接范围为历史记录。本轮曲面开发交付后停止；下一会话得到新指令前，不扩展曲面求解、相场、C5/C7、断裂、CFD、微泡或代理建模。


---

## 2026-10-07 — POLICY_AND_FAST_MATERIAL_SCREENING_V1 (active policy addition)

CONTROL_ORIENTED_SUFFICIENCY_V2 = ACTIVE
DEFAULT_SOLVER = FEM
DEFAULT_MODEL_LEVEL = LEVEL_1_FAST_FEM
PD2_ROLE = SPOT_CHECK_REFERENCE
DEFAULT_SCREENING_TOLERANCE = 30_PERCENT_WITH_DIRECTION_MAGNITUDE_RANKING_GATES
DEFAULT_TIER_A_RELATIVE_UNCERTAINTY = 0.30
FRACTURE_ESCALATION = ONLY_FOR_SELECTED_HIGH_VALUE_OR_NEAR_THRESHOLD_CASES
DYNAMIC_ESCALATION = ONLY_WHEN_T_LOAD_NOT_LARGE_RELATIVE_TO_T_CLOT_OR_IMPACT_IS_EXPLICITLY_MODELED

Highest project policies are the exact user-supplied files in docs/policies/.
30% is project screening policy, not experimental uncertainty. Numerical gates,
application validation, robust rankings and threshold assessments are separate.
This round authorizes six native P2 intact cases only; no fracture, PD2, P3,
dynamic, CFD, surrogate, new mesh or changed historical acceptance.
User explicitly accepted isochoric Ogden adaptation; not exact paper reproduction.
Read clot_clearance/verification/policy_and_fast_material_screening_v1/ for this round.


### 2026-10-07 — Round A completed

POLICY_AND_FAST_MATERIAL_SCREENING_V1 is complete for numerical Level-1 relative
screening: six P2 states passed, 30 new targeted tests and 22 prior read-only tests
passed, four 600 dpi PNG/PDF figures inspected. Risk remains RELATIVE_ONLY; no
absolute threshold or patient validation. Nominal order WALL > G_S > L_M > L_S >
L_P > G_P; L_M/L_S uncertain for all three metrics, WALL/G_S displacement uncertain.
Cases cost 34.40–50.29 wall seconds and 724.56–731.35 MiB peak RSS, CPU only.
Read the round REVIEW and FINAL_STATUS in the directory above. No fracture or PD2
executed this round; WALL/G_S are manual next-round suggestions only. Stop here.


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


### LAMMPS_FRAGMENTATION_PARADIGM_SCREENING_V1 completed

Three equal-budget hemisphere LPS cases completed with true broken bonds, independent components and continued fragment motion. N=2971, initial bonds=140555; s00=.20 and amplitude unchanged. Endpoint FI≈.864 saturates near fixed-base limit, RANKING_UNCERTAIN; predominantly node-scale debris, not clinical fragment sizes. B1/B2 preserved and archived; PD2 fracture stopped. Evidence: clot_clearance/verification/lammps_fragmentation_paradigm_screening_v1/. Await human review; no further simulations.


### COHERENT_FAILURE_AND_FORCE_TIME_MAP_V1

COHERENT_FAILURE_AND_FORCE_TIME_MAP_V1 = ACTIVE
PRIMARY_FRAGMENTATION_SOLVER = LAMMPS_BPM_PERI
PRIMARY_FRAGMENTATION_USE = QUALITATIVE_PATTERN_SCREENING
PHASE_FIELD_FRACTURE_ROLE = ARCHIVED_OPTIONAL_REFERENCE_ONLY
B1_B2_ROUTE = TERMINATED_BY_PROJECT_POLICY_CHANGE
PD2_FRACTURE_DEVELOPMENT = STOPPED
PD2_ROLE = OPTIONAL_INTACT_SPOT_CHECK_REFERENCE
PRODUCTION_BULK_MODEL = LPS
VES = DEFERRED
VISCOELASTIC_RATE_DEPENDENCE = NOT_MODELED_V1
FATIGUE_DAMAGE = NOT_MODELED
SURFACE_EROSION = NOT_MODELED
POROELASTICITY = NOT_MODELED
CHEMICAL_THROMBOLYSIS = NOT_MODELED
SCOPE = EXPLORATORY_ONLY

Liu 2024 only for bulk/interface failure; bounded strength calibration and <=39 amplitude-duration-pattern cases. Preserve prior results. Stop for human review after delivery.


### COHERENT_FAILURE_AND_FORCE_TIME_MAP_V1 completed — 2026-10-07

COHERENT_FAILURE_AND_FORCE_TIME_MAP_V1 = COMPLETE_STOPPED_FOR_HUMAN_REVIEW. The user explicitly authorized changing one anchor layer to three and four additional interface calibrations; exactly four were run, with the final 1957.41 Pa interface peak passing the 1700–3500 Pa gate. Frozen s00_bulk=0.15 (6780.9 Pa), s00_interface=0.140624859375. Liu 2024 SI Table S2 established G=940 Pa, K=33333.333 Pa; the interface inferred ranges belong to porcine-skin peeling, not the separate inner-aorta adhesion observation.

Completed 37 formal cases (27 coarse + 10 log-midpoint refinements), campaign wall 788.09 s. Final categories: 16 NO_FAILURE, 21 PULVERIZATION; whole detachment/resolved fragmentation/hybrid not observed as final states. Transient resolved/hybrid phases exist before pulverization. All 10 Pa cases had zero broken bonds. Equal budgets refer to initial/reference mean traction; actual loads follow current-position spatial decay. Source strength mapping does not validate fracture energy or post-break native contact dynamics. All results EXPLORATORY_ONLY.

Four PNG/PDF figures, three native-state MP4/GIF pairs, native checks, Chinese 12-question review and evidence archive are under clot_clearance/verification/lammps_coherent_failure_force_time_map_v1/ and clot_clearance/results/lammps_coherent_failure_force_time_map_v1/. No automatic further calibration, amplitude grid, VES, fatigue, erosion or bubble/CFD work. PRIMARY_FRAGMENTATION_SOLVER=LAMMPS_BPM_PERI; B1/B2 remain terminated, PD2 fracture remains stopped and phase-field remains archived optional reference. Await human review.
