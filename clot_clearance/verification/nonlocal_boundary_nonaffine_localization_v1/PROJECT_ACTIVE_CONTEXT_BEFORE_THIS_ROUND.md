# clot 项目当前上下文

更新：2026-10-05。后续 Codex **优先读取本文件**，再按当前任务打开具体源码和证据；不默认递归读取历史审计目录。

**当前 `NONLOCAL_BOUNDARY_FORMULATION_V1` 已完成，等待人工审核。** B1 的 D’Elia 风格 surrogate→BAQP 冻结体积力转换已在 N19 完整态法向／剪切仿射基准通过合力、参考力矩、12 个虚功模式、粒子重编号和扰动初值恢复检查。六组 B0/B1/B2 力学平衡均收敛；B1 平衡后的最大虚功误差分别为 1.17e-9／3.61e-9，几何 core F RMS 误差为 1.47e-9／2.72e-9（按施加梯度 0.003 归一化）。B0 虚功误差为 12.0%／17.3%；B2 一层 δ 的物理域体积约束误差为 1.94e-5／2.39e-5，未通过 1e-6 门槛，不能当作严格等效外部 Dirichlet 模型。

主基准：19³ 粒子，L=1.125 mm，h=62.5 µm，δ=251.25 µm，SVK μ=1000 Pa、λ=0。B1 依赖支持覆盖整个 N19 试件；d50/d90/d95/d99 为 1/2/3/4 h，是转换载荷分布，不是物理边界层厚度。N19 几何 core 不是完全闭合 bulk 区。固定 δ 与固定 δ/h 的 N19/N21 小比较通过，但不证明渐近相容性。状态为 `BOUNDARY_DATA_CONVERSION_STATUS=VERIFIED_FOR_INTACT_AFFINE_MECHANICS`、`SURFACE_OPERATOR_CORRECTION_STATUS=SEPARATE_ISSUE`。Scabbia 控制延期；perifit 仅兼容性讨论。

优先读[本轮 REVIEW](clot_clearance/verification/nonlocal_boundary_formulation_v1/REVIEW_NONLOCAL_BOUNDARY_FORMULATION_V1.md)、[边界规范与公式](clot_clearance/verification/nonlocal_boundary_formulation_v1/BOUNDARY_SPECIFICATION.md)、[四篇文献映射](clot_clearance/verification/nonlocal_boundary_formulation_v1/LITERATURE_TO_BAQP_MAPPING.md)、[五组图件索引](clot_clearance/results/nonlocal_boundary_formulation_v1/README.md)。新代码在 `pd_clot/pfpd_candidate/nonlocal_boundary/`，入口为 `scripts/nonlocal_boundary_audit.py`；11 项定向测试通过，56 个受保护输入／源码 SHA256 不变。未运行全量 pytest、activation、stability、phase、fracture、clot-wall 或 CFD。历史 C5/C7 状态保持原样。

**前一轮仿射系数参考审计已完成。** `PFPD_BAQP_AFFINE_COEFFICIENT_REFERENCE_AUDIT_V1` 复现旧 raw Z 误差 1.1362855097e-12 和 exact-P Z 误差 1.7684070888e-15；全部 257 个 Z 张量的 ΔP→ΔZ 因果闭合误差为 2.6930471364e-15（按固定 ||Zbulk||F 归一化），通过 1e-13 门槛。精确 F 高精度路径相对旧参考误差为 1.3529099860e-14，其中包含本构参考精度偏移；同应力高精度 bulk 参考下仅 1.0605625869e-15。误差主要来自 F 重构经本构灵敏度放大，Z 主导分量没有严重抵消，求和顺序影响约 2e-15。新解释为 `Z_MISMATCH_SOURCE=RECONSTRUCTED_STRESS_FLOATING_POINT_PROPAGATION`、`BULK_OPERATOR_CONSISTENCY=SUPPORTED`、`FINITE_DOMAIN_MULTI_HOP_BOUNDARY_INFLUENCE=SUPPORTED`；仅支持进入独立边界模型评审。六项定向测试通过，51 个历史／输入文件 SHA256 不变。优先读[最新系数参考报告](clot_clearance/verification/pfpd_baqp_affine_coefficient_reference_audit_v1/REVIEW_AFFINE_COEFFICIENT_REFERENCE_V1.md)。未运行新力学平衡、activation、fracture、stability 或 CFD。

**前一轮完整力模板审计的原始记录保持不变。** `PFPD_BAQP_FULLY_CLOSED_FORCE_STENCIL_AUDIT_V1` 固定 h=62.5 µm、δ=251.25 µm、ε=0.003，禁用相场。N27 中心实际三跳图及全部所需体积闭合；中心 |f|=1.94594e-22 N，按固定局部尺度 |Pzz|h² 归一化为 8.26546e-15，力门槛通过。但原始重构应力得到的 Z 与解析无限晶格参考值相对差 1.13629e-12，略超预声明 1e-12 系数门槛；因此完整验收 `BULK_AFFINE_FORCE_PATCH_TEST=NOT_ESTABLISHED`，`ROOT_CAUSE_PRIMARY=INCONCLUSIVE`，没有把这个超限改记为通过。选择性计算器与完整 N11/N13 对照通过，七项定向测试通过。仅作说明的恒定解析 P 诊断给出 Z 误差 1.77e-15，未替代原始验收。优先读[最新科学报告](clot_clearance/verification/pfpd_baqp_fully_closed_force_stencil_audit_v1/REVIEW_FULLY_CLOSED_FORCE_STENCIL_V1.md)。

此前[有限域边界审计](clot_clearance/verification/pfpd_baqp_finite_domain_boundary_consistency_v1/REVIEW_BAQP_FINITE_DOMAIN_BOUNDARY_CONSISTENCY_V1.md)及[定向理论来源](clot_clearance/verification/pfpd_baqp_finite_domain_boundary_consistency_v1/BOUNDARY_THEORY_NOTES.md)原样保留：L075 的 27 个自身邻域完整点并未隔离多跳模板，旧 patch FAIL 不改。本轮未运行边界修正、activation、稳定性或 CFD，未改生产公式或 C5/C7 接受态。

上一轮 Phase A 已确认 affine 无热点、真实平衡第一内层应变集中、axial-only 控制不能降低热点；没有选择 Phase B 修正。用户已明确取消旧 cleanup byte-hash 与论文逐式等价的科学开发硬门。[上一轮报告](clot_clearance/verification/pfpd_at1_activation_root_cause_and_fix_v1/REVIEW_ACTIVATION_ROOT_CAUSE_AND_FIX_V1.md)保留原结论。

历史恢复仍有 72 项旧字节依赖未恢复、用户没有备份；[恢复报告](clot_clearance/verification/cleanup_recovery_audit_v1/REVIEW_CLEANUP_RECOVERY.md)原样保留。这不再是本轮科学开发 gate，也不授权继续清理或改旧测试/manifest。

## CURRENT MODEL IDENTITY

**MODEL: `PFPD_AT1_SPECTRAL_CANDIDATE_V1`**。

当前是 BAQP mechanics + AT1 + Green 应变/SVK spectral split + PETSc RSLS、`h_kin(s)=1` early regime 的已验证候选组合。它不是 **exact PFPD 2026 reproduction**：论文打印的完整能量退化、拉伸驱动与演化公式之间仍有映射歧义。

工作区：`/home/lzy/projects/temp_storage/clot_clearance_github_sync_20260930`；Python 子项目为 `clot_clearance/`。2026-10-05 同步分支：`clot-clearance/pfpd-progress-20261005`，基于原归档提交 `f0e891686b9b107510ab3292aa27ce83638b3339`。当前源码、关键日志和科学数据已纳入本次同步；具体范围见 [CLOT_SYNC_20261005.md](CLOT_SYNC_20261005.md)。Git 未纳入的本地文件不等于可删除文件；不以 `??` 或 ignore 状态作为清理依据。

历史 C5/C7 相场数值试件：corrected boundary-node cube，L=0.25 mm，δ=0.25125 mm，μ=1000 Pa，ν=0，Gc_demo=1.26884390625×10⁻⁵ J/m²，quadratic kernel，c0=1/4，Yc=0.1010025 Pa，g(s)=(1−s)²，sc=0.95。材料参数是数值演示值。当前候选未接回真实 BraVa 流场、侵蚀或碎片清除。

## CURRENT ACTIVE SOURCE MODULES

以下路径均相对于 `clot_clearance/`：

| 入口 | 职责 |
|---|---|
| `pd_clot/pfpd_candidate/baqp.py`、`kernel.py`、`normalization.py` | BAQP 重构、核／能量权重及 c0 |
| `pd_clot/pfpd_candidate/svk_spectral_split.py` | SVK/Green 应变正负谱分裂 |
| `pd_clot/pfpd_candidate/at1_evolution.py`、`adr_adapter_at1.py` | AT1 能量、固定 x 更新及力学接口 |
| `pd_clot/pfpd_candidate/petsc_vi_candidate.py` | 已验证的 PETSc VI 候选 |
| `scripts/spatial_audit_v2.py`、`spectral_split_audit.py` | corrected C3/C5/C7 几何和谱分裂试件 |
| `scripts/c5_activation_audit.py`、`c7_early_branch.py` | 真实完整态平衡、激活 bracket 和已存分支 |
| `scripts/c7_reduced_vi.py`、`c7_spatial_evidence.py` | 自动 reduced pool、全键验收及空间指标 |
| `pd_clot/pfpd_candidate/boundary_force_diagnostics.py`、`scripts/finite_domain_boundary_audit.py`、`finite_domain_normalization.py` | s=0 有限域仿射力、精确分解及完整力模板诊断 |
| `scripts/finite_domain_author_oracle.py`、`finite_domain_author_oracle.jl`、`plot_finite_domain_boundary.py` | 固定作者 package 的 L075 对照及六组诊断图 |
| `pd_clot/pfpd_candidate/force_stencil_audit.py`、`scripts/fully_closed_force_stencil_audit.py`、`plot_closed_force_stencil.py` | 选择性完整态力、D1–D3 图闭合、无限晶格参考和四张科学诊断图；不接生产求解器 |

前一轮系数参考审计入口：`pd_clot/pfpd_candidate/affine_reference_audit.py`、`scripts/affine_coefficient_reference_audit.py`、`scripts/plot_affine_coefficient_reference.py`；不接生产求解器。新测试为 `tests/test_affine_coefficient_reference_audit.py`（6 passed），未重跑全量测试。

冻结物理来源：[AT1 配置](clot_clearance/verification/pfpd_at1_evolution_v1/provenance/PREDECLARED_CONFIG.json)、[C7 PHYSICS_FROZEN](clot_clearance/verification/pfpd_at1_c7_early_branch_v1/provenance/PHYSICS_FROZEN.json)。读取配置不代表可以修改它。

## CURRENT SCIENTIFIC STATUS

| 项目 | C5 | C7 |
|---|---:|---:|
| 粒子／自由机械 DOF／键数 | 125 / 225 / 5392 | 343 / 735 / 46059 |
| 完整态激活诊断应变 | 0.004497390736212 | 约 0.003137607455254 |
| 首激活位置 | 上下夹持关联 | 上下夹持关联 |
| 首激活键数 | 2 | 2 |
| 激活距离 | d/h≈0.5 | d/h≈0.5 |
| 正相场分支 | 已有审计窗口 | 连续 20 个接受态 |
| 最后 max(s) | 0.105282593434 | 0.044188650755 |
| 稳定性 | 15 个审计态稳定，窗口内未找到首次失稳 | `NOT_ASSESSED` |

C5/C7 activation location：**qualitatively consistent**；activation strain：**SHIFTING**，相对变化约 −30.2%；mesh objectivity：**NOT ESTABLISHED**。

尚未解决的历史科学问题包括 **activation mesh dependence**；当前先独立验证边界数据转换。共同应变、affine、support 和 axial-only 控制现已完成；热点的直接机制是近似有限的第一层位移增量在 h 缩小时放大局部应变，但有限域支持、BAQP 力分配及轴向边界施加方式尚未因果区分。论文 mapping=`AMBIGUOUS` 仍是模型身份限制，已不作为硬门。不能仅凭 fixed-x Hss 对角宣布缺少相场梯度项。

## CURRENT VALIDATED MILESTONES

1. corrected C3/C5/C7 几何、体积和表面权重守恒，以及 BAQP 固定状态力学核对。
2. SVK 谱分裂核对；s=0 恢复完整 SVK，允许混合拉压状态。
3. C5 激活及同载荷正根复核：旧 projected solver 的障碍是活动集停滞。
4. PETSc RSLS 自动恢复 7 个 C5 根；C7 reduced pool 的两点 C5 回归通过。
5. C5 物理 Hessian／Schur／临界锥窗口审计；C7 20 点早期驻定分支。
6. 有限域 BAQP affine 能量与离散一阶变分核对；四尺寸力分解、制造力及虚功比较完成。L075 作者 package 全场力相对差 5.86e-15；几何 core patch test FAIL，物理根因仍 INCONCLUSIVE。C5/C7 半域制造反力/R_cont=1.37586467/1.44235175，不能声称已恢复目标表面牵引。
7. N27 中心实际三跳图闭合且局部归一化力为 8.27e-15；N27 非中心一跳／两跳闭合观测点残力分别约 1.53e-2／3.53e-7。中心对称性会提前抵消，不能据 N11/N19 中心零力宣称最小充分闭合深度。原始 Z 系数门槛仍略超限，完整验收 NOT_ESTABLISHED；旧 C5/C7 结果不变。

优先读小型报告：[C7 REVIEW](clot_clearance/verification/pfpd_at1_c7_early_branch_v1/REVIEW_C7_EARLY_BRANCH_V1.md)、[最新停止报告](clot_clearance/verification/pfpd_at1_activation_mesh_dependence_audit_v1/REVIEW_ACTIVATION_MESH_DEPENDENCE_V1.md)、[论文映射](clot_clearance/verification/pfpd_at1_activation_mesh_dependence_audit_v1/references/PFPD_PHASE_FORMULATION_MAPPING.md)。C7 `SUMMARY.json` 很大，按字段读取，避免整份输出。

下一轮必须保留并按需读取：C5/C7 `intact/`、`activation/`、C5 `states/`、C7 `branch/`、空间指标、physics/config、输入和源码哈希。C5 `results/pfpd_at1_c5_stability_v1/cache/` 是测试与报告依赖，**不能当临时缓存删除**。

## CURRENT KNOWN FAILURES

清理前全量记录为 **412 passed / 3 known failed / 0 skipped**（当时 415 项）。上次恢复审计的全量结果为 **407 passed / 8 failed / 0 skipped**；清理遗留 5 个缺文件失败仍未修复。上一轮 Phase A 的 5 个直接相关测试通过；本轮有限域边界诊断的 7 个直接相关测试通过，没有重跑全量套件。旧字节失败不再阻塞科学开发。详见[历史测试对照](clot_clearance/verification/cleanup_recovery_audit_v1/TEST_COMPARISON.json)及[本轮定向测试](clot_clearance/verification/pfpd_baqp_finite_domain_boundary_consistency_v1/TARGETED_TESTS.json)。原三项均在 `tests/test_pfpd_at1_staggered.py`：

- `test_saved_coupled_mechanics_phase_and_energy_all_converged`
- `test_coupled_no_bond_heals_on_unload_and_reload`
- `test_extra_stagger_iterations_change_no_accepted_state_materially`

原始 staggered 路线未完成加载—卸载—重载及额外迭代验证。不能删除、xfail 或放宽这些测试来消除失败。PETSc 驻定成功不代表稳定性或空间收敛。

## CURRENT NEXT TASK

**等待人工审核 `REVIEW_NONLOCAL_BOUNDARY_FORMULATION_V1.md`。** 本轮只验证完整态仿射法向／剪切边界转换；没有验证一般非仿射边界、真实 clot-wall、fracture 或 CFD，也没有恢复 C5/C7 的 mesh objectivity。

唯一建议的下一开发任务为 `NONLOCAL_BOUNDARY_NONAFFINE_LOCALIZATION_V1`：在有已知局部解及相容体力／表面数据的非仿射弹性试件中，区分真正的边界转换载荷与 bulk 中的非局部－局部算子差异。必须先人工审核；本轮结束后不自动开始。

旧 raw Z 超限与原始 `N27_COEFFICIENT_GATE=FAIL` 保持原样；后续系数参考审计的因果解释另行保留。不要改历史 manifest／assertion，不继续 cleanup，不自动运行 activation、C7 stability、phase、权重优化、Ogden、adhesion、fracture、侵蚀、CFD 或 microbubble loading。

## DO NOT SEARCH BY DEFAULT

以下内容保留用于回溯；没有明确问题时不递归展开：

- `clot_clearance/verification/` 下旧直管、fragmentation、regularization、streaming_erosion、mechanical_equilibrium、damage_energy 历史实验。
- 早期 `pfpd_candidate_v1/`、`pfpd_at1_evolution_v1/`、C3 coupling/branch、spatial V1/V2 的完整历史树；只按需打开冻结配置或具体公式证据。
- `external/`、第三方源码、Julia/PETSc runtime/depot、`build_*`，以及任何嵌套 `.git/`。
- PRE/POST 大型 manifest、全量数组、历史 test 临时输出、`results/` 全树和动画逐帧目录。
- 本轮 `verification/repository_cleanup_v1/` 的逐文件清单；默认只读 REVIEW 和摘要，争议路径再查 CSV。

历史导航：[HISTORICAL_INDEX.md](HISTORICAL_INDEX.md)。完整公式／来源：[项目总览](clot_clearance/CLOT_PROJECT_WORKFLOW_STATUS_FORMULAS_ZH.md)。原始删除记录：[DELETE_MANIFEST.csv](clot_clearance/verification/repository_cleanup_v1/DELETE_MANIFEST.csv)，须结合恢复日志及最终 disposition 读取。历史 PRE/POST manifest 原样保留；若测试将其逐文件读取，其中的缓存与临时副本就是实际测试依赖，不能因其可再生而删除或改写旧哈希。
