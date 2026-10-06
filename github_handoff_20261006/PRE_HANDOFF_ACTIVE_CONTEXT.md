# clot 项目当前上下文

更新：2026-10-06。后续 Codex **优先读取本文件**，再按当前任务打开具体源码和证据；不默认递归读取历史审计目录。

**当前 `BAQP_NEAR_BOUNDARY_KINEMATIC_ENERGY_AUDIT_V1` 已完成，等待人工审核。** 固定 N19、h=62.5 µm、δ=251.25 µm，复用轴向/剪切二次制造目标态，未重跑任何平衡。连续完整球 M3=0；256 邻居完整晶格 M3 数值零，二次重构通过。截断均匀体积与实际梯形体积的最大 Fbar 偏差/0.004 分别为 2.852686e-2 / 2.975839e-2；独立 M3/K 恒等式、Ftilde 横向投影、轴向/剪切行映射、完整非线性 Ψ⁺ 重建全部通过。因果链为 **支持截断 → 非零奇矩 → Fbar 偏差 → 投影后的 Ftilde 偏差 → Ψ⁺ 偏差**。

`FBAR_BIAS_PRIMARY=TRUNCATED_ODD_MOMENT`，`FTILDE_BIAS_PRIMARY=PROJECTED_NODAL_RECONSTRUCTION_BIAS`，`PSI_PLUS_BIAS_PRIMARY=KINEMATIC_RECONSTRUCTION`，`SURFACE_OPERATOR_CORRECTION_READINESS=READY_FOR_BAQP_SPECIFIC_DESIGN`。体积权重次级贡献最大/RMS=5.829070e-3 / 2.749087e-3（/0.004），不能忽略，也不能把最大范数标量相加。K 连续参考 modifier=0.9963–7.5173，含幅值损失及各向异性；K 本身不是本轮根因。完整晶格 K 有 0.3724% 积分差，但二次重构仍精确到舍入误差。结论限于本轮两种二次场。

剪切 Fbar/Ftilde 恒等式最大残差/0.004=7.230457e-13 / 3.381112e-13；跨场实际行映射差=1.508566e-12。剪切 Ψ⁺ 实际最大偏差约 6.43%，矩预测后的因果残差约 1e-12（/各场共同 Ψ_ref）。剪切完整能量 U0/U1/U2=6.62365454e-12 / 6.63988403e-12 / 6.64456656e-12 J，运动学与权重/积分偏差分别为 −1.62294913e-14 / −4.68252790e-15 J；U1 仅作诊断。

优先读[本轮 REVIEW](clot_clearance/verification/baqp_near_boundary_kinematic_energy_audit_v1/REVIEW_BAQP_NEAR_BOUNDARY_KINEMATIC_ENERGY_V1.md)、[修正设计义务](clot_clearance/verification/baqp_near_boundary_kinematic_energy_audit_v1/BAQP_OPERATOR_CORRECTION_TARGETS.md)、[六组图](clot_clearance/results/baqp_near_boundary_kinematic_energy_audit_v1/README.md)、[数据字典](clot_clearance/verification/baqp_near_boundary_kinematic_energy_audit_v1/DATA_DICTIONARY.md)。10 项指定测试通过；425 个冻结旧文件 SHA256 不变。科学诊断 8 workers、17.397 s、较大内存高水位 1.070 GiB；GPU_NOT_USED_BY_DESIGN。未改物理、旧测试、旧 manifest 或 C5/C7 状态，未实施表面修正。Shojaei 仅作理论参考，BB/OSB 权重未移植。新设计须先处理单侧非负标量权重同时消除奇矩并保持 K 满秩的不可行约束。

**上一轮 `NONLOCAL_BOUNDARY_NONAFFINE_SHEAR_LOCALIZATION_V1` 已完成，以下保留其历史结论。** 固定 N19、γ=0.003、β_s=0.001 的二次剪切场；独立闭式 SVK、六面牵引、19 层完整 bulk oracle、仿射 B1 回归、逐节点／通道分解、N27 完全闭合中心以及两个非零扰动平衡控制均已完成。`SHEAR_BOUNDARY_LOCALIZATION=VERIFIED`，`SHEAR_BULK_MODEL_DIFFERENCE=SMALL`，`BOUNDARY_DATA_CONVERSION_STATUS=NONAFFINE_NORMAL_AND_SHEAR_LOCALIZATION_VERIFIED`。结论限于已验证的完整态制造基准。

仍采用 `f_BC=q_bulk V−f_FD`。剪切 bulk 模型差最大 7.520126e-5 N/m³，相对局部体力 4.229733e-8；逐节点分解误差/F_ref=3.618848e-17。N27 Z=0 的侧面 4h/8h/12h/13h 修正/F_ref 为 3.354435e-2/4.904946e-6/1.738699e-15/1.628360e-15；只有 13h 中心完全闭合，naive 中心仍有 1.964243e-2 F_ref 的内部力。localized d50/d90/d95/d99=62.5/125/187.5/250 µm，naive=62.5/125/187.5/312.5 µm；均为转换支持，不是物理 BC 厚度。

NL／经典平衡力残差分别为 8.709579e-12／8.767056e-12；相对离散目标 Fbar 的最大差/0.004 为 8.748958e-11／2.370509e-8。经典场位移偏差 2.072678e-14 m；NL 位移恢复误差为其 0.3401%，两解力响应与 e_model V 的相对闭合误差为 2.443965e-4。剪切的参考坐标连续体力矩有非零源项 `(0,4.271484e-14,0)` N·m，不能与当前位置零角动量混同；体力梯形积分误差另行保留。

**当前下一限制是近边界重构偏差。** 精确目标节点处 Fbar／Ftilde 最大解析偏差/0.004=0.0297584／0.0294285，`KINEMATIC_BOUNDARY_BIAS=MATERIAL`。Fbar 在节点、Ftilde 在参考键中点比较；Ψ⁺ 最大绝对偏差分别 2.582889e-4／2.578865e-4 Pa，约为统一解析参考的 6.43%／6.42%，加权 RMS 为 6.614026e-5／5.083853e-5 Pa。`PSI_PLUS_BOUNDARY_BIAS=REPORTED_NOT_ACCEPTANCE_RATED`；未修正算子、未作 fracture-safe 判定。

优先读[剪切 REVIEW](clot_clearance/verification/nonlocal_boundary_nonaffine_shear_localization_v1/REVIEW_NONAFFINE_SHEAR_LOCALIZATION_V1.md)、[数据字典](clot_clearance/verification/nonlocal_boundary_nonaffine_shear_localization_v1/DATA_DICTIONARY.md)、[六组图件索引](clot_clearance/results/nonlocal_boundary_nonaffine_shear_localization_v1/README.md)。11 项定向测试通过；336 个受保护旧文件哈希不变；三个允许扩展的诊断模块保留原版及前后哈希。旧 manifest 不改，但其历史源码 hash 不代表这三个扩展后的版本。旧 BAQP／材料／核／求解器、旧测试和科学结果未改。8 个 oracle workers、1 个串行平衡 worker；科学阶段累计 2082.045 s，峰值进程树 RSS 2.199 GiB，最小可用 RAM 53.37%。未运行 C5/C7、phase、fracture、clot 几何或 CFD。

**前一轮 `NONLOCAL_BOUNDARY_NONAFFINE_LOCALIZATION_V1`（轴向）已完成，以下保留该轮历史结论。** 固定 N19、ε=0.003、β=0.001 的轴向二次位移场；19 个 Z 层的完整邻域 bulk oracle、仿射回归、逐节点分解、N27 完全闭合中心及两个扰动初值平衡控制均已完成。`BOUNDARY_LOCALIZATION=VERIFIED`，`BULK_MODEL_DIFFERENCE=SMALL`，`BOUNDARY_DATA_CONVERSION_STATUS=NONAFFINE_LOCALIZATION_VERIFIED`。结论限于此完整态制造基准，不是一般物理牵引定律。

新修正为 `f_BC=q_bulk V−f_FD`，不再直接抵消非仿射总残差。bulk 模型差最大为 7.521586e-5 N/m³，相对局部体力 2.090312e-8；N27 中心修正/F_ref=3.196871e-15，直接负残差仍含 2.786077e-2 F_ref 的内部力。局部化修正 d50/d90/d95/d99=62.5/125/187.5/250 µm；naive 为 62.5/187.5/250/375 µm，均为转换支持而非物理 BC 厚度。N27 检查是 Z=0 的侧面 x 深度，只有 13h 中心实际闭合。

NL／经典体力平衡残差分别为 2.030006e-11／2.029948e-11，相对离散目标 Fbar 的最大梯度差为 1.628608e-10／1.113802e-8（按 0.004 归一化）。目标 Fbar 相对解析 F 本身存在 2.975839e-2 的截断重构偏差，不能与求解恢复误差混同。经典场不是精确 BAQP 解：两解力响应与 e_model V 的相对闭合误差为 2.462991e-4；NL 位移恢复误差为经典偏移的约 2.33%。全局力核算另保留 −1.5625e-14 N 的体力梯形积分误差。

优先读[本轮 REVIEW](clot_clearance/verification/nonlocal_boundary_nonaffine_localization_v1/REVIEW_NONAFFINE_BOUNDARY_LOCALIZATION_V1.md)、[数据字典](clot_clearance/verification/nonlocal_boundary_nonaffine_localization_v1/DATA_DICTIONARY.md)、[五组图件索引](clot_clearance/results/nonlocal_boundary_nonaffine_localization_v1/README.md)。9 项定向测试通过；251 个受保护旧文件哈希不变。只增加既有 nonlocal_boundary 命名空间的独立诊断模块，未改旧 BAQP／材料／求解器、旧结果、旧测试或旧 manifest。该轴向轮当时未运行非仿射剪切；本轮剪切结果见上文。没有启动 C5/C7、phase、fracture、曲面 clot-wall 或 CFD。

**上一轮 `NONLOCAL_BOUNDARY_FORMULATION_V1` 已完成，其记录原样保留。** B1 的 D’Elia 风格 surrogate→BAQP 冻结体积力转换已在 N19 完整态法向／剪切仿射基准通过合力、参考力矩、12 个虚功模式、粒子重编号和扰动初值恢复检查。六组 B0/B1/B2 力学平衡均收敛；B1 平衡后的最大虚功误差分别为 1.17e-9／3.61e-9，几何 core F RMS 误差为 1.47e-9／2.72e-9（按施加梯度 0.003 归一化）。B0 虚功误差为 12.0%／17.3%；B2 一层 δ 的物理域体积约束误差为 1.94e-5／2.39e-5，未通过 1e-6 门槛，不能当作严格等效外部 Dirichlet 模型。

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

本轮近边界因果审计入口：`nonlocal_boundary/kinematic_energy_diagnostics.py`、`scripts/near_boundary_kinematic_energy_audit.py`、`scripts/plot_near_boundary_kinematic_energy.py`、`scripts/report_near_boundary_kinematic_energy.py`；测试 `tests/test_near_boundary_kinematic_energy_audit_v1.py`（10 passed）。仅新增诊断，不接生产求解器。

剪切审计入口：`scripts/nonaffine_shear_boundary_audit.py`、`scripts/plot_nonaffine_shear_boundary.py`、`scripts/report_nonaffine_shear_boundary.py`；既有 `nonlocal_boundary` 扩展 `QuadraticShear`、bulk 通道和可选力矩源项，新增 `field_bias.py` 只作诊断。测试 `tests/test_nonaffine_shear_localization_v1.py`（11 passed）。

轴向非仿射边界诊断入口：`nonlocal_boundary/{nonaffine_manufactured,bulk_reference,localized_conversion,residual_decomposition}.py`，`scripts/nonaffine_boundary_audit.py`、`scripts/plot_nonaffine_boundary.py`、`scripts/report_nonaffine_boundary.py`；测试 `tests/test_nonaffine_boundary_localization_v1.py`（9 passed），未跑全量套件。

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

**等待人工审核 `REVIEW_BAQP_NEAR_BOUNDARY_KINEMATIC_ENERGY_V1.md`。** 轴向和剪切的边界转换已验证，近边界二次重构偏差的因果链已建立；未选择或验证修正方案。

唯一建议下一任务：`BAQP_SURFACE_OPERATOR_CORRECTION_DESIGN_V1`，设计 BAQP 专属多项式一致性与能量—力统一约束，先证明单侧支持下的可行性，不直接复制 BB/OSB/perifit 权重。本轮不执行该任务。

旧 raw Z 超限、旧 C5/C7 科学结论及历史 cleanup 缺失保持原样。不改旧 manifest/assertion，不继续 cleanup，不自动运行算子修正、权重优化、ghost、activation、C7 stability、phase、curved clot、fracture、adhesion、CFD 或 microbubble loading。

## DO NOT SEARCH BY DEFAULT

以下内容保留用于回溯；没有明确问题时不递归展开：

- `clot_clearance/verification/` 下旧直管、fragmentation、regularization、streaming_erosion、mechanical_equilibrium、damage_energy 历史实验。
- 早期 `pfpd_candidate_v1/`、`pfpd_at1_evolution_v1/`、C3 coupling/branch、spatial V1/V2 的完整历史树；只按需打开冻结配置或具体公式证据。
- `external/`、第三方源码、Julia/PETSc runtime/depot、`build_*`，以及任何嵌套 `.git/`。
- PRE/POST 大型 manifest、全量数组、历史 test 临时输出、`results/` 全树和动画逐帧目录。
- 本轮 `verification/repository_cleanup_v1/` 的逐文件清单；默认只读 REVIEW 和摘要，争议路径再查 CSV。

历史导航：[HISTORICAL_INDEX.md](HISTORICAL_INDEX.md)。完整公式／来源：[项目总览](clot_clearance/CLOT_PROJECT_WORKFLOW_STATUS_FORMULAS_ZH.md)。原始删除记录：[DELETE_MANIFEST.csv](clot_clearance/verification/repository_cleanup_v1/DELETE_MANIFEST.csv)，须结合恢复日志及最终 disposition 读取。历史 PRE/POST manifest 原样保留；若测试将其逐文件读取，其中的缓存与临时副本就是实际测试依赖，不能因其可再生而删除或改写旧哈希。
