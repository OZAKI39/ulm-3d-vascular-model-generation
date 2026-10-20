# clot 项目当前上下文

更新：2026-10-04。后续 Codex **优先读取本文件**，再按当前任务打开具体源码和证据；不默认递归读取历史审计目录。

**当前还有清理恢复阻塞，必须先处理。** `REPOSITORY_CLEANUP_AND_CONTEXT_REDUCTION_V1` 误将旧测试通过 PRE manifest 间接读取的缓存和临时输出当作可删除项。已精确恢复 5,919 项，尚缺 75 项；不能将清理标为完成。生产源码与当前 C5/C7 状态保持原哈希。先读[恢复清单](clot_clearance/verification/repository_cleanup_v1/RECOVERY_REQUIRED.csv)和[清理报告](clot_clearance/verification/repository_cleanup_v1/REVIEW_REPOSITORY_CLEANUP.md)，从备份恢复，不修改历史测试或其 manifest。

## CURRENT MODEL IDENTITY

**MODEL: `PFPD_AT1_SPECTRAL_CANDIDATE_V1`**。

当前是 BAQP mechanics + AT1 + Green 应变/SVK spectral split + PETSc RSLS、`h_kin(s)=1` early regime 的已验证候选组合。它不是 **exact PFPD 2026 reproduction**：论文打印的完整能量退化、拉伸驱动与演化公式之间仍有映射歧义。

工作区：`/home/lzy/projects/temp_storage/clot_clearance_github_sync_20260930`；Python 子项目为 `clot_clearance/`。HEAD：`f0e891686b9b107510ab3292aa27ce83638b3339`；branch：`clot-clearance/energy-regularization-20260930`。大量当前源码与结果未 tracked，不能将 HEAD 或 `??` 当作完整性／可删除性判据。

当前数值试件：corrected boundary-node cube，L=0.25 mm，δ=0.25125 mm，μ=1000 Pa，ν=0，Gc_demo=1.26884390625×10⁻⁵ J/m²，quadratic kernel，c0=1/4，Yc=0.1010025 Pa，g(s)=(1−s)²，sc=0.95。材料参数是数值演示值。当前候选未接回真实 BraVa 流场、侵蚀或碎片清除。

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

当前 blocker：**activation mesh dependence**。最新审计停在 `STOPPED_AT_PAPER_MAPPING_GATE`，paper mapping=`AMBIGUOUS`；common-strain、affine、boundary-layer 根因诊断尚未运行，根因=`INCONCLUSIVE`。fixed-x Hss 对角与本候选的局部 phase 能量相符，不能仅凭此宣布缺少论文要求的梯度项。

## CURRENT VALIDATED MILESTONES

1. corrected C3/C5/C7 几何、体积和表面权重守恒，以及 BAQP 固定状态力学核对。
2. SVK 谱分裂核对；s=0 恢复完整 SVK，允许混合拉压状态。
3. C5 激活及同载荷正根复核：旧 projected solver 的障碍是活动集停滞。
4. PETSc RSLS 自动恢复 7 个 C5 根；C7 reduced pool 的两点 C5 回归通过。
5. C5 物理 Hessian／Schur／临界锥窗口审计；C7 20 点早期驻定分支。

优先读小型报告：[C7 REVIEW](clot_clearance/verification/pfpd_at1_c7_early_branch_v1/REVIEW_C7_EARLY_BRANCH_V1.md)、[最新停止报告](clot_clearance/verification/pfpd_at1_activation_mesh_dependence_audit_v1/REVIEW_ACTIVATION_MESH_DEPENDENCE_V1.md)、[论文映射](clot_clearance/verification/pfpd_at1_activation_mesh_dependence_audit_v1/references/PFPD_PHASE_FORMULATION_MAPPING.md)。C7 `SUMMARY.json` 很大，按字段读取，避免整份输出。

下一轮必须保留并按需读取：C5/C7 `intact/`、`activation/`、C5 `states/`、C7 `branch/`、空间指标、physics/config、输入和源码哈希。C5 `results/pfpd_at1_c5_stability_v1/cache/` 是测试与报告依赖，**不能当临时缓存删除**。

## CURRENT KNOWN FAILURES

清理前科学基线为 **412 passed / 3 known failed / 0 skipped**，共 415 项。清理后曾有 7 个新增缺文件失败，恢复后仍需按[清理报告](clot_clearance/verification/repository_cleanup_v1/REVIEW_REPOSITORY_CLEANUP.md)的实际最终结果处理；不要将清理前基线当作当前测试结果。原三项均在 `tests/test_pfpd_at1_staggered.py`：

- `test_saved_coupled_mechanics_phase_and_energy_all_converged`
- `test_coupled_no_bond_heals_on_unload_and_reload`
- `test_extra_stagger_iterations_change_no_accepted_state_materially`

原始 staggered 路线未完成加载—卸载—重载及额外迭代验证。不能删除、xfail 或放宽这些测试来消除失败。PETSc 驻定成功不代表稳定性或空间收敛。

## CURRENT NEXT TASK

**common-strain / affine / boundary-layer root-cause audit**。

先明确记录如何处理既有论文映射门槛，再在冻结物理、边界和验收条件下开展共同应变平衡／仿射对照、d/h 与 d/δ profile、support 关联以及有条件的 clamp 控制。保持结论可为 `INCONCLUSIVE`，不预设原因。

先恢复清理误删造成的测试依赖并重新通过无新增失败门槛，再等待阶段审核；不自动启动上述诊断、C7 稳定性、C9、CFD、Ogden、真实 Gc 标定或 solver 开发。

## DO NOT SEARCH BY DEFAULT

以下内容保留用于回溯；没有明确问题时不递归展开：

- `clot_clearance/verification/` 下旧直管、fragmentation、regularization、streaming_erosion、mechanical_equilibrium、damage_energy 历史实验。
- 早期 `pfpd_candidate_v1/`、`pfpd_at1_evolution_v1/`、C3 coupling/branch、spatial V1/V2 的完整历史树；只按需打开冻结配置或具体公式证据。
- `external/`、第三方源码、Julia/PETSc runtime/depot、`build_*`，以及任何嵌套 `.git/`。
- PRE/POST 大型 manifest、全量数组、历史 test 临时输出、`results/` 全树和动画逐帧目录。
- 本轮 `verification/repository_cleanup_v1/` 的逐文件清单；默认只读 REVIEW 和摘要，争议路径再查 CSV。

历史导航：[HISTORICAL_INDEX.md](HISTORICAL_INDEX.md)。完整公式／来源：[项目总览](clot_clearance/CLOT_PROJECT_WORKFLOW_STATUS_FORMULAS_ZH.md)。原始删除记录：[DELETE_MANIFEST.csv](clot_clearance/verification/repository_cleanup_v1/DELETE_MANIFEST.csv)，须结合恢复日志及最终 disposition 读取。历史 PRE/POST manifest 原样保留；若测试将其逐文件读取，其中的缓存与临时副本就是实际测试依赖，不能因其可再生而删除或改写旧哈希。
