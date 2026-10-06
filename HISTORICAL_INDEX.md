# clot 历史证据索引

更新：2026-10-06，交接导航补充。当前先读 [PROJECT_ACTIVE_CONTEXT.md](PROJECT_ACTIVE_CONTEXT.md) 与 [HANDOFF_CONTEXT_20261006.md](HANDOFF_CONTEXT_20261006.md)。最新终点为 PD2 完整态边界回接 VERIFIED；下方旧索引保留其当时状态，旧“当前/下一步”不代表今天的开发指令。

## 新增的已完成主线

路径均在 `clot_clearance/verification/`；各目录的 REVIEW、FINAL_STATUS/STATUS、配置、定向测试和紧凑数值在本分支中。

| 目录 | 结论 |
|---|---|
| `pfpd_baqp_fully_closed_force_stencil_audit_v1` | N27 中心力通过；原 raw Z 门槛略超，原 NOT_ESTABLISHED 保留 |
| `pfpd_baqp_affine_coefficient_reference_audit_v1` | Z 差异由重构应力浮点传播解释；支持 bulk 一致性 |
| `nonlocal_boundary_formulation_v1` | 完整态仿射 B1 边界转换通过 |
| `nonlocal_boundary_nonaffine_localization_v1` | 轴向非仿射 bulk/边界分离通过 |
| `nonlocal_boundary_nonaffine_shear_localization_v1` | 剪切非仿射局部化通过；C1 重构偏差被量化 |
| `baqp_near_boundary_kinematic_energy_audit_v1` | 奇矩→节点梯度→键梯度→正能量因果链验证 |
| `higher_order_baqp_rk_oracle_v1` | 完整 P2/RK 运动学 VERIFIED；Julia 代表节点对照通过 |
| `higher_order_baqp_energy_force_candidate_v1` | 完整态保守变分力学 VERIFIED |
| `higher_order_baqp_boundary_reintegration_v1` | 四个制造恢复+经典剪切；VERIFIED_FOR_INTACT_BENCHMARKS |

仅建议下一轮 `CURVED_CLOT_WALL_INTACT_MECHANICS_V1`；未启动。

## 交接前索引（历史原文）

# clot 历史证据索引

更新：2026-10-05。默认先读 [PROJECT_ACTIVE_CONTEXT.md](PROJECT_ACTIVE_CONTEXT.md)。下表用于按问题定位历史；**KEEP FOR REFERENCE** 不代表当前主线仍使用该模型，也不代表目录中所有科研结论成功。路径均相对于 Git 根的 `clot_clearance/`。

当前先读[有限域边界审计报告](clot_clearance/verification/pfpd_baqp_finite_domain_boundary_consistency_v1/REVIEW_BAQP_FINITE_DOMAIN_BOUNDARY_CONSISTENCY_V1.md)：affine 全局能量正确，L075 几何 core patch FAIL，作者实现和能量差分复现残差；完整力模板未隔离，ROOT_CAUSE_PRIMARY=INCONCLUSIVE，未运行 BC control，等待审核。上一轮 Phase A 和 C5/C7 接受态原样保留。用户明确不再将旧 byte-hash failures 或 paper exact equivalence 作为科学开发硬门。历史 72 项恢复问题仍保留在[恢复报告](clot_clearance/verification/cleanup_recovery_audit_v1/REVIEW_CLEANUP_RECOVERY.md)，未修改其结论。

| 目录／入口 | 一行用途 | 最终状态或结论 | 后续何时需要 |
|---|---|---|---|
| `REPORT_ZH.md`、`results/straight_pipe_demo/` | 初始流场—NOSB-PD—损伤输出链 | 轻度损伤、未断裂；软件基准完成 | 流体接口、旧模型回归 |
| `verification/fragmentation_extension/`、`runs/fragmentation_pilot_003/` | 强制碎裂、低秩分量、脱离／输运／清除账本 | 演示链跑通，参数未物理标定 | 图拓扑、碎片账本和旧动画溯源 |
| `verification/regularization/`、`results/streaming_regularized_demo/` | coupon 能量尺度与 18 组网格／周期比较 | `NOT CONVERGED`；无新断键或脱落 | 演示 Gc 定义、旧网格比较局限 |
| `verification/streaming_erosion_v1/` | 从循环幅值改为独立慢曝光时间的经验侵蚀 | 初始机械稳态不足，发现损伤增能 | 侵蚀背景和已知失败证据 |
| `verification/mechanical_equilibrium_v1/` | 初期平衡诊断与独立小问题参考 | 旧报告被 ADR completion 补充；保留历史 | 审查平衡门槛演变 |
| `verification/mechanical_equilibrium_v1/adr_completion/` | ADR 方法／阻尼保护及力学验收完成 | 41 平衡态、35 跳变；损伤耗散仍未建立 | ADR 算法来源、旧加载/平衡回归 |
| `verification/damage_energy_audit_v1/` | 固定 x 的损伤储能分解和重构恒等式 | 35 次增能重现；重构项为主要机制 | 理解为何转向完整能量候选 |
| `verification/pfpd_candidate_v1/` | BAQP 与作者 Julia 固定版本对照 | 固定状态力学验证；原文 phase 演化仍歧义 | BAQP/归一化/力式原始来源 |
| `verification/pfpd_at1_evolution_v1/`、对应 `results/` | 用户选定 AT1 一致候选与原始 staggered 路线 | fixed-x 通过；联合加载路线遗留 3 项失败 | 冻结物理、原失败回归；按文件读取 |
| `verification/pfpd_at1_evolution_v1/miehe_fulltext_addendum_20261002/` | Miehe 全文身份及公式补充 | 已获得并读完全文；公式符号歧义保留 | 核对 AT1/AT2 与文献适用范围 |
| `verification/pfpd_at1_coupling_stability_audit_v1/` | C3 联合残量、能量 Hessian 和映射诊断 | 有可行负曲率；约 1.93 仅为条件分支特征值 | 区分算法发散、驻定与稳定 |
| `verification/pfpd_at1_branch_stability_tracking_v1/` | C3 分支追踪与非光滑稳定性定位 | 约 1% 附近临界变化；旧几何适用范围受限 | C3 历史回归，不外推为 C5/C7 结果 |
| `verification/pfpd_at1_spatial_discretization_audit_v1/` | 识别 full-voxel 体积与点包围盒不一致 | 停在几何硬门 | 解释旧 C3 与 corrected C3 的区别 |
| `verification/pfpd_at1_spatial_discretization_audit_v2/` | corrected C3/C5/C7 体积／面积积分 | 几何通过；当时细网格触发纯拉伸域门 | 下一轮几何定义，按文件读取 |
| `verification/pfpd_at1_spectral_split_v1/`、对应 `results/` | Green/SVK 正负谱分裂及完整态验证 | 本构通过；当时 C5 projected solver 停滞 | 当前本构身份与 oracle 来源 |
| `verification/pfpd_at1_c5_activation_stationarity_audit_v1/`、对应 `results/` | C5 激活 bracket、同载荷正根和活动集诊断 | `PROJECTED_SOLVER_ACTIVE_SET_STAGNATION` | 下一轮 C5 完整态／激活参考，必须保留 |
| `verification/pfpd_at1_petsc_vi_solver_v1/`、对应 `results/` | PETSc RSLS、Jacobian、7 个 C5 根回归 | RSLS 验证；SSLS 并非全部成功 | 当前 runtime、solver 配置及回归根 |
| `verification/pfpd_at1_c5_stability_v1/`、对应 `results/` | C5 Hessian、Schur、临界锥及能量核对 | 15 态稳定，至 max(s)≈0.105，未找到首次失稳 | KEEP FOR REFERENCE；科学 cache 被测试直接读取 |
| `verification/pfpd_at1_c7_early_branch_v1/`、对应 `results/` | C7 自动 reduced pool 与早期正相场分支 | 20 态驻定；max(s)≈0.0442；稳定性未评估 | 当前下一轮直接输入，不默认展开大 SUMMARY |
| `verification/pfpd_at1_activation_mesh_dependence_audit_v1/`、对应 `results/` | 论文映射与已有夹持／激活事实核对 | `STOPPED_AT_PAPER_MAPPING_GATE`；根因未决 | 当前阻塞与下一轮边界，优先读 REVIEW |
| `verification/pfpd_at1_activation_root_cause_and_fix_v1/`、对应 `results/` | common-strain、affine、profile、support、axial-only 对照 | affine 无热点；平衡第一层集中；clamp WEAK，support MIXED；DIAGNOSIS_INCONCLUSIVE | 前轮直接机制及共同应变机械参考 |
| `verification/pfpd_baqp_finite_domain_boundary_consistency_v1/`、对应 `results/` | 四尺寸完整态力分解、制造力、虚功、固定作者对照与边界理论 | global energy EXACT；几何 core patch FAIL；ROOT_CAUSE_PRIMARY=INCONCLUSIVE；7 tests passed；无 BC control | 当前优先读 REVIEW；下一步只建议完整力模板 affine oracle，等待审核 |
| `references/`、各 milestone 的 `references/`、`resources/` | 论文、固定版作者源码、文献核对 | 保留全部原始证据 | 只在公式／身份需要时展开 |
| `provenance/`、各 milestone 的 provenance 和 PRE/POST manifests | 命令、环境、源码与数据身份 | 旧清单原样保留；实际测试读取的缓存不豁免 | 文件来源和复现争议时查 |
| `build_*` 中保留项 | 工具、下载原件、未明确的测试/失败副本 | 不确定项列入 REVIEW_REQUIRED，未自动删除 | 人工评估后续空间清理 |
| `verification/repository_cleanup_v1/` | 本轮 inventory、删除清单、哈希／引用／测试验证 | 以 [REVIEW](clot_clearance/verification/repository_cleanup_v1/REVIEW_REPOSITORY_CLEANUP.md) 和 SUMMARY 为准 | 默认只读摘要；追查具体删除项再查 CSV |
| `verification/cleanup_recovery_audit_v1/` | 123 项分类、SHA/Git/归档检索、只读备份、语义暂存、测试对照 | `CRITICAL_UNRECOVERED=72`；407 passed / 8 failed；等待人工审核 | 当前优先读 REVIEW；具体路径查 CSV |

生产 `pd_clot/`、`tests/`、`configs/` 和全部源码／正式科研证据仍可直接搜索。本轮没有压缩历史目录，没有改写 Git 历史，也没有把遗留科学失败改为通过。
