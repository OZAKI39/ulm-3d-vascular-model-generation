# clot 项目当前上下文

更新：2026-10-06，GitHub handoff。新会话先读[交接正文](HANDOFF_CONTEXT_20261006.md)，再按需读具体 REVIEW；不递归展开历史树。

**当前科学终点：`HIGHER_ORDER_BAQP_BOUNDARY_REINTEGRATION_V1` 已完成，`HIGHER_ORDER_BOUNDARY_MECHANICS=VERIFIED_FOR_INTACT_BENCHMARKS`。** 本轮只做仓库交接；没有继续科学开发。唯一建议下一任务是 `CURVED_CLOT_WALL_INTACT_MECHANICS_V1`，尚未启动。

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

三个原定向套件各 12 项；本次 handoff 服务器复核 **36 passed in 4.42 s**。只用同步文件的 20 项可移植子集见 [环境与命令](HANDOFF_ENVIRONMENT.md)。未运行全量 pytest、平衡 campaign 或图件重绘。

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

曲面/不规则点云仍须检查局部满秩、条件数、体积与边界。`CURVED_CLOT_WALL_READINESS=READY` 仅表示下一轮验证可开始，不表示曲面模型已通过。PD2 二次再现和变分一致不证明积分精确、稳定性、mesh objectivity、phase、fracture、adhesion、Ogden、CFD 或微泡加载。它们均未在本轮启动。

历史 cleanup 最终仍有 72 项旧字节依赖未恢复，用户已说明无备份。清理前全量记录 412 passed / 3 failed；恢复审计 407 passed / 8 failed，不是今天的全量测试结果。原 3 项 staggered 失败见 `tests/test_pfpd_at1_staggered.py`。用户已取消历史 byte-hash/论文 exact-equivalence 的科学开发硬门，但未允许修改旧断言、旧 manifest 或继续清理。保留[恢复报告](clot_clearance/verification/cleanup_recovery_audit_v1/REVIEW_CLEANUP_RECOVERY.md)。

不依据 Git `??`、ignore 或 cache 名称删除文件。C5 stability cache、C7 branch 和全部 accepted states 仍是科学证据。旧 scientific manifests 描述完整服务器产物，不是本轮 GitHub 文件清单。

## 交接范围与停止点

这是选择性 ADD/UPDATE/SYNC。大数组、重复图/论文/环境保留本地；精确路径/大小/哈希与替代摘要见 [LOCAL_ARTIFACTS_NOT_SYNCED.md](LOCAL_ARTIFACTS_NOT_SYNCED.md)，提交范围见 [HANDOFF_SYNC_MANIFEST_20261006.md](HANDOFF_SYNC_MANIFEST_20261006.md)。没有删除历史 tracked 大文件、重写 Git 历史或创建 PR。

交接后停止。下一会话得到新开发指令前，不启动曲面血栓、相场、C5/C7、断裂、CFD、微泡或代理建模。
