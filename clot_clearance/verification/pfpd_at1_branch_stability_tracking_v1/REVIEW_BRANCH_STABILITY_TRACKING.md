# PFPD_AT1_BRANCH_STABILITY_TRACKING_V1 人工审核报告

当前候选分支首次丧失局部变分稳定性的已分辨区间为 **[0.0099921875, 0.0100000000] 工程应变**，即 **[0.99921875%, 1.00000000%]**。两套预声明加载网格得到相同结果，区间宽 **7.8125e-6**，小于目标 `1e-5`。这是 **NONSMOOTH_STABILITY_BRACKET**：在上端点，185 根键一起进入弱下界状态，出现可行的单侧负曲率方向。它不是一条已证明光滑的特征值曲线穿过零点，更不是血栓材料的实验断裂应变。

工作区：`/home/lzy/projects/temp_storage/clot_clearance_github_sync_20260930/clot_clearance`。Git 根目录为其父目录；分支 `clot-clearance/energy-regularization-20260930`，HEAD `f0e891686b9b107510ab3292aa27ce83638b3339`。原有未跟踪成果也纳入保护。适用祖先与工作区未发现另行约束的 AGENTS.md，读取的用户级 AGENTS.md 为空。保护清单涵盖 **17,084 个原文件、273 个原链接**；生产、历史结果、既有测试均未修改，详见 [PRESERVATION_CHECK](provenance/PRESERVATION_CHECK.json)。

## 1. 为什么做这一轮

上一轮说明：工程应变 0.009 的 A 状态局部稳定，0.013 的 B 状态虽被旧交替算法接受，却具有负曲率；失败载荷的两个联合驻点也不是局部极小值。本轮要缩小第一次稳定性丧失的位置，并检查它是否依赖加载步长。

本轮保持原 3×3×3 立方体、185 根键、SVK 材料、Gc_demo、`g(s)`、`h=1`、`sc=0.95`、ADR、原能量与纯拉伸门限。没有修复原交替算法，也没有重复研究普通 fixed-point 谱。完整原提示词见 [USER_PROMPT](provenance/USER_PROMPT.txt)，计算前声明的数值配置见 [PREDECLARED_CONFIG](provenance/PREDECLARED_CONFIG.json)。

## 2. “平衡”与“稳定”有什么区别

平衡意味着附近的一阶驱动力接近零；稳定还要求沿允许的小扰动方向不能下坡。球可以停在谷底，也可以暂时停在山脊上。求解器说“收敛”只回答是否找到了平衡，不能单独回答是否稳定。

下界使这个比喻多一层限制：有些方向根本不允许走；有些方向虽然允许，但一开始仍然上坡。这后一类方向即便二阶曲率为负，也不能只凭该曲率宣称失稳。你已明确选择：**采用临界锥判定，完整保留辅助探针**。本轮据此只让一阶能量变化为零的可行方向参与失稳判定。强下界相场方向在临界锥中固定为零；弱下界方向只能向相场增加的一侧移动。授权记录见 [CRITICAL_CONE_DECISION](provenance/CRITICAL_CONE_DECISION.json)。

## 3. 怎么追踪同一分支

新增 [branch_stability_tracking.py](../../pd_clot/pfpd_candidate/branch_stability_tracking.py) 用 SciPy `least_squares` 求残差方程，而不是最小化候选能量。自由变量是 27 个自由位移分量与 185 个相场变量；固定边界不进入未知量。默认相场初值严格等于前一个已接受加载状态，没有预先增加损伤；位移初值仅作简单仿射加载增量预测。

残差沿用上一轮尺度：机械分量是 `force_free / 4.2822421875e-6 N`；相场分量是原 `projected_kkt(s,q_scaled,s_lower,sc)`，其中 `q_scaled` 是上一轮的密度梯度归一化定义，不是另加的优化权重。物理下界 `s_lower` 来自前一**接受加载状态**，不是同一加载步的失败内循环。位移变量尺度为 `1.25e-6 m`，能量尺度为 `5.3263037109375e-12 J`。这些尺度固定保存于 [RESIDUAL_SCALING](provenance/RESIDUAL_SCALING.json)。

采用 `least_squares(method='dogbox', tr_solver='exact', x_scale=1, loss='linear')`，`ftol=xtol=gtol=1e-13`、`max_nfev=300`，均在运行前声明。dogbox 支持边界且可从原下界出发，不需人工抬升相场。使用的是成熟残差求解工具，未替换原生产 ADR 或交替求解器。[SciPy least_squares 文档](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html)

能量、力、相场梯度、投影 KKT、bounds、权重和 Hessian 均直接调用上一轮实现。新增 Jacobian 仅是对原残差作链式求导，原 Hessian 组装函数保持不变。原 Hessian writer 的输出目录在隔离进程内临时指向新目录，并在 `finally` 中恢复；未复制另一套 Hessian，未写入旧目录。测试验证复用结果逐元素完全一致。

接受一个驻点必须同时满足 least_squares 成功、数值有限、力残差 ≤ `1e-9`、投影 KKT ≤ `1e-8`、不可逆和上下界、`s<sc`、原纯拉伸域 `min Green ≥ -1e-9`、固定边界精确不动。负 Hessian 不阻止保存 `STATIONARY_UNSTABLE`。任一求解试探点越纯拉伸域或达到 sc 即停止；仅导数探针可以读取能量代数延拓，不能成为分支状态。

先用 A/B 各自的保存状态和前一加载步预测值做 **4/4 复现**。最大位形差 `5.912e-17 m`，最大相场差 `4.322e-11`，能量相对差 ≤ `9.0e-15`，最小 Hessian 特征值差 ≤ `1.047e-10`，均通过预声明门槛。A 仍稳定，B 仍不稳定。复现通过后才开始新加载网格，见 [REPRODUCTION](comparisons/REPRODUCTION.json)。

每个点记录位形/相场距离、能量/反力变化和活跃集合的键 ID 变化。没有从不稳定上端点反向“愈合”；细化点始终从稳定下端点向前求解。仅执行提示词明确要求的两套网格至 0.013 及区间细化，没有超过该范围或切换到另一分支。

## 4. Coarse tracking 结果

Schedule A 的全部九个预声明点都保留，没有跳点。表中能量为候选内部能，不代表已验证全局外功—能量平衡。

| 工程应变 | 拉伸反力（μN） | 最大 s | 总能量（pJ） | 常规约化最小特征值 | 分类使用的曲率 | 稳定性 |
|---:|---:|---:|---:|---:|---:|---|
| 0.0090 | 3.84828659 | 0 | 4.31001423 | +0.062826346 | +0.062826346 | 稳定 |
| 0.0095 | 4.06510468 | 0 | 4.80459786 | +0.062880881 | +0.062880881 | 稳定 |
| 0.0100 | 4.28224219 | 约 4.53e-14 | 5.32630371 | +0.062935444 | **−0.014702137，单侧见证** | 不稳定 |
| 0.0105 | 3.69823251 | 0.09342161 | 5.82389559 | −0.045201906 | −0.045201906 | 不稳定 |
| 0.0110 | 3.21569546 | 0.17437544 | 6.25508027 | −0.045640643 | −0.045640643 | 不稳定 |
| 0.0115 | 2.81352177 | 0.24498364 | 6.63116100 | −0.045625055 | −0.045625055 | 不稳定 |
| 0.0120 | 2.47565926 | 0.30693547 | 6.96113527 | −0.045211278 | −0.045211278 | 不稳定 |
| 0.0125 | 2.18975411 | 0.36158907 | 7.25223694 | −0.044459563 | −0.044459563 | 不稳定 |
| 0.0130 | 1.94619381 | 0.41004639 | 7.51033534 | −0.043429468 | −0.043429468 | 不稳定 |

Schedule A：9 个驻点，**2 稳定、7 不稳定、0 未分辨、0 marginal**。Schedule B：17 个驻点，**4 稳定、13 不稳定、0 未分辨、0 marginal**。细化另增加 6+5 个稳定点。合计 37 个追踪/细化驻点，外加 4 个复现点。原始表含规定边界位移、加权平均相场、弹性/裂纹/总能量，见 [coarse_tracking/TRACKING](coarse_tracking/TRACKING.json)、[fine_tracking/TRACKING](fine_tracking/TRACKING.json)、[COUNTS](comparisons/COUNTS.json)。

## 5. Stability eigenvalue 怎样变化

对每个驻点都采用原无量纲约定下的 `1e-5、3e-6、1e-6` 三档差分；41 个状态共保存 123 份 Hessian 审计。原始矩阵的最大相对非对称误差约 `1.582e-8`，Hss 与原解析对角结构的最大相对差约 `2.273e-10`，均低于原门槛。误差先于矩阵对称化量测。每点保存第一、第二特征值、Hxx 最小特征值、Schur 最小特征值和 Hxx 条件数；Hxx 条件数最大约 3.4481，未在病态机械块上硬求逆。对称特征分解直接使用 [SciPy eigh](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.eigh.html)。

在稳定下端点，185 根键均为强下界变量，临界锥的相场分量为零，常规机械约化空间的最小特征值约 **+0.06293459093**。在上端点，185 根键均变为弱下界变量。若仅删除这些键，得到的常规约化最小特征值仍约 **+0.06293544367**，这时单看该值会漏掉允许的单侧损伤方向。

本轮保留弱下界变量的单侧约束，另找到归一化方向 v，使一阶能量变化约 `1.54e-16`，在预声明 `1e-8` 临界方向容差内；三档 Hessian 都给出 **vᵀHv ≈ −0.0147021366**，且方向满足物理下界、sc、固定边界与原模型域。这足以证明该驻点不是局部极小值。

必须区分两个量：常规约化特征值是双侧线性子空间上的最小值；弱下界的单侧见证是临界锥内一个负 Rayleigh 商，**不是已经求出了锥内精确最小值**。因此只能说锥内最小曲率不大于约 −0.01470214。所有图件均把这种见证与普通特征值区分标注。

每点还保留确定性随机位移、内部相场、混合方向、下界开启方向，以及从负特征向量投影并利用原机械块求响应的方向。强下界开启而一阶能量上升的辅助方向不参与失稳分类；正的随机探针也不用于证明弱下界锥处处正定。实际上，本轮弱下界点都找到了负的可行临界方向，未靠“没采到负值”判稳定。见每个 [hessian](hessian/) 状态目录内的 `CONE_PROBES.json`、`CONE_AND_MODE.npz` 和三档块矩阵。

## 6. 临界 bracket 是多少

Schedule A 的第一对相邻稳定/不稳定状态是 **[0.0095, 0.0100]**，Schedule B 是 **[0.00975, 0.0100]**。分别完成 6 次与 5 次二分细化，最终都得到：

```text
CRITICAL_LOAD_INTERVAL = [0.0099921875, 0.0100000000]
width = 0.0000078125 = 7.8125e-6
midpoint estimate ≈ 0.009996
method = NONSMOOTH_STABILITY_BRACKET
```

中点只是区间摘要，不能当作更高精度的真实根。区间横跨活跃集合变化，允许方向集合发生扩张；没有证据证明存在一条光滑的 `lambda_min(epsilon)=0` 曲线。因此没有使用 brentq，也没有对切换点作样条或光滑插值。[SciPy brentq 文档](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.brentq.html)

所有细化点从当前稳定下端点作向前预测，父状态与下界逐一保存。最终上下端点的下界也相容：上端点求解所用物理相场下界与最后稳定下端点相场逐元素相同。见 [SCHEDULE_A](critical_bracket/SCHEDULE_A.json)、[SCHEDULE_B](critical_bracket/SCHEDULE_B.json)、[FINAL_BRACKET_BOUND_COHERENCE](comparisons/FINAL_BRACKET_BOUND_COHERENCE.json)。区间仅在规定的采样和细化分辨率内定位第一次失稳，未声称穷尽更细的未采样行为。

## 7. 两种 load schedule 是否一致

两套网格都从 0.009 开始，步长分别为 `5e-4` 与 `2.5e-4`，都追踪到用户要求的 0.013。它们的最终区间完全相同，端点应变差为 0。

在 15 个匹配载荷点上，两套结果的最大相场差约 **6.295e-14**，最大坐标差约 **3.663e-20 m**，反力最大相对差约 **1.042e-14**，最小约化特征值/曲率指标最大绝对差约 **2.989e-11**；活跃键集合一致。均通过预声明门槛，`LOAD_RESOLUTION_CONSISTENCY=PASSED`。

两套网格的下界来源是各自前一个接受点，比较没有强行替换加载历史。实际计算显示它们跟踪到相同的近均匀驻定分支。没有使用可选 Schedule C，因为两套规定网格已在目标分辨率内一致；这不是空间网格收敛或材料标定证明。详见 [SCHEDULE_CONSISTENCY](comparisons/SCHEDULE_CONSISTENCY.json)。

## 8. Active set 是否在临界附近改变

有明确的大规模变化：最后稳定点为 **185 strong lower / 0 weak / 0 interior**；0.01 点变为 **0 strong / 185 weak / 0 interior**；随后采样的损伤驻点为 **0 strong / 0 weak / 185 interior**。

进入弱下界表示相场尚未明显增长，但阻止它增长的一阶能量斜率已接近零，单侧损伤方向可以进入临界锥。这个集合变化正是普通约化 Hessian 仍正、完整约束稳定性却已失败的关键。每根键的 ID、进入和离开集合的事件都已保存；标记为 `BRANCH_CONTINUITY_EVENT / MASS_ACTIVE_SET_SWITCH`。它不等于状态本身发生有限跳跃，也不允许把前后曲率作光滑插值。

位形与相场增量、反力和能量增量没有被忽略，均在各状态的连续性记录中。见 [active_sets](active_sets/)，特别是 `A_grid_02_EVENT.json` 和 `B_grid_04_EVENT.json`。

## 9. 临界模态是什么

本轮分辨了一个在 0.01 状态可行、具有负曲率的单侧方向。它来自原耦合 Hessian 的负方向，在弱下界条件下保留非负相场分量，再用已经检查过的 Hxx 解出相应机械响应。它是**数值/变分临界方向的见证**，并非唯一的临界特征向量，更不是裂纹外形。

保存了普通最小特征向量、选中的锥方向和全部辅助探针。比较普通特征向量时先做符号对齐；相场分量为零时，该分量余弦记为未定义，不补成 1。最终区间两端完整模态的对齐余弦分别约 `3.15e-12`（A）和 `2.06e-13`（B），因此不能报告模态连续，应标 **MODE_SWITCH_NEAR_ONSET**。这也反映了允许相场方向的变化。

从 B 的弱下界临界方向到最近内部不稳定点 0.01025，完整模态余弦约 0.77269，位移分量约 0.999967，相场分量约 0.70710。单侧锥方向与普通双侧特征向量是不同对象，数值接近程度不能证明同一光滑模态。尤其 B 的区间跨越比较中，为数学对齐而翻转锥方向会使它不可行；该标志已记录，展示和能量探针始终使用原可行方向，没有为了提高余弦反转实际数据。见 [CONTINUITY](critical_mode/CONTINUITY.json)。

空间图左侧为按原键权重聚合到粒子的相场开启分量，右侧为同一方向的位移幅值；各自归一化仅用于显示。坐标使用保存位形，不放大几何位移，不画箭头。图上醒目标明 `NUMERICAL / VARIATIONAL CRITICAL MODE` 与 `NOT PHYSICAL CRACK PREDICTION`。

## 10. Energy landscape 如何从谷底变成不稳定方向

按预声明无量纲幅值 `±1e-6、±3e-6、±1e-5、±3e-5` 逐项试探并记录可行性，没有把非法点画成有效能量。另有计算前允许的小幅值 `±1e-7、±3e-7` 用于模型域受限时的补充检查。

| 检查点 | 合法方向范围与结果 |
|---|---|
| 最后稳定点 0.0099921875 | 原要求的较大幅值均越过原纯拉伸门限；补充的 ±1e-7、±3e-7 合法。±3e-7 的能量上升约 1.5375e-26 J，Hessian 二次参考约 1.5084e-26 J。 |
| 弱下界点 0.0100000000 | 负幅值违反相场不可逆下界；正幅值合法。+3e-5 能量下降约 3.5204e-23 J，与负曲率二次参考一致。 |
| 最近内部不稳定点 0.01025 | 所有规定正负幅值均合法；±3e-5 两侧能量均下降约 1.074e-22 J。 |

因此图 5 用三个面板同时呈现稳定的小幅向上趋势、切换点合法的单侧向下趋势，以及最近内部不稳定点的双侧向下趋势。无法满足约束的对称曲线没有伪造。

很小幅值的能量差接近 float64 相消误差，部分不稳定点的极小探针会出现正负噪声。原数据保留，未改符号、未滤波；图中的虚线是 Hessian 二次参考，实点才是原能量调用。可靠性依赖多档 Hessian、一阶中性、可行性和较大合法幅值的交叉证据，而非单个极小能量差。见 [energy_landscape](energy_landscape/)。

## 11. Damage onset 与 stability onset 的先后

相场“已检测到增长”的预声明阈值为 `max s > 1e-8`。按实际采样，Schedule A 第一次超过该阈值是 **0.0105**，Schedule B 是 **0.01025**；对应采样区间分别为 **(0.0100, 0.0105]** 与 **(0.0100, 0.01025]**。

局部 AT1 模型原有的解析损伤阈值为 0.0100。稳定性区间的上端点也为 0.0100，但此时相场只有约 `4.53e-14`，远低于检测阈值，不能将它宣称为已经观测到有限损伤。

本轮采用 **COINCIDENT_WITHIN_RESOLUTION**，含义限于：损伤启动采样区间与稳定性区间在该切换处相接，并与原局部 AT1 阈值相符。实际正相场启动没有另行细化，分辨率较稳定性区间粗；不能以第一次非零采样的应变声称严格先后，更不能声称“损伤导致失稳”。见 [DAMAGE_ONSET](comparisons/DAMAGE_ONSET.json)。

## 12. Onset 时离 sc 还有多远

最终区间下端点 `max s=0`；上端点 `max s≈4.527e-14`，距 `sc=0.95` 约 **0.95**。因此本轮发现发生在当前 `h=1` 的早期损伤设定中，远未到高损伤门限。

最近内部不稳定点 0.01025 的 `max s≈0.04842233`，同样远低于 sc；整个规定网格到 0.013 的最大相场约 0.41004639。没有改变 sc，也没有实现尚未验证的高损伤 `h(s)<1`。

## 13. 当前能否叫 fracture initiation

不能称为真实血栓裂纹启动，也不能当作材料失效应变。本轮定位的是**当前候选能量、当前 3×3×3 几何、夹持条件、SVK 参数和演示 Gc 下的分支局部变分失稳区间**。

它说明过去的“已接受驻点”可能处在不稳定分支上，且失稳可在弱下界切换时发生。它没有证明真实血栓失稳，没有验证裂纹传播，没有建立实验对应关系，也没有修好原交替求解器。原普通 fixed-point 谱仍未分辨，原耦合演化和完整加载—卸载—重载验证状态不变。

## 14. Tests

基线 **231 passed / 3 known failed / 0 skipped**；新增诊断最终独立运行 **24 passed**；完整回归 **255 passed / 相同 3 failed / 0 skipped**，原 62 条 warning 保留。原测试的 ID、失败消息和 traceback 逐项比对未变，未改断言、未加 xfail/skip。三个原失败仍是 `tests/test_pfpd_at1_staggered.py` 中：

1. `test_saved_coupled_mechanics_phase_and_energy_all_converged`
2. `test_coupled_no_bond_heals_on_unload_and_reload`
3. `test_extra_stagger_iterations_change_no_accepted_state_materially`

新增测试覆盖 A/B 实际复现、全部驻点的力/KKT 门限、前一加载步下界、初值相场不预增、不可逆性、固定边界、纯拉伸/sc 停止、残差 Jacobian 独立校验、原 Hessian 逐元素复用、三档差分符号与块一致性、强/弱下界临界锥规则、完整网格、键集合事件、模态符号对齐、稳定/不稳定分类、二分收缩、向前父状态、步长一致性、合法能量探针及历史保护。

初版新增测试中，残差 Jacobian 的独立差分参考使用 `3e-6`，一项近相消分量超过 `atol=1e-8` 的检查。进一步保存了预声明三档 epsilon 的 15 项独立检查：最大向量相对误差约 `1.59e-8、5.59e-8、1.82e-7`，随步长减小呈相消误差放大。最终独立参考改用原已声明的最大差分步长 `1e-5`，**原测试断言的 rtol/atol 未放宽、求解 Jacobian 步长未改、物理和验收门槛未改**。初次失败日志仍在 [DIAGNOSTIC.log](tests/DIAGNOSTIC.log)，三档数据见 [RESIDUAL_JACOBIAN_EPSILON_STUDY](tests/RESIDUAL_JACOBIAN_EPSILON_STUDY.json)，最终通过日志见 [DIAGNOSTIC_FINAL](tests/DIAGNOSTIC_FINAL.log)。

完整回归日志：[FINAL_FULL](tests/FINAL_FULL.log)；精确比较：[COMPARISON](tests/COMPARISON.json)。完整 pytest 返回非零退出码源于保留的三个原失败，不能称为全部测试通过。

机器为 i7-13700HX、12 个物理核/24 个逻辑核，RAM 约 8.18 GB。Python 3.13.11、NumPy 2.5.3、SciPy 1.18.1、Numba 0.67.0。ProcessPoolExecutor 按运行时可用 RAM 扣除至少总 RAM 的 20% 和父进程预算，实际选 **1 个 worker**；没有为了并行占满用户正在使用的内存。每个 worker 的 BLAS/OMP/MKL 均单线程，保留原 1 GiB/worker 预算。见 [HARDWARE](provenance/HARDWARE.json)、[WORKER_BUDGET_TRACK](provenance/WORKER_BUDGET_TRACK.json)。

所有缓存和 pytest 临时文件写入本轮新目录。原文件/链接核对见 [PRE_RUN_MANIFEST](provenance/PRE_RUN_MANIFEST.json)、[POST_RUN_MANIFEST](provenance/POST_RUN_MANIFEST.json)。配置 SHA256 为 `83f80ffbcdd3e14c33cdcdbb38062f118eeae633be8b1ec1fe0a1a4ed5974a36`。原有 17,084 个文件和 273 个链接保持不变，生产/历史/测试改动均为 0；最终检查见 [FINAL_REVIEW](provenance/FINAL_REVIEW.json)。

## 15. 当前科学状态

```text
BRANCH_TRACKER: VERIFIED
STATIONARY_BRANCH: TRACKED
STABILITY_CLASSIFICATION: RESOLVED
FIRST_STABILITY_LOSS: BRACKETED
CRITICAL_STRAIN_INTERVAL: [0.0099921875, 0.0100000000]
LOAD_RESOLUTION_CONSISTENCY: PASSED
ACTIVE_SET_AT_ONSET: NONSMOOTH
CRITICAL_MODE: RESOLVED
DAMAGE_ONSET_RELATIVE_TO_STABILITY_LOSS: COINCIDENT_WITHIN_RESOLUTION
HIGH_DAMAGE_GATE: NOT_REACHED_AT_ONSET
FIXED_POINT_ORDINARY_SPECTRUM: REMAINS_UNRESOLVED
PFPD_AT1_EVOLUTION_V1: REMAINS_NOT_VERIFIED
PRODUCTION MODEL: UNCHANGED
```

`CRITICAL_MODE: RESOLVED` 指分辨了合法负曲率方向；跨区间连续性另为 **MODE_SWITCH_NEAR_ONSET**。它不表示求得唯一的锥最小方向，也不表示普通 fixed-point 特征向量已变成可解释的物理裂纹。

八组图均为黑底、Arial 常规字重、英文、无箭头的 PNG/PDF：

| 图 | 审核问题 | 文件 |
|---|---|---|
| stability_eigenvalue_vs_strain | 首次负曲率在哪里出现；区分普通约化特征值与锥见证 | [PNG](figures/stability_eigenvalue_vs_strain.png) / [PDF](figures/stability_eigenvalue_vs_strain.pdf) |
| critical_bracket_zoom | 两套区间是否一致；不作光滑根插值 | [PNG](figures/critical_bracket_zoom.png) / [PDF](figures/critical_bracket_zoom.pdf) |
| branch_force_phase | 反力、最大/平均相场和两个启动区间的关系 | [PNG](figures/branch_force_phase.png) / [PDF](figures/branch_force_phase.pdf) |
| active_set_evolution | 185 根键何时从强下界变弱下界，再进入内部 | [PNG](figures/active_set_evolution.png) / [PDF](figures/active_set_evolution.pdf) |
| curvature_along_critical_mode | 合法小扰动的上坡、单侧下坡与双侧下坡 | [PNG](figures/curvature_along_critical_mode.png) / [PDF](figures/curvature_along_critical_mode.pdf) |
| critical_mode_map | 数值临界方向的相场分量与位移幅值 | [PNG](figures/critical_mode_map.png) / [PDF](figures/critical_mode_map.pdf) |
| stationarity_quality | 不稳定点是否仍满足平衡/KKT，而不是未求平衡 | [PNG](figures/stationarity_quality.png) / [PDF](figures/stationarity_quality.pdf) |
| schedule_consistency | 同一端点区间、相场与反力是否依赖加载步长 | [PNG](figures/schedule_consistency.png) / [PDF](figures/schedule_consistency.pdf) |

代码入口为 [run_branch_stability_tracking.py](../../scripts/run_branch_stability_tracking.py)、[analyze_branch_stability_tracking.py](../../scripts/analyze_branch_stability_tracking.py)、[visualize_branch_stability_tracking.py](../../scripts/visualize_branch_stability_tracking.py)、[review_branch_stability_tracking.py](../../scripts/review_branch_stability_tracking.py)。本轮所有科学状态保存在新的 verification 子目录内，无需另外复制旧 results；已有状态目录拒绝重复创建，避免误覆盖。复算应在新副本中进行。

## 16. 唯一下一步建议

人工审核后，进行一次**临界锥负曲率的空间离散一致性审计**：固定物理尺寸、材料、Gc_demo 与物理视界，检查加密粒子离散后该弱下界负方向和稳定性区间是否保留。这样可先判断当前 3×3×3 开发算例的发现是否对离散方式敏感，再决定后续模型或算法工作。本轮没有实施该步骤，未启动分支切换、失稳后新分支追踪或新的物理耦合。
