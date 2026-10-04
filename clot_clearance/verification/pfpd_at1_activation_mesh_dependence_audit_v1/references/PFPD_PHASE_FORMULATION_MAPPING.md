# PFPD phase formulation：原文与当前候选的逐项映射

结论：**PFPD_PAPER_PHASE_MAPPING = AMBIGUOUS**。裂纹能、核归一化和低损伤阶段的直接相场耦合结构可核对；完整的谱分裂能量和演化律不能标为原文的精确或等价实现。本轮停止在用户第 8、46A 节规定的门槛，未进入新的网格/边界数值诊断。

本结论不撤销历史 AT1、谱分裂及 PETSc 验证。它们验证的是已由用户选定的候选模型；论文整体等价性是另一个命题。没有把此前明确选择的候选解释为本轮新发现的软件缺陷，也没有修改该选择。

## 来源与版本

本轮重新读取 `/home/lzy/projects/2604.09215v1.pdf`，30 页，封面日期 April 13, 2026。SHA256：`a77d29c88bd45bbd1cba78f1e823c09a9157b3be436b5eb6b3855c97a260a157`。页码均指本地 PDF 页码，公式直接从此次全文提取核对，未沿用旧摘要替代原文。提取工具提示存在旋转文本；本审计依据的第 5–9、14 页正文与公式有完整文本。

[arXiv HTML](https://arxiv.org/html/2604.09215v1)用于交叉核对数学排版。当前 HTML 显示 August 24, 2026，与本地封面日期不同，不能声称两者字节相同；本地 PDF 为权威版本。已读的关键公式 (3.4)、(3.17)、(3.25)、(3.59) 内容相符。论文为 CC BY 4.0。以下是数学转写与分析，非全文复制。

代码路径相对项目根 `/home/lzy/projects/temp_storage/clot_clearance_github_sync_20260930/clot_clearance`。行号取本轮未修改的源文件。

## 逐项映射

设无向键 b=(i,j)，`W_b=2*wbar_b*V_i*V_j`；该因子 2 合并有向双积分的两个端点贡献。所有公式使用参考构形体积。

| PAPER：公式、符号、物理含义 | CURRENT CODE：文件、函数、行号、表达式 | STATUS 与适用范围 |
|---|---|---|
| (2.1)，p.3：H(X) 为半径 δ 的物体内邻域 | `kernel.py:40`，`Family.from_arrays`；检查键长与 horizon；`scripts/c7_early_branch.py:28` 读取历史 corrected geometry | EQUIVALENT：离散邻域；δ 固定 |
| (2.4)–(2.5)、(2.10)–(2.11)，pp.4–5：K、Fbar、Ftilde，非局部重构与键向修正 | `pd_clot/pfpd_candidate/baqp.py:37`，`evaluate`，K 加权二阶矩，Ftilde=Favg+(eta−Favg*xi)⊗xi/r²，至行 53 | EQUIVALENT：位移形式与总位置形式在满秩离散重构下等价；此项未重新开展 affine 数值测试 |
| (2.13)–(2.14)，p.5：wbar=ω/2*(1/ω0_i+1/ω0_j)，ω0 为当前几何的 intact family 核积分 | `pd_clot/pfpd_candidate/kernel.py:49`，`Family.from_arrays`，w0 累加 ωV，wb=.5*w*(1/w0[i]+1/w0[j]) | EQUIVALENT：离散双端点归一化已实现，不是遗漏 |
| (2.15)–(2.16)，p.5：完整球邻域的核幅值归一化；示例 cubic spline | `kernel.py:6`，`radial` 默认 (1−r/δ)²；`normalization.py:8` 对同一核求 c0 | EQUIVALENT（核幅值）、APPROXIMATION（选用不同允许径向核）：常数幅值在 wbar、K⁻¹ 重构及 c0 比值中抵消；不是复现 spline 数值例 |
| (3.1)–(3.4)，p.6：总能量中的动能和准静态驻值 | `coupling_stability_audit.py:86`，`AuditProblem.parts`；当前求静态能量与力平衡 | EQUIVALENT：仅准静态范围，动能为零 |
| §3.1、(3.4)，p.6：s(t,X,X′) 是键变量，0 完整、1 断裂 | `kernel.py:19`，`Family` 唯一 i<j 键；`AuditProblem.__init__:54`，ns=len(pairs)，单份 s 向量 | EXACT：键相场，不是每粒子独立相场 |
| (3.2)–(3.3)，p.6：经典局部相场含 ℓc|∇s|²；转入 PFPD 时取 ℓc=δ，并明确去除此显式梯度项 | `AuditProblem.parts:86`；`at1_evolution.py:70`：无梯度/邻键差项 | EXACT（项的缺席）：不得把经典式 (3.2) 的梯度项当成 PFPD 式 (3.4) 漏项 |
| (3.4)–(3.8)，pp.6–7：Ecrack=∫∫wbar*Gc*s/(c0δ)，且逐键裂纹能不依赖邻键 | `normalization.py:14`，`crack_energy`；`adr_adapter_at1.py:32`，`pair_weights`；`coupling_stability_audit.py:87`：weights@(2*yc*s) | EQUIVALENT：2Yc=Gc/(c0δ)，无向键和有向积分计数一致 |
| (3.8)，p.7：γc(s)=s，AT1 线性裂纹能 | `at1_evolution.py:71`：2*yc*s；`AuditProblem.parts:87` | EXACT |
| (3.9)，p.7：g(s)=(1−s)² | `degradation.py:11`，`energetic`；`svk_spectral_split.py:87` | EXACT |
| (3.4)、(3.10)、(3.15)，pp.6–8：打印能量为 g(s)Ψ0，打印应力为 g(s)P0 | `svk_spectral_split.py:84`，`degraded_svk_split_response`：能量 g*Ψplus+Ψminus，应力 F@(g*Splus+Sminus)；`scripts/spectral_split_audit.py:119` 将其接入 BAQP | APPROXIMATION / 有意候选扩展：混合应变且 s>0 时不等价；s=0 或 Ψminus=0 时能量相等。不是缺少 fracture gradient |
| (3.11)–(3.14)，pp.7–8：s≤sc 时 h=1；大损伤时 h 改变 K、B | `degradation.py:14`，`kinematic` 给出完整 h；当前 `spectral_split_audit.py:130` 在 s<.95 范围显式 h=1，门限禁止越界 | EXACT：本次及历史 early-branch 域内；不宣称实现高损伤耦合演化 |
| (3.16)–(3.18)，p.8：q=−2(1−s)Ψ0+Gc/(c0δ)，Yc=Gc/(2c0δ) | `at1_evolution.py:13`，`critical_driver`；`scripts/c5_activation_audit.py:47`，`q_fields`，q=−2(1−s)Y+2yc | EQUIVALENT：采用 (3.19) 的 tensile driver 后；不是打印的全 Ψ0 混合应变能量导数 |
| (3.19)，p.8：Y=Ψ0+，仅拉伸驱动 | `svk_spectral_split.py:32`，`_parts` 用 Green 应变谱正负部分；`degraded_svk_split_response:89` 返回 Y=Ψplus | EXACT：Y 符号关系；AMBIGUOUS：原文未在此给出与 (3.4)/(3.10) 同时一致的具体 Green/SVK 分裂能量和应力公式 |
| (3.23)，p.9：以最大正主 Cauchy 应力定义另一种驱动 | 当前使用 (3.19)，不用 (3.23) | EQUIVALENT：论文提供的驱动选择，不是漏项 |
| (3.24)，p.9：历史最大驱动 Yhistory；不可逆性 | `at1_evolution.py:66`：max(old,trial)；`AuditProblem:57` 保留前一接受载荷的 lower；`petsc_vi_candidate.py:20` 向 PETSc 传递上下界 | APPROXIMATION：实现相场不可逆约束，不是逐字采用原文历史闭式演化 |
| (3.25)，p.9：s=min(1,Yhistory/(Yhistory+Yc)) | `at1_evolution.py:32`：max(0,1−Yc/Y)；PETSc 求相同 AT1 能量的受约束驻值 | AMBIGUOUS：当前不实现 (3.25)；原文 (3.25) 与其 (3.8)/(3.17) 的 AT1 驻值关系本身不相容。是既知、已记录的候选选择，不是本轮新发现 |
| (3.26)–(3.27)，p.9：粒子损伤是键相场的 wbarV 聚合 | `scripts/c7_spatial_evidence.py`，`region`：按端点聚合用于空间诊断；动力学未知量仍是 bond s | EQUIVALENT：可视化/测度聚合，不产生新的粒子相场自由度或邻键相场方程 |
| (3.28)–(3.61)，pp.10–14：无限域平面完全断裂标定；(3.59) c0=∫ωρ³/(2∫ωρ²) | `normalization.py:8`，`c0`：相同两个积分；quadratic 核 c0=1/4 | EXACT：连续标定常数；有限、强截断小试件的离散裂纹能收敛不由这个等式自动保证 |
| (2.14)–(2.15)，p.5；§2.1，p.3：截断 family 归一化；Dirichlet 施加于边界层 | `kernel.py:49` 使用真实 family；`scripts/c7_early_branch.py:30` 和 `AuditProblem:53` 固定上下表面整粒子 xyz | APPROXIMATION：当前 boundary-node 试件仅约束最外层，未复现某一论文试件的体积边界层；本轮不替换夹持方式 |

## 必答问题

**Q1：论文是否有显式邻相场正则项？** PFPD 式 (3.4)–(3.8) 没有 |∇s|²、(s_b−s_c)² 或邻键 phase-phase 项。第 7 页明确将单键裂纹能描述为局部于该键。式 (3.2) 的梯度项属于前述经典相场，原文随后明确移除。没有找到 MISSING_REQUIRED_PHASE_TERM。

**Q2：Hss 为什么是对角？** 在固定 x、h=1、几何积分权重不依赖 s 时，当前能量为

\[
 E(x,s)=\sum_b W_b[(1-s_b)^2\Psi_b^+(x)+\Psi_b^-(x)+2Y_c s_b].
\]

因此 \(\partial E/\partial s_b=W_b[-2(1-s_b)\Psi_b^++2Y_c]\)，且
\[
 (H_{ss})_{bc}=2W_b\Psi_b^+\delta_{bc}.
\]

打印的原文 AT1 能量同样逐键可分（将 Ψ+ 换成 Ψ0），故此处没有应实现而未实现的非对角梯度块。原 `scripts/c5_stationarity_diagnostics.py:73` 的 Qs=2Y 与该结构一致。加权物理 Hss 与按 D=max(Y,Yc) 缩放的求解器 Jss 必须区分；二者此域内都只有相场对角项。本轮只作代数/源码检查，未构造任何 C7 Hessian 或 Schur 补。

若离开 h=1 区域，K、Ftilde 会随其他键的 s 改变，上述 fixed-x 可分性不能直接推广。即使 h=1，经过位移再平衡也可能间接耦合各键；Hss 对角本身不能证明损伤区宽度必然趋零。

**Q3：c0、δ、kernel 控制什么？** δ 决定非局部力学支持和裂纹能密度的 Gc/(c0δ) 系数；核通过 BAQP 加权及 (3.59) 进入 c0。原文以跨越平面裂纹的完全损伤键，在无限域、宽约 2δ 的端点贡献带上匹配 Griffith 能量。它不是在当前小立方体中证明某个固定相场局部化宽度的定理。当前 δ≈L，没有真正包含完整半径 δ 球的 bulk 点；这一几何事实不能单独给出网格偏移的原因。

**Q4：是否完整实现？** 低损伤域的 bond AT1 线性裂纹能、双端点核权重、c0 和 Yc 已实现；没有已证实缺失的直接相场正则项。完整的打印能量/应力与演化律未被原样实现。谱分裂和 AT1 约束演化是此前明确采用的候选，不能将候选的数值验证改写为整篇论文等价性证明。

## 可直接复核的代数差异

把原文 Ψ0=Ψ+ + Ψ− 代入打印能量：

\[
 e_{candidate}-e_{paper}=[1-(1-s)^2]\Psi^-.
\]

在 s>0、Ψ−>0 时差异严格非零。在 s=0 两个能量值相等，但对 s 的单侧导数差为 2Ψ−。原文 (3.19) 选择拉伸驱动，却没有在 (3.4)/(3.10) 写出当前代码采用的对应能量/应力替换，因此这是需要明确规范解释的地方。

原文 (3.17) 采用 tensile driver 后给出的受约束 AT1 最小值为 max(0,1−Yc/Y)。以 Y=Yc、旧相场为零为例，该值为 0；打印的 (3.25) 给出 1/2，代入 AT1 q 后得到 q=Yc>0，内部点不满足驻值条件。这个矛盾不可能靠调整 PETSc 容差消除，也不能通过重新命名 Y 自动消失。

上述差异并没有解释 C5→C7 的阈值偏移：两级历史结果使用的是同一候选。仍需同应变、真实机械平衡下的 Ψ+ 集中证据才能归因。因此 ROOT_CAUSE_PRIMARY=INCONCLUSIVE，而非 PHASE_FORMULATION_INCOMPLETE。

## 硬门决定

触发用户第 46A 节的“论文 phase energy 与当前实现不一致”，且第 8 节要求的 VERIFIED / EQUIVALENT 整体映射未成立。停止后续 common-strain equilibrium、affine oracle、profile collapse、support-driver 和 clamp-control 数值诊断；不实施修正。缺失项计数为 0，未使用 MISSING_TERM 或 PHASE_FORMULATION_INCOMPLETE 标签。

唯一待人工审核事项：是否将已授权的 AT1+谱分裂候选作为后续网格诊断的规范，并仅要求其裂纹能/核归一化/耦合结构与原文对应。此次未自行放宽用户新设的论文等价性门槛。
