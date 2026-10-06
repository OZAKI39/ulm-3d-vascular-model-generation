# Clot CFD–PD / PFPD：2026-10-06 科学交接

本文件供新会话首先阅读。当前科学终点是 **HIGHER_ORDER_BAQP_BOUNDARY_REINTEGRATION_V1 = VERIFIED_FOR_INTACT_BENCHMARKS**。本轮只整理并同步已经完成的工作；没有运行新的平衡、相场或流体计算。下一建议任务为 **CURVED_CLOT_WALL_INTACT_MECHANICS_V1**，尚未启动，仍需新的开发指令。

## 1. 项目目标与控制导向路线

用户最终目标链是 **PZT 激励 → 微泡响应 → 泡形态变化 → streaming / shear → 血栓响应**。面向控制与代理模型的策略，是先让各模块的输入、输出和适用范围可验证，再以可信的局部高保真结果建立紧凑响应关系，最终支持激励选择与控制。这个目标来自本轮交接要求；本分支没有提供已训练代理、已标定 PZT—微泡模型或端到端闭环控制结果。

目前实际工作集中在最后一个模块的固体力学基础：给定几何、材料和局部替代位移/载荷，验证离散重构、储能、内力与非局部边界转换。这里的 *local surrogate field* 是用于验证的解析连续体替代场，不应与未来统计/机器学习代理模型混同。历史直管、侵蚀、碎片输运演示存在，但当前 PD2 候选尚未接回真实 BraVa 流场，也没有完成生物材料标定。

开发采取逐层证据链：运动学正确 → 同一能量的变分内力 → 有限域边界转换与平衡 → 曲面完整态 → 后续再考虑损伤与流固作用。后续任务应围绕当前验证缺口展开，不应把过去的演示链当作已验证的完整物理系统。

## 2. 分支、位置与模型身份

| 项目 | 值 |
|---|---|
| 来源分支 | `clot-clearance/pfpd-progress-20261005` |
| 交接前 HEAD | `695e75166aa864809460ea9a1845fe72abce229c` |
| 交接分支 | `clot-clearance/pd2-boundary-handoff-20261006` |
| Git 根目录 | `/home/lzy/projects/temp_storage/clot_clearance_github_sync_20260930` |
| Python 子项目 | `clot_clearance/` |
| 最新科学状态 | `HIGHER_ORDER_BOUNDARY_MECHANICS=VERIFIED_FOR_INTACT_BENCHMARKS` |
| 生产 BAQP | 保持原状；PD2 位于独立候选命名空间 |

本文件随交接提交一起保存；新提交自身的 SHA 以 `git rev-parse HEAD` 和远端分支为准，避免文档自引用哈希。

必须区分两个身份：历史 `PFPD_AT1_SPECTRAL_CANDIDATE_V1` 包含 C1 BAQP、AT1、Green/SVK 谱分裂、PETSc RSLS；当前新主线是 **PD2 完整态变分 BAQP + 已验证基准的边界转换**。PD2 尚无相场耦合力，不能称为完整 PFPD 2026 复现，也不能把历史 C5/C7 的相场结论直接转移到 PD2。

## 3. 最短阅读顺序

1. [当前上下文](PROJECT_ACTIVE_CONTEXT.md)。
2. [PD2 边界回接报告](clot_clearance/verification/higher_order_baqp_boundary_reintegration_v1/REVIEW_HIGHER_ORDER_BAQP_BOUNDARY_REINTEGRATION_V1.md)，重点是 §3 闭合图、§7 全局核算、§8/13 恢复、§14 经典剪切及 §19/20 结论。
3. [PD2 变分能量—力报告](clot_clearance/verification/higher_order_baqp_energy_force_candidate_v1/REVIEW_HIGHER_ORDER_BAQP_ENERGY_FORCE_V1.md)，重点是 §3–7 的离散能量及伴随推导。
4. [PD2 RK oracle 报告](clot_clearance/verification/higher_order_baqp_rk_oracle_v1/REVIEW_HIGHER_ORDER_BAQP_RK_ORACLE_V1.md)，重点是完整九项基、矩条件、条件数和 Julia 对照。
5. 需要解释设计动机时，再读 [C1 近边界因果报告](clot_clearance/verification/baqp_near_boundary_kinematic_energy_audit_v1/REVIEW_BAQP_NEAR_BOUNDARY_KINEMATIC_ENERGY_V1.md) 与 [历史索引](HISTORICAL_INDEX.md)。

各报告保留原轮次的 `NOT_RUN` 和“下一任务”记录。它们描述当时状态；当前先后关系以本文件和最新 REVIEW 为准，不回写旧报告。

## 4. 从 C1 到 PD2 的证据链

| 阶段（均在 `clot_clearance/verification/`） | 结论与保留的限制 |
|---|---|
| `pfpd_baqp_finite_domain_boundary_consistency_v1` | 有限域仿射能量正确；几何 core 力 patch FAIL；作者代码复现该力，不等于证明物理边界正确。 |
| `pfpd_baqp_fully_closed_force_stencil_audit_v1` | N27 中心三跳图闭合，力通过；原始 Z 误差 1.13629e−12 略超 1e−12 门槛，原 NOT_ESTABLISHED 不改。 |
| `pfpd_baqp_affine_coefficient_reference_audit_v1` | 后续独立审计解释 Z 差异为重构应力浮点传播，支持 bulk 一致性；没有修改原失败门槛。 |
| `nonlocal_boundary_formulation_v1` | 仿射法向/剪切 B1 数据转换通过；直接表面载荷 B0 和单层体积约束 B2 的局限保留。 |
| `nonlocal_boundary_nonaffine_localization_v1` | 二次轴向分离 bulk 模型差与边界截断力，恢复制造目标。 |
| `nonlocal_boundary_nonaffine_shear_localization_v1` | 二次剪切验证同一局部化；目标附近存在 C1 重构/正能量偏差。 |
| `baqp_near_boundary_kinematic_energy_audit_v1` | 确认截断支持 → 非零三阶矩 → Fbar 偏差 → 键投影偏差 → Ψ⁺ 偏差。 |
| `higher_order_baqp_rk_oracle_v1` | 完整 P2 重构消除已测场偏差，运动学 VERIFIED。 |
| `higher_order_baqp_energy_force_candidate_v1` | 同一离散能量的一阶变分内力 VERIFIED。 |
| `higher_order_baqp_boundary_reintegration_v1` | 四个制造恢复与一个经典剪切控制完成，完整态边界力学 VERIFIED。 |

C1 的近边界误差并非孤立力公式实现错误。对二次位移，节点梯度误差含

$$
(\bar F_{C1}-F_{exact})_{ak}=\tfrac12 u_{a,bc}M_{3,bcd}^{(i)}(K_i^{-1})_{dk},\qquad
M_{3,bcd}^{(i)}=\sum_j\omega_{ij}\xi_b\xi_c\xi_dV_j.
$$

完整对称邻域的奇矩消失，单侧截断邻域通常不消失。既有轴向/剪切目标的 max Fbar 偏差 / 0.004 为 0.0297584，Ftilde 为 0.0294285；剪切节点 Ψ⁺ 最大偏差约 6.43%。物理梯形体积亦有次级影响，不能把不同最大范数直接相加。

不能通过单侧支持上的非负标量权重任意消除奇矩并同时保持满秩。已接受方向是完整二次 RK 导数权重，未移植 BB/OSB 表面修正，也未实施 ghost、GMLS 或重新调参。

## 5. 当前离散定义及公式来源

主验证域每轴 19 点，共 6859 点，690020 条无向键、1380040 个有向邻接。$h=62.5\,\mu m$，$\delta=251.25\,\mu m$，$\delta/h=4.02$，$L=18h=1.125\,mm$。核为 $(1-r/\delta)^2$；节点体积按内部/面/棱/角取 $h^3,h^3/2,h^3/4,h^3/8$。SVK 参数 $\mu=1000\,Pa,\lambda=0$ 为演示参数。

令 $q=\xi/\delta$，完整基为

$$Q=[q_x,q_y,q_z,q_x^2,q_y^2,q_z^2,q_xq_y,q_xq_z,q_yq_z]^T,$$
$$M_i=\sum_j\omega_{ij}Q_{ij}Q_{ij}^TV_j,\quad M_iC_i=\delta^{-1}[I_3;0],\quad
g_{ij}=\omega_{ij}V_jC_i^TQ_{ij},$$
$$F_i=I+\sum_j(u_j-u_i)\otimes g_{ij}.$$

$g$ 已含邻点体积，单位 m⁻¹，不能再乘一次 $V_j$。参考权重随变形固定；使用直接矩阵求解，无伪逆/正则化。矩条件 $\sum_jQ_{ij}g_{ij}^T=\delta^{-1}[I_3;0]$ 在局部满秩时给出二次多项式梯度再现，这是代数结论，不是稳定性定理。

来源由既有 [RK 来源清单](clot_clearance/verification/higher_order_baqp_rk_oracle_v1/references/SOURCE_PROVENANCE.json) 定位：Hillman–Pasetto–Zhou，§4.1–4.2，Eqs. 63–70、77–81；Behzadinasab–Trask–Bazilevs Part I，§2.2，Eqs. 3、5 与 Example 1。Part II 仅为未来边界/FSI 文献，本项目未采用其完整动力学算法。原文 PDF 未新增到本分支。

全部节点矩秩为 9，最大条件数 218.069；九个物理多项式最大归一化残差不超过 1.658e−14。轴向/剪切/XZ 的 PD2 max Fbar 误差 /0.004 分别约 3.332e−13、2.276e−13、3.628e−13。省略交叉项的 RK2_DIAG 在 XZ 场失败（约 0.06055），是预期负对照。Peridynamics.jl 固定版本的代表节点及权重对照通过。

当前键关联投影和完整态能量为

$$R_b=I-\xi_b\otimes\xi_b/r_b^2,\qquad
\widetilde F_b=\tfrac12(F_i+F_j)R_b+(x_j-x_i)\otimes\xi_b/r_b^2,$$
$$E=\tfrac12(F^TF-I),\quad\Psi_0=\tfrac\lambda2(\operatorname{tr}E)^2+\mu E:E,$$
$$U=\sum_{i<j}m_b\Psi_0(\widetilde F_b),\qquad m_b=2\bar w_{ij}V_iV_j,\qquad f_{int}=-\nabla_xU.$$

双端能量约定中的因子 2 与原 BAQP 一致。力来自该能量的解析伴随：逐键累积 $A_b=\tfrac12m_bP_bR_b$ 和 $t_b=m_bP_b\xi_b/r_b^2$，再把节点 $G_i=\sum_{b\ni i}A_b$ 通过 $G_ig_{ik}$ 回传到有向邻接的两端。详细符号、单位与推导见能量—力 REVIEW §3–7；不能仅替换 F 后继续使用旧 C1 Z 力式。

该实现通过 C1 回归、平移与 37° 转动客观性、当前位置总力/力矩、六组方向能量差分和 12 个选定自由度检验。C1 力回归最大归一化差 2.182e−13；方向差分固定中间步误差最大 2.306e−10。离散能量—力一致不等于连续体积分精确或稳定性已证明。

谱分裂沿用原 `svk_spectral_split.py`：$E^+=\sum_a\max(e_a,0)n_a\otimes n_a$，$\Psi^+=\lambda\langle\operatorname{tr}E\rangle_+^2/2+\mu E^+:E^+$。在当前完整态 PD2 轮次中，它只用于诊断，未形成相场退化力。

## 6. 边界回接究竟验证了什么

边界转换沿用同一替代场定义：

$$f_{BC,i}=q_{bulk}(Z_i)V_i-f_{FD,i}(x_L),\quad e_{model}=q_{bulk}+b_{local},$$
$$f_{FD}+f_{BC}+b_{local}V=e_{model}V.$$

$q_{bulk}$、$b_{local}$、$e_{model}$ 的单位为 N/m³；乘物理体积后才是 N。非仿射场的边界转换不能简单设为 $-f_{FD}$，否则会抵消真实体域响应。该分解借鉴已记录的 local-to-nonlocal 数据转换思想；BAQP 实现与本项目推导的对应关系见 [边界文献映射](clot_clearance/verification/nonlocal_boundary_formulation_v1/LITERATURE_TO_BAQP_MAPPING.md)，不宣称照搬外部论文的完整模型或定理。

实际力图必须追踪到 D3：累计节点数 D0/D1/D2/D3 = 1/257/1797/5821，最大坐标偏移 0/4/8/12h。D1 的 257 包括观察点；旧记法的 256 只数直接邻点。N19 不闭合；N25 包含所需节点却在 D3 读到半体积边界；N27 边界位于 13h，才使所需 12h 范围处于完整体积区。不能把对称中心的提前数值抵消当作闭合证明。

独立 PD2 bulk oracle 在 N27 上计算全部 19 个 Z 层，并以平移 N29 与完整 N19 选点交叉核对；没有复用 C1 的 bulk 力。扩大模板时制造场的分母 L 仍是主域的 1.125 mm。仿射 bulk 残差 /Fref 约 1e−15；轴向/剪切 max |e_model| 分别为 7.520789e−5 / 7.520065e−5 N/m³，相对体力约 2.090091e−8 / 4.229698e−8。

仿射法向、剪切合力/参考力矩/12 种虚功检查均通过 1e−6 门槛，最大全通道虚功误差约 1.53e−13。剪切六面牵引及法向小分量保留。参考坐标力矩允许非零；不能与客观性要求的当前位置全局内力矩零混同。

四个制造载荷使用 $b_{NL}=-q_{bulk}$，从非零扰动出发，复用原 TRF/LSMR 求解器：

| 场 | 机械残差 | max ΔFbar/尺度 | max ΔFtilde/尺度 | max ΔΨ⁺/Ψref |
|---|---:|---:|---:|---:|
| normal_affine | 5.191387e−10 | 1.186512e−7 | 1.295678e−7 | 2.400850e−7 |
| shear_affine | 8.658278e−11 | 2.438900e−8 | 2.435422e−8 | 2.168188e−8 |
| nonaffine_axial | 1.150552e−12 | 6.089131e−13 | 7.658356e−13 | 1.311080e−12 |
| nonaffine_shear | 5.778522e−11 | 4.489059e−10 | 5.076200e−10 | 8.738645e−10 |

机械/梯度/正能量门槛分别为 1e−8/1e−5/1e−4；仿射梯度尺度 0.003，二次场 0.004。所有门槛通过，输入外载在每次求解前后逐位一致。仿射也走同一 NL 接口，因此实际外载与存储 fBC 相差舍入量级的 −qV；原始数组未事后改写，单独核算见 `PD2_AFFINE_LOAD_ROUNDOFF_ACCOUNTING.json`。

唯一经典剪切控制改用 $b_{local}V+f_{BC}$，未加 −e_model 补偿。机械残差 3.780054e−12，目标位移偏移 1.693642e−14 m；两解响应相对模型载荷闭合误差约 1.40%，NL 位移恢复误差为经典偏移的约 2.37%。状态为 `CONSISTENT_WITH_BULK_MODEL_DIFFERENCE`，不赋予这些微小位移每位小数物理意义。解析目标自身在宽松力门槛下亦可能“收敛”，所以不能仅凭 solver success 断言精确局部解。

PD2 载荷幅值分位深度 d50/d90/d95/d99 = 0/2/2/4h，历史 C1 为 1/2/3/4h。**这是转换支持分布，不是物理边界层厚度。** C1−PD2 载荷差只说明离散表示变化，不是可唯一识别并加回模型的物理补偿力。C1 与 PD2 场质量图比较的是历史 C1 目标态和实际求解的 PD2 状态，不是两个新求解结果。

## 7. 必须保留的科学限制

- 多项式再现只对当前满秩支持及已测场成立；不规则、薄层、曲面或任意单侧点云仍需逐节点秩、条件数和解析场检查。
- PD2 消除二次重构偏差，不消除键/体积积分误差。保留 $U_{PD2}\simeq U_{exact\,kin/current\,weights}$ 与其相对连续体积分的差异，不通过改权重强迫一致。
- N19 几何 core 不是完全闭合 bulk；制造场转换目前也不是任意真实流体牵引的通用物理边界定律。
- `CURVED_CLOT_WALL_READINESS=READY` 只表示可以开展下一轮完整态曲面验证，不是曲面血栓已验证。
- PD2 phase force 未实现；PD2 稳定性、网格客观性、损伤/断裂、黏附、Ogden、真实 clot 几何、CFD、微泡、代理训练均未验证。
- 历史 C5/C7 激活应变随网格变化约 −30.2%，mesh objectivity 未建立；C5 15 个审计态稳定，C7 20 个早期接受态仅证驻定，稳定性未评估。未把它们换成 PD2 结果。

## 8. 代码导航与复现边界

以下代码路径相对于 `clot_clearance/`。依赖的生产模块、旧脚本和配置随祖先提交保留；新候选源全部纳入交接。

| 路径 | 职责 |
|---|---|
| `pd_clot/pfpd_candidate/higher_order_rk/{basis,moment_matrix,rk_weights,rk_kinematics}.py` | P2 空间、直接矩阵解、冻结权重和节点梯度 |
| 同目录 `{mechanics_candidate,variational_energy,adjoint_accumulation,variational_force}.py` | 参考数据、原键能量、统一 C1/PD2 伴随力 |
| `pd_clot/pfpd_candidate/nonlocal_boundary/mechanics_backend.py` | PD2/C1 backend、原求解器所需方向导数 callback |
| 同目录 `pd2_bulk_reference.py` | 根据真实伴随图构造选择性闭合 PD2 bulk oracle |
| 同目录 `{mechanics_adapter,localized_conversion}.py` | 可选 backend 与既有转换代数；扩展前源归档于本轮 `source_before/` |
| 同目录 `{geometry,local_surrogate,nonaffine_manufactured,diagnostics,field_bias}.py` | 物理几何、制造场、全局核算和读出 |
| `pd_clot/pfpd_candidate/{baqp,kernel,energy,svk_spectral_split}.py` | 保留的 C1 与材料/谱分裂基础 |
| `pd_clot/regularization/equilibrium.py` | 复用的 TRF/LSMR 求解器 |
| 三个 `verification/higher_order_baqp_*/` | 各轮配置、`run_audit.py`、定向 tests、图表脚本、最终 REVIEW |
| `scripts/*boundary*`、`*force_stencil*`、`*kinematic_energy*` | 先前 C1 边界与因果诊断入口 |

交接保留六张 PNG 及相应紧凑数字。九项多项式热图的 36 个类内最大值另存 `github_handoff_20261006/PD2_POLYNOMIAL_FIGURE_DATA.csv`；大型逐节点 CSV 的描述统计另存 `LARGE_TABLE_AGGREGATES.json`。这些是现有数据的摘录，不是新科学计算，也不能替代原门槛和逐节点证据。

新分支是选择性快照。未上传的 N19 梯度缓存、大型键诊断、完整逐节点表、重复图件和论文原文仍保存在服务器，路径、大小、SHA256 与替代摘要见 [本地数据清单](LOCAL_ARTIFACTS_NOT_SYNCED.md)。小 N7 缓存作为已存在测试的固定输入保留。不要为让远端测试全部通过而修改旧断言或旧清单。

服务器本轮已重新运行三个现有定向套件：**36 passed in 4.42 s**，日志在 [SERVER_TARGETED_TESTS.log](github_handoff_20261006/SERVER_TARGETED_TESTS.log)。完整套件的一些测试读取大缓存、被保护的本地证据或逐节点数组，精简分支不承诺脱离服务器运行全部 36 项。仅使用同步文件的复核命令与结果在 [HANDOFF_ENVIRONMENT.md](HANDOFF_ENVIRONMENT.md)。全量 pytest、平衡 campaign 和图件重绘均未执行。

原科学 OUTPUT/PRE/POST manifests 描述当时完整服务器产物；不能将它们理解为本次 GitHub 文件清单。当前同步范围以 [HANDOFF_SYNC_MANIFEST_20261006.md](HANDOFF_SYNC_MANIFEST_20261006.md) 为准。旧报告中的本地绝对路径及未上传图/数组链接保留原义；不要为消除这些链接而更改历史报告。

## 9. 不应重新打开的工作与历史缺失

没有新的相反证据时，不重新推翻已完成的完整九项 RK oracle、C1 变分回归或本轮完整态边界恢复；不重新把省略交叉项的空间当成完整 P2。不要仅为重复获得相同结论运行 C5/C7、Hessian、稳定性或长平衡。

旧 cleanup 造成的 123 项缺失中，恢复审计最终仍有 72 项旧字节依赖未恢复，用户已明确无备份。该事实和失败记录保留；用户后来取消其作为后续科学开发硬门，未授权篡改历史 tests/manifests，也未授权继续清理。旧全量测试记录为清理前 412 passed / 3 failed、恢复后 407 passed / 8 failed；这不是本轮重跑结果。原 3 个 staggered 路线失败仍不能改成通过。

本地未跟踪、ignore 或名为 cache 的内容，不等于可删除；已发生过 manifest 的间接依赖破坏。此次只 ADD/UPDATE/SYNC，没有清理、覆盖科学态、改变 physics 或重写历史。

## 10. 下一会话如何继续

先确认本分支 HEAD，阅读上述三个 REVIEW 和各自 FINAL_STATUS，再定位新任务真正依赖的源与输入。唯一建议下一任务是 **CURVED_CLOT_WALL_INTACT_MECHANICS_V1**：在新的明确范围下，把已验证完整态 PD2/边界接口用于曲面血栓—壁面几何，并先验证几何、支持可辨识性及力学边界。这里是建议，不是已执行结果或自动启动许可。

本轮终点为已核实的 GitHub handoff branch；不创建 PR，不继续科学开发。
