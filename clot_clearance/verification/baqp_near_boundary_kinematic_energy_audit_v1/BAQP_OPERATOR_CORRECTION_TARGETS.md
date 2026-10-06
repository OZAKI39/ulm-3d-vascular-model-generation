# BAQP operator correction targets — design only

状态：READY_FOR_BAQP_SPECIFIC_DESIGN。没有选择或实现算子修正。

本轮证据及适用范围见 [REVIEW](REVIEW_BAQP_NEAR_BOUNDARY_KINEMATIC_ENERGY_V1.md)。以当前 quadratic kernel、固定 h/δ、两种 Z-Z 二次场为证据；任意 Hessian、非规则点云及曲面边界尚须另行设计验证。

## Target algebra and invariants

当前 $\Delta\bar F=\tfrac12HM_3K^{-1}$。保持当前线性重构形式时，对任意二次 H 消除误差需要合适的三阶矩目标；仅将 K 调成连续 K 不会消除 M3。二阶一致性与可逆性需同时满足，不能通过削弱法向支持“制造”零误差。

| Obligation | Future design acceptance evidence |
|---|---|
| Affine exactness | 对任意常数 F 验证重构精确，不限于一次轴向拉伸 |
| Shape tensor invertibility | 保持 K 满秩，分析条件数与接近退化点云的行为 |
| Quadratic / odd-moment consistency | 明确目标是完整有限 horizon 算子，推导一般 Hessian 的误差；区分逐行 G 与完整 M3 条件 |
| Pair reciprocity | 端点重编号及键反向不改变标量总能量与一致的成对力分配；不能假设节点权重各自拟合自然互易 |
| One energy and first variation | 修正后的依赖关系须进入同一总能量，完整变分得到所有力项；若权重依赖构形须包含其导数 |
| Physical horizon | 不通过扩大/改变 δ 或 kernel 来掩盖原问题；未来设计必须明确支持和物理长度 |
| Boundary-data conversion | 载荷转换与表面算子分别定义；若改变算子，重新验证转换与该能量的匹配，不直接沿用旧数值结果 |
| Irregular clouds | 验证体积、邻域各向异性、曲率和不规则布点的可解性；正则化不得未经说明改变物理 |
| No activation fitting | 不以 C5/C7 首激活位置/应变或旧热点为权重优化目标 |
| Intact regression and energy quadrature | 仿射、二次与运动学—权重能量分解分别验收；U1 不是替代模型 |

## Feasibility constraint for nonnegative scalar weights

设平直边界法向 n 的保留键都满足 $\xi_n\ge0$，有效标量权重 $a_j=\omega_j V_j c_j\ge0$。若 $K_{nn}=\sum a_j\xi_n^2>0$，则至少一项满足 $a_j>0,\xi_n>0$，因此

$$M_{3,nnn}=\sum_j a_j\xi_n^3>0.$$

反向单侧家族同理为负。若强制该三阶矩为零，则在这些条件下只能令所有非零法向键权重为零，导致法向 K 退化。因此“现有单侧家族 + 非负标量权重 + K 可逆 + 所有二次奇矩消失”不能同时满足。对于当前纯法向二阶导数，M3_nn: K⁻¹=0 等价于 M3_nn:=0（K 可逆），也受此约束。

这是本轮从有限和直接推出的必要性结果，未声称 Shojaei 论文给出此 BAQP 定理。未来可比较更丰富的重构空间或显式高阶项等数学途径的可行性，但不得由此默认采用有符号权重、ghost、扩大 horizon 或任何未审核方案。对称完整 lattice 外延仅是本轮参考，不是 ghost 实现。

## Relation to Shojaei PD-to-PD matching

Shojaei 等的完整 horizon 目标匹配思想提供设计方向，BB/OSB 权重与距离幂次服务于各自算子。BAQP 的目标需从 K、M3、Ftilde、能量及其一阶变分重新推导。原始来源和具体位置见 [LITERATURE_PROVENANCE.json](LITERATURE_PROVENANCE.json) 及 [原论文](https://doi.org/10.1002/nme.70385)。未调用 perifit，没有验证其权重适用于 BAQP。

## One next task

BAQP_SURFACE_OPERATOR_CORRECTION_DESIGN_V1：给出满足上述义务且数学可行的 BAQP 专属候选设计和独立验收计划，交人工审核。当前不执行，不启动 activation、fracture、clot-wall 或 CFD。
