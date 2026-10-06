# 有限域 BAQP 边界理论与实现对照

任务：`PFPD_BAQP_FINITE_DOMAIN_BOUNDARY_CONSISTENCY_V1`。仅讨论 s=0 力学；下文区分原始来源、本轮推导与尚待验证的物理判断。没有据此修改模型。

## 1. PFPD 原文对边界和内力的实际陈述

来源：[用户本地 2604.09215v1.pdf](/home/lzy/projects/2604.09215v1.pdf)，PDF 第 3–5、20–21 页；公开版本索引：[arXiv:2604.09215v1](https://arxiv.org/abs/2604.09215v1)。直接读取本地 PDF，未重新审计 AT1。

- §2.1、式 (2.1)：family 是物体内与当前点相距不超过 δ 的点，不自动包含物体外虚点。
- 式 (2.2) 后的正文：位移或外力密度施加在有限厚度的材料点层；这与单纯的经典连续表面条件有区别。
- 式 (2.4)–(2.10)：K 逆保证仿射运动学完整性；该陈述不等于有限域离散内力逐点为零。
- 式 (2.12)–(2.19)：端点对称的 ω̄、应力积分 Z 与力状态。打印版 (2.19) 的 direct 分母写成 |ξ|；固定作者代码使用 |ξ|²，本地既有实现也使用后者，本轮保持不变。量纲和独立能量差分支持当前代码形式。
- §5.1 的动态裂纹示例给顶底速度条件，并提到一层不失效点；没有由此给出任意静态 BAQP cube 的通用 collar 厚度。
- §5.2、式 (5.2)–(5.3) 明确在一层点上施加 ±σ/Δx 的体力密度。这里用的是其点云体积约定；本项目端面节点为半权重，不能原样套 σ/h 后假定总力相同。

因此，“原文要求体积处理”成立；“所有问题必须约束恰好一个 δ 厚度”并不是该文提供的结论。

## 2. 固定作者实现：边界 benchmark 与 BAQP 力

仓库：[kaipartmann/Peridynamics.jl](https://github.com/kaipartmann/Peridynamics.jl/tree/b6e86563845a83b5cb06c805af19a817bde6b106)，commit `b6e86563845a83b5cb06c805af19a817bde6b106`；本地 checkout 的实际 HEAD 已核对。

[uniform-tension test](https://github.com/kaipartmann/Peridynamics.jl/blob/b6e86563845a83b5cb06c805af19a817bde6b106/test/simulation/test_analytic_solutions.jl#L18-L71)：长度 l=1、截面 0.1×0.1；盒子长度扩展为 l+3Δx，中心移至 −1.5Δx。`x < −l/2` 选左侧三个额外点层；`x > l/2−Δx` 选右侧一个加载层。horizon=3.015Δx。左端按求解器采用三方向零速度或零位移；右端施加 `F / sum(V_right)` 的轴向力密度。Δx=1/30 的普通测试使用 CMaterial 或带表面修正的 BBMaterial，**不是 RKCMaterial 的 BAQP 边界验证**。

[boundary_conditions.jl](https://github.com/kaipartmann/Peridynamics.jl/blob/b6e86563845a83b5cb06c805af19a817bde6b106/src/physics/boundary_conditions.jl) 的实现按 point set 写入分量值；它不会为任意几何自动推导所需层厚。上面的额外左层与右层总力分配是具体 benchmark 的定义。

[rk_correspondence.jl](https://github.com/kaipartmann/Peridynamics.jl/blob/b6e86563845a83b5cb06c805af19a817bde6b106/src/physics/rk_correspondence.jl#L414-L465)：`:C1` 的 gradient weight 包含邻点体积，相当于本项目 BᵢⱼVⱼ。第 583–615 行用对称端点权重计算 Z；第 618–639 行的有向 direct 使用 ω/ω0ᵢ。两个方向都装配后，direct 的系数是 ω(1/ω0ᵢ+1/ω0ⱼ)=2ω̄，与本地无向 bond 形式一致。不要逐条把两种 force-state 中间变量当成相同对象。

本轮直接调用未修改的该 package，对 L=0.75 mm、同一坐标/节点体积、quadratic kernel、SVK μ=1000、λ=0 做一次仿射力计算。全场相对差 5.86×10⁻¹⁵，完整 family core 内相对差 2.32×10⁻¹⁴；没有时间推进。证据：[author_L075/COMPARISON.json](author_L075/COMPARISON.json)。这反驳本地移植失配，不能证明其有限域结果就是正确的物理牵引。

## 3. Peridigm 的点集层厚是显式输入

固定：[peridigm/peridigm@f25c4928e850b2aa0012bc5d3b8b1ec5c78c8903](https://github.com/peridigm/peridigm/tree/f25c4928e850b2aa0012bc5d3b8b1ec5c78c8903)。本轮只读两个固定版本源文件，没有编译或运行 Peridigm。

[MeshConverter 第 147–205 行](https://github.com/peridigm/peridigm/blob/f25c4928e850b2aa0012bc5d3b8b1ec5c78c8903/src/mesh_converter/Peridigm_MeshConverter_Main.cpp#L147-L205) 读取 `Boundary Layer Thickness`，以包围盒的 `x<x_min+t`、`x>x_max−t` 等条件生成六个 node set。厚度 t 来自输入，不由此自动等于 horizon。

[BoundaryCondition.cpp](https://github.com/peridigm/peridigm/blob/f25c4928e850b2aa0012bc5d3b8b1ec5c78c8903/src/core/Peridigm_BoundaryCondition.cpp#L157-L189) 在指定 node set 上施加 Dirichlet 分量；其 Neumann 类也使用指定点集。源码支持显式体积层表达，但不是 BAQP 牵引正确性的证明，也不能据此确定本项目应选的 t。

## 4. 非局部 volume constraint 的数学边界

D’Elia–Yu：[On the prescription of boundary conditions for nonlocal Poisson’s and peridynamics models](https://arxiv.org/pdf/2107.04450)，§2–4。方法先求已知局部问题，再将该解的位移或非局部通量延拓到 interaction domain，构造 DtD 或 DtN 体积约束。

关键细节：Poisson 的外部 interaction layer 厚 δ；文中 LPS 因双重积分，式 (2.13) 的厚度是 **2δ**。因此作用半径要从完整算子确定。其 Theorem 4.1/4.2 的 O(δ²) 结论依赖平滑局部解及对应 Poisson/LPS 算子，不是当前非线性 BAQP 的收敛定理。

## 5. 本轮从既有离散式推导的三个区别

**全局归一化与逐点平衡。** 对均匀 W、对称 ω 和唯一无向键：

\[
\sum_{i<j}2\bar\omega_{ij}V_iV_j
=\sum_i\frac{V_i}{\omega_{0i}}\sum_{j\in H_i}\omega_{ij}V_j
=\sum_iV_i.
\]

它保证 U_aff=WΣV，只沿仿射族约束能量；不保证对任意单粒子虚位移的导数为零。

**一个 family 与完整内力模板。** 仿射 P 恒定时，节点内力仍读取邻居的 ω0ⱼ、Kⱼ⁻¹ 和 Zⱼ；Zⱼ 又读取其邻居的 ω0。由现有代码可得保守几何系数依赖链

\[
i\longrightarrow j\longrightarrow k\longrightarrow \ell .
\]

这是最长 3 个 δ 邻接的上界，**不是已验证的最小 collar 厚度或建议 BC**。dᵢ≥δ 只保证自身 family；dᵢ≥3δ 在规则无限晶格条件下是消除这条边界系数依赖的充分几何条件。当前 L/2=1.493δ 的大试件未隔离这样的中心区；27 个自身 family 完整的点，没有一个满足“所有邻居 family 也完整”。

**半域合力与仿射虚功。** 本轮严格分开比较

\[
R_+=\sum_{Z_i>0}f^{ext,req}_{iz},\qquad
-\sum_i f^{int}_{iz}(Z_i+L/2)=P_{zz}L^3.
\]

对称情况下定义有符号力的重心 z̄_F，则 2R₊z̄_F=R_cont L。若力分布在域内多层，R₊通常不由仿射能量恒等式强制等于 R_cont。相应的非单调累计曲线也不能当成正概率分布；本轮保留 overshoot，不裁剪负层力。

## 6. 本轮可接受的结论边界

文献/实现标签为 `MIXED_BY_METHOD`：体积点层、外部延拓区与具体一层力密度 benchmark 各有定义。没有一个来源证明“直接给当前 C5/C7 加 δ 厚约束”就是正确修复。几何 full-support patch test 已按用户原阈值标 FAIL；作者实现相同、能量差分一致，使问题进一步定位到有限域算子及其边界表达，而非本地程序差错。由于力模板尚未与边界系数完全隔离，本轮不把它判为无限域 bulk 缺陷，也不运行 BC fix。
