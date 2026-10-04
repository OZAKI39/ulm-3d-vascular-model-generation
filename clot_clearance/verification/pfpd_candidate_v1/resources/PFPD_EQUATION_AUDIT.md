# PFPD v1 原图与公式审计（编写候选力学代码之前）

原文：Partmann, Wieners, Ortiz, Weinberg, *Phase-Field Peridynamics*, arXiv:2604.09215v1。
本地 PDF 身份见 PAPER_IDENTITY.json。使用本地解包 poppler-utils 24.02.0 的 pdftotext -layout 提取，并渲染、逐页目视核对 pages/page_04、05、06、07、08、09、13、14、19.png；不是只据 HTML。完整提取文本保存在 paper_pdftotext.txt。公开资源检索单独保存原始 API 响应及全部 refs 搜索证据。

## 逐项核对

记 ξ=Xj−Xi，η=xj−xi，r=|ξ|；所有积分的离散项必须带邻居体积 Vj。

|公式|PDF 页|原图核对结果与实现边界|
|---|---:|---|
|2.10|4|Ftilde=Fav+(η−Fav ξ)⊗ξ/r²，沿键方向精确重现 η。|
|2.11|4|Fav=(Fbar_i+Fbar_j)/2。|
|2.12|5|Ψbar_i=Σj wbar_ij Ψ(Ftilde_ij)Vj。全局能量再乘 Vi。|
|2.13|5|wbar_ij=(ωij/ω0i+ωji/ω0j)/2；端点归一化后对称。|
|2.14|5|ω0i=Σj ωij Vj。边界处随 family 截断而不同；不假设为 1。|
|2.17|5|Ptilde=P(Ftilde)。固定规定 s 时为 g(s)P0。|
|2.18|5|Zi=Σj wbar_ij Ptilde_ij(I−ξ⊗ξ/r²)Vj。|
|2.19|5|原图确实写成 wbar Ptilde ξ/r+Zi Bij。此处 r 未平方，量纲不一致；算法 2 第 20 行重复。见下方作者源码核对，不能照印刷式声称通过能量一致性。|
|3.4|6|动能加有向双积分 wbar[gΨ0+Gc γ/(c0δ)]，没有 1/2。|
|3.6|7|键裂纹项 Gc γ(s)/(c0δ)，量纲为 Pa。|
|3.8|7|γ(s)=s，明示 AT1。|
|3.9|7|g(s)=(1−s)²。|
|3.11|7|s≤sc 时 h=1；否则 h=((1−s)/(1−sc))²。sc=1 明示禁用运动学退化（含 s=1），sc=0 时 h=g。|
|3.13|8|Ki=Σj ωij h(sij) ξ⊗ξ Vj。|
|3.14|8|Bij=ωij h(sij) Ki^-1 ξ。|
|3.15|8|Eb=gΨ0(Ftilde)+裂纹项。|
|3.16|8|写成 g'Ψ0+Gc γ'/(c0δ)=0；在 Ftilde 固定时才是此局部导数。h 改变重构时还需要审计对邻接键能量的依赖，不能自行扩展热力学结论。|
|3.17|8|−2(1−s)Ψ0+Gc/(c0δ)=0。|
|3.18|8|Yc=Gc/(2c0δ)，并声明阈值以下 s=0。|
|3.19|8|Y=Ψ0^+(Ftilde)，正能量分裂的具体选择不能与应力驱动混同。|
|3.23|9|Y=<σ1(Ftilde)>_+²/(2E)，最大主 Cauchy 应力的另一驱动定义。|
|3.24|9|历史 Yhistory(t)=max_{t'≤t}Y(t')。|
|3.25|9|s=min(1,Yhistory/(Yhistory+Yc))；与上面的 AT1 驻值及阈值声明不相容。|
|3.47|13|c0=2∫0¹ fω(ξ)dξ，2 来自裂纹上下两侧。|
|3.59|14|c0=∫ω(ρδ)ρ³dρ/[2∫ω(ρδ)ρ²dρ]。|
|Algorithm 2|19|先按旧 s 更新形状矩阵/梯度，再算 Fbar、Ftilde、驱动、历史、s、gP、Z 和力。第16行采用 Eq3.25；第20行重复 r 而非 r²。伪代码逐点初始化力再散加邻点可能覆盖已累加力，实现必须统一清零后累加。尚不能作为无歧义完整演化算法。|

## 两个不同的问题，分开处理

1. **损伤演化硬门**：AT1 固定 F 的约束极小值为 max(0,1−Yc/Ψ0)，而 Eq3.25 在任意正驱动下均给出正损伤。后者与 γ=s² 的驻值代数一致。这只是公式关系核对，不是选择模型。完整 evolution 的状态由 FORMULA_CONSISTENCY_GATE.json 决定；未解决则不实现历史推进/损伤演化 API。
2. **力公式排印问题**：固定官方仓库 commit b6e86563845a83b5cb06c805af19a817bde6b106，src/physics/rk_correspondence.jl:618–638 使用 ξ/(L*L)；651–658 给出与 Eq2.10 一致的键梯度与横向投影。候选采用这一作者公开实现中的平方分母，并独立有限差分验证。此处依据明确，不把它当作损伤演化歧义已解决。

## 作者实现与本轮离散化的对应

官方公开 BAQP 是 RKCMaterial(monomial=:C1)，不是 Chen–Spencer 相交邻域的 BACMaterial。源码 rkc_weights!:412–460、rkc_defgrad!:536–552、rkc_stress_integral!:583–615、rkc_force_density!:618–641、bond_avg:651–655、strain_energy_density_point!:737–764 是外部基准。

官方点能量使用 ω/ω0i，候选使用论文对称 wbar；全局总能量对于对称键本构和 s 完全等价。官方直接力项也是单端归一化，两个方向累加后等价于 2*wbar。逐点能量因此不可直接当成同一定义比较。本轮比较全局能量、节点力、Fbar 与键 Ftilde。

官方当前代码的 weighted_volume 在运动学权重改变后也改变；论文 Eq2.14 的 ω0 未包含 h。本轮只对完整无损状态使用官方 oracle。规定损伤的 B/C/D 使用固定完好 ω0，不声称与作者未公开 PFPD 演化实现等价。

一个无向物理键只保存一次 s。另建端点有向索引（物理键 ID、起终点、符号）来解释 family，不复制独立物理状态。因为 Ψ、wbar、s 在交换端点时对称：

E = Σ_{i<j} 2 Vi Vj wbar_ij g(sij) Ψ0(Ftilde_ij)。

裂纹能量同理：EΓ=Σ_{i<j} 2 Vi Vj wbar_ij Gc γ(sij)/(c0δ)。没有再除以 2。

候选节点力使用有向状态的端点差乘 ViVj，直接项为 2wbar Ptilde ξ/r²，非均匀项为 Zi Bij−Zj Bji。这是上述离散能量在固定 s 下的负梯度。秩不足时明确报错，不补造失效/脱离算法。

Saint-Venant–Kirchhoff 的能量和 P 使用官方 src/physics/constitutive_models.jl:34–48：E=(FᵀF−I)/2，Ψ=λ(trE)²/2+μE:E，P=F[λ tr(E)I+2μE]。通过该方法核对后，另以当前项目的可压缩 Neo-Hookean 做兼容性诊断；不更改原材料。

## 归一化手算

二次核：∫0¹(1−ρ)²ρ³dρ=1/4−2/5+1/6=1/60；2∫0¹(1−ρ)²ρ²dρ=2(1/3−1/2+1/5)=1/15；c0=1/4。
常数核 c0=(1/4)/(2/3)=3/8；线性核 c0=(1/20)/(1/6)=3/10。
有限 coupon 存在自由边界偏差，固定 δ/Δx 加密不会消除无量纲晶格求积误差。必须分别报告 spacing 与 horizon-ratio 的作用，不能把有限网格偏差改写为材料 Gc。
