# Foster 2011 与当前 LPS 的能量失效审查

已读取用户全文 `sources/FOSTER2011.pdf`，而非只读摘要。[原论文 DOI](https://doi.org/10.1615/IntJMultCompEng.2011002407)。以下公式按本项目记号重写。

Section 2、Eq.8–12 使用 shape tensor / deformation-gradient correspondence：Kshape=∫ω ξ⊗ξ dV，Fhat=(∫ωY⊗ξ dV)Kshape⁻¹，T=ω P Kshape⁻¹ξ，后续例子将率相关 J2 塑性应力接入。它与当前 ordinary、central-force LPS 的 dilatation/deviatoric extension 本构不同。原文声称 Kshape 正定依赖充分邻域；不能把这句话当成任意断裂点云都可逆的保证。

Section 3 Eq.13 定义路径功密度 w_ij=∫(T_i⟨ξ⟩−T_j⟨−ξ⟩)·dη。状态力/dual density 单位 N/m⁶，w 单位 J/m⁶。实际等体积离散相互作用功是 W_ij=V_iV_j w_ij，单位 J。以无向 i<j 键统计一次；若遍历两个有向状态则要 1/2，不能双计数。原文 Eq.15 的连续体极限被作者明确称为假说、未正式证明，不能升级成当前离散 LPS 的既有定理。

Eq.16 在平面裂纹两侧球冠上积分一常数临界功密度，几何积分等于 πδ⁴/4，从而 Eq.17 得到 w_c=4G_c/(πδ⁴)。这里 G_c 为 J/m²；不是把宏观 J/m² 直接存为每键 J 或 J/m⁶。每条离散键跨裂纹面的体积权重、方向分布和边界截断必须重新核对。δ 的解释及原文的过程区讨论属于其模型假定，不能用当前数值 horizon 自动宣称生物过程区已匹配。

Eq.18 以累计功超阈值令 influence=0；正文要求不可逆断键，并在后续步从 Fhat/Kshape 中去掉相互作用。图4及脚注4说明率相关路径可以在不同位移达到同一功：超弹性时功可恢复，有耗散模型时部分做功为永久耗散。当前 LPS 不含该塑性/粘弹耗散，数值阻尼不是其替代品。

当前锁定均匀等体积 LPS：r_ij=|X_j−X_i|，h_ij 为原生边缘体积 taper，e_ij=|x_j−x_i|−r_ij，m_i=Σ_reference h_ij r_ij V，初始缺口键已排除，后续 m 固定。

θ_i=(3V/m_i)Σ_alive h_ij e_ij。

U=(V/2)(K−5G/3)Σ_i θ_i² +(15G V²/2)Σ_{i<j,alive}(h_ij/r_ij)(1/m_i+1/m_j)e_ij²。

相互作用力向量（作用于较小 ID 端、向另一端）为 V² h_ij[(3K−5G)(θ_i/m_i+θ_j/m_j)+15G(1/m_i+1/m_j)e_ij/r_ij] Y_ij/|Y_ij|。无损伤、固定参考拓扑下，Σ_{i<j} f_ij·d(x_j−x_i)=dU。该结论来自同一势能的导数，并已与锁定独立 C++ evaluator 比较；不能用旧原生 PE 变量代替此完整势能。

θ_i² 含多键交叉项。单键路径功受其它存活键状态影响；沿闭合的多分量变形路径，个别键的功可以不归零，而总和为零。即使强行把总 U 拆成逐键条目，分配 θ² 交叉项的规则也不唯一。断一条键会改变相邻节点 θ，进而改变其它键贡献；能量释放 U(x,old)−U(x,new) 是全局重评值，不等同那条键的历史累计功。

因此本轮证明的是“有可计算、可核验的逐相互作用路径功”，不是“每条键具有独立唯一可恢复的弹簧能”。闭合路径证据和原 M1/M2 的数值阻尼账本见 JSON 与图 A3。Foster 原始论文不能替当前 LPS 完成这一步移植推导。

`ENERGY_FAILURE_PORTABILITY=NOT_ESTABLISHED`。本轮不设置 5.90、6–8 或 18.27 J/m² 的生产阈值，不改 s00，不启用能量断键或软化。下一次独立研究若继续，须先定义功变量的历史更新、负功/卸载、破裂释放、交叉项转移和材料耗散，再通过单裂面能量/尺寸/离散与实验验证；本轮到此停止。
