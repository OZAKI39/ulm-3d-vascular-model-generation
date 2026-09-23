# Shear-lift 文献及公式审计

日期：2026-09-23。此文件在候选公式实现之前完成。角色均相对于本项目的 lipid-shelled SonoVue；没有一篇在当前血管、壳层和间隙下构成已验证的直接模型。

## A. 小刚性球

**Saffman, P. G. (1965), The lift on a small sphere in a slow shear flow.** DOI [10.1017/S0022112065000824](https://doi.org/10.1017/S0022112065000824)。原文从 [Caltech](https://authors.library.caltech.edu/records/v662k-mg234) 下载并检查 pp.385–387。对象：无滑移刚性球，半径 a，均匀简单剪切，滑移平行流线。原摘要 81.2 不是应实现的系数；必须连同 1968 勘误理解。半径定义 Re_p=a|u−v|/ν、Re_G=a²|G|/ν，要求 Re_p≪√Re_G≪1，旋转惯性也小。无限域；不覆盖近壁。刚性球适用，洁净/污染/壳层气泡不能直接套用。项目角色：**ORDER_OF_MAGNITUDE_ONLY**。

**Saffman (1968), The lift on a small sphere in a slow shear flow – Corrigendum.** DOI [10.1017/S0022112068999990](https://doi.org/10.1017/S0022112068999990)，JFM 31:624。出版社条目核对成功，全文 PDF 请求返回网页，未声称直接读到勘误正文。修正后的系数由下列可读权威原研究 Eq.(1.1) 和独立系数转换交叉核对。对象、a、Re、无限域及界面限制同上。角色：**ORDER_OF_MAGNITUDE_ONLY**。

**Joseph, D. D. & Ocando, D. (2002), Slip velocity and lift.** DOI [10.1017/S0022112001007145](https://doi.org/10.1017/S0022112001007145)。[作者稿](https://dept.aem.umn.edu/~./faculty/joseph/PL-correlations/docs-ln/slip_velocity_and_lift.pdf) p.2 Eq.(1.1) 明确写出 6.46 μ a² U_s√(γ/ν)，a 是半径、U_s=U_f−U_p。此处只采用其重述的经典领先项，不采用该文圆柱 DNS 作为本项目球体结果。简单剪切的 γ=|du_t/dn|，并非任意三维应变不变量。经典项不含壁面，适用低 Re 的刚性球假设；不构成气泡/壳层定律。角色：**ORDER_OF_MAGNITUDE_ONLY**。

## B. 壁面与方向

**Cherukat, P. & McLaughlin, J. B. (1994), The inertial lift on a rigid sphere in a linear shear flow field near a flat wall.** DOI [10.1017/S0022112094004015](https://doi.org/10.1017/S0022112094004015)，另有 1995 勘误 [10.1017/S0022112095000590](https://doi.org/10.1017/S0022112095000590)。对象：刚性球/单一平壁/线性剪切；结构为 F=ρa²U_s² I(a/l,aG/U_s)，I 随间距和旋转条件变化，不是通用 6.46。a 半径，l 球心到壁面距离；壁面位于扰动内区，l≪min(ν/|U_s|,√(ν/|G|))，可研究 h≪a。与无限域经典式不同。此处仅核对摘要及 Shi 等重述，未实现其有勘误的拟合系数。不能直接覆盖曲壁多分叉或壳层。角色：**NOT_DIRECTLY_APPLICABLE**。

**Shi, P., Rzehak, R., Lucas, D. & Magnaudet, J. (2021), Drag and lift forces on a rigid sphere immersed in a wall-bounded linear shear flow.** DOI [10.1103/PhysRevFluids.6.104309](https://doi.org/10.1103/PhysRevFluids.6.104309)，[作者原文 arXiv:2108.07196](https://arxiv.org/html/2108.07196)。全文 Eq.(1),(5),(7),(14)–(20)。d 是直径，Re=|U_s|d/ν，Sr=Gd/U_s，L_R=2l/d。无限域 Eq.(7) 的 C_L=(18/π²)sgn(Sr)εJ(ε)，ε=√(|Sr|/Re)，J(∞)≈2.254；力的归一化为 πd²ρU_s²/8，转成半径系数 (9/π)2.254≈6.457，与 6.46 舍入一致。方向 U_s×curl(u)。壁面扰动长度 L_G=√(ν/|G|)，经典无限域需要 l/L_G≫1；小 h/a 以外也可能受壁影响。DNS 参数 Re=0.1–250，L_R=1.5–8，|Sr|≤0.5，不能外推为本项目已验证模型。刚性无滑移球，非气泡/壳层。角色：**ORDER_OF_MAGNITUDE_ONLY**（尺度和方向审计；不执行壁面拟合）。

## C. 洁净气泡、污染界面和脂质壳层

**Legendre, D. & Magnaudet, J. (1997), A note on the lift force on a spherical bubble or drop in a low-Reynolds-number shear flow.** DOI [10.1063/1.869466](https://doi.org/10.1063/1.869466)。出版社全文受访问限制，结果与定义由作者 1998 原文交叉核对：同尺寸/同滑移/同简单剪切下，洁净气泡低 Re 领先升力是刚性球的 4/9，不能误用 2/3。等效半径系数 (4/9)6.46；低 Re、无限域、球形、可移动剪切自由界面。不含壁/脂质膜/吸附输运。角色：**ORDER_OF_MAGNITUDE_ONLY**，不在计算中替换成 SonoVue 系数。

**Legendre & Magnaudet (1998), The lift force on a spherical bubble in a viscous linear shear flow.** DOI [10.1017/S0022112098001621](https://doi.org/10.1017/S0022112098001621)。[作者上传全文](https://www.researchgate.net/publication/231943238_The_lift_force_on_a_spherical_bubble_in_a_viscous_linear_shear_flow) pp.109–111，Eq.(16)–(18)：C_L=6J(ε)/(π²√(Re Sr))，J≈2.255/(1+0.2ε⁻²)^(3/2)，体积型升力归一化，Re 用直径，Sr=|ω|d/|U_s|。洁净球形气泡，外域无壁；数值范围 Re=0.1–500、Sr≤1，低 Re 理论另有渐近限制。不能用于已污染或壳层气泡；也不能把其系数和刚性球面积型 C_L 混用。角色：**ORDER_OF_MAGNITUDE_ONLY**。

**Fukuta, M., Takagi, S. & Matsumoto, Y. (2008), Effects of Surfactant on the Lift Force Acting on a Bubble in a Shear Flow.** DOI [10.1299/kikaib.74.1679](https://doi.org/10.1299/kikaib.74.1679)，[出版社原文](https://www.jstage.jst.go.jp/article/kikaib1979/74/744/74_744_1679/_article/-char/en)。并列研究 Numerical study on the shear-induced lift force acting on a spherical bubble in aqueous surfactant solutions，DOI [10.1063/1.2911040](https://doi.org/10.1063/1.2911040)（后者全文未取得）。对象：污染球形气泡；求解三维界面/液相浓度及 Langmuir 吸脱附，切向应力跳跃由 ∇_sσ 决定，而非固定 Saffman 系数。不存在可从摘要获得的通用粒径/系数/Re 公式，不能杜撰；不实施数值闭合。剪切流中表面浓度分布会改变升力；未验证曲壁、低 Re SonoVue 脂质膜。角色：**NOT_DIRECTLY_APPLICABLE**，用于说明界面条件必须单独建模。

**Marmottant, P. et al. (2005), A model for large amplitude oscillations of coated bubbles accounting for buckling and rupture.** DOI [10.1121/1.2109427](https://doi.org/10.1121/1.2109427)，[作者论文](https://liphy-annuaire.univ-grenoble-alpes.fr/pages_personnelles/philippe_marmottant/Publications_files/Marmottant2005_JASACoatedBubble.pdf)。对象：脂质单层包覆微泡径向振荡；σ(R) 分为屈曲、弹性和破裂区，R 是瞬时半径。没有横向剪切升力系数；没有可移用的升力 Re 或近壁适用范围。此文不能补全 lateral lift。角色：**NOT_DIRECTLY_APPLICABLE**（壳层力学背景）。SonoVue 的磷脂壳层另由 [EMA 产品评估](https://www.ema.europa.eu/en/documents/variation-report/sonovue-h-c-303-x-0034-g-epar-assessment-report-variation_en.pdf) 核实。

## 本轮实现约定（只读后处理）

采用 F_c=6.46 μ a² |u−v| √(γ_E/ν)，字段 `CANDIDATE_SAFFMAN_MAGNITUDE`，γ_E=√(2E:E)。这是将简单剪切标量推广到三维局部应变的**量级代理**，不是三维严格 Saffman 定律，也不是经证明的上界。额外保存 curl 与简单剪切偏离指标 |γ_E−|curl||/γ_E，候选方向仅以 (u−v)×curl(u) 示意，零叉积方向为未定义，不虚构箭头。该方向不使标量推广获得理论有效性。

半径 convention：Re_p=a|u−v|/ν，Re_G=a²γ_E/ν；若换直径，Re_p,d=2Re_p，Re_G,d=4Re_G。6.46 的直径形式系数为 1.615。记录连续裕度 m_slip=Re_p/√Re_G、m_shear=√Re_G、m_wall=(a+h)/√(ν/γ_E)，而不把任意 0.1 阈值称作验证边界。文献的 ≪ / ≫ 是渐近关系，没有通用有限截断值。若 m_wall<1，则连无限域的宽松必要次序都失败；反之也不能宣布有效。给出描述性 m_slip,m_shear<0.1 且 m_wall>10 的筛查灵敏度（0.03/0.1/0.3），明确不代表科学有效阈值。

壁面分层使用现有 V1 模型状态，并逐状态保留 h/a；V1 的 0.01/0.05 只决定既有 lubrication 开关，不是 Saffman 有效性界限。等价润滑贡献为 |ζ_n(v·n)|，ζ_n=w(h/a)6πμa²/max(h,h_lower)，不是另外积分的惯性力。所有系数、速度、间隙必须在同一个正式求解采样位置求值。
