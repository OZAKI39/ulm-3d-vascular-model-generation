# AT1 能量与开源交叉检查（候选编码之前）

本轮方法名 PFPD_AT1_CONSISTENT_CANDIDATE。用户明确选择 γ=s、g=(1−s)² 和不可逆约束；它不是原文 Eq3.25 的精确复现。最新作者资源检查见 AUTHOR_UPDATE_CHECK.json，未找到改写该选择的新版 PFPD 实现。

## 主要依据及代数

本地 PFPD PDF 与上一轮 SHA 相同；已复核上一轮原图/逐式审计对应的 Eq3.8、3.9、3.17、3.18。固定 h=1 和 x 时，BAQP 键梯度与相场无关。纯拉伸取 Y=Ψ0，单键密度 e=(1−s)²Y+2Yc*s，q=−2(1−s)Y+2Yc，二阶导数 2Y≥0。Y≤Yc 的处女键最小值在 0；Y>Yc 的驻值为 1−Yc/Y。限制 s≥s_old 后取 max(s_old,s_trial)，即在 [s_old,1] 上最小化。

该结论只用于规定纯拉伸：不从有限变形中自行构造谱分裂、体积/偏差分裂或主应力驱动。每次耦合相场更新前检查所有键 Green 应变的最小特征值 ≥ 预声明 −1e-9（舍入/平衡误差容限）；压缩或混合状态则停止。此特征值只是适用域检查，不参与 Y 的能量分裂。

全局积分使用已有 BAQP 的无向键系数 `2*wbar*Vi*Vj`。c0 调用既有 normalization.c0('quadratic')，Yc=Gc/(2*c0*delta)，不在新模块复制 0.25。

## PhaseFieldX

官方：https://github.com/CastillonMiguel/phasefieldx
固定提交：`a9714c122497f8829860efc68d776e4324d48eca`。MIT，完整许可证和 SHA 已保存。

读取 `src/phasefieldx/Element/Phase_Field/geometric_crack.py:7–89`：AT1 几何函数为 phi、导数 1，FEM 归一化系数 cw=8/3。读取 `Element/Phase_Field_Fracture/g_degradation_functions.py:62–95`：二次 g=(1−phi)²。只通过文件级 import 调用这些无 FEM 依赖的公开纯函数，构造零梯度标量能量并交给 SciPy bounded minimizer，作为外部函数交叉检查；不复制函数、不安装 FEM solver。

换算关系：FEM 局部裂纹系数为 Gc_F/(cw*ell)，PD 局部系数为 Gc_PD/(c0*delta)。本轮取 ell=delta、Gc_F=Gc_PD*cw/c0，使局部系数恰好相同；这不是两个完整空间模型的等价性。不能将 cw=8/3 直接当成 PD 的 c0，也不把 Gc_F 称为另一种血栓实测值。

检查了 fracture solver 的 history/penalty/energy 版本：存在默认 AT2 的路径，不将这些默认 solver 当作本轮 AT1 数值 oracle。不采用其罚函数、粘性、残余刚度、空间离散或张压分裂。

## NewFrac

官方：https://github.com/newfrac/fenicsx-fracture
固定提交：`f5e112dca19f664bf7ff345209ff4678d09df58d`，MIT。
`notebooks/phase-field/phase-field.py:227–233` 定义 w(alpha)=alpha 和二次退化（加微小残余刚度）；273–296 定义 cw 与能量；386–404 使用 PETSc SNESVI 和显式上下界；432–469 展示位移/损伤交替迭代；522 行用已接受相场作为不可逆下界。

用于交叉核对 AT1 约定、变分不等式、不可逆界与交替最小化思路。没有运行其 FEM notebook，也没有复制 solver；不采用残余刚度。示例迭代上限处容许 pass，本轮按照用户要求改为明确失败退出，不继承这一教学示例的宽松行为。

## Miehe–Hofacker–Welschinger (2010)

期刊原始来源：https://doi.org/10.1016/j.cma.2010.04.011
已读取 ScienceDirect 期刊页摘要及公开章节概要，其内容支持历史场与算子分裂的背景说明。Scuola Normale Superiore 的检索命中 PDF 链接实际返回 404，全文未取得；获取失败完整记录在 MIEHE_FETCH_ERROR.json 和 EXTRA_FETCH_LOG.json。

不把该论文的历史场方法冒充本轮严格 AT1 下界约束，不凭摘要补全文公式。当前核心能量已由用户、本地 PFPD 原图、直接代数和上述公开 AT1 函数明确给定；不需要借缺失全文发明新的核心公式，也没有实现一般张压分裂。

## 数值 oracle 与机械复用

原生 SciPy minimize_scalar(method='bounded') 在 [s_old,1] 上求解。额外显式比较上下端点，因为 bounded 一般不会精确返回边界；不将解析 s 输入 oracle。预声明相场绝对容差 5e-8，64 组 Y/Yc×s_old 全覆盖，保存每个原始优化解和端点选择。

机械部分只读复用已有 `regularization/equilibrium_adr.py:equilibrate` 的 force_callback，不复制 ADR。边界适配器仅给该框架提供移动夹持点坐标，真正 BAQP family 的参考坐标保持原样。内部力不含流体载荷，归一化残差使用事先计算的解析参考反力作为 force floor，防止用人为外力充当收敛分母。

预声明 3×3×3 SVK 拉伸体、nu=0，目的为避免泊松压缩并得到有独立仿射平衡解的小验证算例。自由中间层在每次夹持位移变化后从旧位置开始，仍通过 ADR 求平衡。载荷、Gc、迭代及门限在任何本轮候选计算前固定于 PREDECLARED_CONFIG.json，不根据结果调节。

## 作者分支关键词命中的人工复核

`baqp_equiv_stress_dmgmodel` (`f02c376956b03ce1083e163a026f6dc4528cab60`, 2025-12-18) 的 EquivalentStress 草稿确实引用 phase-field 与 Miehe。实际接口采用临界等效应力和应力软化，g(d) 也进入重构；没有本任务的 AT1 裂纹能、Yc 约束最小化或早期 h=1 分离。它早于 PFPD v1，不能据关键词称为新发布的 PFPD 演化实现，也不解决 Eq3.25 歧义。此次没有采用该应力驱动草稿。
