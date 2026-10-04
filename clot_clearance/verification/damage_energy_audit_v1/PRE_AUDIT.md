# 只读能量机制审计（新增计算代码之前）

已读 mechanics.py、fragment_mechanics.py、fast_mechanics.py、energy.py、erosion_energy.py、
erosion.py、erosion_equilibrium_runner.py；本轮承接 adr_completion 下的完成版机械平衡报告。
上级旧报告仍是缺论文阶段记录。未发现适用的 AGENTS.md。工作分支沿用已有 energy-regularization 基线。

## A. 总能量

potential 返回 sum_i V_i E_i 与 sum_i V_i H_i，分别是本构和非仿射稳定化储能。
本次不含物理/伪动能、外力势能或损伤耗散假设，不把流体做功混入固定坐标比较。

## B. 粒子本构能量

full-rank 时 E_i=r_i W(F_i)，W=μ/2(F:F−3)−μ log J+λ/2(log J)^2，
λ=bulk−2μ/3，J=det F。支持维数下降时使用原有 intrinsic 子空间形式，不当作三维F。
本轮通过小幅损伤保持满秩，但API也报告真实rank，不掩盖退化。

## C. remaining

r_i=sum_j(g_ij V_j)/sum_j(V_j)，这里没有 influence 权重，分母为原始邻居体积和。
本构能量乘r，first Piola应力进入力态前也乘r。稳定化系数是alpha*μ*r_i/trace(K_i)。

## D. integrity 如何进入

ξ_ij=X_j−X_i，η_ij=x_j−x_i，q_ij=w_ij*g_ij。
K_i=sum_j q_ij V_j ξ⊗ξ；A_i=sum_j q_ij V_j η⊗ξ；F_i=A_i K_i^{-1}。
双端力含q*Vi*Vj、两端r*P*Kinv以及稳定化系数乘(η−Fξ)。
H_i=0.5*(alpha*μ*r_i/trace K_i)*sum_j q_ij V_j |η−F_i ξ|²。
现有fast kernel只是等价标量循环路径；审计直接复用它，不写新求解器。

## E. 更新时重新计算的量

损伤更新D后g=1−D，重新 prepare_supported_shape 得到K逆、trace、remaining、rank和basis，
再调用evaluate_supported重新计算F、J、本构能量、稳定化能量及力。旧runner在同一个x下比较前后储能。
本次读取这些保存的x/D，绝不重新平衡或运行侵蚀。

## F. x 不动而 F 改变

F是邻域键变形的加权最小二乘重构，权重随损伤改变。如果变形不完全仿射，改变各键相对权重
会改变最佳拟合F。若严格仿射且满秩，A=F_macro K，因此任何保持满秩的权重变化应保持F。

## 实现前的代数核对

全键g同时乘c，K、A、trace、remaining同时乘c，F不变；稳定化系数r/trace不变，
稳定化求和中的q乘c。因此本构与稳定化能量都应乘c，包括非仿射满秩状态。
是否满足由后续数据检验，不据此预先给出PASS。

精确分解：Vi(r0W0−r1W1)=Vi(r0−r1)W0+Vi*r1(W0−W1)，稳定化另列。
W从原内核输出E/r提取，不额外替换本构表达式。保存逐粒子及全局闭合残差。

## 预声明与范围

provenance/PREDECLARED_PLAN.json 声明重放、分解、零态、F不变性和有限差分平台判据。
重放超容差立即停止后续批量任务，保留失败证据。35个历史跳变每个先采样最大ΔD的10条键诊断Y；
更完整基线抽样和小算例全键敏感性另列，不将抽样负Y比例冒充全模型比例。
有限差分按D边界选中央、前向或后向；用逐粒子能量差相加降低无关大数相减误差，保留全部epsilon。
只有连续三档形成预声明平台才报告稳定导数；否则标记NUMERICALLY_UNRESOLVED。
本轮只作诊断，所有counterfactual均为DIAGNOSTIC ONLY / NOT A CANDIDATE MODEL / NOT USED IN PRODUCTION。
冻结旧K但更新A的情形可以通过原kernel传入旧inverse和trace、新remaining来定义；它失去严格仿射再现，
只用于定位，不称为替代模型。冻结F反事实仅改变本构r，保留旧稳定化能量以明确隔离项。

全部新输出仅在本目录；新增审计helper、脚本及测试，production源码不修改。
先完成本地机制审计，再进行外部文献与开源代码对照。
