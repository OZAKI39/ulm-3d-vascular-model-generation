# HIGHER_ORDER_BAQP_ENERGY_FORCE_CANDIDATE_V1

范围：完整态高阶变分力学候选。branch `clot-clearance/pfpd-progress-20261005`；HEAD `695e75166aa864809460ea9a1845fe72abce229c`。

[原始任务](TASK_REQUEST.txt) · [起始 Git 状态](STARTING_GIT_STATE.json) · [预声明配置](CANDIDATE_CONFIGURATION.json)

## 1. Why PD2 kinematics are already accepted

**HIGHER_ORDER_VARIATIONAL_MECHANICS_CANDIDATE = VERIFIED。** 本轮已将接受的完整 PD2 运动学嵌入完整态、保守的 BAQP 能量—力候选，并通过同一离散能量的一阶变分检验。结论限于本报告列出的完整态基准，不表示已可替换生产模型。

前轮 [RK oracle 报告](../higher_order_baqp_rk_oracle_v1/REVIEW_HIGHER_ORDER_BAQP_RK_ORACLE_V1.md) 已验证九项基、全部物理节点多项式再现、轴向/剪切/XZ 二次场、诊断键梯度和独立 Julia 对照。本轮复用该 N19 缓存，三个场的节点 F 与目标位置再次逐位匹配。

主域为 $19^3=6859$ 节点、690020 条无向键、1380040 个有向邻接；$L=1.125\,mm$、$h=62.5\,\mu m$、$\delta=251.25\,\mu m$、$\delta/h=4.02$。二次径向核、梯形体积、完整九项基、无量纲局部坐标、直接矩阵求解均不变。无伪逆、无正则化、无 GMLS。

小型实现检查另用 N7、343 节点、21356 条键、$L=0.375\,mm$；保持相同 h、δ、核和体积规则，调用旧 `block_family` 及已接受的 RK 构造函数一次后缓存。它用于 C1 回归、解析变化量与选定自由度检查，不把 N19 权重错误地裁剪成小网格权重。


N19 权重缓存 SHA256：`ffc7064e7b35b8ffabe52063ae6d14340ba76f4308c9881f23a815ba69634420`；数组/配置身份、小缓存秩与条件数见 [PD2_CACHE_IDENTITY.json](PD2_CACHE_IDENTITY.json)。计算前门槛见 [CANDIDATE_CONFIGURATION.json](CANDIDATE_CONFIGURATION.json)。

## 2. Why old C1 force cannot simply be reused

梯度精度提升本身不定义正确的力算子。[Behzadinasab, Trask, Bazilevs, Part I](https://arxiv.org/abs/2004.11477) 区分了高阶梯度重构与完整对应模型的稳定性问题；[Hillman, Pasetto, Zhou](https://www.mcharleshillman.com/uploads/7/1/8/9/71894909/2019_generalized_reproducing_kernel_peridynamics.pdf) 的梯度构造已在前轮引用和缓存。本轮没有直接移植文献的应力散度力或宣称稳定性已验证。

当前 C1 的 Z/力公式包含原节点梯度的导数结构。只替换 F 而保留不匹配的力映射，不能保证 $f=-\partial U/\partial x$。因此从冻结梯度算子、原键投影与原能量图重新解析求导，并让 **C1 与 PD2 共用一个伴随后端**。没有实现“PD2 F + 旧 C1 Z”的控制分支。

## 3. PD2 energy definition

从旧 [baqp.py](../../pd_clot/pfpd_candidate/baqp.py) 的 `factor_pair` 和 `bond_energy` 读取精确计数约定：

$$m_b=2\bar w_{ij}V_iV_j,\qquad U(x)=\sum_{b=(i,j),\ i<j}m_b\Psi_0(\widetilde F_b(x)).$$

每条物理键只出现一次；系数 2 来自原双端能量约定。$\bar w$ 的单位是 m⁻³，$m_b$ 是 m³，$\Psi_0$ 是 Pa，U 是 J。N19 的 m 数组与前轮保存值逐位一致。

完整态 SVK 直接调用旧 [energy.py](../../pd_clot/pfpd_candidate/energy.py)：

$$E=\tfrac12(F^TF-I),\quad \Psi_0(F)=\tfrac\lambda2(\operatorname{tr}E)^2+\mu E:E,$$
$$P=\frac{\partial\Psi_0}{\partial F}=F[\lambda\operatorname{tr}(E)I+2\mu E],\qquad \lambda=0,\ \mu=1000\,Pa.$$

没有相场变量、退化因子、AT1、Gc、损伤或谱退化力。检查过原谱分裂实现，但本轮力只使用完整态 $\Psi_0$。四个固定张量上对每个 F 分量做微小中心差分验证原 P；步长为 10⁻⁶。


本构最大归一化差 5.780383267e-11，门槛 1.000000000e-07，PASS。见 [CONSTITUTIVE_STRESS_CHECK.json](CONSTITUTIVE_STRESS_CHECK.json)。

## 4. Exact variation of nodal F

保持前轮张量约定：第一下标是位移输出分量，第二下标是参考坐标导数。物理梯度权重 g 已包含邻点体积，单位 m⁻¹。

$$F_i=I+\sum_{k\in H_i}(u_k-u_i)\otimes g_{ik},\qquad u=x-X.$$

因为 g 仅依赖参考配置，变形期间 $\delta g=0$，故

$$\delta F_i=\sum_{k\in H_i}(\delta x_k-\delta x_i)\otimes g_{ik}.$$

这是解析线性变化量，没有在力函数内有限差分 F。`GradientOperator` 的 source、target、g 与 `Reference` 数组均设为只读。节点 F 继续调用前轮 `rk_kinematics.deformation`；每次能量/力求值不求解任何矩矩阵。

## 5. Exact variation of Ftilde

对固定参考键 $\xi=X_j-X_i$，$r^2=\xi\cdot\xi$，$\eta=x_j-x_i$，保留

$$F_{avg}=\tfrac12(F_i+F_j),\qquad
\widetilde F=F_{avg}+\frac{(\eta-F_{avg}\xi)\otimes\xi}{r^2}.$$

令 $R=I-\xi\otimes\xi/r^2$，则 $R^T=R$ 且

$$\widetilde F=F_{avg}R+\eta\otimes\xi/r^2,$$
$$\boxed{\delta\widetilde F=\tfrac12(\delta F_i+\delta F_j)R+(\delta x_j-\delta x_i)\otimes\xi/r^2.}$$

在 N7 二次剪切状态用固定平滑单位方向直接扰动全部位置，对整批键做中心差分，与解析式独立比较。因为投影对 x 为仿射，理论上没有高阶截断项；所见残差来自差分相消。


步长 6.250000000e-10 m，最大绝对导数残差 9.962151997e-07 m⁻¹，导数尺度 7.550589361e+02 m⁻¹，归一化残差 **1.319387338e-09** < 1.000000000e-08。见 [DELTA_FTILDE_IDENTITY.json](DELTA_FTILDE_IDENTITY.json)。

## 6. Matrix-free adjoint derivation

逐键变分为 $\delta U_b=m_bP_b:\delta\widetilde F_b$。利用 Frobenius 内积恒等式

$$P:(\delta F R)=(PR^T):\delta F=(PR):\delta F,$$
$$P:((\delta x_j-\delta x_i)\otimes\xi)=(P\xi)\cdot(\delta x_j-\delta x_i),$$

得到

$$A_b=\tfrac12m_bP_bR_b,\qquad t_b=m_bP_b\xi_b/r_b^2,$$
$$\boxed{\delta U_b=A_b:\delta F_i+A_b:\delta F_j+t_b\cdot(\delta x_j-\delta x_i).}$$

其中 A 的单位是 J，t 的单位是 N。累积 $G_i=\sum_{b\ni i}A_b$ 后，

$$G_i:[(\delta x_k-\delta x_i)\otimes g_{ik}]=(G_ig_{ik})\cdot(\delta x_k-\delta x_i).$$

因此节点梯度项等价于有向边上成对的能量梯度增量 $q_{ik}=G_ig_{ik}$。直接项与节点项相加就是 $\nabla_xU$，内力取负号。上述推导完全来自本轮使用的同一个离散 U，无需显式构造逐键四阶 Jacobian。

对 C1，原公式中的 $Z_i$ 满足 $G_i=V_iZ_i$；代入 C1 的含体积梯度权重可还原原双端力表达式。这只是回归的代数解释；PD2 实现本身不调用旧 Z。

## 7. Efficient two-stage force accumulation

第一阶段逐无向键计算同一 $\widetilde F,\Psi_0,P$：

$$G_i\mathrel{+}=A_b,\quad G_j\mathrel{+}=A_b,\quad
(\nabla U)_j\mathrel{+}=t_b,\quad(\nabla U)_i\mathrel{-}=t_b.$$

第二阶段逐有向邻接 $i\to k$：

$$q_{ik}=G_ig_{ik},\quad(\nabla U)_k\mathrel{+}=q_{ik},\quad
(\nabla U)_i\mathrel{-}=q_{ik},\qquad f=-\nabla U.$$

两个阶段均采用分块、直接索引累加；没有在一条键内循环端点的所有邻居，没有每块长度 N 的散射零缓冲。时间为 $O(B+D+N)$，参考缓存为 $O(B+D+N)$，当前状态中间张量按块保存。每次力调用只有一次节点重构、一次键应力/能量计算、两个伴随累加阶段。没有全局稠密矩阵、Hessian、Schur 补或稳定性特征值计算。

代码职责：

| 文件 | 职责 |
| --- | --- |
| `mechanics_candidate.py` | 冻结参考数据、C1/PD2 算子和统一 `VariationalBAQP` 接口 |
| `variational_energy.py` | 精确原键权重、原键投影、原 SVK 本构 |
| `adjoint_accumulation.py` | 解析变化量和两阶段索引累加 |
| `variational_force.py` | 同一能量图的负梯度 |
| `validation_states.py` | 固定制造态、有限转动和平滑扰动 |

全部位于新候选模块；旧 `baqp.py` 未改。历史 `block_family`、`scripts.verify_pfpd_baqp.deformation`、制造场和连续体 Gauss 工具直接复用。历史 FD 程序绑定旧求值器及不同方向单位，因此保留其中心差分方式，并按本轮明示的 h 尺度和稳健分母实现适配，没有运行历史全量测试。

## 8. C1 energy regression

在 N7 上，先用历史 BAQP 求值器独立比较有限仿射、历史二次剪切、固定 seed=581 的平滑小扰动。随机场最大节点位移为 0.002h；三种状态没有相场。门槛为真正相对能量差 ≤10⁻¹²。所有 C1 门槛通过后才调用 PD2 力。

| 状态 | 历史 U [J] | 通用 C1 U [J] | 相对差 |
| --- | --- | --- | --- |
| affine | 3.051517055e-11 | 3.051517055e-11 | 1.249472170e-14 |
| quadratic_shear | 2.378140639e-13 | 2.378140639e-13 | 6.368895571e-16 |
| seeded_smooth | 1.279654981e-14 | 1.279654981e-14 | 4.315246232e-15 |

[C1_ENERGY_REGRESSION.json](C1_ENERGY_REGRESSION.json)。

## 9. C1 force regression

同一 N7 三态比较全部节点与全部分量。预先定义

$$F_{ref}=\max(4\,Pa\,h^2,\max_i\|f_{historical,i}\|),$$
$$e_{max}=\max_i\|\Delta f_i\|/F_{ref},\quad
e_{RMS}=\sqrt{\sum_iV_i\|\Delta f_i\|^2/\sum_iV_i}/F_{ref}.$$

门槛分别为 10⁻¹⁰、10⁻¹¹。


| 状态 | max absolute [N] | max normalized | weighted RMS |
| --- | --- | --- | --- |
| affine | 5.704138904e-21 | 3.811438357e-14 | 9.974372781e-15 |
| quadratic_shear | 1.896864368e-21 | 1.213993196e-13 | 4.044007147e-14 |
| seeded_smooth | 3.409691753e-21 | 2.182202722e-13 | 4.374548189e-14 |

结果不是逐位相等。节点求和顺序、分块能量求和、K 的线性求解与原显式逆、以及 r² 的点积与范数平方之间的舍入顺序不同，产生约 10⁻²¹ N 的差异；没有放宽门槛。[C1_FORCE_REGRESSION.npz](C1_FORCE_REGRESSION.npz) 保存两套完整力数组。

[600 dpi PNG](../../results/higher_order_baqp_energy_force_candidate_v1/figures/c1_variational_force_regression.png) · [Vector PDF](../../results/higher_order_baqp_energy_force_candidate_v1/figures/c1_variational_force_regression.pdf)

## 10. Rigid translation / rotation objectivity

非零平移为 $c=(0.0007,-0.0005,0.001)$ m；有限转动 Q 是绕归一化轴 (1,2,3) 的 37° 正交转动，检查 $x=QX+c$。节点/键 F 分别应为 I 或 Q，Green 应变、能量与内力应为数值零。

节点/键 F 与应变门槛均为 5×10⁻¹²，能量绝对门槛为 10⁻²⁵ J。刚体力门槛在计算前声明为 $\max(10^{-12}(4\,Pa\,h^2),10^{-15}\,N)$；后一项直接继承历史 `test_translation_zero_energy_force` 的门槛。实际残差远小于该历史绝对门槛。


| 状态 | U [J] | max |f| [N] | max node F error | max bond F error | max Green strain |
| --- | --- | --- | --- | --- | --- |
| translation | 4.262497135e-37 | 6.163512675e-21 | 9.887173992e-16 | 1.925351325e-15 | 1.852352740e-15 |
| rotation | 1.917299899e-36 | 1.753551053e-20 | 9.179839318e-15 | 8.946524484e-15 | 6.335387498e-15 |

还对非零有限剪切态叠加同一 Q 与 c，检验 U 不变和 $f(Qx+c)=Qf(x)$：能量相对差 6.774764149e-15，最大归一化力协变误差 1.317360143e-13，均通过 10⁻¹⁰。

理论上，冻结 g 的线性完备性使 $F(Qx+c)=QF(x)$，当前键投影也随 Q 左乘；SVK 仅依赖 FᵀF，故离散能量客观。见 [OBJECTIVITY_RESULTS.json](OBJECTIVITY_RESULTS.json)。

## 11. Global force balance

两阶段中每个能量梯度增量以等大反号分配给两个节点，所以平移不变性导出的全局内力和在精确算术下为零。

对每个非刚体验证态计算 $\|\sum_i f_i\|/F_{ref}$；此处 $F_{ref}=\max(4\,Pa\,h^2,\max_i\|f_i\|)$。N7 C1 三态、N19 PD2 六态、叠加转动态以及 N7 选定 DOF 态都列入结果。


| 模型/域 | 状态 | global force residual |
| --- | --- | --- |
| C1 N7 | affine | 4.117001977e-14 |
| C1 N7 | quadratic_shear | 1.592549051e-14 |
| C1 N7 | seeded_smooth | 4.081538296e-15 |
| PD2 N19 | finite_affine_shear | 6.438459058e-13 |
| PD2 N19 | axial | 1.765678029e-14 |
| PD2 N19 | shear | 6.939636356e-13 |
| PD2 N19 | cross | 1.523270249e-13 |
| PD2 N19 | affine_axial | 4.586310021e-15 |
| PD2 N19 | affine_shear | 5.949789728e-13 |
| PD2 N19 | superposed_finite_shear | 2.323310313e-13 |
| PD2 N7 | cross_selected_DOF | 1.059406428e-14 |

全部小于 10⁻¹⁰；按节点向量和检验，未用所有节点力范数之和扩大分母。完整向量与尺度见 [GLOBAL_BALANCE_RESULTS.csv](GLOBAL_BALANCE_RESULTS.csv)。

## 12. Global moment balance

主验收明确使用当前坐标及当前算术质心：

$$M_{sum}=\sum_i(x_i-x_c)\times f_i,\quad x_c=\frac1N\sum_i x_i,\quad M_{ref}=F_{ref}L.$$

客观性能量对无穷小叠加转动的变化为零，代入 $f=-\nabla U$ 即给出当前配置角动量平衡。节点间单项贡献无需都是中心力。参考坐标力矩不是此处的主判据。


| 模型/域 | 状态 | current moment residual |
| --- | --- | --- |
| C1 N7 | affine | 3.690402164e-15 |
| C1 N7 | quadratic_shear | 1.314158047e-14 |
| C1 N7 | seeded_smooth | 1.237695758e-15 |
| PD2 N19 | finite_affine_shear | 1.547397807e-13 |
| PD2 N19 | axial | 5.848692511e-15 |
| PD2 N19 | shear | 7.672242667e-14 |
| PD2 N19 | cross | 9.397078701e-14 |
| PD2 N19 | affine_axial | 2.094125522e-15 |
| PD2 N19 | affine_shear | 7.528075880e-13 |
| PD2 N19 | superposed_finite_shear | 3.021314573e-13 |
| PD2 N7 | cross_selected_DOF | 4.559235617e-15 |

全部小于 10⁻¹⁰。

[600 dpi PNG](../../results/higher_order_baqp_energy_force_candidate_v1/figures/pd2_objectivity_and_global_balance.png) · [Vector PDF](../../results/higher_order_baqp_energy_force_candidate_v1/figures/pd2_objectivity_and_global_balance.pdf)

## 13. Energy-force directional finite differences

严格使用三个力学态：有限仿射剪切 γ=0.03、旧二次剪切、旧 XZ 场；每态两个平滑方向，seed=581/582。方向由低频三角/线性函数的固定随机系数组合得到，移除均值平移并归一化 $\|v\|_F=1$。

$$D_{FD}(\epsilon)=\frac{U(x+\epsilon v)-U(x-\epsilon v)}{2\epsilon},\qquad D_f=-\sum_i f_i\cdot v_i,$$
$$\epsilon=h\,[10^{-5},3\times10^{-6},10^{-6}],$$
$$E_{ref}=\max(|U(x)|,\mu\,0.004^2\sum_iV_i),\quad
D_{ref}=\max(|D_{FD}|,|D_f|,E_{ref}\|v\|_F/h).$$

主门槛是固定中间步长 $|D_{FD}-D_f|/D_{ref}\le5\times10^{-6}$。全部三步还须小于预先声明的 5×10⁻⁵ 粗异常门槛；全部低于 10⁻⁷ 的序列按预声明判为数值噪声区。没有挑选最优步长替代中间步长。


| 状态 | 方向 | error @1e-5 | error @3e-6 | error @1e-6 | 中间步不加稳健下限的相对误差 |
| --- | --- | --- | --- | --- | --- |
| finite_affine_shear | 1 | 2.235825988e-11 | 3.311185721e-11 | 8.687982320e-11 | 4.831288180e-10 |
| finite_affine_shear | 2 | 1.974694317e-11 | 3.939782664e-11 | 1.251382816e-11 | 4.660608638e-10 |
| shear | 1 | 3.734471971e-11 | 2.305946258e-10 | 3.901588144e-10 | 1.182101101e-09 |
| shear | 2 | 1.752219528e-11 | 4.334856680e-11 | 1.398546963e-10 | 2.097002052e-10 |
| cross | 1 | 2.105443124e-11 | 7.217406147e-11 | 2.967458457e-10 | 2.501334439e-09 |
| cross | 2 | 3.837602898e-11 | 1.296396035e-11 | 1.318257431e-10 | 4.358139649e-10 |

全部六组 PASS。更小步长的残差不严格单调；它们已进入舍入/相消区，不能将这些点宣称为已观测到的二阶收敛斜率。中间步没有显示影响验收的相消污染，且不加稳健分母下限的误差也一并公开。原始 U±、Dfd、−f·v、尺度与全部 18 行数据见 [ENERGY_FORCE_DIRECTIONAL_DERIVATIVES.csv](ENERGY_FORCE_DIRECTIONAL_DERIVATIVES.csv)。没有安装或运行 AD/复步框架。

[600 dpi PNG](../../results/higher_order_baqp_energy_force_candidate_v1/figures/pd2_energy_force_directional_derivative.png) · [Vector PDF](../../results/higher_order_baqp_energy_force_candidate_v1/figures/pd2_energy_force_directional_derivative.pdf)

## 14. Selected-DOF derivative oracle

在 N7 的 XZ 目标态，仅取几何中心、面、棱、角四点的 xyz 分量，共 12 个 DOF；中心点不等于视界完全未截断。直接逐坐标对总能量做中心差分，步长固定为 3×10⁻⁶h，与快速伴随内力比较。对于近零真力，仍使用共同 F_ref，未逐分量除以零。

| 节点（零基） | 分量 | slow force [N] | adjoint force [N] | normalized error |
| --- | --- | --- | --- | --- |
| 171 | 0 | 6.731613058e-20 | 1.444408508e-24 | 4.308139915e-12 |
| 171 | 1 | -0.000000000e+00 | -2.359405139e-26 | 1.510019289e-18 |
| 171 | 2 | -1.446576239e-10 | -1.446576244e-10 | 3.080611845e-11 |
| 168 | 0 | -3.702387182e-19 | 1.939068068e-23 | 2.369651897e-11 |
| 168 | 1 | 3.365806529e-20 | 7.094387734e-26 | 2.154111638e-12 |
| 168 | 2 | 1.071156531e-09 | 1.071156522e-09 | 5.740525587e-10 |
| 21 | 0 | -7.750273084e-09 | -7.750273087e-09 | 2.181721820e-10 |
| 21 | 1 | -0.000000000e+00 | 3.059441495e-25 | 1.958042556e-17 |
| 21 | 2 | -2.504816339e-09 | -2.504816343e-09 | 2.353715566e-10 |
| 0 | 0 | -2.963589168e-09 | -2.963589169e-09 | 8.258497120e-11 |
| 0 | 1 | 1.312664546e-18 | 6.733387182e-25 | 8.401048787e-11 |
| 0 | 2 | -9.427969830e-10 | -9.427969842e-10 | 7.609067695e-11 |

全部小于 10⁻⁶。[SELECTED_DOF_FORCE_ORACLE.csv](SELECTED_DOF_FORCE_ORACLE.csv) 同时保存 U±，可独立复算符号和尺度。未形成全 Jacobian 或刚度矩阵。

## 15. Affine C1 vs PD2 equivalence

N19 仿射轴向与剪切均采用前轮 0.003 梯度。节点 F 完备、同一键投影、同一本构和权重意味着 U 应一致；同时逐键比较 Ftilde 与能量密度，覆盖全部物理键。

| 场 | U_C1 [J] | U_PD2 [J] | 相对能量差 | max ΔFtilde / 0.004 | max density relative difference |
| --- | --- | --- | --- | --- | --- |
| affine_axial | 1.285292532e-11 | 1.285292532e-11 | 8.170370644e-15 | 1.252866636e-13 | 2.955576294e-13 |
| affine_shear | 6.407255395e-12 | 6.407255395e-12 | 1.134673375e-15 | 1.237059364e-13 | 2.068166571e-13 |

能量门槛 10⁻¹²；逐键 F 与密度门槛分别按预声明尺度取 10⁻¹⁰。全部通过，见 [AFFINE_ENERGY_COMPARISON.json](AFFINE_ENERGY_COMPARISON.json)。

## 16. Quadratic energy decomposition

完整二次场中解析 F 是坐标的一次函数，解析键中点 F 同时等于端点平均且满足键割线关系。所以验证目标是

$$U_{PD2}\simeq U_1=\sum_bm_b\Psi_0(F_{exact}(X_{mid,b})).$$

$U_2=\int_\Omega\Psi_0(F_{exact}(X))\,dV$ 是另一项连续体积分参考，不能通过调整 m 强迫匹配。轴向/剪切 U0/U1/U2 与此前保存的分解再次核对；XZ 的 U1 直接用解析中点 F，U2 用 XZ 张量 Gauss 并乘精确 Y 长度，32 与 64 阶相符。没有新平衡模拟。


| 场 | U_C1 [J] | U_PD2 [J] | U_exactkin/currentweights [J] | U_continuum [J] |
| --- | --- | --- | --- | --- |
| axial | 1.328909535e-11 | 1.332236345e-11 | 1.332236345e-11 | 1.333181265e-11 |
| shear | 6.623654542e-12 | 6.639884033e-12 | 6.639884033e-12 | 6.644566561e-12 |
| cross | 2.684446957e-12 | 2.791470854e-12 | 2.791470854e-12 | 2.847659794e-12 |

| 场 | U_PD2−U_C1 [J] | U_PD2−U_continuum [J] | PD2 vs U1 normalized error |
| --- | --- | --- | --- |
| axial | 3.326809706e-14 | -9.449208931e-15 | 4.623814810e-14 |
| shear | 1.622949134e-14 | -4.682527896e-15 | 1.914769936e-15 |
| cross | 1.070238968e-13 | -5.618893952e-14 | 1.843852531e-15 |

所有 U1 比较通过 10⁻¹⁰ 门槛。残留 U1−U2 属于既定权重/积分差异，本轮原样保留；不能将能量—力自洽误读为积分已精确或物理边界已经正确。

[600 dpi PNG](../../results/higher_order_baqp_energy_force_candidate_v1/figures/c1_pd2_energy_decomposition.png) · [Vector PDF](../../results/higher_order_baqp_energy_force_candidate_v1/figures/c1_pd2_energy_decomposition.pdf)

## 17. Optional manufactured recovery

**NOT_RUN。** 本轮第一变分、客观性、平衡和独立导数门槛已经提供所需证据；可选求解不是强制门槛。没有调用力学平衡求解器，没有制造外载、求解历史或可选恢复图。

## 18. What is now verified

**已经接受的完整 PD2 键关联运动学，可以构成一个完整态、保守的 BAQP 候选；其内力已按本轮测试验证为同一离散能量的负一阶变分。**

| 状态 | 结果 |
| --- | --- |
| PD2_KINEMATICS | REUSED_VERIFIED |
| DELTA_FTILDE_IDENTITY | PASS |
| C1_VARIATIONAL_ENERGY_REGRESSION | PASS |
| C1_VARIATIONAL_FORCE_REGRESSION | PASS |
| PD2_RIGID_TRANSLATION | PASS |
| PD2_RIGID_ROTATION | PASS |
| PD2_GLOBAL_FORCE_BALANCE | PASS |
| PD2_GLOBAL_MOMENT_BALANCE | PASS |
| PD2_ENERGY_FORCE_FD | PASS |
| PD2_SELECTED_DOF_FORCE_ORACLE | PASS |
| PD2_AFFINE_ENERGY_REGRESSION | PASS |
| PD2_QUADRATIC_ENERGY_KINEMATIC_REFERENCE | PASS |
| PD2_OPTIONAL_MANUFACTURED_RECOVERY | NOT_RUN |
| HIGHER_ORDER_VARIATIONAL_MECHANICS_CANDIDATE | VERIFIED |
| PRODUCTION_BAQP | UNCHANGED |
| BOUNDARY_DATA_CONVERSION | NOT_REINTEGRATED_THIS_ROUND |
| PHASE_FORCE | NOT_IMPLEMENTED |
| C5_C7 | NOT_RUN |
| FRACTURE | NOT_RUN |
| CLOT_GEOMETRY | NOT_RUN |

12 项指定测试全部通过，耗时 0.82 s；[测试日志](TARGETED_TESTS.log)。只运行本轮套件。最终计算验收 25.045 s，3 个独立状态/方向线程、BLAS=1，峰值 RAM 高水位 0.635 GiB，最小可用内存占比 61.40%，GPU_NOT_USED_BY_DESIGN。资源值见 [PERFORMANCE_SUMMARY.json](PERFORMANCE_SUMMARY.json)。

数值先通过后，代码审查将分块全节点缓冲改成直接索引累加以满足线性复杂度，随后完整复验同一门槛。初次通过的日志/摘要保留在 [assembly_implementation_check](assembly_implementation_check/)，不作为最终结果。最终时间只对应最终实现的一次验收，初次计算时间单独保留；没有改物理、权重、步长或门槛。

所有 600 个受保护旧文件在计算前后哈希一致，含前轮缓存、生产源码和旧 manifest。最终文档/图件/哈希验收见 [FINAL_ARTIFACT_QA.json](FINAL_ARTIFACT_QA.json)。输出字段见 [DATA_DICTIONARY.md](DATA_DICTIONARY.md)。

## 19. What remains unverified

未验证 PD2 的物理边界条件、边界数据转换重接、非线性平衡求解兼容性、稳定性谱、相场耦合力、C5/C7、断裂、网格客观性、不规则血栓/壁面力学或生产可用性。未运行 Ogden、黏附、侵蚀、CFD 或微泡加载。

本轮通过不能替代这些未来验证。尤其当前内力在精确制造态上非零，反映有限域能量及其离散梯度；尚未按 PD2 重新建立边界等效载荷。旧 C1 生产实现与科学结论保持不变。

## 20. One next development task only

**HIGHER_ORDER_BAQP_BOUNDARY_REINTEGRATION_V1**。人工审核本报告后，才以已验证的 PD2 变分内力重建既有边界数据转换中的等效载荷，并复核高价值法向/剪切/非仿射基准的力、力矩、虚功、F、Ftilde、Psi_plus。

本轮到此停止；没有自动开始下一阶段。
