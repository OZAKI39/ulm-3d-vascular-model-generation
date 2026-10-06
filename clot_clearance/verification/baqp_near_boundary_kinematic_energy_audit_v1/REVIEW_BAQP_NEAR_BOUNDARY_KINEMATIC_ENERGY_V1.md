# BAQP near-boundary kinematic / energy audit V1

2026-10-06 · 固定态诊断 · 原物理公式未改

## 1. Why boundary-force conversion is no longer the main blocker

已有轴向和剪切非仿射制造基准的边界数据局部化验证保持有效，范围仍是这两个完整态基准。当前问题出现在**精确制造目标节点坐标本身**，无需新的平衡求解即可复现，因此本轮检查的是运动学算子。旧边界转换、载荷、残差和 C5/C7 接受态均未改。

读取来源：[剪切审计](../nonlocal_boundary_nonaffine_shear_localization_v1/REVIEW_NONAFFINE_SHEAR_LOCALIZATION_V1.md)、[轴向审计](../nonlocal_boundary_nonaffine_localization_v1/REVIEW_NONAFFINE_BOUNDARY_LOCALIZATION_V1.md)、[本轮冻结配置](PREDECLARED_CONFIG.json)、[目标态复现](SAVED_TARGET_REPRODUCTION.json)。两个保存的 X、V、目标坐标和解析 F 均逐字节匹配；用原 `block_family` 恢复相同拓扑，仅各调用一次固定态 `baqp.evaluate()`，不求平衡。

N19=6,859 个物理粒子、690,020 个无向键，L=1.125 mm，h=62.5 µm，δ=251.25 µm，δ/h=4.02；quadratic 核。SVK μ=1000 Pa、λ=0、bulk=666.6666666666666 Pa。PHASE=DISABLED；equilibrium、fracture、clot、C5/C7 均未运行。

## 2. Observed Fbar/Ftilde/Psi_plus bias

F 偏差使用 Frobenius 范数并除以固定 0.004；这不是除以局部 F 的相对误差。Ψ⁺ 列为最大绝对能量密度偏差，单位 Pa，并除以各场同一个解析节点最大值。剪切 Ψ_ref=0.004016032032 Pa，轴向 Ψ_ref=0.016064064 Pa；节点和键共用各自场的参考。

| Field | Channel | Actual max (F dimensionless / Ψ Pa) | Actual / scale | Predicted / scale |
| --- | --- | --- | --- | --- |
| shear | Fbar | 1.19033543e-04 | 2.97583857e-02 | 2.97583857e-02 |
| shear | Ftilde | 1.17714149e-04 | 2.94285372e-02 | 2.94285372e-02 |
| shear | Psi_Fbar | 2.58288901e-04 | 6.43144523e-02 | 6.43144523e-02 |
| shear | Psi_Ftilde | 2.57886483e-04 | 6.42142495e-02 | 6.42142495e-02 |
| axial | Fbar | 1.19033543e-04 | 2.97583857e-02 | 2.97583857e-02 |
| axial | Ftilde | 1.17714149e-04 | 2.94285372e-02 | 2.94285372e-02 |
| axial | Psi_Fbar | 9.35209793e-04 | 5.82175092e-02 | 5.82175092e-02 |
| axial | Psi_Ftilde | 9.27745038e-04 | 5.77528226e-02 | 5.77528226e-02 |

剪切近边界 Fbar 约 2.976%、Ftilde 约 2.943%、Ψ⁺ 约 6.43%，描述的是目标态的重构偏差，不能与上一轮平衡求解的目标恢复精度混合。[全量因果摘要](CAUSAL_IDENTITY_SUMMARY.json)。

## 3. Quadratic BAQP reconstruction identity

对任意二次映射 y，Taylor 展开没有余项：

$$\eta_a=F^{exact}_{ac}(X_i)\xi_c+\tfrac12H_{acd}\xi_c\xi_d.$$

原 BAQP 定义 $A_{ab}=\sum_j\omega\eta_a\xi_bV_j$，故

$$A_{ab}=F^{exact}_{ac}K_{cb}+\tfrac12H_{acd}M_{3,cdb},\qquad
\boxed{\Delta\bar F_{an}=\tfrac12H_{acd}M_{3,cdb}(K^{-1})_{bn}}.$$

这是由当前 [baqp.py](../../pd_clot/pfpd_candidate/baqp.py) 的 $\bar F=AK^{-1}$ 推导的有限和恒等式，条件是 K 可逆；不是从 BB/OSB 论文移植的公式。实际实现的位移形式 $I+A_uK^{-1}$ 在精确算术下等价。

剪切 $u_x=\gamma Z+\beta_s Z^2/L$，$\gamma=0.003,\beta_s=0.001$；轴向 $u_z=\epsilon Z+\beta Z^2/L$，$\epsilon=0.003,\beta=0.001$。令 $a=x$ 或 $z$：

$$H_{azz}=2\beta/L,\qquad \Delta\bar F_{an}=(\beta/L)M_{3,zzb}(K^{-1})_{bn}.$$

其余输出行只有浮点误差。通用 Hessian 收缩与此专门形式分别计算并核对。

## 4. Second- and third-moment definitions

$$K_{cb}=\sum_{j\in H_i}\omega_{ij}\xi_c\xi_bV_j,\qquad
M_{3,cdb}=\sum_{j\in H_i}\omega_{ij}\xi_c\xi_d\xi_bV_j.$$

ω 无量纲、ξ 单位 m、V 单位 m³，所以 K、M3、G=M3_zz: K⁻¹ 单位分别为 m⁵、m⁶、m；H 单位 m⁻¹。独立累加器按标量分量与端点 scatter 求和，反向键的 K 项同号、M3 项异号；不借用 evaluate 的 F 输出生成预测。

[MOMENT_TENSORS.npz](MOMENT_TENSORS.npz) 保存全部物理节点的 K、invK、M3、均匀体积对照和 G，并包含最近面、棱、角距离，family size、ω0、condition(K) 和几何类。分类只数“距某轴最近面小于 δ”的轴数，0/1/2/3 分别为 bulk-like/face-near/edge-near/corner-near；它不宣称多跳力模板完全闭合。

| Geometry | Nodes | Family min | Family max | max cond(K) | max \|\|M3\|\| [m6] | max predicted / 0.004 |
| --- | --- | --- | --- | --- | --- | --- |
| bulk-like | 729 | 256 | 256 | 1.00000000e+00 | 8.53828316e-39 | 1.79371510e-17 |
| face-near | 2430 | 152 | 256 | 1.36269417e+00 | 3.61065714e-24 | 2.63086546e-02 |
| edge-near | 2700 | 90 | 256 | 3.61669692e+00 | 3.61065714e-24 | 2.87318798e-02 |
| corner-near | 1000 | 53 | 256 | 4.92504538e+00 | 3.61065714e-24 | 2.97583857e-02 |

## 5. Continuous full-ball reference

对当前 $\omega(r)=(1-r/\delta)^2$，$r<\delta$：

$$K_{full}=\frac{4\pi}{3}I\int_0^\delta r^4(1-r/\delta)^2\,dr
=\frac{4\pi\delta^5}{315}I,\qquad \omega_{0,full}=\frac{2\pi\delta^3}{15}.$$

球域和径向核在 $\xi\mapsto-\xi$ 下不变，而每个三阶单项式变号，因此所有 M3 分量积分为零，包括 zzx、zzy、zzz。二次 Fbar 重构由第 3 节恒等式精确成立。

K 标量系数 **3.99419867e-20 m⁵**；解析 ||M3||=0。独立 8×12×24 极坐标 Gauss/方位积分的 K 相对误差 2.43922254e-15，M3 数值范数 1.38690803e-39 m⁶，按 ||K||δ 缩放为 7.97905212e-17。

[连续球参考](CONTINUOUS_FULL_BALL_REFERENCE.json)。

## 6. Discrete full-lattice reference

原有 inversion-symmetric offset 集含 256 个邻居，所有体积为 h³，使用同一核和 δ。K 相对连续球参考误差为 **3.72424068e-03**；||M3||=2.54392944e-40 m⁶（按 ||K||δ 缩放 1.45812343e-17）。

在全部 6,859 个物理观测位置外延同一多项式，以完整 offset 家族直接重构，shear/axial 最大 Fbar 误差/0.004 分别为 **6.18228203e-14 / 1.11022303e-13**，低于 1e-10。

K 的 0.3724% 积分误差不妨碍此二次场重构，因为同一个离散 K 在 A 与 K⁻¹ 中配对，奇矩才决定二次误差。H3/H4 作为本轮偏差来源不获支持；这不证明 BAQP 对所有非仿射场都精确。

[完整晶格参考](FULL_LATTICE_REFERENCE.json)、[全观测点直接重构](FULL_LATTICE_RECONSTRUCTION.npz)。

## 7. Truncated-support reference

R1 为完整对称晶格；R2 保留实际物理邻居坐标，诊断体积统一设 h³；R3 为实际梯形体积。R2/R3 不改变已有生产状态，只建立独立诊断有限和。

| Variant | max error / 0.004 | V-weighted RMS / 0.004 |
| --- | --- | --- |
| R1_full | 5.11811737e-19 | 5.11811737e-19 |
| R2_support_uniform | 2.85268581e-02 | 9.55378380e-03 |
| R3_physical | 2.97583857e-02 | 1.11038010e-02 |

R2 单独已产生约 2.85% 最大偏差，R3 约 2.98%；支持截断足以产生主体偏差。连续与晶格 full references 均为零或数值零，符合用户 CASE A 的判据。[逐节点三变体表](FAMILY_VARIANT_COMPARISON.csv)。

## 8. Effect of physical boundary volumes

$$A=E_{support}-E_{full},\quad B=E_{physical}-E_{support},\quad C=E_{physical}-E_{full},\quad A+B=C.$$

这些是逐节点向量贡献。各自的最大范数和 RMS **不具有标量加法关系**，不能用两列最大值相减定义体积贡献。

| Contribution | max / 0.004 | V-weighted RMS / 0.004 |
| --- | --- | --- |
| support_contribution | 2.85268581e-02 | 9.55378380e-03 |
| volume_contribution | 5.82906974e-03 | 2.74908711e-03 |
| total_contribution | 2.97583857e-02 | 1.11038010e-02 |

向量闭合残差/0.004=8.47136339e-19。体积贡献最大约 0.583%、RMS 约 0.275%，可测且次于支持贡献，不能忽略。距边界 4h 的 602 个节点均有 256 邻居，R2 最大 1.63910138e-17，R3 仍为 2.61687769e-07（均 /0.004）：边界半体积可破坏完整坐标集合的加权对称性。

已检查相对两端面、侧面、棱、角的有符号 M3_zz: 和 G，见 [GEOMETRY_CAUSAL_DETAILS.json](GEOMETRY_CAUSAL_DETAILS.json)。Z 轴中心线上上下端偏差反号；侧面可产生输入 x/y 分量，不能把所有误差都看成 zz 分量或赋予统一符号。[向量分解](SOURCE_CONTRIBUTIONS.npz)、[摘要](SOURCE_DECOMPOSITION.json)。

## 9. Exact Fbar error reconstruction

| Field | max residual / 0.004 | V-weighted RMS / 0.004 |
| --- | --- | --- |
| shear | 7.23045700e-13 | 1.42773173e-13 |
| axial | 1.49879654e-12 | 4.45801840e-13 |

全部节点通过预声明 max≤1e-10、weighted RMS≤1e-11；未改门槛。实际残差是完整张量之差的范数，不是两个误差范数之差。通用/专门预测之间的误差也接近舍入尺度。独立 K 与原 evaluate 的 K 对照及保存目标态 bitwise 检查通过。

[逐节点闭合](FBAR_MOMENT_CLOSURE.csv)、[深度统计](FBAR_DEPTH_PROFILES.csv)、[因果摘要](CAUSAL_IDENTITY_SUMMARY.json)。

## 10. Role of K anisotropy

保留物理 M3，分别使用物理 K、连续完整球 K 及完整晶格 K：

$$R_K=\frac{\|\tfrac12HM_3K_{physical}^{-1}\|}{\|\tfrac12HM_3K_{full}^{-1}\|}.$$

连续球参考的 R_K 范围 **[0.996294270729267, 7.517341817553924]**；晶格参考 **[1.0000047103798373, 7.545338207743268]**；有效节点 5770。安全条件为分母/0.004>1e-10，否则存 NaN，不能将 0/0 当放大率。condition(K) 范围 [1.0, 4.925045384359479]。

这个按用户定义命名的 “K anisotropy modifier” 同时包含 K 幅值损失与各向异性，不是纯各向异性量。它会放大或转向已有的奇矩误差；对任意可逆 K，若 M3=0，二次误差仍为零，因此 K 本身不是本轮根因。

## 11. Exact Ftilde projection identity

二次 y 的 secant 等于中点解析梯度乘 ξ，且解析 F 沿位置线性：

$$\eta=F^{exact}_{mid}\xi,\qquad F^{exact}_{mid}=\tfrac12(F^{exact}_i+F^{exact}_j).$$

代入原 $\widetilde F=F_{av}+(\eta-F_{av}\xi)\otimes\xi/|\xi|^2$，令 $D=(\Delta\bar F_i+\Delta\bar F_j)/2$、$n=\xi/|\xi|$，得到

$$\boxed{\Delta\widetilde F=D(I-n\otimes n)}.$$

先用实际节点 ΔF 独立核验投影恒等式；再用**独立矩预测**替换两个节点 ΔF，验证从 M3 一直到键的完整因果链。

| Field | Prediction | max residual / 0.004 | Weighted RMS / 0.004 |
| --- | --- | --- | --- |
| shear | Ftilde | 3.38111190e-13 | 4.86566554e-14 |
| shear | Ftilde_moment_chain | 6.15768096e-13 | 9.85955805e-14 |
| axial | Ftilde | 2.77555756e-13 | 4.71230804e-14 |
| axial | Ftilde_moment_chain | 1.49879636e-12 | 3.07574938e-13 |

全部 690,020 个键满足相同 1e-10/1e-11 门槛；RMS 权重为原 BAQP 能量权重 $2\bar\omega_{ij}V_iV_j$。结论是节点重构偏差经横向投影保留，不是本轮证据所支持的额外键修正缺陷。

## 12. Bond-orientation dependence

逐键保存 |nx|、|ny|、|nz|、r/δ、中点最近面深度/δ、实际和预测 ΔFtilde 范数、D 范数及削减比。投影的 Frobenius 范数满足 $\|D(I-nn^T)\|_F\le\|D\|_F$；仅在 D 安全非零时报告比值。

| \|nz\| bin | Bonds | Actual max / 0.004 | Weighted RMS / 0.004 | Weighted mean reduction ratio |
| --- | --- | --- | --- | --- |
| [0,0.25) | 136762 | 2.94285372e-02 | 1.10207689e-02 | 7.90610011e-01 |
| [0.25,0.5) | 190656 | 2.49527451e-02 | 9.22642610e-03 | 8.17062223e-01 |
| [0.5,0.75) | 191024 | 2.09908205e-02 | 7.80387925e-03 | 7.96017572e-01 |
| [0.75,0.95) | 147752 | 1.34938291e-02 | 5.62744763e-03 | 7.24816132e-01 |
| [0.95,1.0] | 23826 | 1.23601515e-02 | 5.03644234e-03 | 5.61061083e-01 |

接近横向的键组保留较大偏差；但削减比不是单纯的 $\sqrt{1-n_z^2}$，因为实际 G 在棱角含 x/y 分量且 D 随空间变化。比值均值也不必按 |nz| 单调。接近 z 的键仍可能保留横向输入分量。[键方向表](BOND_ORIENTATION_PROFILES.csv)、[逐键几何](BOND_GEOMETRY.npz)。

## 13. Axial/shear geometry cross-check

$$G_{i,n}=M_{3,i,zzb}(K_i^{-1})_{bn},\qquad
\Delta\bar F^{shear}_{xn}=(\beta_s/L_s)G_{i,n},\quad
\Delta\bar F^{axial}_{zn}=(\beta/L_a)G_{i,n}.$$

映射输出行后的实际差最大/0.004=1.50856571e-12，体积加权 RMS=4.67775687e-13；均低于 1e-10。矩预测映射差仅 1.51758428e-17。保留历史输入 L_axial=0.001125、L_shear=18h=0.0011250000000000001 的一 ulp 差异，未为求完全重合修改输入。

几何与核相同，G 相同；这个交叉结果将根因定位到运动学算子，不能归因于剪切本构或牵引专属现象。[逐节点映射](AXIAL_SHEAR_GEOMETRY_CROSSCHECK.csv)。

## 14. Psi_plus causal reconstruction

沿用 [svk_spectral_split.py](../../pd_clot/pfpd_candidate/svk_spectral_split.py)：

$$E=\tfrac12(F^TF-I),\quad E_+=\sum_\alpha\max(e_\alpha,0)q_\alpha\otimes q_\alpha,$$
$$\Psi^+(F)=\tfrac\lambda2\langle\operatorname{tr}E\rangle_+^2+\mu\,E_+:E_+.$$

当前候选使用 Green 应变/SVK 谱分裂；公式来源是已冻结实现及其原谱分裂审计，本轮不重新宣称它等同论文的完整 PFPD 模型。节点取 $F_{pred}=F_{exact}+\Delta\bar F_{M3}$，键取 $F_{pred}=F_{exact,mid}+\Delta\widetilde F_{M3}$，完整非线性调用同一谱分裂。预测不是对实际 ΔΨ 拟合，也不采用线性灵敏度作为证据。

| Field | Sample | Actual max bias [Pa] | Predicted max bias [Pa] | max causal residual / Psi_ref |
| --- | --- | --- | --- | --- |
| shear | Psi_Fbar | 2.58288901e-04 | 2.58288901e-04 | 1.04261837e-12 |
| shear | Psi_Ftilde | 2.57886483e-04 | 2.57886483e-04 | 1.01886115e-12 |
| axial | Psi_Fbar | 9.35209793e-04 | 9.35209793e-04 | 2.57680569e-12 |
| axial | Psi_Ftilde | 9.27745038e-04 | 9.27745038e-04 | 2.57669771e-12 |

四组因果残差均远低于 1e-8；因而本轮约 6.4% 的剪切正能量偏差可由重构 F 解释，没有独立谱映射缺陷的证据。这是这些状态上的因果核验，不是对所有谱分裂分支或一般材料的证明。[深度数据](PSI_BIAS_RECONSTRUCTION.csv)。

## 15. Global energy source decomposition

以无向键计数，保留原能量权重 $b_{ij}=2\bar\omega_{ij}V_iV_j$：

$$U_0=\sum_{i<j}b_{ij}W(\widetilde F_{ij}),\quad
U_1=\sum_{i<j}b_{ij}W(F^{exact}_{mid}),\quad
U_2=\int_\Omega W(F^{exact}(X))\,dV.$$

$$\Delta U_{kin}=U_0-U_1,\quad\Delta U_{weight/quad}=U_1-U_2,\quad
U_0-U_2=\Delta U_{kin}+\Delta U_{weight/quad}.$$

U1 只在本构输入中替换 F，未创建新生产模型。正能量类比将 W 换成 Ψ⁺，没有相场演化。下表全为 J。

| Field | Channel | U0 [J] | U1 [J] | U2 [J] | U0-U1 [J] | U1-U2 [J] | U0-U2 [J] |
| --- | --- | --- | --- | --- | --- | --- | --- |
| shear | intact | 6.62365454e-12 | 6.63988403e-12 | 6.64456656e-12 | -1.62294913e-14 | -4.68252790e-15 | -2.09120192e-14 |
| shear | positive | 3.32241093e-12 | 3.33059967e-12 | 3.33296201e-12 | -8.18873672e-15 | -2.36233482e-15 | -1.05510715e-14 |
| axial | intact | 1.32890953e-11 | 1.33223634e-11 | 1.33318127e-11 | -3.32680971e-14 | -9.44920893e-15 | -4.27173060e-14 |
| axial | positive | 1.32890953e-11 | 1.33223634e-11 | 1.33318127e-11 | -3.32681024e-14 | -9.44920893e-15 | -4.27173114e-14 |

剪切 intact 的运动学、权重/积分及总相对偏差分别为 -0.24425207%、-0.07047153%、-0.31472360%。全局负号只表示净能量差，并不意味着各节点正能量偏差都为负。

32 点 Gauss 对照 64 点的最大相对差 3.87785139e-15；四组分解闭合均为 0 J，重复求和确定性检查通过。深度壳层按键中点半 h 分组，连续参照用嵌套中心立方体能量差，详见 [ENERGY_DEPTH_PROFILES.csv](ENERGY_DEPTH_PROFILES.csv)。不能把壳层积分或整体能量吻合等同于局部无偏。

## 16. Primary root cause

在这两个固定二次场上，完整连续球与完整晶格的奇矩为零；截断均匀体积家族产生显著奇矩；物理体积改变其幅度和分布；独立矩恒等式、投影恒等式、跨场行映射和非线性能量重建全部通过。因此：

**HORIZON TRUNCATION → NONZERO ODD MOMENTS → FBAR BIAS → PROJECTED FTILDE BIAS → PSI_PLUS BIAS**。

H1 是主体来源；H2 为可测的次级贡献。H3 的“全晶格仍有显著二次重构误差”及 H4 的“完整 BAQP 对此二次场天然不精确”不成立于本基准；H5/H6 的额外键或谱映射误差解释不获支持；H7 无需作为此次结论。这个结论不扩展为所有非仿射场、任意边界或非规则点云的普遍误差定理。

完整状态标签保存为 [FINAL_STATUS.json](FINAL_STATUS.json)：

| Status | Value |
| --- | --- |
| CONTINUOUS_FULL_BALL_M3 | ZERO |
| DISCRETE_FULL_LATTICE_M3 | NUMERICAL_ZERO |
| FULL_LATTICE_QUADRATIC_REPRODUCTION | PASS |
| FBAR_MOMENT_IDENTITY | PASS |
| FBAR_BIAS_PRIMARY | TRUNCATED_ODD_MOMENT |
| BOUNDARY_VOLUME_WEIGHTING_EFFECT | {"max_scaled": 0.005829069741766859, "weighted_RMS_scaled": 0.00274908710997336} |
| K_ANISOTROPY_MODIFIER | {"continuous_reference_range": [0.996294270729267, 7.517341817553924], "lattice_reference_range": [1.0000047103798373, 7.545338207743268]} |
| FTILDE_PROJECTION_IDENTITY | PASS |
| FTILDE_BIAS_PRIMARY | PROJECTED_NODAL_RECONSTRUCTION_BIAS |
| AXIAL_SHEAR_GEOMETRY_INVARIANCE | PASS |
| PSI_PLUS_CAUSAL_RECONSTRUCTION | PASS |
| PSI_PLUS_BIAS_PRIMARY | KINEMATIC_RECONSTRUCTION |
| GLOBAL_ENERGY_KINEMATIC_BIAS | {"shear": {"intact": -1.6229491340158054e-14, "positive": -8.188736722203713e-15}, "axial": {"intact": -3.326809706507113e-14, "positive": -3.326810243932658e-14}} |
| GLOBAL_ENERGY_WEIGHT_QUADRATURE_BIAS | {"shear": {"intact": -4.68252789645226e-15, "positive": -2.362334824442443e-15}, "axial": {"intact": -9.449208932394328e-15, "positive": -9.449208932394328e-15}} |
| SURFACE_OPERATOR_CORRECTION_READINESS | READY_FOR_BAQP_SPECIFIC_DESIGN |
| BOUNDARY_DATA_CONVERSION | UNCHANGED_VERIFIED |
| C5_C7_ACTIVATION | NOT_RUN |
| CLOT_GEOMETRY | NOT_RUN |
| FRACTURE | NOT_RUN |

## 17. What is and is not justified about surface correction

证据足以进入 BAQP 专属算子修正的设计，但尚未选取或验证任何修正。

Shojaei 等（2026）讨论截断非局部算子向完整 horizon PD 算子的匹配；其 Eq.33/40 分别展示完整/半域奇矩行为。这为本轮诊断提供动机。其 BB/OSB 算子与权重目标具有模型依赖性，不能直接作为 BAQP 权重。来源：[原始出版 PDF](https://tore.tuhh.de/bitstreams/8fda7c78-cb40-401e-b8ce-206cd986fafb/download)，[DOI](https://doi.org/10.1002/nme.70385)，[缓存哈希及阅读位置](LITERATURE_PROVENANCE.json)。本轮的三维 M3/K 及 Ftilde 恒等式由本项目 BAQP 定义直接推导。

公开原始全文的既有缓存经 SHA256 复核，本轮网页打开受工具 10 MiB 限制；未重新下载 50 MB 文件，也未声称与用户上传 PDF 逐字节核对。没有应用 perifit、BB/OSB 权重、ghost/mirror、修改 K/M3 或改变物理 horizon。

## 18. BAQP-specific future correction requirements

设计义务完整列在 [BAQP_OPERATOR_CORRECTION_TARGETS.md](BAQP_OPERATOR_CORRECTION_TARGETS.md)：仿射精确、K 可逆与条件性、适当的二次奇矩目标、互易性、共同能量一阶变分、保持物理 horizon、区分边界数据转换、非规则点云可用性、不针对 C5/C7 首激活调参。

还需正视一个直接代数约束：平直单侧边界上，若所有权重非负且法向 ξ 全同号，则法向三阶矩不可能在法向二阶矩为正时变为零。仅靠现有家族内正标量权重强制“所有奇矩为零”可能不可行。此处是 BAQP 设计约束推导，不能误引为论文结论；见设计说明中的条件和证明。未实现任何替代构造。

## 19. Implication for future clot-wall damage initiation

当前候选的损伤驱动力依赖 Ψ⁺，因此真实 clot-wall 附近的首激活可能同时受真实应变和截断重构偏差影响；这是基于本轮算子证据的后续风险推断。不能据此断言现存 C5/C7 热点已完全解释、首激活应变已修正、mesh objectivity 已恢复，或真实血栓将在哪里断裂。

本轮没有 clot 几何、相场、断裂、CFD 或 activation 计算。未来须在完成算子设计并单独验证后，再区分表面误差与物理驱动；旧 C5/C7 科学状态保持原样。

## 20. One next task only

唯一建议下一任务：**BAQP_SURFACE_OPERATOR_CORRECTION_DESIGN_V1**。先研究 BAQP 多项式一致性目标、单侧正权重可行性及统一能量—力约束，形成可审阅设计。本轮不执行该任务。

10 项指定定向测试全部通过（0.67 s）；未运行全量 pytest。主诊断使用最多 8 CPU workers，BLAS=1，墙钟 17.397 s；采样峰值进程树 RSS=1134133248 B，主进程高水位=1148493824 B，较大者 1.070 GiB。最小可用系统 RAM 比例 57.80%（要求 ≥20%）；GPU_NOT_USED_BY_DESIGN。资源值对应科学诊断阶段，绘图和报告不包含在这段墙钟中。

425 个冻结输入/旧源码/历史结果的 SHA256 在计算前后保持不变，最终检查另见 FINAL_ARTIFACT_QA.json；新模块仅作诊断，不接求解器。旧 manifest、旧 assertions、kernel、volumes、材料、谱分裂、能量、边界转换均未改。分支 `clot-clearance/pfpd-progress-20261005`，HEAD `695e75166aa864809460ea9a1845fe72abce229c`；不提交、不推送。

交付索引：[数据字典](DATA_DICTIONARY.md)、[六组图](../../results/baqp_near_boundary_kinematic_energy_audit_v1/README.md)、[定向测试](TARGETED_TESTS.json)、[历史完整性](HISTORICAL_INTEGRITY.json)、[资源记录](RESOURCE_USAGE.json)。所有图为英文 600 dpi PNG 与矢量 PDF；散点抽样规则见 FIGURE_SAMPLING.json，全部 gate 使用全量数据。

**停止并等待人工审核；不自动开始修正或后续科学计算。**
