# HIGHER_ORDER_BAQP_BOUNDARY_REINTEGRATION_V1

**VERIFIED_FOR_INTACT_BENCHMARKS**。分支 `clot-clearance/pfpd-progress-20261005`；HEAD `695e75166aa864809460ea9a1845fe72abce229c`。

[任务原文](TASK_REQUEST.txt) · [冻结配置](PD2_BOUNDARY_CONFIGURATION.json) · [最终状态](FINAL_STATUS.json) · [图件](../../results/higher_order_baqp_boundary_reintegration_v1/README.md)

## 1. Why boundary reintegration is now the correct task

上一轮已验证完整 PD2 重构、离散能量与其负梯度内力。本轮检验把该内力接回既有边界数据转换后，边界积分、闭合体域局部化和扰动起步平衡恢复是否同时成立。只在完整态 SVK 候选中进行，未引入新材料、边界模型或多项式场。

来源：[已接受的能量—力报告](../higher_order_baqp_energy_force_candidate_v1/REVIEW_HIGHER_ORDER_BAQP_ENERGY_FORCE_V1.md)、[仿射边界报告](../nonlocal_boundary_formulation_v1/REVIEW_NONLOCAL_BOUNDARY_FORMULATION_V1.md)、[非仿射轴向报告](../nonlocal_boundary_nonaffine_localization_v1/REVIEW_NONAFFINE_BOUNDARY_LOCALIZATION_V1.md)、[非仿射剪切报告](../nonlocal_boundary_nonaffine_shear_localization_v1/REVIEW_NONAFFINE_SHEAR_LOCALIZATION_V1.md)。

## 2. What is reused unchanged: PD2 mechanics and boundary architecture

主网格 $N=19^3=6859$，690020 条无向键，1380040 条有向邻接；$L=18h=1.125\,mm$，$h=62.5\,\mu m$，$\delta=251.25\,\mu m$，$\delta/h=4.02$。内部/面/棱/角体积为 $h^3,h^3/2,h^3/4,h^3/8$。核为 $(1-r/\delta)^2$；$\mu=1000\,Pa,\lambda=0$。

接受的梯度缓存 SHA256 为 `ffc7064e7b35b8ffabe52063ae6d14340ba76f4308c9881f23a815ba69634420`。完整九项 P2 权重与参考键权重逐字节复用：
$$F_i=I+\sum_j(u_j-u_i)\otimes g_{ij},\quad \widetilde F_b=\tfrac12(F_i+F_j)R_b+(x_j-x_i)\otimes\xi_b/r_b^2,$$
$$R_b=I-\xi_b\otimes\xi_b/r_b^2,\quad U=\sum_{i<j}2\bar w_{ij}V_iV_j\Psi_0(\widetilde F_b),\quad f_{int}=-\partial U/\partial x.$$
这些公式和一阶伴随推导的直接来源是上一轮报告第 3–7 节及未修改的 `higher_order_rk` 源码。PD2/RK 文献归属见已接受 RK 报告及其缓存文献；本轮未宣称新的外部定理或稳定性证明。

本轮仅在既有 `nonlocal_boundary` 内新增 backend 和选择性 PD2 oracle；`localized_conversion.py` 添加调用接口，`mechanics_adapter.solve` 添加可选 backend 参数。默认 C1 路径不变。边界函数只接收 `energy/force/node_F/bond_Ftilde/psi_plus_diagnostic` 接口；仍调用原 `localize` 代数。`source_before` 保存两份扩展前源码，详见 [扩展来源记录](SOURCE_EXTENSION_PROVENANCE.json)。旧 manifest 未改；其中的旧源码 hash 对应归档版本，不声称扩展后仍匹配旧 hash。

三个完全相同目标数组的 PD2 力/能量与上一轮逐位相同。法向仿射采用原边界基准的 `X @ F.T` 表达式，与上轮制造场表达式有浮点差异；其接受源码和缓存仍受 hash 保护。没有通过改变目标表达式去制造逐位一致。定向测试1另在保存的旧法向目标坐标上求值，也验证了力和能量逐位相同。

## 3. PD2 force-dependency closure

从实际伴随图反向追踪：观测点力需要该点直接键项以及所有 source/target 指向它的 $G_i g_{ik}$。因此 D1 是观测点与 256 个邻点；每个 $G_i$ 又需要所有邻接键的应力，端点梯度位于 D2；D2 梯度和键权重的 $\omega_{0,i}$ 必须读取 D3 的位置、体积和完整邻域。

$$f_p\leftarrow G_{D1},P_{incident(p)},g_{D1}\leftarrow F_{D2},m_{D2}\leftarrow (X,V,x)_{D3}.$$
本报告 D1 包含观测点，因此计数为 257；旧报告仅列 256 个直接邻点，差异是集合定义，不是新增邻接。

| 集合 | 节点数 | 最大坐标偏移/h | 最大半径 [m] |
| --- | --- | --- | --- |
| D0 | 1 | 0 | 0.000000e+00 |
| D1 | 257 | 4 | 2.500000e-04 |
| D2 | 1797 | 8 | 5.000000e-04 |
| D3 | 5821 | 12 | 7.500000e-04 |

N19 不闭合；N25 虽含所有所需节点，但 D3 外层仍遇到半体积边界点。N27 的边界在 13h，而 D3 最大坐标偏移为 12h，故完整闭合。依据：[PD2_FORCE_STENCIL_CLOSURE.json](PD2_FORCE_STENCIL_CLOSURE.json)。

## 4. Fully closed PD2 bulk oracle

每个物理 Z 层在独立 N27 模板中心计算 PD2 体域力密度。模板扩大时制造场分母L仍固定为N19的1.125 mm，不随N27/N29域大小重新缩放。选择性计算调用已接受的 `local_system`、`deformation`、`energy_force` 和伴随函数；只保留对观察点必要的键与有向回传，不使用 C1 的 K/Z/力或 q_bulk 数值。选择性模板的其余节点力与能量是部分和，不作为全域结果。

19 层分别计算两个非仿射场；各三个选定层另用偏移 N29 模板核对。N19 的三个选定点再与完整接受求值器比较。[独立验证](PD2_ORACLE_VALIDATION.json) 全部通过。仿射体域先行门槛：

| 场 | \|q\|h³/F_ref | 门槛 |
| --- | --- | --- |
| normal_affine | 1.838261e-15 | 1e-10 |
| shear_affine | 1.043280e-15 | 1e-10 |

非仿射体域差异：

| 场 | max \|e_model\| [N/m³] | RMS [N/m³] | max/b_ref | RMS/b_ref |
| --- | --- | --- | --- | --- |
| nonaffine_axial | 7.520789e-05 | 7.520053e-05 | 2.090091e-08 | 2.089886e-08 |
| nonaffine_shear | 7.520065e-05 | 7.520032e-05 | 4.229698e-08 | 4.229679e-08 |

这里 RMS 是 19 个 Z 层的等权 RMS。非零差异是固定 h、δ 下的体域模型/离散差异，不是边界污染。

## 5. Affine normal boundary conversion

沿用 $F=I+0.003e_z\otimes e_z$、六个参考法向、z± 两面位移约束和原虚位移集合。体域 q 在数值精度内为零，故 $f_{BC}=qV-f_{FD}\simeq-f_{FD}$。为保留可核查的体域数据，本轮使用测得的微小 q，不人为截零。两个仿射平衡也使用共同 NL 载荷接口，因此实际外载为 $f_{BC}-qV$；这里 $-qV$ 仅为仿射 oracle 的舍入量级，不是非零物理体力。该实现与解析 $q=0$ 的仿射边界载荷在数值精度内一致，但不声称实际数组逐位等于保存的 $f_{BC}$。实际外载及其独立力/力矩/虚功核算见 [舍入核算](PD2_AFFINE_LOAD_ROUNDOFF_ACCOUNTING.json)，原始载荷未事后改写。法向面牵引仍为 $t=PN$。所有面法向、牵引及积分保存在 [目标摘要](PD2_TARGET_SUMMARIES.json)。

## 6. Affine shear boundary conversion

沿用 $F=I+0.003e_x\otimes e_z$。令剪切量为 k，完整 SVK 的
$$P=\mu\begin{pmatrix}k^2&0&k+k^3\\0&0&0\\k&0&k^2\end{pmatrix}.$$
六个面的 $PN$ 全部保留，包括 x 面切向牵引与小的法向分量。对固定应力拆分 $P_n=\operatorname{diag}(P)$、$P_t=P-P_n$，经接受的力伴随映射计算各载荷通道。没有给分布式节点力规定任意物理面归属；该分解也不意味着非线性解可叠加。

通道和最大误差/F_ref = 1.670052e-13。

## 7. Force / moment / virtual-work consistency

复用旧 `load_checks`：合力误差除以 $\|P\|_F L^2$；参考坐标力矩误差除以 $\|P\|_F L^3$。参考力矩不强制为零，剪切 P 一般非对称。九个仿射虚位移和三个平移的非零目标虚功使用真正相对误差；零目标使用既有稳定应力—体积尺度。
$$\sum_i f_{BC,i}\cdot(AX_i+c)\quad\leftrightarrow\quad |\Omega|P:A.$$
各项门槛均为 $10^{-6}$，无事后收紧。

| 场 | 通道 | 合力 | 参考力矩 | 最大虚功误差 |
| --- | --- | --- | --- | --- |
| normal_affine | combined | 3.327079e-14 | 1.530414e-17 | 3.324610e-14 |
| shear_affine | combined | 1.983462e-14 | 1.604017e-15 | 1.529294e-13 |
| shear_affine | normal | 1.833347e-15 | 3.266828e-17 | 3.565950e-16 |
| shear_affine | tangential | 1.125286e-15 | 8.798725e-16 | 2.738638e-16 |

非仿射另沿用 [全局核算](PD2_NONAFFINE_ACCOUNTING.json) 和 [虚功分解](PD2_NONAFFINE_VIRTUAL_WORK.csv)，保留体积积分误差、内部参考力矩和连续体力矩源项，未将转换力强制等同于局部表面积分。

## 8. Affine equilibrium recovery

复用原 `regularization.equilibrium._scipy` 的 TRF/LSMR，不编写新 Newton 或 ADR。新增 callback 用接受的 PD2 函数计算每次内力残差；解析方向导数只供原求解器线性步骤使用，采用稀疏线性映射求和。未组装全局 Hessian，也未做稳定性分析。方向导数独立中心差分相对误差 4.389478e-10，互易性误差 7.058115e-16。

五个算例使用同一冻结配置：max_nfev=80，LSMR maxiter=120，xtol=ftol=gtol=1e-12，LSMR atol=btol=1e-10，单例墙钟上限3600 s。前三个求解串行；最后两个剪切算例在第42条允许的实测内存证明后并行，求解配置、载荷及停止条件均不变。原两 z 面固定，其余 18411 个自由度。非零扰动为 0.02×0.003×L×Πsin(πq)×(1,−0.4,0.7)，固定面恢复精确目标。

| 场 | 机械残差 | max ΔFbar/尺度 | max ΔFtilde/尺度 | max ΔΨ⁺/Ψ_ref | max Δu [m] |
| --- | --- | --- | --- | --- | --- |
| normal_affine | 5.191387e-10 | 1.186512e-07 | 1.295678e-07 | 2.400850e-07 | 6.376315e-14 |
| shear_affine | 8.658278e-11 | 2.438900e-08 | 2.435422e-08 | 2.168188e-08 | 1.197591e-14 |

机械门槛 1e-8；两种梯度门槛 1e-5；Ψ⁺门槛 1e-4。仿射梯度误差尺度0.003。固定节点外载照样保存，反力单独保存，未混入自由残差。

## 9. Nonaffine axial localized conversion

复用 $u_z=0.003Z+0.001Z^2/L$，其他位移分量为零。令 $a=1.003+0.002Z/L$，
$$F=\operatorname{diag}(1,1,a),\quad P_{zz}=\mu(a^3-a),\quad b_z=-\mu(3a^2-1)(0.002/L).$$
闭式参考来自已接受 `QuadraticAxial`，没有通用本构自证参考。转换仍是
$$f_{BC,i}=q_{bulk}(Z_i)V_i-f_{FD,i}(x_L).$$
保留完整物理节点体积；$-f_{FD}$ 仅保存在 naive 诊断通道，不作为该场的边界转换。

## 10. Nonaffine shear localized conversion

复用 $u_x=0.003Z+0.001Z^2/L$；$k=0.003+0.002Z/L$，$k'=0.002/L$。应力矩阵同第6节但 k 随 Z 变化，
$$b_{local}=\big[-\mu(1+3k^2)k',\ 0,\ -2\mu kk'\big].$$
横向位移并不意味着 z 体力为零。使用既有 `QuadraticShear` 的全部面牵引和精确 Gauss 参考。N27 oracle 在每个绝对 Z 上计算，不把 N19 几何中心误当作闭合体域。没有增加 XZ 或其他制造场。

## 11. Nodewise bulk/boundary decomposition

逐节点保留单位明确的五个量：$f_{FD}$、$q_{bulk}V$、$f_{BC}$、$b_{local}V$、$e_{model}V$，其中 $e_{model}=q_{bulk}+b_{local}$ 是密度，单位 N/m³。力数组单位 N。
$$f_{FD}+f_{BC}+b_{local}V=e_{model}V.$$
该代数检查不是边界物理定律的独立验证；必须与闭合体域局部化和实际扰动恢复合看。

| 场 | 最大节点闭合误差/F_ref |
| --- | --- |
| normal_affine | 6.464299e-17 |
| shear_affine | 4.904983e-17 |
| nonaffine_axial | 5.168965e-17 |
| nonaffine_shear | 6.951574e-17 |

全部低于预声明1e-10。逐节点源数据见 [NPZ](PD2_NODEWISE_BOUNDARY_DECOMPOSITION.npz)。

## 12. Fully closed center localization

N27 从侧面沿 X 方向取4h、8h、12h和13h深度点，保持Z=0。独立偏移N29供给体域参考。

| 场 | 深度/h | \|f_BC\|/F_ref | 实际图闭合 |
| --- | --- | --- | --- |
| nonaffine_axial | 4 | 3.101696e-04 | False |
| nonaffine_axial | 8 | 2.506949e-07 | False |
| nonaffine_axial | 12 | 1.999689e-15 | False |
| nonaffine_axial | 13 | 2.058890e-15 | True |
| nonaffine_shear | 4 | 2.187671e-02 | False |
| nonaffine_shear | 8 | 4.906529e-06 | False |
| nonaffine_shear | 12 | 2.263162e-15 | False |
| nonaffine_shear | 13 | 3.571868e-15 | True |

12h 处数值很小仍不能称完全闭合，因为依赖图读到边界体积；只有13h中心满足完整图条件。两种中心结果均通过1e-10门槛。这里“约1δ/2δ/3δ”是测点深度，不是指定物理边界厚度。

## 13. Nonlocal-manufactured recovery

使用 $b_{NL}=-q_{bulk}$，冻结 $f_{ext}=b_{NL}V+f_{BC}$；保存的两项分开核算。虽然该组合在目标上等价于 $-f_{FD}(x_L)$，恢复结果来自非零扰动的完整平衡求解，不是目标代入。PD2 再现二次场，故两种梯度均直接与解析 F 比较；键梯度的解析参考位于参考键中点。

| 场 | 机械残差 | max ΔFbar/尺度 | max ΔFtilde/尺度 | max ΔΨ⁺/Ψ_ref | max Δu [m] |
| --- | --- | --- | --- | --- | --- |
| nonaffine_axial | 1.150552e-12 | 6.089131e-13 | 7.658356e-13 | 1.311080e-12 | 1.626303e-19 |
| nonaffine_shear | 5.778522e-11 | 4.489059e-10 | 5.076200e-10 | 8.738645e-10 | 4.022389e-16 |

梯度尺度0.004；Ψ_ref取解析节点场最大Ψ⁺，节点和键诊断共用同一尺度。原谱分裂仅作读出，没有改谱分裂或引入相场。

## 14. Classical-shear result

第五个且唯一经典载荷求解使用 $f_{ext}=b_{local}V+f_{BC}$，没有加入 $-e_{model}V$。求解前独立验证目标残差等于 $e_{model}V$。这类平衡不要求精确恢复局部目标。

| 场 | 机械残差 | max ΔFbar/尺度 | max ΔFtilde/尺度 | max ΔΨ⁺/Ψ_ref | max Δu [m] |
| --- | --- | --- | --- | --- | --- |
| classical_shear | 3.780054e-12 | 2.223610e-08 | 2.434597e-08 | 5.534395e-08 | 1.693642e-14 |

| 交叉核算 | 结果 |
| --- | --- |
| NL_to_classical_displacement_max_m | 1.672652e-14 |
| NL_to_classical_Fbar_max | 8.794387e-11 |
| external_difference_equals_model_scaled | 1.278888e-16 |
| force_response_closure_relative_to_model | 1.400220e-02 |
| NL_recovery_over_classical_shift | 2.374993e-02 |
| classical_target_residual_on_solver_scale | 4.093433e-09 |

解析目标自身的机械残差可能已低于宽松的1e-8接受门槛；因此不能从“converged”推出精确局部解。表中同时量化 NL 数值恢复误差与经典偏移、两解内力差对模型载荷的响应闭合。剩余偏移按固定视野体域差异、积分与求解精度共同解释，不称作边界失败，也不赋予每一位小数物理意义。

## 15. PD2 boundary-load support

**CONVERSION SUPPORT / NOT PHYSICAL BC THICKNESS**。分位数按Σ|f_BC|累积，深度是到全部六个物理面的最近距离。

| 场 | d50 [µm / h / δ] | d90 | d95 | d99 |
| --- | --- | --- | --- | --- |
| normal_affine | 0 / 0 / 0.000000 | 125 / 2 / 0.497512 | 125 / 2 / 0.497512 | 250 / 4 / 0.995025 |
| shear_affine | 0 / 0 / 0.000000 | 125 / 2 / 0.497512 | 125 / 2 / 0.497512 | 250 / 4 / 0.995025 |
| nonaffine_axial | 0 / 0 / 0.000000 | 125 / 2 / 0.497512 | 125 / 2 / 0.497512 | 250 / 4 / 0.995025 |
| nonaffine_shear | 0 / 0 / 0.000000 | 125 / 2 / 0.497512 | 125 / 2 / 0.497512 | 250 / 4 / 0.995025 |

每层 mean/max、节点数、累积幅值分数保存在 [PD2_BOUNDARY_SUPPORT.csv](PD2_BOUNDARY_SUPPORT.csv)。N19 的较小尾部不代替 N27 的闭合证明。

## 16. How PD2 boundary load differs from historical C1

C1仅从已接受仿射与非仿射NPZ读取，核对坐标完全一致；未重跑C1大算例。12项定向测试中的一个小N7兼容性检查用于验证接口保持旧转换。

| 场 | 后端 | Σ\|f\| [N] | L2 [N] | max \|f\| [N] | C1−PD2 L2 [N] |
| --- | --- | --- | --- | --- | --- |
| normal_affine | C1 | 2.126686e-05 | 4.981792e-07 | 1.810829e-08 | 4.840524e-07 |
| normal_affine | PD2 | 2.168743e-05 | 6.286481e-07 | 2.741451e-08 | 4.840524e-07 |
| shear_affine | C1 | 1.963689e-05 | 3.506894e-07 | 1.060685e-08 | 3.407449e-07 |
| shear_affine | PD2 | 2.031398e-05 | 4.425319e-07 | 1.378381e-08 | 3.407449e-07 |
| nonaffine_axial | C1 | 2.135717e-05 | 5.214625e-07 | 2.342710e-08 | 5.064632e-07 |
| nonaffine_axial | PD2 | 2.177613e-05 | 6.602726e-07 | 3.618409e-08 | 5.064632e-07 |
| nonaffine_shear | C1 | 1.964397e-05 | 3.615995e-07 | 1.361848e-08 | 3.512865e-07 |
| nonaffine_shear | PD2 | 2.032496e-05 | 4.572109e-07 | 1.817004e-08 | 3.512865e-07 |

C1的分位数为1/2/3/4h，PD2为0/2/2/4h。PD2最大载荷更集中于表面，但非局部转换依然存在。两种载荷的差只表示同一替代场下离散力学表示不同；不能唯一拆成“C1重构补偿”和“真实边界作用”，不能称为物理修正力，更未把它添加为新的模型项。

## 17. C1 vs PD2 F/Ftilde/Psi_plus near-boundary quality

按本轮要求，比较历史 C1 **目标状态诊断** 与 PD2 **制造载荷求解状态**，而非重跑C1。C1标量数组直接来自近边界运动学审计；两者使用匹配的节点/键中点深度、0.004梯度尺度和相同Ψ_ref。

| 场 | 数据来源 | 采样 | max ΔF/0.004 | max \|ΔΨ⁺\|/Ψ_ref |
| --- | --- | --- | --- | --- |
| nonaffine_axial | C1_TARGET | Fbar | 2.975839e-02 | 5.821751e-02 |
| nonaffine_axial | C1_TARGET | Ftilde | 2.942854e-02 | 5.775282e-02 |
| nonaffine_axial | PD2_SOLVED | Fbar | 6.089131e-13 | 9.975916e-13 |
| nonaffine_axial | PD2_SOLVED | Ftilde | 7.658356e-13 | 1.311080e-12 |
| nonaffine_shear | C1_TARGET | Fbar | 2.975839e-02 | 6.431445e-02 |
| nonaffine_shear | C1_TARGET | Ftilde | 2.942854e-02 | 6.421425e-02 |
| nonaffine_shear | PD2_SOLVED | Fbar | 4.489059e-10 | 8.738645e-10 |
| nonaffine_shear | PD2_SOLVED | Ftilde | 5.076200e-10 | 8.315667e-10 |

结论据实际求解值和冻结门槛判定；目标PD2接近机器精度不被当成求解结果。图5为主要近边界比较，图6另列目标与求解后的误差曲线。[深度数据](C1_PD2_FIELD_QUALITY.csv)。

## 18. What remains quadrature rather than boundary error

PD2恢复二次梯度不使有限域键积分自动等于连续体积分。上一轮已把重构贡献与键/体积积分误差分开；本轮沿用原积分规则。局部体力梯形积分、键能量积分和有限视野体域差异均保留，不通过边界补偿调零。

| 场 | U_target [J] | U_solved [J] | Σf_ext·u_target [J] | Σf_ext·u_solved [J] | ΔU_solved-target [J] |
| --- | --- | --- | --- | --- | --- |
| normal_affine | 1.285293e-11 | 1.285293e-11 | 2.574435e-11 | 2.574435e-11 | -1.740957e-23 |
| shear_affine | 6.407255e-12 | 6.407255e-12 | 1.281457e-11 | 1.281457e-11 | -6.988222e-24 |
| nonaffine_axial | 1.332236e-11 | 1.332236e-11 | 2.668743e-11 | 2.668743e-11 | 5.977672e-26 |
| nonaffine_shear | 6.639884e-12 | 6.639884e-12 | 1.327984e-11 | 1.327984e-11 | -3.106330e-22 |
| classical_shear | 6.639884e-12 | 6.639884e-12 | 1.327984e-11 | 1.327984e-11 | -2.001059e-20 |

表中外载功是固定节点载荷在对应位移的功，未包含未知加载路径的夹持反力功；不要求其等于储能。额外的曲面/视野/网格极限未在本轮验证。

## 19. Whether the intact higher-order boundary mechanics are ready for a curved clot geometry

状态为 **READY**，指具备进入下一轮曲面完整态验证的资格，不代表曲面血栓已验证。仅这四个完整态恢复基准及经典剪切控制被计算。

| 最终状态 | 值 |
| --- | --- |
| PD2_VARIATIONAL_MECHANICS | REUSED_VERIFIED |
| PD2_BOUNDARY_BACKEND | IMPLEMENTED |
| PD2_FORCE_STENCIL_CLOSURE | VERIFIED |
| PD2_AFFINE_BULK_FORCE | PASS |
| PD2_NONAFFINE_BULK_ORACLE | VERIFIED |
| PD2_NODEWISE_BOUNDARY_DECOMPOSITION | PASS |
| PD2_NORMAL_RESULTANT | PASS |
| PD2_NORMAL_MOMENT | PASS |
| PD2_NORMAL_VIRTUAL_WORK | PASS |
| PD2_NORMAL_EQUILIBRIUM | PASS |
| PD2_SHEAR_RESULTANT | PASS |
| PD2_SHEAR_MOMENT | PASS |
| PD2_SHEAR_VIRTUAL_WORK | PASS |
| PD2_SHEAR_EQUILIBRIUM | PASS |
| PD2_NONAFFINE_AXIAL_RECOVERY | PASS |
| PD2_NONAFFINE_SHEAR_RECOVERY | PASS |
| PD2_CLASSICAL_SHEAR | CONSISTENT_WITH_BULK_MODEL_DIFFERENCE |
| PD2_CLOSED_CENTER_BOUNDARY_CORRECTION | PASS |
| C1_PD2_BOUNDARY_LOAD_COMPARISON | COMPLETE |
| PD2_NEAR_BOUNDARY_F_BIAS | REMOVED_TO_SOLVER_TOLERANCE |
| PD2_NEAR_BOUNDARY_PSI_PLUS_BIAS | REMOVED_TO_SOLVER_TOLERANCE |
| HIGHER_ORDER_BOUNDARY_MECHANICS | VERIFIED_FOR_INTACT_BENCHMARKS |
| CURVED_CLOT_WALL_READINESS | READY |
| PRODUCTION_BAQP | UNCHANGED |
| PHASE_FORCE | NOT_IMPLEMENTED |
| C5_C7 | NOT_RUN |
| FRACTURE | NOT_RUN |
| CLOT_GEOMETRY | NOT_RUN |

保护的旧文件 675 项，变化 0 项。生产BAQP、接受的PD2、谱分裂、旧结果、旧断言及manifest未修改；两项允许扩展的边界源码已单独归档。仅运行12项本轮测试，见 [TARGETED_TESTS.log](TARGETED_TESTS.log)。

资源：8个CPU oracle workers；前三个平衡串行，最后两个在内存证明后并行，最多2个equilibrium workers；BLAS=1，GPU_NOT_USED_BY_DESIGN；扣除并行重叠后的科学阶段墙钟 4311.65 s，峰值进程树RSS 2.436 GiB，最低可用RAM 43.79%。进程树RSS含重复共享页；同时监测系统可用内存，始终要求≥20%。并行前按实测单例峰值的1.5倍预留，证明记录见 [并行决策](PARALLEL_EXECUTION_DECISION.json)，实际两个工作进程合计RSS与可用内存见 [并行记录](PARALLEL_EXECUTION_SUMMARY.json)。原串行调度父进程被暂停以防重复提交，并仅在两个科学工作进程均结束后退出；未终止科学计算。五个算例各运行一次。未安装新依赖。

## 20. One next task only

人工审核通过后，唯一建议任务为 **CURVED_CLOT_WALL_INTACT_MECHANICS_V1**。

本轮到此停止，等待人工审核；未开始 curved clot、phase、AT1、Gc、C5/C7、stability、fracture、Ogden、adhesion、CFD 或 microbubble。
