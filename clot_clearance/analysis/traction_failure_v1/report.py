from analysis.traction_failure_v1.common import *
import html,subprocess,importlib.metadata as md
def load(p):return json.loads(Path(p).read_text())
def write(p,s):Path(p).write_text(s)
def main():
 q=load(B/'LOCAL_MAPPING_QUALITY.json');au=load(B/'CONSERVATIVE_MAPPING_AUDIT.json');sf=load(B/'CFD_SURFACE_FORCE_AUDIT.json');ga=load(B/'GEOMETRY_AND_UNIT_AUDIT.json');iv=load(B/'INDEPENDENT_VERIFICATION.json');a=load(A/'SURFACE_CORRECTION_FEASIBILITY.json');me=load(A/'FORCE_ENERGY_COMPATIBILITY.json');pw=load(A/'LPS_INTERACTION_WORK_AUDIT.json');nd=load(A/'SURFACE_NEIGHBORHOOD_AUDIT.json')
 assert iv['passed'];nextstep='只做一个现有 V3 完整态、禁止断键的原生静态受载接口验证：施加本轮含栓 CFD 节点力，核对单位、总载荷、约束反力/力矩与小变形响应，并明确排除既有脱离碎片阻力的重复作用。先完成这一项，不开展动态断裂、材料升级或 AI 扫描。'
 write(B/'REVIEW_Y_CLOT_CFD_TO_LAMMPS_CONSERVATIVE_TRACTION_MAP_V1.md',f'''# 含栓 CFD 到 LAMMPS 离线载荷映射审查

Run `{RUN}`。**B 线通过冻结的离线守恒与局部工程筛查，可以提供未来接口验证用的节点力；没有执行生产耦合，也没有认证逐点载荷收敛或材料行为。**

1. **真实输入。** 重新读取上一轮 `20261011T154854Z/refined` 的原生体网格、`result_131.vtu` 和全精度 `CLOT_SURFACE_TRACTIONS.npz`。速度/压力 float64；原生坐标按精确 float32 舍入双射对应冻结 double 几何。文件 SHA、3,611 个 tag=5 面、相邻四面体、梯度、分项牵引均重新核对，非报告均值反推。`INPUT_SHA.json`、`NATIVE_SURFACE_IDENTITY.json` 保存证据。

2. **配准。** 原 `body.data` 的 3,488 个 ID/坐标与几何数组完全一致；材料 ID 1–1927，锚点 ID 1928–3488。原世界坐标 m、x≈−6 mm、母管 +z 内壁、单位变换矩阵 I；不作平移、旋转、缩放。连续、CFD 三角化、材料积分体积分别为 0.271574388、0.271576593、0.240875 mm³；离散比连续小约 11.30%。没有按体积比放大任何力。原附着内部面排除，暴露三角面坐标集合与原 tag=5 完全相同。

3. **单位与符号。** 长度 m、面积 m²、牵引 Pa、面元/节点力 N、力矩 N·m。相邻流体四面体确认 n_fluid 外法向；n_clot=−n_fluid。σ=−pI+µ(∇u+∇uᵀ)，µ=0.001 Pa·s，t=σn_clot。保存压力、切向黏性、法向黏性与总项；法向黏性未忽略。p 是原出口零参考牵引条件下的表压，开放面合力会随压力参考偏移变化。例如统一 +1 Pa 会增加约 {np.array(sf['open_surface_pressure_offset_1Pa_delta_force_N'])*1e6} µN；它不是绝对血压或壁面反力。

4. **原始积分。** 总力 {sf['total_force_N']} N，模长 {np.linalg.norm(sf['total_force_N'])*1e6:.9f} µN；关于实际材料质心 {sf['reference_point_m']} m 的力矩 {sf['total_moment_Nm']} N·m，恢复旧审计。面元中心力适合合力；线性压力与位置相乘后是二次量，力矩采用原三点二阶三角积分。压力的“精确积分−中心近似”为 {sf['centroid_to_exact_moment_difference_Nm']['pressure']} N·m。独立解析三角公式再核对，未改法向或表面以凑合力。

5. **哪些点受载。** 距实际暴露三角面 ≤50 µm 的 580 个材料点；1,347 个内部材料点外载为零，1,561 个锚点完全排除。所有面均在 ≤150.5 µm 支持内找到邻点，Wendland C2 权重非负、逐面归一化。表面到材料点最近距离中位/最大约 34.44/131.48 µm；仅 1.185% 表面积最近距离超过 100 µm。有效面积 A_i=Σ_f w_if A_f，总和等于真实暴露面积 {ga['exposed_area_m2']*1e6:.9f} mm²。A_i/dx²范围 {ga['effective_area_over_dx2']['min']:.3f}–{ga['effective_area_over_dx2']['max']:.3f}，不是统一 dx²。

初次支持检查的切平面条件排除了 10 个接触边缘面，占面积 0.0232621%；尚未生成节点力即停止。唯一一次局部修正对这十个面使用实际三角面邻接图路径 ≤150.5 µm，替代把曲面局部切平面误当全局内部判据的条件。这里的邻接距离按相邻面中心间弦长累加，不是严格曲面测地距离，也不保证其上界；冻结文件中过强的上界措辞另以 `REPAIR_GRAPH_METRIC_CLARIFICATION.json` 更正留证，算法与数值不变。距离上限、法向一致性、点筛选和全部质量门槛不变。没有远距离补点、删面、半径扫描或重新计算 CFD。原拒绝结果和修正前冻结文件保留。

6. **保力误差。** 归一化尺度 Σ_f|F_f|，压力/剪切/法向黏性分别为 {iv['components']['pressure']['force_error_scaled']:.3e} / {iv['components']['shear']['force_error_scaled']:.3e} / {iv['components']['viscous_normal']['force_error_scaled']:.3e}，总项为 {au['components']['total']['corrected']['force_scaled_error']:.3e}。这些是浮点代数误差，不是物理精度。

7. **保力矩误差。** 尺度为上述 Σ_f|F_f| 乘 L={float(np.linalg.norm(np.ptp(np.load(B/'CLOT_NODAL_TRACTION_FORCES.npz')['reference_coordinates_m'],axis=0))):.9g} m。总项为 {au['components']['total']['corrected']['moment_scaled_error']:.3e}；独立解析积分核查各分项最大 {max(v['moment_error_scaled'] for v in iv['components'].values()):.3e}，远低于冻结 1e−8。还通过移动参考点的力矩平移恒等式。

8. **分项守恒。** 压力、切向黏性、法向黏性分别做一次相同算法，最后相加，不以总项抵消分项误差。初始非负 W 已保力，但有力矩偏差。六约束均满秩；长度归一化 Gram 条件数为 26.29、41.29、87.20，未用有量纲原始矩阵条件数作判断。

9. **修正大小。** 约束最小二乘使 C f=(F,M/L)，最小化 Σ_i |Δf_i|²/D_i；D_i=|f_i^0|²+(0.01 Σ|F_f| A_i/A)²，后项是固定正则尺度。压力/剪切/法向黏性/总项 L2 修正分别 {q['metrics']['pressure']['correction_L2_fraction']*100:.4f}% / {q['metrics']['shear']['correction_L2_fraction']*100:.4f}% / {q['metrics']['viscous_normal']['correction_L2_fraction']*100:.4f}% / {q['metrics']['total']['correction_L2_fraction']*100:.4f}%。剪切接近但低于预先冻结 10% 警示尺度；未放宽门槛。单点最大相对修正压力 1.81%、剪切 18.53%；尾部分布完整披露。

10. **局部合理性。** 压力/剪切反向力份额均为 0；压力向外作用份额为 0。最大单点总力份额 {q['metrics']['total']['single_node_force_share_max']*100:.4f}%。压力最大转向 0.414°，剪切最大转向 8.563°；剪切相对面积平均法向的绝对投影约 5.284%，来自跨曲面平均及小幅全局修正，不能称作每点严格切向。剪切前 1% 节点力份额由 4.498% 变为 5.131%；区域变化通过冻结门槛。原面接触边缘占总牵引模积分约 21.08%，初始离散边缘节点约 23.83%；离散支持跨越区域边界，二者不要求逐点相同。四个共同 x 区域的原面/节点向量另在 `INDEPENDENT_VERIFICATION.json`。仍有 50 µm 点距和薄边缘限制，未证明网格客观性或局部收敛。

11. **未来接口。** `CLOT_NODAL_TRACTION_FORCES.npz` 含全部 1,927 个真实 ID/坐标、分项初始/修正节点力和支持标志；`CFD_TO_LAMMPS_MAPPING_MATRIX.npz` 含 CSR 权重与六约束低秩修正数据；元数据记录 SHA、法向、表压与单位。原生 BPM 先累加力密度 N/m³，现有 `force_conversion()` 乘本点体积一次，再加外力 N。本轮力已乘面元面积，禁止再乘 V_i。乱序 ID、锚点拒绝、密度先转 N 后叠加的离线测试通过；生产输入未修改。输出适合下一轮原生完整态接口验证，并非已获认证的动态耦合。

12. **范围。** 本轮新 CFD=0、新 LAMMPS 动态断裂=0，仅一个接受加密场、一次边缘支持修正和一次全局最小守恒修正。现有脱离碎片已有背景流体阻力，未来须明确对象/阶段，避免重复计力。当前牵引仅对应固定完整几何，破裂后不能无限沿用。未含微泡 streaming、PZT、几何反馈或双向 FSI。

下一轮唯一优先建议：{nextstep}

图 B1–B4 与完整审核 JSON 见 `OPEN_RESULTS.html`。原生接口依据为锁定本地源码及 [LAMMPS addforce](https://docs.lammps.org/fix_addforce.html)、[external](https://docs.lammps.org/fix_external.html) 的接口文档；本轮没有调用新的生产动力学加载。
''')
 write(A/'FOSTER_FAILURE_CRITERION_AUDIT.md',r'''# Foster 2011 与当前 LPS 的能量失效审查

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
''')
 write(A/'FAILURE_PARAMETER_REQUIREMENTS.md',r'''# 独立失效参数与实验缺口

现有 G=940 Pa、K=33333.333 Pa 来自 Liu 2024 Supplement Table S2 的未损伤小应变映射：C1=0.47 kPa，G0=2C1；D1=0.06 kPa⁻¹，K0=2/D1。C2=0.17 kPa、C3=−0.012 kPa 和 Mullins r=1、β=0.26、m=0.87 J/m³ 未进入当前线弹性 LPS。ρ=1000 kg/m³、s00_bulk=0.15、alpha=0、dx=50 µm、δ=150.5 µm 保持原值；它们不是完整 Yeoh/Mullins/cohesive-zone 材料卡。

| 要独立控制的量 | 还需参数 | 需要的独立实验与验证 |
|---|---|---|
| 本体起裂强度 | 多轴强度/临界变形，缺陷与速率依赖 | 同组成、同成熟时间、同速率的无缺口拉伸/剪切；区分夹持热点与材料强度 |
| 裂纹传播阻力 | 本征 Gc、软化长度/位移或相场/非局部长尺度、R 曲线 | 已预缺口、稳定扩展、实际裂面积，几何/裂长/厚度与网格控制；起裂位置本身不能代替 Gc |
| 壁面脱粘 | 接口强度、I/II 模式粘附能、混合模态规则 | 指定内动脉壁、表面处理和血栓组成的剥离/剪切；不能用猪皮剥离参数冒充内动脉壁 |
| 体耗散 | 粘弹松弛谱、Mullins/塑性历史参数与不可逆能量 | 循环加载、保持/松弛和卸载；独立扣除机器、接触及数值阻尼 |
| 载荷速率 | 强度、Gc、耗散的速率函数及时间尺度 | 覆盖拟用脉冲时间尺度的多速率实验；当前参数不能由静态文献外推至超声 |

Liu 2021 全文的人体全血凝块 Gc=5.90±1.18 J/m²，剪切模量302.92±80.42 Pa；正文长度尺度约0.2/7 mm，图5与用户给出的近似值约0.17/6.5 mm相应。l=Gc/E，按 E≈3G 得约6.49 mm。它们描述该研究的人体样本；当前 Liu 2024 是牛血，不能拼成同一张认证卡。

Taylor 2018 Section5指出：裂纹、剩余韧带或试样尺寸过小可使试验落在强度控制区，即便裂纹从缺口尖端长出也不证明测得材料韧性。存在 R-curve 时必须分清起裂与传播阻力、裂纹增长量和约束。有限试样的 PVCE G_drive、无自然扩展的阈值、低载 G 均不能直接作为材料 Gc。

来源：[Liu 2021](https://doi.org/10.1016/j.eml.2021.101444)、[Taylor 2018](https://doi.org/10.1016/j.jmbbm.2017.07.007)、[Liu 2024](https://doi.org/10.1016/j.jmps.2024.105858)，本目录 sources 保存用户全文与 SHA。禁止将所列缺口通过同一条缺口曲线同时拟合到看似正确；强度和耗散应分别识别。
''')
 write(A/'REVIEW_LPS_SURFACE_CORRECTION_AND_ENERGY_FAILURE_FEASIBILITY_V1.md',f'''# LPS 表面权重与能量失效可行性

Run `{RUN}`。**M5C 的权重预处理及一个明确标注的变分诊断通过；V3 权重有负值和极大值，当前生产 LPS 的力学表面修正与能量失效移植均未建立。** 没有修改原生模型、s00、材料卡、断键规则或旧结果。

13. **哪里有邻域偏差。** 真实 M5C 2,880 点、132,610 条参考键，2,570 条初始缺口键始终不在计算图中；V3 1,927 材料点和实际锚层/界面拓扑分别分析。完整相同立方格点的 δ 邻域有 {nd['M5C']['reference_full_stencil_neighbors']} 个邻居。表面、夹具角点与缺口附近的邻居数、m/full-m、方向二阶矩与奇方向偏置均偏离内部，数值见 `SURFACE_NEIGHBORHOOD_AUDIT.json`。二阶矩用 (3/m)Σ(h/r)V ξ⊗ξ；物理自由表面确实缺邻居，不等同必须修掉的数值误差。任意仿射形变未必满足无牵引边界，不能用所有边界节点力为零作为验收。

14. **Shojaei 适用性。** 已读全文 36 页中相关方法/算例与锁定 perifit commit `66b2192...`。论文 Eq.50/64–69 是局部多项式 collocation、导数与能量算子的正则化最小二乘，通过全局稀疏 K_w w=F_w 耦合；对称键因子 (w_i+w_j)/2，保持参考拓扑。论文 Eq.46–47 明确采用完整球解析 m，不随权重改变；常数核 ω=1。当前 native 则是 ω=1/r、带 h taper、初始截断离散 m。不能把权重文件直接交给 LAMMPS 后宣称论文复现。源码未施加非负约束；本次冻结源码实际 OSB 为54算子行，而 README 广告为109，按实读源码计算并披露。使用无量纲 x/δ、V/δ³ 避免源码绝对 1e−20 截断在 SI 小量上吞掉高阶行。

原文静态3D OSB示例支持在其定义下研究；动态裂纹 Section5.4 是 PMB Kalthoff–Winkler，不能替当前 LPS 失效机制背书。原文权重一次预处理，不意味着动态断键后可以任意重算而能量不变。若未来更新权重或 m，会额外改变势能，必须计入能量账本；本轮没有动态更新、恢复裂面键或删除坏节点。参考图 A1、A2 不能被当作物理表面“消失”。

15. **实际权重。** M5C 范围 {a['M5C']['weights']['min']:.6f}–{a['M5C']['weights']['max']:.6f}，均正，全局残差 {a['M5C']['global_relative_residual']:.3e}。V3 范围 {a['V3']['weights']['min']:.6f}–{a['V3']['weights']['max']:.6f}，15 个负节点、59 条非正平均键，局部正则矩阵最大条件数 {a['V3']['local_regularized_gram_condition']['max']:.3e}。V3 只对实际材料 bulk 图拟合，未把真实锚层当作自由材料扩充，故不能等同粘附边界修正。没有截断负值、清理点云或再次调参。实际算子残差并非普遍改善：M5C 中位 {a['M5C']['operator_residual_before']['median']:.3f}→{a['M5C']['operator_residual_after']['median']:.3f}，V3 {a['V3']['operator_residual_before']['median']:.3f}→{a['V3']['operator_residual_after']['median']:.3f}。全局权重方程收敛不等于每个匹配算子精确。

16. **相容势能。** 仅对正权重的 M5C 定义 kernel-preserving 诊断：保持原 h/r 与初始 m；把对称 w 同时放入 θ 的每项和 e² 储能每项。使用锁定原 LPS 的展开势能，并重新求其精确导数。它是“保留核与参考归一化的变分适配”，不是 Shojaei 逐式复现或已认证表面修正。改写成标准偏伸长平方时还需 (5G/6)(1−m_weighted/m)θ² 项；只给原力乘权重、或者改 θ 却不改储能都不一致。正权重、K−5G/3>0 时该展开能量非负；这不构成任意点云无零模/动力稳定的证明。

17. **力—能量数值检查。** 只用一个实际 M5C 试样的确定性小形变场及刚体变换。原 evaluator 与锁定独立 C++ 储能相对差 {me['original_cpp_reference']['relative_energy_error']:.3e}，力差 {me['original_cpp_reference']['relative_force_error']:.3e}。加权诊断方向/分量差分误差为 {me['weighted_diagnostic']['directional_FD_relative_error']:.3e}/{me['weighted_diagnostic']['component_FD_relative_error']:.3e}，h=1 nm；刚体参考力/能量为0，旋转后的能量/力客观性、净力/当前位置力矩为舍入量级。U 从 {me['original']['U_J']:.9e} 增到 {me['weighted_diagnostic']['U_J']:.9e} J，最大节点力从 {me['original']['max_node_force_N']:.9e} 增到 {me['weighted_diagnostic']['max_node_force_N']:.9e} N。没有连续体边值解作为物理真值，不能称这种变硬是改善。V3 非正权重不进入力学应用。

18. **Foster 不能直接移植什么。** Section2是 correspondence/率相关材料示例；Eq.13是双向状态力差的路径积分，单位J/m⁶；Eq.16–17由理想平面裂纹球冠积分导出，Eq.18永久去掉影响。当前 LPS 多键 θ² 交叉项、原生冻结 m、有限自由/初裂面、破裂瞬间能量转移及材料耗散均需单独处理。`FOSTER_FAILURE_CRITERION_AUDIT.md` 逐项说明，不只做 4Gc/(πδ⁴) 数值代入。

19. **逐相互作用指标是否唯一。** 可唯一计算“指定完整形变路径和拓扑历史下的状态力路径功”，不能从本次证据定义路径无关、独立的逐键弹性势能。四角闭合弹性形变路径的总功与 U 最大归一化误差 {pw['max_total_work_energy_error_scaled']:.3e}，闭合总功 {pw['closed_total_work_J']:.3e} J；单键闭合功最大 {pw['closed_max_abs_pair_work_J']:.3e} J，各键绝对闭合功之和 {pw['closed_sum_abs_pair_work_J']:.3e} J，正负在总量抵消。证据是实际计算，不是示意假曲线。原 M1/M2记录另复用，外功差由动能与数值阻尼解释；未把阻尼当生物耗散。

20. **能否加新断键/软化。** 否。数学上一种加权势能可导、总功闭合，不证明 Foster w_c 对当前 LPS 是可校准传播阻力。原 s00=0.15 仍只保留为已有定性工程基线；本轮不把其瞬断行为升级为正确材料韧性，也没有足够证据生产启用独立能量耗散机制。`ENERGY_FAILURE_PORTABILITY=NOT_ESTABLISHED`、`MATERIAL_GC_UNAVAILABLE=true`。

21. **材料缺口。** 同组成/成熟时间/速率的无缺口强度、稳定扩展本体 Gc/R-curve、指定血管壁界面粘附、独立循环/松弛体耗散，以及速率关系。Liu2021人体全血与Liu2024牛血不可混卡；Taylor2018要求排除小试样/短缺口导致的强度控制。当前只采用Liu2024 TableS2未损伤小应变 G/K，未标定裂纹能、Mullins 或界面。见 `FAILURE_PARAMETER_REQUIREMENTS.md`。

来源：[Shojaei2026](https://doi.org/10.1002/nme.70385)、[perifit固定源码](https://github.com/Hereon-InstituteMS/perifit/tree/66b2192a9f95db29d16adefa162acb579759b925)、[Foster2011](https://doi.org/10.1615/IntJMultCompEng.2011002407)、[Liu2021](https://doi.org/10.1016/j.eml.2021.101444)、[Taylor2018](https://doi.org/10.1016/j.jmbbm.2017.07.007)、[Liu2024](https://doi.org/10.1016/j.jmps.2024.105858)。用户全文与补充材料、本次读取源码、SHA见 sources 和共同清单。
''')
 # Clarify evidence-based partial statuses and the local rim amendment.
 ast=load(A/'FINAL_STATUS.json');ast['states']['SURFACE_CORRECTION_WEIGHT_VALIDITY']='M5C_POSITIVE_PREPROCESSING_ONLY; V3_FAIL_NEGATIVE_AND_LARGE_WEIGHTS';save(A/'FINAL_STATUS.json',ast)
 meta=load(B/'CLOT_NODAL_TRACTION_METADATA.json');meta.update(contact_rim_support_repair='One10-face actual-surface graph repair; see frozen amendment; all other gates unchanged',repair_freeze_sha256=sha(B/'CONTACT_RIM_SUPPORT_REPAIR_FREEZE.json'),formal_field_count=1,support_screen_attempts=2,actual_force_mapping_count=1,global_constraint_corrections=1,independent_verification_sha256=sha(B/'INDEPENDENT_VERIFICATION.json'));save(B/'CLOT_NODAL_TRACTION_METADATA.json',meta)
 bcost=load(B/'RESOURCES.json');acost=load(A/'RESOURCES.json');fcost=load(COMMON/'FIGURE_RESOURCES.json')
 write(COMMON/'COMBINED_REVIEW.md',f'''# 含栓载荷映射与 LPS 失效可行性：共同审查

Run `{RUN}`。**含栓 CFD 已能产生通过本轮离线工程筛查的 LAMMPS 节点载荷；仍未接入生产动态。表面权重/能量失效研究没有充分依据升级当前生产 LPS。** 两线独立验收，无总科学 PASS。B问题1–12和A问题13–21见各自完整报告。

22. **新增代码与计算。** 代码仅在 `analysis/traction_failure_v1/`，数据仅在本 run 三个 verification 和两个 results 目录。B 一个原接受加密场，初次支持检查拒绝十个面，唯一一次真实表面邻接修正后实际生成一次载荷，再做六约束最小修正；独立重载/解析积分/接口检查。A 两个已存拓扑的权重预处理，一个 M5C 受控静态力—能量诊断、少量解析变形路径积分和原 M1/M2只读账本。原生新CFD、LAMMPS动力学/新平衡、高载、失效参数扫描均为0。没有安装或升级服务器生产依赖；临时断连期间仅创建WSL独立分析环境。

23. **资源与成本。** 维护结束后恢复原服务器，HEAD `cc5ee77149f55e58d715bb5a4660c6467231a826`，CPU cgroup v1 配额1536000/100000=15.36核，亲和0–31，RTX4090 24GB；Python3.13.11等版本快照保存。正式离线进程 BLAS/OMP/MKL单线程、GPU计算0。B实际映射 {bcost['wall_s']:.3f}s，A {acost['wall_s']:.3f}s，七组绘图 {fcost['wall_s']:.3f}s；CPU/RSS见 `RESOURCE_AND_COST_SUMMARY.json`。这些相加是计算任务成本，文献阅读/开发/传输与用户维护等待不伪算成求解成本。实际受控作业由本轮专属Supervisor托管，禁止自动重启；未触碰其它任务。

24. **审查图。** B1真实面/材料点/锚点及支持距离；B2原生表压/剪切/面积/法向与力；B3初始/修正分布、守恒、修正与集中；B4真实全1927点的压力/剪切/总载荷箭头。A1真实M5C/V3邻域；A2真实权重、残差/条件数及同场力—能量适配；A3全局功/单键闭合路径与原M1/M2。共7组PNG/SVG，各线与共同静态索引。图上向量缩放明确；没有新医学外壳或断裂动画。

25. **保护。** 最终 `PROTECTION_AUDIT.json` 重核任务前注册的1611个源/材料/原生库/旧结果 SHA；正式完成以该最终文件为准。Git HEAD和环境前后比较另存。未改生产force provider、s00、原生break、材料卡或历史轨迹；未commit/push/清理。初次支持拒绝与后续修正原样留证，不反改旧结论。

26. **实际WSL文档。** 由独立串行写入脚本更新 `/home/lzy/projects/clot_rupture.md`，执行结果见 `document_edit/AUDIT.json`；没有把服务器同名路径当WSL文件。保留①–⑨、旧无栓轨迹、先前高载失败与Gc未识别结论。

27. **证据位置。** `document_edit/before.md`、`after.md`、`change.diff`、`AUDIT.json`记录最新SHA/并发保护/读回。`SOURCE_MANIFEST.json`、`CODE_ADDITIONS.diff`、`SERVER_ARTIFACT_MANIFEST.json`、`COMPACT_MANIFEST.json`、`LOCAL_RETRIEVAL_AUDIT.json`分别记录源码/原生库、差异、服务器产物和WSL回收。实际命令见 `REPRODUCE.txt`。完整旧CFD体网格/原生全场留服务器，不重新下载或舍入。

28. **唯一下一优先项。** {nextstep} 本轮停止，不自动执行。

要点：580个材料点受载、1347个材料点零外载、1561锚点排除；压力/剪切/总项修正0.618%/9.015%/0.772%，代数误差~1e−16；剪切局部支持与曲面平均尚有离散限制。M5C变分适配可导不代表表面修正已验证；V3有15负节点/59非正键。总弹性能可恢复不意味着单键路径功有唯一弹簧势能。当前独立Gc仍不可用。
''')
 save(COMMON/'RESOURCE_AND_COST_SUMMARY.json',dict(B=bcost,A=acost,figures=fcost,geometry_preflight=load(B/'GEOMETRY_PREFLIGHT.json')['runtime'],job_count=dict(CFD=0,native_LAMMPS_dynamic=0,native_new_equilibrium=0,accepted_CFD_fields_mapped=1,support_screen_attempts=2,single_local_support_repair=1,global_minimum_conservative_correction=1,weight_geometries=2,controlled_small_specimen_force_energy=1,figure_groups=7),GPU_compute=False,CPU_threads_per_analysis=1,counts_scope='Offline evaluations only; path quadrature nodes are evaluations, not equilibrium solves. Reading, development and maintenance wait not included in compute totals.'))
 def index(folder,title,figdir,report):
  figs=sorted(figdir.glob('*.png'));links=''.join(f'<a href="{os.path.relpath(folder/f,folder)}">{html.escape(f)}</a><br>' for f in report);cards=''
  for p in figs:
   src=os.path.relpath(p,folder);cards+=f'<article><h2>{html.escape(p.stem)}</h2><a href="{src}"><img src="{src}" loading="lazy"></a><p><a href="{os.path.relpath(p.with_suffix(".svg"),folder)}">SVG</a></p></article>'
  write(folder/'OPEN_RESULTS.html',f'<!doctype html><meta charset="utf-8"><title>{title}</title><style>body{{font:17px system-ui;max-width:1400px;margin:35px auto;padding:0 20px;color:#17324b}}a{{color:#14699b}}img{{max-width:100%}}article{{margin:40px 0;border-top:1px solid #ccd5dd}}p{{line-height:1.6}}</style><h1>{title}</h1><p>Run {RUN} · 真实数组绘图 · 离线工程与研究结论分项验收 · 未运行新断裂/CFD</p>'+links+cards)
 index(B,'B · 含栓 CFD 到真实材料点的守恒载荷',FB,['REVIEW_Y_CLOT_CFD_TO_LAMMPS_CONSERVATIVE_TRACTION_MAP_V1.md','FINAL_STATUS.json','LOCAL_MAPPING_QUALITY.json','CONSERVATIVE_MAPPING_AUDIT.json','INDEPENDENT_VERIFICATION.json','CLOT_NODAL_TRACTION_FORCES.npz','CLOT_NODAL_TRACTION_METADATA.json','CFD_TO_LAMMPS_MAPPING_MATRIX.npz'])
 index(A,'A · LPS 表面权重与能量失效可行性',FA,['REVIEW_LPS_SURFACE_CORRECTION_AND_ENERGY_FAILURE_FEASIBILITY_V1.md','FINAL_STATUS.json','FOSTER_FAILURE_CRITERION_AUDIT.md','FAILURE_PARAMETER_REQUIREMENTS.md','SURFACE_CORRECTION_FEASIBILITY.json','FORCE_ENERGY_COMPATIBILITY.json','LPS_INTERACTION_WORK_AUDIT.json'])
 write(COMMON/'OPEN_RESULTS.html',f'<!doctype html><meta charset="utf-8"><title>Y clot traction / failure audit</title><style>body{{font:20px system-ui;max-width:1050px;margin:60px auto;line-height:1.7}}a{{color:#176699}}</style><h1>含栓载荷映射与 LPS 失效可行性</h1><p>Run {RUN}</p><p><b>B：通过离线守恒与局部工程筛查。</b> 580点受载，锚点排除；生产耦合未运行。</p><p><b>A：M5C变分诊断通过，V3权重不通过；能量失效移植未建立。</b> 生产模型和材料保持不变。</p><p><a href="{os.path.relpath(B/"OPEN_RESULTS.html",COMMON)}">B线4组图、报告与节点力</a></p><p><a href="{os.path.relpath(A/"OPEN_RESULTS.html",COMMON)}">A线3组图、报告与数学审查</a></p><p><a href="COMBINED_REVIEW.md">共同中文评审与下一步</a></p><p>完整原始CFD留服务器；未运行新断裂、streaming或双向FSI。</p>')
 save(COMMON/'FINAL_STATUS.json',dict(run_id=RUN,B=load(B/'FINAL_STATUS.json')['states'],A=ast['states'],overall_scientific_PASS_not_assigned=True,stopped_after_delivery=True,next_unique_priority=nextstep))
 print('reports and indexes generated')
if __name__=='__main__':main()
