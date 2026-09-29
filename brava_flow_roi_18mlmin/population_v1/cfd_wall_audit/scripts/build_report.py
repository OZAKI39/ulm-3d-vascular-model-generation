"""Build a standalone Chinese audit report from machine-readable executed evidence."""
from pathlib import Path
import argparse,json,html
import numpy as np
import pandas as pd
ap=argparse.ArgumentParser();ap.add_argument('--population-root',type=Path,required=True);root=ap.parse_args().population_root;out=root/'cfd_wall_audit'
def read(name):return json.loads((out/name).read_text())
def table(df):
 def value(v):
  if isinstance(v,float):return f'{v:.6g}'
  return str(v).replace('|','/').replace('\n',' ')
 return '| '+' | '.join(map(str,df.columns))+' |\n| '+' | '.join(['---']*len(df.columns))+' |\n'+'\n'.join('| '+' | '.join(value(v) for v in row)+' |' for row in df.itertuples(index=False,name=None))+'\n'
i=read('AUDIT_IDENTITY.json');s=read('SUMMARY.json');sp=read('SUPPLEMENT_SUMMARY.json');part=read('PARTICLE_OVERLAP_SUMMARY.json');val=read('validation/independent_recomputation.json');a=read('known_regression_388769.json');b=read('known_regression_138584.json')
tiers=pd.read_csv(out/'layer_tiers.csv');dist=pd.read_csv(out/'layer_distributions.csv');cor=pd.read_csv(out/'correlations.csv');clusters=pd.read_csv(out/'anomaly_clusters.csv');ex=pd.read_csv(out/'particle_layer_exposure.csv');assoc=pd.read_csv(out/'outcome_associations.csv');bc=pd.read_csv(out/'boundary_identity_flux.csv');ws=pd.read_csv(out/'wall_structure_distributions.csv');near=tiers[(tiers.region.isin(['W0','W1','W2']))&(tiers.metric=='eta_div')&(tiers.upper==np.inf)].copy();near['数量占比_%']=near.count_fraction*100;near['体积占比_%']=near.volume_fraction*100
ncluster=len(clusters[clusters.family=='eta_ge_0p5']);ncells=int(clusters.loc[clusters.family=='eta_ge_0p5','cell_count'].sum())
body=[r'''# BraVa 近壁 CFD/FEM 全域审计报告

本轮范围：只读审计当前冻结流场，关联已有历史轨迹；不重算 CFD，不改速度、压力、网格、壁面、WSS 或微泡物理。下述“高值”“异常筛查”都不是自动的 CFD 不合格判据。

## 1. 通俗结论
''',f'''当前证据支持“强异常局部集中”，不支持“整个 wall layer 普遍失效”。W0/W1/W2 分别有 {s['layer_counts']['W0']:,}/{s['layer_counts']['W1']:,}/{s['layer_counts']['W2']:,} 个单元；η≥0.5 仅有 18/14/0 个，η≥0.1 有 660/3893/4757 个。全部正式壁面节点速度严格为零。388769 的强朝壁面速度被独立复现，且周围形成小簇；138584 是唯一的四节点全零单元，四角都在壁面上，体现 P1 近壁离散的局部局限。微泡访问会显著放大少数单元的影响，但并非所有加速都经过这些单元。没有发现节点错配、WSS 公式或单位错误的证据；这不能证明原 CFD 的物理精度。\n''',r'''
## 2. 输入身份、单位和实际流程

`labelled_core_m.vtp → sv.meshing.TetGen → canonical mesh → svMultiPhysics fluid → frozen steady_flow.vtu → 原始节点 P1 插值 → population_v1`。

网格生成脚本 `brava_flow_roi_18mlmin/scripts/generate_mesh.py:7–17` 调用 SimVascular 2023-05-31 的 TetGen，目标边长 0.0002 m，**未设置边界层网格**。材料及边界来源 `cases/balanced_pressure_final/run/solver.xml:44–110`：恒定牛顿动力黏度 0.00345312 Pa·s，密度 1056 kg/m³；WALL 为无滑移 Dirichlet 0；入口为 flat、周缘置零且归一化到 −3×10⁻⁷ m³/s（外法向符号），即 18,000 μL/min。出口自然牵引压力参数依次为 O1=556.97065554884477 Pa、O2=0 Pa、O3=17.137427153562019 Pa。自然牵引参数不等于“所有出口节点压力强制等于该值”。

CFD 时间步为 0.01 s，非线性容差 10⁻¹⁰，GMRES/PETSc-Jacobi；这是已有结果的配置记录，本轮没有再次运行求解器或宣称重验求解器收敛。微泡生产使用的名义时间步为 0.0005 s，二者不能混淆。四节点速度/压力都是连续 P1；稳定化不可压缩有限元并不要求每个四面体的 trace(G) 严格为零。本审计不把非零散度自动判为错误。

本轮独立求解每个四面体的 3×3 梯度系统；原始冻结场没有经验流速替代或修改。已有 WSS 入口为 `vendor/flow_solver_support/wss_case.py:51–71` 的 `recover`。原粒子插值为 `particle_3d/field.py:49–57`；独立审计不读取旧梯度缓存。
''',f"\n网格节点 **{i['node_count']:,}**；四面体 **{i['tetra_count']:,}**；WALL 三角形 **{i['wall_triangle_count']:,}**。审计 NumPy {i['numpy']}、PyVista {i['pyvista']}、VTK {i['vtk']}。Git 工作树基线 `{i['source_code_commit']}`；历史 production 的实际版本与文件哈希保留在旧 `IDENTITY.json`。本轮源码逐文件 SHA 见交付清单；不把工作树基线当作所有旧生产文件的版本。\n",
 table(pd.DataFrame([dict(input=k,SHA256=v) for k,v in i['sha256'].items()])),
 r'''
单位全链路：坐标 m、速度 m/s、时间 s、压力 Pa、动力黏度 Pa·s、梯度及应变率 s⁻¹、通量 m³/s、WSS Pa。图中几何坐标换成 mm，只用于显示。

### 对应关系实际验证

- canonical GlobalNodeID=1…135852、GlobalElementID=1…692313；审计 tetra_id 与历史轨迹 tetra_id 都是 **0-based canonical row**。
- **冻结流场没有 GlobalNodeID 字段**，不能声称已经读到了该字段。实际核查了逐节点坐标精确相等、逐单元四个节点的无序集合精确相等；全部坐标唯一。局部顶点置换不影响 canonical 复算。
- WALL 及所有端口的 GlobalNodeID 坐标与体网格精确一致，GlobalElementID 确实指向拥有该三角面的唯一四面体。
- 外表面 101456 个面完全由官方 WALL/INLET/OUTLET_01/02/03 无重叠覆盖。端口身份来自官方文件及编号，不按空间位置或标签大小猜测。
- 全部单元体积正、canonical signed det 正，没有发现翻转/零体积或共面邻接超过两个体单元。此项不是近壁分辨率足够的证明。

端口位置与本轮独立通量积分：
''',table(bc),f"\n独立边界净流量相对入口的绝对偏差为 {sp['boundary_relative_imbalance']:.6g}。整个冻结场节点压力范围 {sp['pressure_range_Pa']} Pa，速度模长范围 {sp['velocity_node_speed_range_m_s']} m/s。整体守恒不能替代局部精度评估。\n",
 r'''
## 3. 区域、公式及可复算数据

W0 是拥有正式 WALL 面的四面体；W1 是 W0 的共面邻居，排除 W0；W2 是 W1 的下一圈共面邻居，排除 W0/W1。三组互不重叠。另对**全部 692313 单元**计算重心到有限 WALL 三角形的最近欧氏距离，按 d/h_max 分成 [0,0.5)、[0.5,1)、[1,2)、[2,4)、[4,∞)。这些是审计集合，不是新增物理边界。

设 X 的行是 x₁−x₀、x₂−x₀、x₃−x₀，U 的行是相应节点速度差，独立求解：

\[
G=(X^{-1}U)^T,\quad \nabla p=X^{-1}\Delta p,\quad
\eta=\frac{|\operatorname{tr}G|}{\|G\|_F}.
\]

仅在分母恰好为零时将比值记录为 0，同时另存 `all_four_u_zero`；没有用大 epsilon 掩盖低速区域。应变率张量 D=(G+Gᵀ)/2，curl 按完整 G 的反对称部分计算。

每个面使用三个节点速度的算术均值，乘以向外的面积向量：

\[
Q_f=\bar{u}_f\cdot(n_f A_f),\qquad \sum_f Q_f=(\operatorname{tr}G)V.
\]

这对 P1 是精确面积分。两边相符只是散度定理的**离散实现恒等式**，绝不是“局部物理质量守恒通过”。

保存 h_max、h_min=3V/max(A_f)、等体积球直径、六条边长、最小内部二面角。形状指标 mean-ratio=12(3V)^(2/3)/Σl²，正四面体=1；scaled Jacobian 为各角 |det|/三边乘积的最小值，再乘 √2 归一化。另存有符号 det。未套用统一好坏阈值。

W0 每个 WALL 面沿面重心至第四节点的 10%、25%、50%、75% 取样，并单列四面体重心（恰好是该线段的 25%）。使用**朝管腔的法向**；u_n<0 表示朝壁面。多 WALL 面的单元逐面保存；层统计选其面积最大的 WALL 面，避免把一个单元重复计数。W1/W2 重心法向来自最近的有限壁面三角形，是局部诊断约定，边/角处可能不唯一。

核心计算实码（完整版本在 scripts/audit_field.py）：

```python
G = np.linalg.solve(X[:,1:]-X[:,:1], U[:,1:]-U[:,:1]).transpose(0,2,1)
div = np.trace(G, axis1=1, axis2=2)
eta = np.divide(abs(div), norm_G, out=np.zeros_like(div), where=norm_G!=0)
traction = mu * np.einsum('nij,nj->ni', G_owner + G_owner.transpose(0,2,1), n_out)
tangent = traction - np.einsum('ij,ij->i', traction, n_out)[:,None]*n_out
WSS = np.linalg.norm(tangent, axis=1)
```

## 4. Q1/Q2：单元比例和体积比例

以下阈值仅用于筛查优先级，**不是物理合格/不合格标准**。
''',table(near[['region','lower','total_count','flag_count','数量占比_%','体积占比_%','total_volume_m3']]),
 '\n距离分组（分母 h_max）：\n',table(tiers[(~tiers.region.isin(['W0','W1','W2']))&(tiers.metric=='eta_div')&(tiers.upper==np.inf)][['region','lower','total_count','flag_count','count_fraction','volume_fraction']]),
 '\n主要分布（min/median/P90/P95/P99/max；SI 单位在 metric 列）：\n',table(dist[(dist.region.isin(['W0','W1','W2'])) & (dist.metric.isin(['abs_div_s_inv','eta_div','G_norm_s_inv','speed_m_s','abs_u_normal_m_s','normal_fraction','h_max_m','h_min_altitude_m','mean_ratio_quality']))]),
 r'''
## 5. Q3/Q4：位置和空间聚集性
''',f"\nW0–W2 的 η≥0.5 共 {ncells} 单元，形成 {ncluster} 个共面连通簇。最大簇有 {int(clusters[clusters.family=='eta_ge_0p5'].iloc[0].cell_count)} 个单元，包含已知 388769 周边；强筛查单元对应的 WALL 面积占全部 WALL 的 {s['W0_eta_ge_0p5_wall_area_fraction']*100:.6f}%。中等阈值 η≥0.1 共 9310 单元，形成 {sp['moderate_eta_cluster_count']} 个簇，最大 135 单元。这说明强现象是局部小簇，中等偏差是分散于多个位置的区域性现象，不能称整个壁层普遍严重失效。\n",
 table(clusters[clusters.family=='eta_ge_0p5'].head(12)),
 '\n全部簇、三维包围盒、单元体积、壁面覆盖面积、关联单元/三角面 ID 均另存。`spatial_octants.csv` 按实际包围盒八分区统计空间分布；不是未经核实的解剖分区。\n',
 r'''
## 6. Q5/Q6：尺寸和网格质量是否解释异常
''',table(cor[cor.region.isin(['W0','W1','W2'])]),
 f"\nW0 中 η 与 h_max 的 Spearman 系数 {float(cor[(cor.region=='W0')&(cor.y=='h_max_m')].spearman.iloc[0]):.6f}，与 mean-ratio 为 {float(cor[(cor.region=='W0')&(cor.y=='mean_ratio_quality')].spearman.iloc[0]):.6f}，没有强单调关联。W1/W2 有弱相关，也不足以给出因果结论。388769 的 mean-ratio={a['metrics']['mean_ratio_quality']:.6f}、最小二面角={a['metrics']['min_dihedral_deg']:.3f}°、h_max={a['metrics']['h_max_m']*1e6:.3f} μm，并非明显扁塌坏单元；但约 0.2 mm 的无边界层 P1 网格仍不能自动保证微米级微泡附近的梯度可靠。壁面相邻面最大法向转角与 η 的整体 Spearman 为 {sp['eta_vs_local_normal_angle_spearman']:.6f}，不支持简单地以角度越大越异常作诊断。\n",
 table(ws),
 r'''
## 7. Q7/Q8/Q9：朝壁面流动、无滑移和 P1 特殊结构
''',f"\n全部 {s['wall_nodes']['count']} WALL 节点速度 max/median/P95/P99 均为 0 m/s，非零个数为 0；空的 `validation/nonzero_wall_nodes.json` 是实际检查结果。所有 {s['wall_face_zero3_count']} WALL 面的三个节点严格零速；其中 {s['wall_face_zero3_opposite_gt_0p1_m_s_count']} 个面的第四节点速度 >0.1 m/s。0.1 m/s 只是描述性筛选阈值，不是‘过快’判据。\n",
 r'''
对三壁面节点全零的 P1 单元，令 n 为 inward、H 为第四节点到壁面平面的高度，则：

\[
u(x)=\frac{d_{plane}(x)}{H}u_{opp},\quad G=\frac{u_{opp}\otimes n}{H},\quad
\operatorname{div}u=\frac{u_{opp}\cdot n}{H},\quad \eta=\frac{|u_{opp}\cdot n|}{|u_{opp}|}.
\]

因此 W0 中 η 和 normal_fraction 几乎完全相等是**同一 P1 结构的代数关系**，不能当作两个独立证据。第四节点高速可以是正常切向剪切；真正值得关注的是它的方向、局部面通量和邻域。
''',f"\nW0 中法向占比≥0.5 的单元共 {s['W0_normal_ge_0p5_count']} 个，全部 {s['W0_toward_normal_ge_0p5_count']} 个朝向壁面。并非大量壁面单元均主要朝壁面。η 与 normal_fraction 的最大差 {sp['W0_identity_eta_equals_normal_fraction_max_abs']:.3g}。原始每面/每取样点速度、符号、高度与法向均已保存，W1/W2 不套用上述秩一恒等式。\n",
 r'''
## 8. Q10：388769 独立复算
''',table(pd.DataFrame([{k:a['metrics'][k] for k in ['tetra_id','layer','div_s_inv','G_norm_s_inv','eta_div','u_normal_m_s','normal_fraction','volume_m3','Q_net_m3_s','d_wall_m']}])) ,
 f"\nGlobalNodeID={a['GlobalNodeIDs']}；节点速度模长={np.linalg.norm(a['velocity_m_s'],axis=1).tolist()} m/s；节点压力={a['pressure_Pa']} Pa。散度 −827.429163466 s⁻¹、η=0.994887343 符合指定回归值。重心位置（mm）={np.array([a['metrics']['centroid_x_m'],a['metrics']['centroid_y_m'],a['metrics']['centroid_z_m']])*1000}。它是极端筛查值，不代表普通 W0 单元；但并非孤立单体尖峰，邻居也有高值：\n",table(pd.DataFrame(a['neighbors'])[['tetra_id','layer','eta_div','div_s_inv','speed_m_s','pressure_Pa','d_wall_m']]),
 r'''
## 9. Q11：138584 全零单元
''',f"\n全 692313 单元只有 **{s['zero_cells']}** 个四节点全零单元，零速度连通分量大小为 **{b['zero_component_cell_count']}**。其四个 GlobalNodeID={b['GlobalNodeIDs']}，全部属于正式 WALL；拥有两个 WALL 面、另两个面连到有流速的体单元。不是拓扑孤岛、连续大零速斑块或本次发现的节点映射错误。\n\n节点压力={b['pressure_Pa']} Pa，|∇p|={b['metrics']['grad_p_norm_Pa_m']:.6g} Pa/m，体积={b['metrics']['volume_m3']:.6g} m³、重心到 WALL={b['metrics']['d_wall_m']*1e6:.6f} μm。连续 P1 无滑移把四角速度全部约束为零，单元内部于是处处零；这里 η=0 不代表可信的低风险。该结构是独立的近壁离散风险，原样保留并单独标记，未填充邻居速度。\n",table(pd.DataFrame(b['neighbors'])[['tetra_id','layer','eta_div','speed_m_s','pressure_Pa','d_wall_m']]),
 r'''
## 10. 压力、WSS 和数值恒等式
''',f"\n逐单元四面通量与 div×V 的最大绝对差为 {s['flux_identity']['max_absolute_error_m3_s']:.6g} m³/s；用 Σ|Q_face| 归一化后最大差 {s['flux_identity']['max_error_scaled_by_sum_abs_face_flux']:.6g}。\n\n由原始 G、μ、外法向计算切向黏性牵引后，与已有 `WSS_raw_Pa` 按**canonical 三节点集合+owner**对应，最大差 {s['WSS_comparison']['max_absolute_difference_Pa']:.6g} Pa，RMS {s['WSS_comparison']['rms_difference_Pa']:.6g} Pa。对比的是原始三角面值，不是显示用面积平均节点值。两种程序一致只排除了相应实现/映射分歧，不构成独立物理 WSS 验证。\n\nP1 压力在共面上共享三个节点，面上的压力连续；允许单元压力梯度跳变。邻单元 |∇p_a−∇p_b| 的 median/P95/P99/max={s['pressure_neighbor_gradient_jump_quantiles_Pa_m']} Pa/m。保存了最大的 2000 个邻接压力变化，但未定义没有物理依据的压力 FAIL 阈值。\n",
 r'''
## 11. Q12：真实微泡访问和加速/失败关联
''',f"\n完整读取历史 {part['track_count']} 条轨迹、{part['state_count']:,} 个去重后接受时刻状态，访问 {part['visited_cell_count']:,} 个体单元，其中 W0–W2 {part['near_wall_visited_cells']:,} 个。所有状态位置的重心坐标检查均落在所记录 canonical 单元内（1e−8 容差），无无效 ID。η>0.1 的状态 {part['state_eta_gt_0p1']:,}（{part['state_eta_gt_0p1_fraction']*100:.4f}%），η>0.5 的状态 {part['state_eta_gt_0p5']:,}（{part['state_eta_gt_0p5_fraction']*100:.4f}%）。注意这里按用户问题用严格 >；全域筛查表用 ≥。\n",
 table(ex),
 '\n直接回答Q12：W0被访问单元中高η比例0.03010%，高于全部W0的0.01781%，且18个高η W0均被访问。重复状态加权则升到2.1795%，驻留时间近似加权为1.6563%，说明少数区域会被反复采样；W1的访问单元比例没有同样升高，不能把该偏向推广到所有近壁层。\n',
 '\n“访问单元集合”与“重复状态计数/停留时间”必须分别解读。若一个单元因慢运动或自适应子步被反复记录，状态权重会很大。`residence_s_start_cell_assignment` 把每段已接受时间间隔归给其起点所在单元；路径可能跨多个单元，因此这是明确的端点归属近似，不是沿线精确分段驻留积分。\n',
 f"\n201 个既有加速候选中，±0.06908 s 时间窗内有 {part['acceleration_window_eta_gt_0p5']} 个经过 η>0.5 单元，加速判定时刻有 {part['acceleration_onset_eta_gt_0p5']} 个位于该类单元。这里的加速关联采用全体体单元的η，不仅限W0–W2；具体onset单元ID已保存，可与all_cells.csv按ID联接。其余案例不能沿用 388769 的诊断。两个 numerical failure（440、1505）都终止在 138584，该单元 η=0，必须通过全零场检查捕获。\n",
 table(assoc),f"\n观察到驻留时长与 η>0.5 区域时间占比的 Spearman={part['residence_vs_exposure_fraction_spearman']:.6f}，是负相关，不支持笼统的‘高 η 必然导致更长驻留’。出口分组、轨迹路线、删失以及旧实验的 outcome-dependent 停灌均构成混杂因素；这些表不能给出因果或校正后的出口概率。\n",
 r'''
## 12. 已执行验证、限制及交付索引
''',table(pd.DataFrame([
 ['原始身份与编号','已运行','全部输入 SHA；坐标、连接、官方面 owner 一致'],
 ['P1 梯度独立复算','已运行',f"692313 单元，与实际粒子采样器重建 G 最大差 {val['original_node_gradient_recompute_vs_actual_particle_sampler_max_abs_s_inv']:.3g} s^-1"],
 ['散度定理','已运行','精确 P1 面积分；仅算法恒等式'],
 ['有限 WALL 距离','已运行',f"全域 VTK CPU 与 RTX4090 CUDA 交叉核对，最大差 {val['vtk_locator_vs_cuda_finite_triangle_distance_max_error_m']:.3g} m"],
 ['质量指标基准','已运行','正四面体 mean-ratio 和 normalized scaled Jacobian 均为 1'],
 ['全部近壁样本/连通簇','已运行','每个 WALL 面的五类取样、W0/W1/W2 与距离带'],
 ['两项指定回归','已运行','388769 强法向；138584 单一全零单元'],
 ['已有 WSS 对应','已运行','原始面片场同离散输入的一致性'],
 ['历史轨迹关联','已运行','全部1565条；不重积分'],
 ['网格收敛/重新求解CFD','未运行','用户明确禁止；无新增物理精度证明'],
 ],columns=['检查','执行状态','证据/限制'])),
 f"\n最近三角形 ID 有 {val['nearest_triangle_id_tie_or_roundoff_differences']} 个 CPU/GPU 差异，距离仍一致；共享边/顶点允许等距多面。W1/W2 并列最近面的法向选择最大可改变0.82/0.87的法向占比，分别有9/23个单元跨过0.5筛查线。因此W1/W2的法向簇需要保留这种方向不确定性，不能像W0自有壁面法向那样解释。完整敏感性：`{json.dumps(val.get('nearest_triangle_normal_sensitivity',{}),ensure_ascii=False)}`。W0 使用自己的正式 WALL 面，不受最近三角形并列选择影响。\n",
 r'''
本轮结论边界：未发现输入映射、无滑移节点或现有 WSS 实现不一致的证据；确认存在少量强朝壁面离散单元以及一个全零速度单元。它们可影响部分微泡轨迹，不能自动解释所有慢后变快现象，也不能由 CPU/GPU parity 推导为真实血流行为。没有修改 CFD、微泡物理、MB–MB 模型或 SonoVue 分布。

主要附件：

- `wall_adjacent_cells.csv`：每个 W0/W1/W2 单元四节点 ID、坐标、速度、压力、六边长和全部指标。
- `all_cells.csv`：整个体网格的指标；`recomputed_arrays.npz`：完整 G、D、curl、grad p、连接与邻接原始精度数组。
- `wall_nodes.csv`、`wall_face_structure.csv`、`wall_sample_metrics.csv`：无滑移、第四节点及每面五类样本。
- `layer_distributions.csv`、`layer_tiers.csv`、`correlations.csv`、`spatial_octants.csv`。
- `anomaly_clusters.csv`、`cluster_cells.csv`、`cluster_wall_triangles.csv`；中等阈值有相应 `moderate_eta_*` 文件。
- `worst_cells.csv` 与 `flagged_cell_neighbors.csv`：严重筛查单元及其所有共面邻居；`pressure_neighbor_extremes.csv`。
- `known_regression_388769.json`、`known_regression_138584.json`。
- `particle_overlap.csv`、`particle_exposure.csv`、`particle_layer_exposure.csv`、`acceleration_overlap.csv`、`outcome_associations.csv`。
- `wall_layers_audit.vtu`、`wall_audit.vtp`：审计派生数据，可在 ParaView 独立查看，原始输入未被改动。
- `figures/01…12*.png/.pdf`：12 张 300 dpi 科学图，原始单元值/面值，无平滑。图号和说明见下。
- `validation/`：原始失败日志、修正后的完整执行日志、输入不变性、独立核对和轨迹哈希。早期诊断脚本的依赖、整数类型和日志序列化错误保留；以 `audit_field_completed.log`、`particle_overlap_final.log` 为成功执行记录，不能把失败尝试当通过。

复现命令见 `README.md`。全部新文件限定在本审计目录及相邻固定协议目录，旧 production 未被覆盖。
''']
figindex=read('figures/FIGURE_INDEX.json')
for k,item in enumerate(figindex,1):body.append(f"\n### 图 {k}\n\n{item['caption']}\n\n![图{k}](figures/{item['file']}.png)\n")
(out/'WALL_CFD_AUDIT_REPORT_ZH.md').write_text('\n'.join(body))
(out/'OPEN_RESULTS.html').write_text('<!doctype html><meta charset="utf-8"><title>近壁 CFD 审计</title><style>body{font:17px system-ui;max-width:1400px;margin:30px auto;background:#f6f8fa;color:#172b3d}img{max-width:100%;background:white}section{margin:35px 0}a{color:#12649b}</style><h1>近壁 CFD / FEM 全域审计</h1><p>只读检查；筛查阈值不代表物理合格判据。<a href="WALL_CFD_AUDIT_REPORT_ZH.md">中文主报告</a> · <a href="../fixed_window_protocol/README.md">固定时间实验协议</a></p>'+''.join(f'<section><h2>{k}. {html.escape(item["file"])}</h2><p>{html.escape(item["caption"])}</p><a href="figures/{item["file"]}.pdf">PDF</a><img src="figures/{item["file"]}.png"></section>' for k,item in enumerate(figindex,1)))
(out/'README.md').write_text('''# 只读 CFD 壁面审计\n\n入口：WALL_CFD_AUDIT_REPORT_ZH.md / OPEN_RESULTS.html。所有 CSV 使用 SI，列名包含单位；layer 0/1/2=W0/W1/W2，3=其它体单元；tetra_id 和 wall_triangle_id 均从0开始。全部原始坐标/节点值在 wall_adjacent_cells.csv 与 recomputed_arrays.npz 中保留。\n\n服务器复现（不调用 CFD 或微泡积分器）：\n\n```bash\ncd /workspace/brava_flow_roi_18mlmin_20260928/population_v1\nexport OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1\nexport PYTHONPATH="$PWD/source/particle_3d/src"\nenv/bin/python cfd_wall_audit/scripts/audit_field.py --environment data/environment.json --output cfd_wall_audit --commit 51257fb4f\nenv/bin/python cfd_wall_audit/scripts/validate_recomputation.py --population-root "$PWD"\nenv/bin/python cfd_wall_audit/scripts/supplement.py --population-root "$PWD"\nenv/bin/python cfd_wall_audit/scripts/particle_overlap.py --population-root "$PWD"\nenv/bin/python cfd_wall_audit/scripts/plot_audit.py --population-root "$PWD"\nenv/bin/python cfd_wall_audit/scripts/build_report.py --population-root "$PWD"\n```\n\n再次运行只覆盖本目录派生结果，建议复制到新的审计目录后执行。`inputs/visible_speedup_candidates.csv` 是上一轮加速审查的只读参考副本；不从运动图像臆造新事件。环境依赖 NumPy/SciPy/Pandas/PyVista/VTK/Matplotlib，GPU 距离核对用已安装 PyTorch CUDA 和项目 GPUWall。\n\n包括体积及面积等比重的分母见 AUDIT_CONFIG.json；粒子逐段驻留近似归给起始单元。不要把状态次数当驻留时间或把 WSS 代码一致性当物理精度验证。\n''')
print(out/'WALL_CFD_AUDIT_REPORT_ZH.md')
