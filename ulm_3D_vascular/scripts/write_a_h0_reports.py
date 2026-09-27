"""Write Chinese review documents from actual saved numerical evidence."""
from pathlib import Path
import json,csv,xml.etree.ElementTree as ET
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
REPORT=ROOT/'reports/a_network_1d0d_boundary_v2_idealized';DATA=REPORT/'data'
def read(name):return json.loads((DATA/name).read_text())
def table(headers,rows):
    return '\n'.join(['| '+' | '.join(headers)+' |','|'+'|'.join(['---']*len(headers))+'|']+['| '+' | '.join(map(str,row))+' |' for row in rows])
def write(name,text): (REPORT/name).write_text(text.strip()+'\n')
def pct(values,digits=4):return ' / '.join(f'{100*x:.{digits}f}%' for x in values)
ASSUMPTION='本轮采用理想化模型：2410 为主动指定的水力源点，其余 122 个结构叶节点为 p_ref=0 Pa 的参考压力终端；0 Pa 是表压基准。这些定义均不是生理标注，结果不代表真实小鼠脑循环 ground truth。'

def run():
    s=read('final_summary.json');op=read('roi_operating_point_H0.json')['ports'];bc=read('network_to_fem_pressure_transfer.json');mass=s['network_mass_residual']
    sensitivity=read('terminal_sensitivity_summary.json');down=read('downstream_resistance_audit.json')
    loo=list(csv.DictReader((DATA/'terminal_leave_one_out_sensitivity.csv').open()))
    radius=list(csv.DictReader((DATA/'radius_hydraulic_sensitivity.csv').open()))
    diff=read('3d_case_config_diff.json');remote=read('3d_remote_run.json')
    done=s['new_3D_case_status']=='PASS'
    v=read('3d_physics_validation.json') if done else None
    optable=table(['端口','真实切点表压 Pa','有符号流量 m³/s','流量 pL/s','Q/Q_ROI','圆管平均速度估计 mm/s'],[
        [r['port'],f"{r['pressure_realcut_Pa']:.8f}",f"{r['signed_Q_m3s']:.12e}",f"{r['signed_Q_m3s']*1e15:.9f}",f"{100*r['signed_fraction']:.6f}%",f"{r['mean_velocity_estimate_m_s']*1e3:.7f}"] for r in op])
    bctable=table(['端口','R_extension Pa·s/m³','H0 Q pL/s','真实切点 Pa','延长段压降 Pa','原始端盖 Pa','平移后端盖 Pa'],[
        [r['port'],f"{r['R_extension_Pa_s_m3']:.9e}",f"{r['signed_Q_m3s']*1e15:.9f}",f"{r['pressure_realcut_Pa']:.6f}",f"{r['extension_deltaP_Pa']:.6f}",f"{r['pressure_cap_raw_Pa']:.6f}",f"{r['pressure_cap_shifted_Pa']:.6f}"] for r in bc['ports']])
    compare=table(['出口','当前 3D %','H0 %','变化 百分点','相对变化 %'],[[n,f"{100*a:.6f}",f"{100*b:.6f}",f"{100*(b-a):+.6f}",f"{100*(b-a)/a:+.6f}"] for n,a,b in zip(['O1','O2','O3'],s['current_3D_split'],s['H0_split'])])
    sentable=table(['出口','H0 %','最小 %','最大 %','P05 / P50 / P95 %','影响最大终端','最大绝对变化 百分点'],[
        [n,f"{100*r['baseline']:.6f}",f"{100*r['min']:.9g}",f"{100*r['max']:.6f}",' / '.join(f'{100*x:.6f}' for x in r['P05_P50_P95']),r['most_influential_terminal_id'],f"{r['maximum_absolute_difference_pp']:.6f}"] for n,r in sensitivity['ports'].items()])
    influential=sorted(loo,key=lambda r:max(abs(float(r[n+'_difference_pp'])) for n in ['O1','O2','O3']),reverse=True)[:3]
    worsttable=table(['移除终端','O1 %','O2 %','O3 %','重新匹配所需源压力 Pa'],[[r['removed_terminal_id'],*[f"{100*float(r[n+'_fraction']):.9g}" for n in ['O1','O2','O3']],f"{float(r['source_pressure_Pa']):.6f}"] for r in influential])
    radtable=table(['作用区域','半径倍数','源压力 Pa','源压力/H0','O1 %','O2 %','O3 %','O2 变化 百分点'],[
        [r['scope'],r['radius_factor'],f"{float(r['source_pressure_Pa']):.6f}",f"{float(r['source_pressure_ratio']):.9f}",*[f"{100*float(r[n+'_fraction']):.6f}" for n in ['O1','O2','O3']],f"{float(r['O2_difference_pp']):+.6g}"] for r in radius])
    summary3d='新 case 已建立并启动；等待实际 FEM 收敛、全局质量守恒及独立验证。此时不报告新 3D 分流。'
    causal='H0 相比原 3D 的 O2 分流下降 35.6150 个百分点，但这同时涉及模型维度差异，不能仅凭两组数字把下降量全部归因于出口边界。独立 3D 验证结果出来后再量化相同 FEM 下的边界效应。'
    if done:
        m=v['measurements'];c=v['comparison'];change=c['new_minus_current_3D_pp'][1]
        summary3d=f"独立新 3D 验证通过。O1/O2/O3 = {pct(s['new_3D_split'])}；全局相对质量残差 {m['epsilon_mass']:.3e}，低于原阈值 1e-6。第 {v['steady']['first_qualifying_step']} 步满足连续 5 个稳态区间，第 {v['steady']['final_step']} 步结束。"
        causal=f"在保持同一 geometry、mesh、入口、μ、ρ、WALL、P1/P1+VMS 和数值设置的受控比较中，O2 从 {100*c['current_3D_remeasured_split'][1]:.6f}% 变为 {100*s['new_3D_split'][1]:.6f}%，改变 {change:+.6f} 个百分点（相对变化 {100*c['new_vs_current_relative_change'][1]:+.4f}%）。这是本次理想化压力替换在当前 3D 模型中的出口边界效应。它说明原 O2≈85% 对等出口压力条件敏感；不能把这一变化比例推广为真实脑血流中“边界造成的比例”。"
    massparagraph=f"内部节点最大绝对不平衡 {mass['max_absolute_nodal_imbalance_m3s']:.6e} m³/s，RMS {mass['rms_nodal_imbalance_m3s']:.6e} m³/s；相对于全 A 源流量的最大节点残差 {mass['max_relative_nodal_imbalance_to_A_source']:.6e}。全 A 源入流 {s['A_total_source_inflow']:.12e} m³/s，终端总出流 {s['A_total_terminal_outflow']:.12e} m³/s，全局相对残差 {mass['global_relative_residual']:.6e}。ROI 相对残差 {mass['ROI_relative_residual']:.6e}，目标流量相对误差 {mass['ROI_target_relative_error']:.1f}。小的非零残差全部保留。"
    write('H0_MODEL_ASSUMPTIONS_ZH.md',f'''# H0 模型假设

{ASSUMPTION}

1. structural root node **2410** 被本模型正式指定为 `IDEALIZED_HYDRAULIC_SOURCE`。机器摘要使用 `source_definition=IDEALIZED_STRUCTURAL_ROOT`。无需再以缺少生理入口标签为由停止。
2. 唯一水力域是已经有 provenance 的 analysis A：**7419 个原节点、7418 条原边、1 个连通分量、0 个 SWC 图环**。精确插入 3 个虚拟切点后求解图为 7422 节点、7421 边。完整发布文件其余 42 个分量仅作溯源背景。
3. H0 = `ALL_STRUCTURAL_LEAVES_REFERENCE_PRESSURE`。除 2410 外全部 degree=1 结构叶节点共 **122 个**，均定义为 `IDEALIZED_REFERENCE_PRESSURE_TERMINAL`，p_ref=0 Pa。统一设零是可重复的参考边界假设，不是测量结论。所有末端同压是一项额外的理想化假设；共同平移压力基准并不能把本来不同的末端压力变成相同。
4. O3 原节点 **4484** 本身就是 H0 终端；p_O3(real)=0 Pa 是模型定义。`geometry-derived R3 = NOT_AVAILABLE`；不存在的 A 下游几何不生成有限 R3。
5. 模型为 geometry-informed、steady、rigid、Newtonian、Poiseuille。复用 v1 的 SI 单位、沿边线性半径和稳定解析阻力积分。μ={s['mu_Pa_s']} Pa·s，ρ={s['rho_kg_m3']} kg/m³；1D 稳态阻力不依赖密度，3D 保持该密度。
6. 不加入 Fahraeus–Lindqvist、hematocrit、phase separation、非牛顿黏度或血管顺应性。没有为让分流更平均而调节半径或终端压力。
7. 用途是考察 ROI 截断和等出口压力对分流的影响。不能称为 true physiological mouse cerebral flow。

# 工作点与流向

先取源点 1 Pa、所有终端 0 Pa。单位源压力下 ROI 入口 Q={s['unit_ROI_inlet_Q']:.12e} m³/s，方向与 FEM 入口一致，因此 λ={s['scaling_lambda']:.12f}。所有压力、流量乘此 λ，使 ROI 入口达到 {s['ROI_target_Q_m3s']:.12e} m³/s。源点表压 {s['scaled_source_pressure_Pa']:.9f} Pa 仅是该工作点的模型压力尺度。

全 A 源入流 {s['A_total_source_inflow']*1e15:.9f} pL/s 与 ROI 入口 {s['ROI_inlet_Q']*1e15:.9f} pL/s 分开输出。ROI 流向由真实切点几何、保留的 ROI 侧相邻边和 outward normal 定义；SWC parent-child 只作代数存储方向。入口正号表示入 ROI，出口正号表示出 ROI，绝不对入口流量取绝对值以强行匹配。

原始来源、精确映射和 40 项输入 SHA 沿用 v1 并在运行前重验；[v1 provenance](../a_network_1d0d_boundary_v1/A_NETWORK_PROVENANCE_ZH.md) 保持只读。求解数据见 [graph](data/analysis_A_H0_graph_si.npz)、[方向定义](data/roi_port_flow_sign_convention.json)、[unit solve](data/unit_solve_summary.json)。
''')
    write('ROI_H0_OPERATING_POINT_ZH.md',f'''# ROI H0 工作点

{ASSUMPTION}

{optable}

三出口均为正向出流，无端口 backflow。Q1+Q2+Q3={sum(r['signed_Q_m3s'] for r in op[1:]):.12e} m³/s，与入口 {s['ROI_inlet_Q']:.12e} m³/s 一致至数值求解精度。分流分母统一为 Q_ROI_in；机器数据同时保留 positive-outflow normalized fraction。

表中速度由各真实 cut 的 SWC 半径以 Q/(πr²) 估计。入口得到 1.9597241 mm/s，与 FEM 端面实际面积定义的 2.0 mm/s 不同，是面积定义和几何位置的差别；匹配的是已经验证的 ROI 流量，不修改 FEM 入口。

# 单位压力与缩放

单位源压力 1 Pa 时 ROI 入口 Q={s['unit_ROI_inlet_Q']:.12e} m³/s；λ={s['scaling_lambda']:.12f}；缩放源压力={s['scaled_source_pressure_Pa']:.9f} Pa。它是理想化 H0 在指定工作点下的 gauge pressure scale，不能解读为真实小鼠动脉压力。全 A 入流 {s['A_total_source_inflow']*1e15:.9f} pL/s，远大于穿过该 ROI 的 {s['ROI_inlet_Q']*1e15:.9f} pL/s。

# 与当前 3D 参考比较

{compare}

O2 仍为 H0 中最大的单一出口（49.5900%），但 O3 已达 44.1077%，原先 85.2050% 的强优势显著改变。描述上属于 C（优势显著改变），并保留“仍略高于 O3”的事实；没有为 A/B/C/D 强设阈值。

# 数值守恒

{massparagraph}

[审计 JSON](data/network_mass_balance_audit.json) 保存稀疏求解诊断。[独立重加和审计](data/independent_scaled_flow_mass_audit.json) 直接从已缩放边流量重新累加，最大相对节点残差约 1.386e-13，确认没有靠人工归零得到结果。

![质量守恒](figures/01_network_mass_balance.png)

![当前 3D 与 H0](figures/04_roi_split_current_vs_H0.png)
''')
    refinement=read('extension_station_refinement.json')['variants']
    refine_rows=[]
    for r in bc['ports']:
        fine=next(x['R_Pa_s_m3'] for x in refinement if x['port']==r['port'] and x['station_count']==160)
        refine_rows.append([r['port'],f'{fine:.9e}',f"{100*(fine/r['R_extension_Pa_s_m3']-1):+.6f}%"])
    write('NETWORK_TO_FEM_PRESSURE_TRANSFER_ZH.md',f'''# 从 analysis A 真实切点到 FEM 人工端盖

{ASSUMPTION}

网络压力属于 **REAL_A_NETWORK_CUT**，FEM 出口面位于 **FEM_ARTIFICIAL_EXTENSION** 的远端。沿出口 outward 方向使用有符号 Q：`p_cap_raw = p_real_cut - R_extension * Q_outward`。正向出流沿延长段压降，反向流则反号；不能直接复制真实 cut 压力。

{bctable}

共同减去最小原始端盖压力 −295.572318375502 Pa，等价于给三端盖都加 295.572318375502 Pa，得到 **585.286433259832 / 2932.015271078201 / 0 Pa**。所有两两压差的保持误差最大 {bc['max_pressure_difference_preservation_error_Pa']:.6e} Pa，仅为浮点舍入。

# 延长段估计与量级

在当前冻结的 FEM exterior surface 上，沿 v1 确认的真实 cut→cap 轴向长度，取 40 个等距半格中心截面。通过 VTK 交线的闭合多边形计算面积，换算等面积圆半径，复用 v1 解析线性半径阻力积分；两端半格采用最近截面半径闭合，覆盖完整长度。每个截面必须唯一、闭合，不使用缺失几何的猜测替代。[截面与网格 SHA](data/current_FEM_extension_sections.json) 保留全部数值。

40→160 截面加密核查：

{table(['端口','160 截面 R Pa·s/m³','相对 40 截面的改变'],refine_rows)}

这验证了本几何下的轴向积分稳定性，不消除等效圆管 Poiseuille 对非圆截面、三维发展和接头效应的模型误差。实际 case 使用声明的 40 截面估计。

真实切点 O2−O1 压差为 2723.266598 Pa，延长段修正使它减少 376.537760 Pa，约 13.827%。O1−O3 压差由 342.341401 Pa 变为 585.286433 Pa，修正达 242.945032 Pa，约原差值的 70.97%。因此不能声称这些延长段修正可以忽略。

# O1/O2 外部阻力审计与 O3 闭合

O1 外部保存了 61 条边、1 个参考终端，独立单位端口压力求解得到 R={down['O1']['R_Pa_s_m3']:.12e} Pa·s/m³；与 full-A 的 p/Q 相对差 {down['O1']['relative_equivalence_error']:.3e}。O2 外部保存了 76 条边、1 个参考终端，R={down['O2']['R_Pa_s_m3']:.12e}，相对差 {down['O2']['relative_equivalence_error']:.3e}。这些是对已存几何的等效阻力审计；本次 FEM 使用 fixed pressure，未改成 coupled resistance。

O3：`geometry-derived downstream R3 = NOT_AVAILABLE`。H0 closure 为 `Dirichlet p_ref terminal`，真实 cut 压力 0 Pa。人工延长段 R_extension 仍可估计，必须与不存在的 A 下游 R3 区分。

# 写入 FEM 的边界语义

保持原 `Type=Neu`、`Time_dependence=Steady`，仅替换三个出口 `Value`；O3 平移后仍为零，实际数值变动只有 O1、O2。源代码的稳态 Neu 分支将正的压力值转换为负外法向 traction（`set_bc_neu_l` 的 `hg=-h*gx`），故无需再人为反号。它是给定压力牵引，不是每个端面节点的 pressure Dirichlet；求解出的面积平均压力可因黏性法向牵引而略不同。

FEM 保留原入口 XML 的 {diff['inlet_Q_preserved_m3s']:.12e} m³/s，网络匹配的是实测 {s['ROI_target_Q_m3s']:.12e} m³/s；既有差别仅约 2.87e-10 相对量，不能为了抹去这一差别而额外修改入口。

此次是单次 H0 工作点压力传递。若新 3D 分流不同，它不会迭代更新 R_extension*Q，也不等价于一个与全 A 实时耦合的多端口阻力边界。

[压力传递 JSON](data/network_to_fem_pressure_transfer.json) · [可写入 BC JSON](data/roi_fixed_pressure_bc_H0.json) · [case 配置审计](data/3d_case_config_diff.json)

![通俗压力边界解释](figures/07_pressure_boundary_explanation.png)
''')
    write('TERMINAL_RADIUS_SENSITIVITY_ZH.md',f'''# H0 的终端与半径敏感性

{ASSUMPTION}

# 逐一移除终端

每次只取消一个终端的 p=0 Dirichlet 条件，把它作为没有源汇的自由叶节点，其连续性条件就是 zero-flow；其余终端维持原条件。每个变体重新做单位压力稀疏求解，再重新匹配同一个 ROI inlet Q。122 个变体全部有效，无 INVALID_VARIANT；未运行任何敏感性 3D case。

{sentable}

P05/P50/P95 是这个有限、确定性的“逐一删除”集合的分位数，不是生理置信区间。119 个非局部终端在固定 ROI 入口流量后几乎不影响下游分流；由树结构、单入口 ROI 与三个下游终端的拓扑可解释这一现象。分位数因此几乎都落在 baseline，不能用它掩盖最坏情形。

{worsttable}

移除 O1 对应的 3368、O2 对应的 4461 或 O3 本身 4484 都会改变路由。O2 在移除 4461 后数值上接近零；移除 4484 后 O2 升至 69.1305%，O1 升至 30.8695%。O1 最小有符号 fraction 为 −7.026233400379352e-15，这是未被归零的求解舍入量，不能解释为已解析的物理 backflow。所有精确符号和值保留在 [逐终端 CSV](data/terminal_leave_one_out_sensitivity.csv)。

结论：H0 对多数终端不敏感，但对决定 ROI 三条出路的少数局部终端很敏感，不能称为“对任意单终端假设都稳健”。

![终端敏感性](figures/05_terminal_sensitivity.png)

# 半径敏感性

{radtable}

全局半径统一乘 a 时所有阻力按 a⁻⁴ 缩放。实际重算 0.90/0.95/1.00/1.05/1.10 后，分流最大变化约 3.3e-12 个百分点，源压力比例与 a⁻⁴ 相符。它验证线性同质缩放，不表示对各支路独立的半径误差同样不敏感。

O1 和 O2 的“external branch”定义为去掉 ROI 内部边后，与各 cut 相连的保存外部分量。改变该分量内所有真实节点半径，**保持虚拟 real-cut 半径不变**，避免把 ROI 内部边一起改变；第一段外部边的线性半径因此连续过渡。O1 外部 ±5% 使其自身分流变化约 −0.6741 / +0.6483 个百分点；O2 外部 ±5% 使 O2 分流变化 −4.9730 / +4.7127 个百分点。

O3 downstream radius perturbation = N/A，原因是未保存任何 A 下游几何，而非假定某个人造 R3。人工 FEM 延长段不混入 A 几何半径敏感性。

[半径 CSV](data/radius_hydraulic_sensitivity.csv) 保存每次真实求解的源压力、ROI 入口压力、各出口压力、流量比例及变化；没有随机采样。
''')
    write('A_NETWORK_H0_REVIEW_ZH.md',f'''# 一句话结论

H0 下 O1/O2/O3 = **{pct(s['H0_split'])}**，相比原 3D 的 4.257917% / 85.205026% / 10.537057%，O2 优势明显改变：下降 35.6150 个百分点，仍是最大单一出口，但与 O3 已接近。{summary3d}

当前阶段：**{s['final_stage_status']}**。

{ASSUMPTION}

# 这一轮改变了什么

采用透明的 idealized hydraulic baseline 来回答当前几何中的模型问题；不再追求不存在的真实边界标注。复用 v1 provenance、analysis A、精确虚拟切点、SI、阻力积分和稀疏求解基础。未修改 Particle、原 production flow 或 v1 报告。

# 2410 为什么可以使用

2410 是 structural root，本轮主动定义它为理想化水力源；没有证明它是真实动脉。源 1 Pa 时 ROI 入口为 {s['unit_ROI_inlet_Q']:.12e} m³/s，方向正确。λ={s['scaling_lambda']:.9f} 使 ROI 入口匹配 15.5135916044 pL/s，源表压得到 5226.146255 Pa，这是该模型的压力尺度。

# terminal 为什么全部设 0

除 2410 外所有 degree=1 节点共 122 个，统一作为参考压力端口，便于可重复地闭合管网。0 是 gauge/reference，不是真实组织压力。O3 原节点 4484 就是其中一个端点，其 real-cut 压力为 0 是定义结果；没有可从 A 推导的有限 downstream R3。

# 完整 A 如何决定 ROI pressure

把 A 看成相互连通的细水管。每段的长度和半径决定阻力；流量在分叉处连续，流经管段会损失压力。ROI 外仍保存的 O1、O2 管段将真实切点与参考端相连，所以二者切点有不同压力。O3 切点就是参考端。只解 analysis A 的 7419 原节点，不把其它 42 个分量混入。

![完整 A 压力](figures/02_analysis_A_pressure.png)

![完整 A 流量](figures/03_analysis_A_flow.png)

# ROI operating point

{optable}

全 A 源入流为 {s['A_total_source_inflow']*1e15:.6f} pL/s；ROI 入口只有 {s['ROI_inlet_Q']*1e15:.6f} pL/s。三出口均向外，无 backflow。{massparagraph}

# O2=85% 到底有多少是 ROI 边界造成的

{compare}

{causal}

# terminal sensitivity

全部 122 个 leave-one-terminal-out 变体均重算并匹配同一 ROI 流量。大多数几乎不改变 ROI，但三个局部终端强烈影响分流：移除 4461 时 O2 接近零；移除 4484 时 O3 为零，O2 变为 69.1305%。O2 最坏变化 49.5900 个百分点。应同时看最坏情形与分位数，不能只看几乎不变的中位数。详见 [敏感性报告](TERMINAL_RADIUS_SENSITIVITY_ZH.md)。

# radius sensitivity

全局半径 0.90–1.10 倍实际数值验证了分流不变、压力按 a⁻⁴ 改变。O2 外部半径 −5% / +5% 则使 O2 分流改变 −4.9730 / +4.7127 个百分点；局部半径误差仍重要。O3 下游为 N/A。

# network pressure 如何传到 FEM cap

真实切点压力不能直接复制到人工延长段端盖。分别扣除 52.6273 / 429.1650 / 295.5723 Pa 压降，再共同加 295.5723 Pa。最终端盖压力 **585.286433 / 2932.015271 / 0 Pa**。40→160 截面的阻力估计变化均小于 0.04%，但圆管近似仍有模型误差。详见 [压力传递报告](NETWORK_TO_FEM_PRESSURE_TRANSFER_ZH.md)。

# 新 3D FEM

{summary3d}

新 case：`{diff['new_case']}`。服务器独立目录：`{remote['remote_root']}`。原 case 的网格、geometry、入口、μ、ρ、WALL、求解器和数值参数保持不变，仍是 P1/P1 + VMS。配置审计只允许三个 outlet Value；其中 O3 数值仍为零。全局 mass_limit 仍为 1e-6，未增加任意内部截面 1e-6 门槛。

{('[三代分流主图](figures/06_three_way_flow_split.png) · [3D 验证](ROI_3D_A_H0_PRESSURE_VALIDATION_ZH.md)' if done else '3D 正在独立验证，未把新场接入 Particle。')}

# 1D/0D 与 3D 为什么不同

1D/0D 是 Poiseuille 圆管网络近似；3D 含真实接头、曲率、非圆截面、三维发展和人工延长段。固定压力传递基于 H0 工作点，并不随新 3D 流量迭代更新。两者分流不要求严格相等，主要受控比较是旧 3D 与新 3D。

# 当前最主要限制

结构 root 和 leaves 是模型边界假设；没有真实压力或入口标注。SWC 树不等于完整生理连通图，几何截断和局部半径不确定性仍存在。少数局部终端删除会改变路由。人工段压力损失用等面积圆管近似。既有 P1/P1+VMS 内部局部质量守恒问题属于独立限制，本轮没有用速度重建来混合改变变量。

# 下一步

1. 审核固定压力下的新 3D 分流及其与 H0 的差异。
2. 如需进一步降低边界假设依赖，再单独研究多端口耦合与局部半径/终端不确定性。
3. 审核完成后再另行决定 Particle 接入；本轮未运行 Particle。
''')
    if done:
        m=v['measurements'];c=v['comparison'];ex=json.loads((Path(diff['new_case'])/'reports/execution.json').read_text())
        split_table=table(['出口','当前 3D %','Full-A H0 %','新 3D %','新−旧 3D 百分点','新 3D−H0 百分点'],[
            [n,f'{100*a:.6f}',f'{100*b:.6f}',f'{100*cval:.6f}',f'{100*(cval-a):+.6f}',f'{100*(cval-b):+.6f}'] for n,a,b,cval in zip(['O1','O2','O3'],c['current_3D_remeasured_split'],s['H0_split'],s['new_3D_split'])])
        flows=table(['端口','有符号约定','流量 m³/s','流量 pL/s'],[['INLET','进入 ROI 为正',f"{m['Q_in_m3_s']:.12e}",f"{m['Q_in_m3_s']*1e15:.9f}"]]+[[n,'离开 ROI 为正',f'{q:.12e}',f'{q*1e15:.9f}'] for n,q in m['outlet_flows_m3_s'].items()])
        write('ROI_3D_A_H0_PRESSURE_VALIDATION_ZH.md',f'''# 新 3D FEM 独立验证

**PASS / {s['final_stage_status']}**。{ASSUMPTION}

{summary3d}

{split_table}

{causal}

![三代结果](figures/06_three_way_flow_split.png)

# 数值、几何与边界验收

{flows}

出口总流量 {m['Q_out_total_m3_s']:.12e} m³/s；原生产双精度求和方法的全局相对质量残差 **{m['epsilon_mass']:.9e}**，入口目标相对误差 {m['epsilon_Q']:.9e}。这两项均用原入口指定流量作分母并满足原阈值 1e-6。对同一组四个有符号边界通量作补偿求和，残差为 {m['global_signed_flux_residual_compensated_m3s']:.12e} m³/s，相对量 {m['epsilon_mass_compensated_boundary_sum']:.9e}；普通求和显示的零来自浮点舍入，没有人工归零，也不意味着任意内部截面严格守恒。端口反流判定：{m['port_backflow']}。有符号 fraction 使用 actual Q_in 作分母，另保存正出流归一化比例。

压力范围 {m['pressure_range_pa'][0]:.8f} 到 {m['pressure_range_pa'][1]:.8f} Pa，速度、压力均有限；WALL 最大速度 {m['wall_velocity_max_m_s']:.9e} m/s，小于 {m['wall_noslip_tolerance_m_s']:.9e} m/s，no-slip PASS。给定的是 Neu 压力牵引，面积平均 p 不需逐点等于指定压力。

原生 solver exit code=0；native log 独立重解析后线性与非线性 convergence PASS。每 10 步保存，最后连续 5 个区间满足 E_u≤1e-5、E_Q≤1e-6；第 {v['steady']['first_qualifying_step']} 步首次达标，最终第 {v['steady']['final_step']} 步。停止后的速度相对变化 {v['steady']['final_relative_velocity_change']:.6e}，归一化通量变化 {v['steady']['final_normalized_flow_change']:.6e}。本次失败的最终线性求解数为 {ex['history']['failed_linear_solves']}，已恢复的 ILU 重试数为 {ex['history']['recovered_attempts']}。原始求解日志完整保留。

网格节点 {v['export']['node_count']}、四面体 {v['export']['tetra_count']}。验证坐标顺序、单元连通、mesh/exterior/WALL 的 SHA 与原 case 一致：

| 输入 | SHA256 |
|---|---|
| mesh | {v['export']['mesh_sha256']} |
| exterior geometry | {v['export']['geometry_sha256']} |
| WALL | {v['export']['wall_sha256']} |

配置比较 [3d_case_config_diff.json](data/3d_case_config_diff.json) 证明唯一允许字段是 O1/O2/O3 outlet Value；O3 数值维持零。未切换 FEM formulation，未做 conservative velocity reconstruction，未将内部任意截面误差作为本轮 gate。

# 1D/0D 与 3D 的差异

新 3D−H0 的分流差为 {' / '.join(f'{x:+.6f}' for x in c['new_minus_H0_pp'])} 个百分点。圆管网络近似和三维接头、曲率、非圆截面、发展段的差异允许造成不一致。H0 端盖压力在其工作点计算后固定，未强迫 FEM 重现网络流量；严格相等不作为验收标准。

# 运行与保留结果

1 MPI rank，OMP_NUM_THREADS=1；复用原 GPU/PETSc 配置。服务器求解耗时 {v['wall_time_s']:.3f} s（{v['wall_time_s']/60:.2f} min）。5 秒采样的进程树峰值 RSS {v['peak_tree_RSS_MiB_sampled']:.2f} MiB，最大单进程 VmHWM {v['peak_single_process_HWM_MiB']:.2f} MiB，wait4 子进程 maxrss {v['child_rusage_maxrss_MiB']:.2f} MiB；采样峰值不宣称捕获采样间瞬时总峰值。完整环境见新 case 的 reports/host.json。

服务器：`{remote['remote_root']}`。

本地独立 CFD 场：`{diff['new_case']}/frozen_flow/steady_flow_mean_2p0_mmps_A_H0.vtu`。该文件是实际最终 native VTU 的逐字节副本，附带 SI 数组、checkpoint 和 SHA manifest；没有由旧场缩放生成。结果未晋升为 Particle production，未执行任何 Particle。

[独立验收 JSON](data/3d_physics_validation.json) · [主审阅报告](A_NETWORK_H0_REVIEW_ZH.md)
''')
    tests=ET.parse(REPORT/'logs/pytest.xml').getroot()
    suites=list(tests.iter('testsuite'));nt=sum(int(x.get('tests','0')) for x in suites);nf=sum(int(x.get('failures','0'))+int(x.get('errors','0')) for x in suites)
    write('README.md',f'''# A-Network Idealized Hydraulic Baseline Model v2

{ASSUMPTION}

状态：**{s['final_stage_status']}**。

- [主审阅报告](A_NETWORK_H0_REVIEW_ZH.md)
- [明确的 H0 假设](H0_MODEL_ASSUMPTIONS_ZH.md)
- [ROI 工作点与守恒](ROI_H0_OPERATING_POINT_ZH.md)
- [真实切点到 FEM cap](NETWORK_TO_FEM_PRESSURE_TRANSFER_ZH.md)
- [终端和半径敏感性](TERMINAL_RADIUS_SENSITIVITY_ZH.md)
{('- [新 3D 独立验证](ROI_3D_A_H0_PRESSURE_VALIDATION_ZH.md)' if done else '')}
- [机器摘要](data/final_summary.json)
- [测试日志](logs/pytest.log)：{nt} 项，{nf} 个失败/错误（v1 30 + v2 19）。

所有图有 PNG 与 PDF；Figure 06 在真实新 FEM 验证通过后生成。网络计算本地 1 worker，实测 {s['network_runtime_seconds']:.3f} s，峰值 RSS {s['peak_RSS_MiB']:.2f} MiB。新 FEM 使用独立服务器目录，不覆盖原 production。

# 复现入口

从 `ulm_3D_vascular` 根目录使用项目已有科学 Python 环境：

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 python scripts/run_a_network_h0_v2.py
python -m pytest tests/network_1d0d tests/network_h0 -q
python scripts/render_a_h0_v2.py
python scripts/write_a_h0_reports.py
```

网络脚本只写本 v2 目录。若完整重跑，网络脚本会将摘要重置为 3D pending；必须再用已保留的实际 FEM 输出执行 `scripts/validate_a_h0_fem.py --case <新case> --flow-root <FEM_SimVascular根目录>`，再生成图表和报告。

创建新 FEM 输入使用 `network_1d0d.fem_h0_case.create_case(original, fresh_target, bc_json)`，目标存在就拒绝覆盖。服务器启动命令与环境参数保存在 `data/3d_remote_run.json`，runner 为 `scripts/solve_a_h0_fem_remote.py --case <remote_case> --reference-root /workspace/flow_mean_2p0_mmps_20260922`。实际 executable/PETSc SHA、命令、MPI/OMP、完整求解历史记录在新 case 的 reports/execution.json。长期保存本地结果，不依赖服务器永久保留目录。

当前分支 `dev/a-network-1d0d-idealized-baseline-v2`，本轮仅新增代码、测试、报告和独立 flow case；不提交或清理已有大量历史修改。Git 快照与最终 SHA 复验见 logs/。

最终 tracked diff-stat：`{s.get('git_diff_stat_summary','待最终快照')}`。并发工作区变化记录：{s.get('concurrent_workspace_change','见最终工作区审计')}。保留已有修改，未回退其它工作的 `.gitignore` 规则。文件保护结果见 [最终输入审核](logs/final_input_protection_audit.json)，工作区差异见 [Git 审核](logs/git_workspace_final_audit.json)。
''')
    print('Wrote Chinese reports; stage:',s['final_stage_status'])

if __name__=='__main__':run()
