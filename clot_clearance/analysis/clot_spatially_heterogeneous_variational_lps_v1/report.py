from common import *
from literature import REFERENCES
import html,datetime

def main():
 data=json.loads((O/'COMPARISON.json').read_text());cs=data['cases'];sp=data['spatial'];meta=json.loads((O/'HETEROGENEITY_METADATA.json').read_text());audit=json.loads((O/'VARIATIONAL_ENERGY_AUDIT.json').read_text());path=json.loads((O/'PRESCRIBED_PATH_ENERGY.json').read_text());videos=json.loads((O/'VIDEO_GENERATION.json').read_text());validation=json.loads((O/'DELIVERY_VALIDATION.json').read_text());resources=json.loads((O/'RUN_RESOURCES.json').read_text());drive=json.loads((O/'SHEAR_DRIVE_FREEZE.json').read_text());freeze=json.loads((O/'SOURCE_FREEZE.json').read_text());U=cs['uniform'];H=cs['heterogeneous'];diff=sp['difference'];full={m:bool(q['objective']['full_target']) for m,q in cs.items()}
 def paired(fn):return '；'.join(('U' if m=='uniform' else 'H')+'：'+fn(cs[m]) for m in CASES)
 def event(q,key):
  e=q['native']['events'].get(key)
  return 'NOT_OBSERVED' if e is None else f"{e['time_s']*1e3:.6f} ms，原始 bond ID {e.get('original_bond_ids',[])}，参考中点 {e.get('reference_midpoints_m',[])} m"
 def percent(x):return '未定义（没有脱落）' if x is None else f'{100*x:.4f}%'
 conclusion=('OBSERVED_EXPLORATORY_TARGET_SINGLE_REALIZATION' if full['heterogeneous'] else 'TARGET_NOT_MET_SINGLE_REALIZATION')
 status=dict(task=TASK,run=RUN,status='COMPLETED',label=HET_LABEL,
 INITIAL_KAPPA_ARRAY_IDENTITY=False,INITIAL_KAPPA_DIFFERENCE='Only already-precut complete-failure sentinel kappa=sc is remapped; initial d and all uncut kappa/d unchanged',HETEROGENEITY_FIELD_STATUS='FROZEN_REPRODUCIBLE_REFERENCE_FIELD',ORIGINAL_TOPOLOGY_PRESERVED=True,
 VARIATIONAL_ENERGY_PARITY='PASS',FORCE_GRADIENT_PARITY='PASS',IRREVERSIBLE_HISTORY_STATUS='PASS',REFERENCE_ELASTIC_RESPONSE='BITWISE_UNIFORM_FALLBACK_AND_COMMON_PRE_ONSET',FROZEN_SHEAR_PROTOCOL_STATUS='SAME_FROZEN_OPEN_LOOP',UNIFORM_NATIVE_STATUS=U['native']['status'],HETEROGENEOUS_NATIVE_STATUS=H['native']['status'],CRACK_PATH_DIFFERENCE=diff,MULTIFRAGMENT_FORMATION_STATUS=full,TARGET_FRAGMENT_20_200_COUNT={m:q['objective']['n_fragments_20_200'] for m,q in cs.items()},TARGET_FRAGMENT_MASS={m:q['objective']['mass_fraction_target_initial']*1792*V*RHO for m,q in cs.items()},LARGEST_FREE_COMPONENT={m:max(q['final_free_sizes'],default=0) for m,q in cs.items()},FINE_DEBRIS_FRACTION={m:q['objective']['fine_debris_mass_fraction_detached'] for m,q in cs.items()},SECONDARY_FRAGMENTATION_STATUS={m:dict(all_events=sp[m]['secondary_all_count'],parent_at_least20=sp[m]['secondary_parent_at_least20_count']) for m in CASES},ENERGY_ADMISSION_STATUS={m:q['native']['max_relative_energy_residual']<=.01 for m,q in cs.items()},VISUALIZATION_STATUS='PASS_PNG_PDF_ACTUAL_DATA',VIDEO_GENERATION_STATUS=videos['status'],SCIENTIFIC_CONCLUSION=conclusion,formal_full_protocol_cases=2,additional_seed_trials=0,optimization_trials=0,acoustic_realizability='UNESTABLISHED',surface_status='UNCHANGED_PARTIAL_GEOMETRIC_FLUID_ACCESSIBILITY_PROXY')
 assert all(q['native']['status']=='COMPLETED' for q in cs.values());assert validation['status']=='PASS';assert videos['status']=='PASS';dump(O/'FINAL_STATUS.json',status)
 lines=[f'# CLOT_SPATIALLY_HETEROGENEOUS_VARIATIONAL_LPS_V1\n\n运行：`{RUN}`。标签 `{HET_LABEL}`。\n',
 '## 本轮结论\n',
 ('在本次单一冻结材料场与同一开环剪切驱动下，H达到预先冻结的探索性碎块目标。仍需重复实现及尺度研究，不能称实验标定。' if full['heterogeneous'] else '在本次单一冻结材料场与同一开环剪切驱动下，H未达到至少5个20–200节点凝聚性游离块的完整探索目标。必须区分损伤路径改变、失效键数量和真实中等块体形成，不能用前两者代替后者。'),
 '\n'+paired(lambda q:f"末态20–200节点块 {q['objective']['n_fragments_20_200']}；历史最多 {q['native']['max_medium_components']}；最大自由块 {max(q['final_free_sizes'],default=0)} 节点；细屑/脱落质量 {percent(q['objective']['fine_debris_mass_fraction_detached'])}；实际净外功 {q['work']['actual_net_J']*1e9:.6f} nJ。"),
 '\n## 冻结模型、场与能量实现\n',
 '1792材料点、336固定锚点、95666原始键、2736原始预切；参考X、m、r0、taper、质量及邻域不变。dx=50 µm、horizon=150.5 µm、G=940 Pa、K=100000/3 Pa、ρ=1000 kg/m³。bulk s0=0.15；interface s0/sc=0.140625/0.16875；固定支架卡不变。所有材料属性以原始atom/bond ID对应。需要明确的初始化细节：未预切键初始κ=0、d=0；原始预切d=1保持。现有初始化用κ=sc作为预切完全失效标记，因此H预切键的κ标记随其局部sc变化，数值上不与U逐键相同。这不是新增变形或新增预切，已完全失效键仍无材料力；但不能声称两组所有初始κ数组逐值完全相同。',
 f"\nGSTools {meta['gstools_version']}；Gaussian z方差1，各向同性，显式rescale=1，C_z(r)=exp[-(r/300 µm)²]；RandMeth 4096模式，seed=20261017。q=tanh(z−c)，c使材料参考体积加权q均值为0。实际q范围[{meta['q_min']:.9f},{meta['q_max']:.9f}]，标准差{meta['q_std']:.9f}。单一冻结实现，没有结果后换种子。", 
 f"\n对全部原始bulk键使用qij=(qi+qj)/2、Δsij=(0.18−0.15)exp(0.30qij)/Z、scij=0.15+Δsij。Z={meta['Z']:.12g}；参考权重w=h*r0*V²*(1/mi+1/mj)，与完好参考水静pair刚度乘r0²成比例，包括原始bulk预切。sc范围{meta['sc_min']:.12f}–{meta['sc_max']:.12f}，加权平均{meta['sc_weighted_mean']:.12f}；排除预切后的均值{meta['sc_uncut_weighted_mean']:.12f}。归一化软化区间不等于归一化所有开裂路径耗能。",
 '\n完整能量仍为每个原始端点star的正定局部弹性/体积能与串联柔度罚能之和，对内部弹性伸长精确消元。固定d时，U(x,d)、−∂U/∂x及力装配公式未改；sc只进入原生不可逆演化 d(κ)=0 (κ≤s0)、sc(κ−s0)/[κ(sc−s0)] (s0<κ<sc)、1 (κ≥sc)。κ=max历史正伸长。D按固定当前x的U(x,d_old)−U(x,d_new)累加，数值R单列，不混入D。没有在节点力上乘随机因子。κ驱动本构没有显式速率变量；本轮惯性/拖曳响应不能视为Liu实验内禀断裂速率依赖的复现。',
 '\n新增complete_by_bond与原始bond ID冻结接口、restart V4；读旧V2/V3时回退typed sc。原有onset字段保留。设置接口必须在初始history导入/commit前调用，拒绝缺失/重复ID、非法sc及界面属性改动。输入逆序与atom排序后仍保持对应。隔离库及其工厂、协调器依赖重编译；旧库保持原哈希。',
 f"\n新库SHA `{freeze['library_sha256']}`；冻结SHA `{sha(O/'SOURCE_FREEZE.json')}`。原surface实现逐字节复用；积分仍为原生VV、速度相关末半步隐式回调、BASE=0.25 µs、最小预测子步0.03125 µs、Δd≤0.1及能量门槛1%。预测拒绝发生于native commit之前；提交后异常会保留failed_unaccepted而非声称已回滚。",
 '\n## 实际材料场结构与尺度\n',
 f"线性梯度R²={meta['linear_gradient_R2']:.6f}；26邻接局部极值统计{meta['extrema_26_neighbor']}；q>0与q<0分别{meta['positive_q_components_26']}、{meta['negative_q_components_26']}个连通区域。试样材料中心跨度650×350×750 µm；lc/dx={meta['lc_over_dx']:.4g}，lc/horizon={meta['lc_over_horizon']:.6f}。正负区域连通不等于只有两个组织，也不证明各局部极值统计独立。没有测得足够多的独立材料块；有限域和非局部作用可能削弱局部属性对分离路径的辨识。",
 '\n## 规定路径耗能（离线，不是native动力学，也不是Gc）\n',
 '完整原图在参考z=0两侧作±opening/2刚体平移；没有新增初始预切。0–180 µm，共721规定状态。初始刚度相同；第一非零损伤网格点7.75 µm并非精确起裂点。该点已经进入软化，反力不要求相同；严格起裂前的弹性与统一κ=s0阈值保持。',
 '\n|离线路径材料卡|完整分离损伤耗能 nJ|末端功积分余额 nJ|\n|---|---:|---:|',
 *[f"|{m}|{q['path_dissipation_J']*1e9:.8f}|{q['final_quadrature_residual_J']*1e9:.8f}|" for m,q in path['summaries'].items()],
 '\n3552条未预切跨面原始键最终全部失效；损伤能量增量非负。耗能仍有有限规定步长误差：外功梯形积分与储能+耗能约1.2%差异，不能将上述值称精确连续极限或GC_EFFECTIVE。最小/中位/最大sc三张卡仅为离线诊断，不是新增正式剪切工况。',
 '\n## 同一开环加载与实际成本\n',
 f"从上一轮checkerboard接受记录提取{drive['knots']}个幅值结点，含t=0时3200 Pa，按时间分段线性插值；单独乘原2 ms升、2.4 ms降、2.8 ms停驱包络。与旧FP64回调幅值最大误差{drive['max_native_key_amplitude_error_Pa']:.3g} Pa。两个正式case共享同一只读文件，完全取消等功反馈。250 µm棋盘空间模板、当前位置/法向、被动拖曳与开放裂面代理不变。实际两组功允许不同。这里是数千Pa机械需求剪切代理，不是已验证的真实微泡声流。",
 '\n|工况|净外功 nJ|驱动功 nJ|拖曳功 nJ|实际峰值剪切 Pa|最大能量残差 %|wall s|CPU s|峰值RSS MiB|\n|---|---:|---:|---:|---:|---:|---:|---:|---:|',
 *[f"|{m}|{q['work']['actual_net_J']*1e9:.6f}|{q['work']['drive_J']*1e9:.6f}|{q['work']['drag_J']*1e9:.6f}|{q['work']['actual_peak_shear_Pa']:.6f}|{q['native']['max_relative_energy_residual']*100:.6f}|{q['native']['wall_s']:.3f}|{q['native']['cpu_s']:.3f}|{q['native']['peak_RSS_bytes']/2**20:.3f}|" for m,q in cs.items()],
 '\n'+paired(lambda q:f"最大面片速度{q['work']['max_surface_speed_m_s']:.6f} m/s；规定外缘速度代理峰值{q['work']['max_outer_speed_m_s']:.6f} m/s；Re_h峰值{q['work']['max_Re_h']:.3f}；初始受力面积{q['work']['initial_area_m2']*1e6:.6f} mm²，按首次受力面积累计的独特面片面积{q['work']['cumulative_unique_area_m2']*1e6:.6f} mm²")+'。数千量级Re_h不支持把该薄层滑移代理解释为低Re Stokes预测；规定外缘速度不是CFD求解或实验测得的微泡声流。各时刻有效受力面积保存在LOADING_HISTORY.csv。',
 '\n完整账本：W_drive+W_drag=ΔU+ΔK+D+ΔUc−Jc+R。实际接受端点力平均点乘原生位移；Jc单列新接触排除解除产生的能量跃变。残差分母是历史最大能量尺度，不是简单末态W。所有接受行、每次损伤变化与真实分量成员保留；完整FP64坐标/速度/κ/d每25 µs及事件保存，未声称每接受步均保存完整x/v。',
 '\n## 必须回答的28项科研问题\n']
 answers=[
 '已实际接入原生不可逆更新及完整变分能量路径。sc参与κ→d；力和U通过实际d变化响应，未改最终力乘数。',
 '是。初态参考文件逐字节沿用；原始节点、95666键、2736预切及m不变；本轮检查和正式末态原生图核验通过。',
 f"不是简单线性梯度（R²={meta['linear_gradient_R2']:.4f}），有平滑三维变化；极值和正负连通域如上，经验相关函数在Fig2。", 
 '尚不能证明足够多统计独立材料区。大多数极值位于边界，正负场各连通，颜色数量不是组织区域数量。',
 f"bulk sc最小{meta['sc_min']:.12f}、最大{meta['sc_max']:.12f}、参考加权均值{meta['sc_weighted_mean']:.12f}。", 
 '保持。q=0新旧库三种状态逐位一致；H和U完好能量/弹性一致，bulk s0全部0.15，界面不改。',
 f"保持固定历史下的力—能量一致性；三个代表状态最大相对有限差分误差{max(q['relative_gradient_error'] for q in audit['checks']['variational']):.4g}，且异质低/高q节点均覆盖。力矩、总力及独立完整能量核对通过。", 
 '保持原生不可逆规则。未预切键初始κ/d一致；预切d=1相同，κ=sc初始化标记的数值差异如上披露。卸载不恢复、κ不下降、未提交试探不改历史；排序改变2126个atom槽位后ID/参数保持，restart V4保留sc、κ/d和D。未重跑完整MPI矩阵。',
 '是，同一冻结开环幅值、空间函数、包络和表面算法。不同当前表面和速度导致实际牵引及功不同，属于材料响应差异。',
 paired(lambda q:f"净外功{q['work']['actual_net_J']*1e9:.6f} nJ，驱动{q['work']['drive_J']*1e9:.6f}，拖曳{q['work']['drag_J']*1e9:.6f} nJ"),
 paired(lambda q:event(q,'first_softening'))+f"。首次原始键相同={diff['first_softening_same_original_bonds']}。",
 paired(lambda q:event(q,'first_failure'))+f"。首次原始键相同={diff['first_failure_same_original_bonds']}。",
 paired(lambda q:f"失效键共享端点网络出生{q['damage_regions']['birth_count']}次，最多同时{q['damage_regions']['max_concurrent_regions']}区")+'。这是损伤区域代理，不能自动解释为独立裂纹起源。',
 '真实三维连续裂面贯通尚未由完整连续表面证明。自由连通分量证明其与锚/主体之间的所有原始存活材料路径已切断，但小点链分量不能当作中等三维裂块。辅助开放裂面图和failed-edge区数均不充分。',
 paired(lambda q:f"自然bulk失效{int(q['native']['last_row']['bulk_failed'])}，interface失效{int(q['native']['last_row']['interface_failed'])}"),
 paired(lambda q:f"末态20–200块{q['objective']['n_fragments_20_200']}，历史最多{q['native']['max_medium_components']}，目标质量/初始{percent(q['objective']['mass_fraction_target_initial'])}"),
 paired(lambda q:f"末态最大{max(q['final_free_sizes'],default=0)}节点，历史最大{q['native']['history_max_free_nodes']}，最大/脱落质量{percent(q['objective']['largest_free_fragment_fraction_detached'])}"),
 paired(lambda q:f"1–2节点细屑/脱落质量{percent(q['objective']['fine_debris_mass_fraction_detached'])}，总脱落/初始{percent(q['objective']['total_detached_mass_fraction'])}"),
 '；'.join(f"{m}：全部二次分裂{sp[m]['secondary_all_count']}次，母块≥20节点{sp[m]['secondary_parent_at_least20_count']}次" for m in CASES)+'。细屑分裂不称凝聚性中等块二次断裂；内部存活键、birth、寿命和体积等效直径保存在真实拓扑历史。',
 '；'.join(f"{m}：q与首次软化时间的条件Spearman={sp[m]['q_vs_first_softening_spearman']}，q与末态d={sp[m]['q_vs_final_damage_spearman']}" for m in CASES)+'。低/高q四分位以全部合格原始bulk键为分母；U中q只是参考协变量。',
 '只能描述关联。只有一个随机实现和一对对照；位置、原始切口、暴露表面、加载及非局部耦合均会混杂。没有独立样本数、显著性或组织组成因果结论。',
 f"存在：lc/δ={meta['lc_over_horizon']:.4f}，lc/dx=6。试样Y方向仅350 µm，约1.17lc；非局部作用与有限域可限制有效局部对比。没有本轮尺度收敛证据。", 
 '本轮不能决定必须改变G/K、界面或损伤模式。单一sc场若无目标，不能据此否定真实组织异质性，也不能同时改多个机制追求成功。',
 '不应立即把不规则删键网络当修复。原始完整非局部图是本轮受控条件；未来如研究不规则离散，须先证明参考m、能量和离散收敛，且不能人工设计目标碎块。',
 '可支持：局部sc确实进入native；均匀回退/历史/能量检查通过；同一开环驱动下的实际功、真实断键、图分量、空间差异及视频。不能超出这一对真实工况。',
 '未标定：300 µm、seed、alpha、q与实际组织组成对应、sc分布及其宏观强度/韧性；当前剪切不是确定的微泡声流；细屑水动力和完整湿润裂面仍未解析。',
 f"动态结果图与三个MP4使用真实已接受FP64状态，位移×1。Fig1/2来自冻结初始材料场，Fig3含明确标记的离线规定路径和原生验证，不能把这些输入/离线图称全时程native动力学。每视频{videos['videos'][0]['native_timeline_frames']}个真实时序帧+60帧末态定格，24 FPS、1920×1080 H.264。各面板显示实际时间，最大两组配对时差{videos['max_time_alignment_difference_s']*1e6:.6g} µs；无轨迹插值。", 
 paired(lambda q:f"CPU {q['native']['cpu_s']:.3f}s、wall {q['native']['wall_s']:.3f}s、峰值RSS {q['native']['peak_RSS_bytes']/2**20:.3f}MiB")+'。构建、离线路径、静态图及视频成本另见RUN_RESOURCES和各资源JSON；并行wall不相加作为总历时。资源计量覆盖所列计算/后处理进程；依赖安装、交互检查、预览修图及网络传输未全量计时，不能解释为完整服务账单。']
 for i,text in enumerate(answers,1):lines.append(f'\n{i}. {text}\n')
 lines += ['\n## 表面与视频限制\n',paired(lambda q:f"末态{q['surface']['final']['free_components_without_surface']}/{q['surface']['final']['true_free_components']}个自由分量没有可表示的受支持三维湿润表面；未解析节点{q['surface']['final']['free_nodes_in_components_without_surface']}")+'。这些分量的流体受力未解析，不能解释为物理上无阻力。当前surface按支持四面体、张开和外部连通及有限遮挡探针近似；未知连续湿润面积无法由点数据确定。该限制两组相同，但会限制后续输运及二次碎裂解释。',
 '\n视频1比较损伤与真正自由分量；视频2展示近邻存活/软化键与透明失效轨迹，右侧q随原始ID运动；视频3展示真实分量和内部全部存活材料键。颜色由每组精确成员集合注册ID稳定生成，跨U/H同色不表示同一块。近邻显示阈值1.01dx不改变使用全部原始存活键的拓扑；无虚构连续壳、放大位移或视觉聚类。',
 '\n## 下一轮唯一优先方向\n',
 ('优先固定本轮成功设置，检验第二个独立冻结空间实现的可重复性；不立即进入多参数优化或声学可实现性结论。' if full['heterogeneous'] else '优先验证材料相关尺度在当前非局部离散中的可分辨性：围绕lc/δ≈2及窄试样边界占比，开展受控尺度研究；保持统一起裂与单独sc异质性，先建立足够空间采样与有效局部参数对比证据。暂不同时引入弹性随机场、界面变化或不规则删键网络。'),
 '\n## 文献支持与适用边界\n']
 for r in REFERENCES:lines.append(f"\n- [{r['citation']} — {r['title']}]({r['url']})：{r['support']} {r['limit']}\n")
 lines+=['\n历史10%异质性实验采用250 µm场、变动s0及球形压入，未产生目标碎块。它说明旧设置的结果，不能作为本轮sc-only剪切的直接均匀对照。', '\n## 交付与复现\n\n`SOURCE_FREEZE.json`、`SOURCE_CHANGES.diff`、`BUILD_AUDIT.json`记录代码/库；`HETEROGENEITY_FIELD.npz`、两个属性CSV和场元数据记录完整实现；两个NATIVE_RESULT、ENERGY_LEDGER、NATIVE_EVENT_HISTORY与TRUE_FRAGMENT_HISTORY记录接受证据；`REPRODUCE.md`给出环境和执行命令；`OPEN_RESULTS.html`链接图、视频与源数据。历史数据未删除；没有优化、微泡CFD或完整V3生产运行。\n']
 lines.append('\n绘图曾因空失效键集合报错，VTK视频曾在重建坐标控件后崩溃；修复仅限后处理并保留失败日志，没有重跑原生工况。失败渲染及短帧渲染检查未计入成功后处理资源总计。')
 report='\n'.join(lines)+'\n';name='REVIEW_CLOT_SPATIALLY_HETEROGENEOUS_VARIATIONAL_LPS_V1.md';(O/name).write_text(report)
 resources['build']=json.loads((O/'BUILD_AUDIT.json').read_text())['wall_s'];resources['path_audit']=path['wall_s'];resources['figures']=json.loads((O/'FIGURE_RESOURCES.json').read_text());resources['videos']=videos['videos'];dump(O/'RUN_RESOURCES.json',resources)
 # One static local index, with source links and explicit exploratory outcome.
 import os
 def href(p):return os.path.relpath(p,F)
 cards=[]
 for p in sorted((F/'figures').glob('*.png')):cards.append(f'<figure><a href="figures/{p.stem}.pdf"><img src="figures/{p.name}" loading="lazy"></a><figcaption>{html.escape(p.stem)} · <a href="figures/{p.stem}.pdf">PDF</a></figcaption></figure>')
 videohtml=''.join(f'<section><h3>{html.escape(v["file"])}</h3><video controls preload="metadata" poster="videos/{Path(v["file"]).stem}_COVER.png" src="videos/{v["file"]}"></video><p><a href="videos/{Path(v["file"]).stem}_FRAME_TIMES.csv">实际帧时间与源状态</a></p></section>' for v in videos['videos'])
 rows=''.join(f'<tr><td>{m}</td><td>{q["objective"]["n_fragments_20_200"]}</td><td>{max(q["final_free_sizes"],default=0)}</td><td>{percent(q["objective"]["fine_debris_mass_fraction_detached"])}</td><td>{q["work"]["actual_net_J"]*1e9:.6f}</td><td>{q["native"]["wall_s"]:.2f}</td></tr>' for m,q in cs.items())
 links=''.join(f'<li><a href="{href(O/f)}">{f}</a></li>' for f in [name,'FINAL_STATUS.json','SOURCE_FREEZE.json','HETEROGENEITY_METADATA.json','NODE_LOCAL_PROPERTIES.csv','BOND_LOCAL_PROPERTIES.csv','NATIVE_EVENT_HISTORY.csv','TRUE_FRAGMENT_HISTORY.json','FRAGMENT_SIZE_DISTRIBUTION.csv','FRAGMENT_MASS_DISTRIBUTION.csv','ENERGY_LEDGER.json','RUN_RESOURCES.json','REPRODUCE.md'])
 page=f'''<!doctype html><html lang="zh"><meta charset="utf-8"><title>Spatially heterogeneous variational LPS</title><style>body{{max-width:1400px;margin:35px auto;padding:0 24px;font:17px/1.6 system-ui;color:#233b4d;background:#f8fafb}}h1{{font-size:30px}}table{{border-collapse:collapse;width:100%;background:white}}td,th{{padding:10px;border:1px solid #b8c6cf}}figure,section{{background:white;padding:15px;border:1px solid #c8d3da;margin:25px 0}}img,video{{width:100%;height:auto}}a{{color:#156b89}}.note{{background:#f4eadc;padding:16px}}</style><h1>空间相关软化能否产生多个凝聚性游离块？</h1><p>1792材料点，原始非局部图和界面不变；只改变bulk软化区间。一个固定300 µm场，同一冻结棋盘涡旋开环剪切，两个真实3.6 ms原生计算。</p><p class="note">{conclusion}。本轮参数未标定；实际机械剪切不是已验证微泡声流。全部动态画面来自实际接受态，位移×1。</p><table><tr><th>工况</th><th>20–200节点块</th><th>最大自由块</th><th>细屑/脱落</th><th>净外功 nJ</th><th>wall s</th></tr>{rows}</table><h2>三个真实三维视频</h2>{videohtml}<h2>静态图（点击下载PDF）</h2>{''.join(cards)}<h2>证据与复现</h2><ul>{links}</ul><p>表面仍是受支持几何代理；无三维表面的细屑水动力未解析。单一随机实现只能支持本次对照和空间关联，不能推出一般血栓组织规律或已标定Gc。</p></html>'''
 (F/'OPEN_RESULTS.html').write_text(page);(O/'OPEN_RESULTS.html').write_text(f'<!doctype html><meta charset="utf-8"><a href="{os.path.relpath(F/"OPEN_RESULTS.html",O)}">打开本轮图、视频与科研报告</a>')
 print('REPORT_COMPLETE',status,flush=True)
if __name__=='__main__':main()
