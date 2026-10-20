"""Report numerical evidence, including a failed/incomplete formal case if present."""
from common import *
import html,csv

def markdown_blocks(lines):
 blocks=[]
 for line in lines:
  if line.startswith('|') and blocks and blocks[-1].startswith('|'):
   blocks[-1]+='\n'+line
  else:blocks.append(line)
 return '\n\n'.join(blocks)+'\n'

def mechanism_tables(cases,spatial):
 """Display measured phase retention and exact final targets, with no new gates."""
 lines=['## 相区保留、质量分布与逐块证据',
 '下表U中的强/弱标签仅为H的参考空间区域，U始终使用均匀材料参数。相区内部连接与整个游离块连接是不同问题：一个相区可以在内部断开而仍经周围材料连成同一个大块。',
 '|案例|10个参考强域中内部仍为单分量的个数|各域最大内部连接点数|存活跨相键|存活跨强域strong–strong键|',
 '|---|---:|---|---:|---:|']
 for mode in CASES:
  s=spatial[mode];domains=s['domain_retention']
  lines.append(f"|{mode}|{sum(r['final_internal_components']==1 for r in domains)}|"+', '.join(str(r['largest_internal_connected_nodes']) for r in domains)+f"|{s['remaining_live_cross_phase_bonds']}|{s['remaining_direct_strong_to_other_strong_domain_bonds']}|")
 lines+=['下表逐强域追踪原始节点最终去向；目标档和大块档依据整个真实游离分量的节点数，而不是域内节点数。',
 '|案例|强域|仍附着节点|已游离节点|位于20–200点分量的节点|位于>200点分量的节点|位于1–2点分量的节点|','|---|---:|---:|---:|---:|---:|---:|']
 for mode in CASES:
  for r in spatial[mode]['domain_retention']:
   lines.append(f"|{mode}|{r['strong_domain_id']}|{r['attached_nodes']}|{r['free_nodes']}|{r['nodes_in_target_components']}|{r['nodes_in_gt200_components']}|{r['nodes_in_1_2_components']}|")
 lines+=['|案例|材料对／参考标签对|自然完全失效／初始可失效键|仍存活键|平均损伤|','|---|---|---:|---:|---:|']
 for mode in CASES:
  for r in spatial[mode]['phase_damage']:
   lines.append(f"|{mode}|{r['kind']}|{r['failed']}/{r['eligible_bonds']}|{r['live']}|{r['mean_damage']:.6f}|")
 lines+=['|案例|节点数区间|游离分量数|节点数|占初始材料质量|占游离材料质量|','|---|---|---:|---:|---:|---:|']
 with (O/'FRAGMENT_SIZE_DISTRIBUTION.csv').open() as f:
  for r in csv.DictReader(f):
   detached=f"{float(r['fraction_detached']):.5%}" if r['fraction_detached'] else 'UNDEFINED_NO_DETACHED_MASS'
   lines.append(f"|{r['case']}|{r['bin']}|{r['count']}|{r['nodes']}|{float(r['fraction_initial']):.5%}|{detached}|")
 support=json.loads((O/'TARGET_SURFACE_SUPPORT.json').read_text())
 lines+=['以下仅列最后接受态仍存在的20–200点目标分量；全部历史目标块及其原始成员ID见TARGET_FRAGMENT_EVIDENCE.csv。参考跨度使用原始坐标，当前跨度使用实际变形坐标；点云三维秩本身不认证连续实体形貌。',
 '|案例|原始成员集合ID|节点|内部存活键|参考X/Y/Z跨度 µm|当前X/Y/Z跨度 µm|受支持加载面|表面节点覆盖|','|---|---:|---:|---:|---|---|---:|---:|']
 count=0
 with (O/'TARGET_FRAGMENT_EVIDENCE.csv').open() as f:
  for r in csv.DictReader(f):
   if r.get('survives_to_last_accepted')!='True':continue
   mode=r['case'];key='final' if cases[mode]['native']['status']=='COMPLETED' else 'last_accepted';cid=int(r['component_id'])
   ss=next(v for v in support[mode][key] if v['component_id']==cid)
   spans=[' / '.join(f'{v*1e6:.1f}' for v in json.loads(r[k])) for k in ['reference_span_m','last_actual_span_m']]
   lines.append(f"|{mode}|{cid}|{r['nodes']}|{r['last_observed_internal_alive_bonds']}|{spans[0]}|{spans[1]}|{ss['active_supported_faces']}|{ss['node_fraction_on_active_faces']:.2%}|");count+=1
 if not count:lines.append('TARGET_FRAGMENT_NOT_OBSERVED')
 lines.append('受支持加载面数来自沿用的当前四面体表面代理。节点覆盖率不是连续表面积的完备性估计；零激活面需要结合辅助四面体支持、开放/可达性和遮挡解释，不能当成真实零阻力的证据。每个目标块完整原始边界的存活键数、自然失效和初始预切分别保留在GRAPH_SEPARATION_BOUNDARY_AUDIT.json。')
 unload=json.loads((O/'UNLOAD_TARGET_LINEAGE.json').read_text())
 lines+=['逐块检查2.8 ms驱动关闭时已存在的目标块；后代由原始成员集合确定，质量保留率没有额外人为准入阈值。',
 '|案例|撤载时成员集合ID|撤载时节点数|末态最大后代节点数|最大后代质量保留率|末态目标档后代个数|原始成员完全不变|','|---|---:|---:|---:|---:|---:|---|']
 for mode in CASES:
  for r in unload[mode]['targets']:
   lines.append(f"|{mode}|{r['drive_off_component_id']}|{r['drive_off_nodes']}|{r['largest_descendant_nodes']}|{r['largest_descendant_mass_fraction']:.5%}|{r['target_descendant_count']}|{r['unchanged_member_set_survives']}|")
 splits=json.loads((O/'COHESIVE_SECONDARY_SPLITS.json').read_text())
 lines+=['以下单列至少产生两个≥20节点子块的原生二次分裂，以区分凝聚块之间的分裂和单点掉落。该描述性筛选包含>200节点大块，不改变20–200目标定义；完整事件保留所有小后代和原始成员ID。',
 '|案例|实际时间 ms|父块ID／节点数|≥20节点子块：ID／节点／出生时内部存活键|其余小后代节点总数|','|---|---:|---|---|---:|']
 for mode in CASES:
  for r in splits[mode]['events']:
   children='; '.join(f"{c['component_id']} / {c['nodes']} / {c['internal_alive_bonds_at_birth']}" for c in r['children'] if c['nodes']>=20)
   lines.append(f"|{mode}|{r['time_s']*1e3:.8f}|{r['parent_component_id']} / {r['parent_nodes']}|{children}|{sum(c['nodes'] for c in r['children'] if c['nodes']<20)}|")
 if not any(r['event_count'] for r in splits.values()):lines.append('NOT_OBSERVED：未观察到上述描述性事件。')
 ending=json.loads((O/'END_STATE_CASE_COMPARISON_READBACK.json').read_text())
 lines+=['末态全部≥20节点游离块的组成如下。大块可能是单强域核心、跨域复合体或纯弱区材料，不能把>200节点一概等同于完整强域释放。U仍只使用参考空间标签。',
 '|案例|成员集合ID|节点数|强区标签节点|弱区标签节点|逐强域节点数|内部存活键|激活加载面|','|---|---:|---:|---:|---:|---|---:|---:|']
 for mode in CASES:
  for r in ending[mode]['components']:
   if r['nodes']<20:continue
   lines.append(f"|{mode}|{r['component_id']}|{r['nodes']}|{r['strong_nodes']}|{r['nodes']-r['strong_nodes']}|{r['strong_domain_node_counts']}|{r['internal_alive_bonds']}|{r['active_facets']}|")
 return lines

def main():
 data=json.loads((O/'COMPARISON.json').read_text());q=data['cases'];sp=data['spatial'];p0=json.loads((O/'PREFLIGHT_STRUCTURE_AUDIT.json').read_text());path=json.loads((O/'PRESCRIBED_PATH_ENERGY.json').read_text());checks=json.loads((O/'VARIATIONAL_ENERGY_AUDIT.json').read_text());vid=json.loads((O/'VIDEO_GENERATION.json').read_text());assert vid['status']=='PASS'
 complete=all(q[m]['native']['status']=='COMPLETED' for m in CASES);h=q['biphasic'];u=q['uniform'];obj=h['objective'];maximum=h['native']['max_medium_components'];largest=max(h['final_free_sizes'],default=0);fines=obj['fine_debris_mass_fraction_detached']
 boundary_audit=json.loads((O/'GRAPH_SEPARATION_BOUNDARY_AUDIT.json').read_text());assert boundary_audit['status']=='PASS'
 large_lineage=json.loads((O/'FIRST_LARGE_COMPONENT_DESCENDANTS.json').read_text())
 if not complete:conclusion='INCOMPLETE_NATIVE_PROTOCOL_NO_FULL_WINDOW_CONCLUSION';nextstep='先定位本轮实际数值停止原因并恢复完整时程的可信求解；本轮不追加正式重跑，也不修改材料参数。'
 elif obj['full_target']:conclusion='MULTIPLE_COHESIVE_TARGET_FRAGMENTS_OBSERVED_SINGLE_REALIZATION';nextstep='优先研究强区排列与独立空间实现的重复性；本轮不能宣称实验标定完成。'
 elif maximum>=5:conclusion='TRANSIENT_TARGET_MULTIFRAGMENTATION_NOT_RETAINED';nextstep='按原始ID追踪目标块出生后的二次碎裂、残余动能及卸载应力重分布，优先验证凝聚块形成后的剪切破坏和混合模式稳定性；不重新调大强弱阈值差异。'
 elif fines is not None and fines>.5 and largest<20:conclusion='FINE_DEBRIS_DOMINATED_TARGET_NOT_ACHIEVED';nextstep='优先检验剪切主导损伤及三维材料分离机制，使用独立混合模式路径验证候选断裂律；不再通过扩大阈值比或更换区域排列追逐目标块。'
 elif sum(size>200 for size in h['final_free_sizes'])>=5:conclusion='MULTIPLE_COHESIVE_BLOCKS_OUTSIDE_SPECIFIED_TARGET_BAND';nextstep='按末态逐块组成区分单强域保留、多强域桥接及纯弱区凝聚体，再核对270–271点初始强域与20–200点目标之间的尺度差异、弱间隔相对horizon的非局部连接及dx可分辨性；不能把所有大块一概解释为完整强域脱落。当前块数本身不足以支持优先更换断裂律或扩大阈值比，候选剪切/混合模式修改应先经独立加载路径验证。'
 elif largest>200:conclusion='TARGET_COUNT_BELOW_FIVE_WITH_LARGE_COHESIVE_COMPONENTS';nextstep='优先区分单强域核心释放、多个强域仍被连接的大块及主分离路径，量化剩余跨相/跨强域连接与实际剪切状态，再用独立混合模式路径检验候选损伤律；仅有存活桥接计数不能证明其因果作用，不立即提高参数比。'
 else:conclusion='BIPHASIC_TARGET_COUNT_BELOW_FIVE';nextstep='优先量化非局部跨区连接和内部剪切分离状态，用独立混合模式断裂数据检验物理律；剩余跨区键计数不能单独证明其因果作用。'
 loading_representation={}
 for mode in CASES:
  key='final' if q[mode]['native']['status']=='COMPLETED' else 'last_accepted'
  targets=q[mode]['surface']['target_component_surface_support'][key]
  loading_representation[mode]=dict(target_count=len(targets),zero_active_facet_targets=[r['component_id'] for r in targets if r['active_supported_faces']==0],zero_supported_tetra_targets=[r['component_id'] for r in targets if r['supported_current_tetra']==0],active_facet_counts=[r['active_supported_faces'] for r in targets],surface_node_coverage=[r['node_fraction_on_active_faces'] for r in targets],interpretation='Connectivity remains an exact original-bond-graph result. Supported-facet coverage does not certify complete physical fluid loading or hydrodynamic survival.')
 if complete and any(r['zero_active_facet_targets'] for r in loading_representation.values()):
  nextstep='优先验证并修复存活非局部材料图与当前表面加载表示之间的一致性：末态存在没有激活加载面的凝聚目标块，需区分辅助四面体支持缺失、开放/可达性未满足及遮挡，并用独立的脱落凝聚块受力路径验证；当前结果不能认证完整流体阻力下的卸载稳定性。在此之后，'+nextstep
 status=dict(STRUCTURE_PREFLIGHT_STATUS=p0['status'],STRONG_DOMAIN_COUNT=p0['strong_domain_count'],STRONG_VOLUME_FRACTION=p0['strong_volume_fraction'],WEAK_S0=.1,WEAK_SC=.12,STRONG_S0=.25,STRONG_SC=.30,CROSS_PHASE_S0=.175,CROSS_PHASE_SC=.210,ORIGINAL_TOPOLOGY_PRESERVED=True,VARIATIONAL_FORCE_PARITY='PASS',IRREVERSIBLE_HISTORY_STATUS='PASS',SHEAR_DRIVE_IDENTICAL=True,UNIFORM_NATIVE_STATUS=u['native']['status'],BIPHASIC_NATIVE_STATUS=h['native']['status'],FIRST_SOFTENING_LOCATION={m:sp[m]['first_events']['softening'] for m in CASES},FIRST_FAILURE_LOCATION={m:sp[m]['first_events']['failure'] for m in CASES},MULTISITE_CRACK_INITIATION={m:dict(separated_early_softening_sites=len(sp[m]['early_spatially_separated_softening_sites']),scope='Spatial event evidence only, independent continuum cracks not certified') for m in CASES},TARGET_FRAGMENT_COUNT_20_200={m:q[m]['objective']['n_fragments_20_200'] for m in CASES},TARGET_FRAGMENT_MASS_FRACTION={m:q[m]['objective']['mass_fraction_target_initial'] for m in CASES},LARGEST_FREE_COMPONENT={m:max(q[m]['final_free_sizes'],default=0) for m in CASES},FINE_DEBRIS_FRACTION={m:q[m]['objective']['fine_debris_mass_fraction_detached'] for m in CASES},SECONDARY_FRAGMENTATION={m:q[m]['secondary_count'] for m in CASES},ENERGY_BALANCE_STATUS={m:dict(status='PASS_ACCEPTED_STATES',maximum_relative_R=q[m]['native']['max_relative_energy_residual']) for m in CASES},VISUALIZATION_STATUS='ACTUAL_ACCEPTED_STATE_PNG_AND_VECTOR_PDF_GENERATED',VIDEO_STATUS='TWO_1080P_REAL_ACCEPTED_STATE_MP4_GENERATED_PENDING_FINAL_DECODE',SCIENTIFIC_CONCLUSION=conclusion,NEXT_PHYSICAL_PRIORITY=nextstep,full_3p6ms_comparison_complete=complete,historical_maximum_target_count={m:q[m]['native']['max_medium_components'] for m in CASES},formal_runs=2,material_label=HET_LABEL,load_label=LABEL)
 status['TARGET_LOADING_REPRESENTATION']=loading_representation
 status['TOTAL_DETACHED_MASS_FRACTION']={m:q[m]['objective']['total_detached_mass_fraction'] for m in CASES}
 status['ATTACHED_BODY_FRACTION']={m:q[m]['objective']['attached_body_fraction'] for m in CASES}
 dump(O/'FINAL_STATUS.json',status)
 lines=[f'# {TASK.upper()}',f'运行标识：`{RUN}`。科研状态：**{conclusion}**。',
 '## 实施、冻结与必要检查',
 f"新材料包络 1.2 × 1.2 × 0.6 mm，dx=50 µm，horizon=150.5 µm。{p0['material_nodes']} 材料节点 + {p0['anchor_nodes']} 固定锚点；{p0['original_bulk_bonds']} bulk、{p0['original_interface_bonds']} interface、{p0['anchor_scaffold_bonds']} anchor scaffold 原始键。共同单边缺口 {p0['initial_precrack']} 键；新奇数网格将 z=0 点归入下侧，以 z>0 的跨侧线段交点 x≤−50 µm 定义缺口。无按材料相区删键。",
 f"P0 已先于材料开发与动力学展示。seed=20261018；受约束面邻接区域生长得到 {p0['strong_domain_count']} 个强区，每区270–271点，强相参考体积分数 {p0['strong_volume_fraction']:.8%}；弱基质面邻接分量为1。区域坐标中心跨度250–350 µm；等参考体积球直径约400.95–401.45 µm。后者由 V=dx³ 计算，不能与中心包络跨度混为一谈。离散节点体积总和1.015625 mm³，中心外包络0.864 mm³：沿用端节点完整V的既有积分约定，没有以包络体积重新缩放质量。",
 f"强区最小节点中心距离 {p0['minimum_interdomain_node_distance_m']*1e6:.1f} µm，小于horizon。跨相原始bulk键 {p0['cross_phase_material_bonds']}，占bulk {p0['cross_phase_fraction_of_bulk']:.5%}；其中初始缺口外存活 {p0['cross_phase_initially_live']}。强区之间存在普通弱节点及跨越弱间隔的非局部连接，强区标签不是预设碎块；每区270–271点本身高于目标200点上限，完整强域脱落必须计入>200档，不能称为目标块。区域生长使用格点单元约束，形态有离散/单元边界痕迹，未声称来自组织图像或统计独立组织块。",
 'WEAK s0/sc=.100/.120，MIXED=.175/.210，STRONG=.250/.300；UNIFORM=.15/.18。强/弱起裂、完成和软化区间比均为2.5，无耗能归一化。混合键使用端点算术平均，端点互换不变。G=940 Pa、K=100000/3 Pa（33333.33 Pa为显示舍入）、rho=1000 kg/m³；interface=.140625/.16875；anchor scaffold=1/1.2。',
 f"节点体积均值 s0/sc={p0['mean_node_s0']:.9f}/{p0['mean_node_sc']:.9f}；原始bulk参考权重 h*r0*V²*(1/mi+1/mj) 均值={p0['weighted_bulk_s0']:.9f}/{p0['weighted_bulk_sc']:.9f}。均值接近均匀卡不意味着整体有效断裂能相等。",
 '完整endpoint-star最小化串联柔度储能 U(x,d)、−∂U/∂x、原始m/r0/taper、节点及键ID、接触和Velocity-Verlet保持。原生仅增加同时冻结 s0/sc 的原始ID接口，并将 reference/restart 的旧4096节点读取上限放宽；没有重写力核。均匀旧实现比较库只修改两处容量读取上限，其他计算代码逐字保持。原production库未覆盖。restart V4已包含两阈值字段，沿用原格式。',
 'U/H全部初始κ/d一致：预切bulk使用共同 κ=.30、d=1；其余κ=d=0。这是失效历史标记，超过所有本轮bulk完成阈值，不代表物理初应变或初始损伤耗能；预切材料能量/力为零。',
 f"Test A：新几何均匀旧/新力、完整储能、κ/d及原始ID逐位一致。Test B：六个代表性状态（intact、weak near-onset/partial、mixed near-onset/partial、strong partial），包含真实跨相端点；最大相对固定损伤能量梯度误差 {max(r['relative_gradient_error'] for r in checks['checks']['variational']):.6g}。Test C：10000个节点数组槽位改变后属性仍按ID保持，restart两阈值与历史保持，卸载不恢复，试算不提交。仅两个新增非正式单步动力学检查，用于当前加载/拖曳和原生接口，不属于额外全时程扫描。",
 'Test D：同一新几何完整endpoint-star模型，以 z=0 两侧节点作 ±opening/2 刚性位移（原始预切保留），0–180 µm、0.5 µm间隔；四组分别将所有bulk设为对应卡，interface不变。反力是完整模型对该位移路径的共轭反力。D为每个路径点固定x下的完整能量下降，非原生剪切时程、非Gc标定。',
 '|路径卡|初始刚度 N/m|路径D nJ|积分残差 nJ|完成失效跨面键|','|---|---:|---:|---:|---:|']
 for m,r in path['summaries'].items():lines.append(f"|{m}|{r['initial_stiffness_N_m']:.7g}|{r['path_dissipation_J']*1e9:.7g}|{r['final_quadrature_residual_J']*1e9:.7g}|{r['complete_crossing_bond_count']}/{r['crossing_bonds']}|")
 lines.extend([f"路径STRONG/WEAK耗能比 {path['strong_weak_path_dissipation_ratio']:.7g}，并非2.5。表中反力梯形积分残差保留；路径采样有限，未宣称连续极限或误差小于原生动力学1%准入阈值。",'## 同一加载与两组正式原生结果',
 'UNIFORM_SURFACE_SHEAR_X：+X投影 P_t eX，不归一化投影、不加空间窗。tau_nominal=1000g(t) Pa；g为指定三次smoothstep，0–2ms升、2–2.4ms恒、2.4–2.8ms降、2.8–3.6ms驱动为零但被动拖曳保留。牵引=β[(1000/β)gP_t eX−P_t v]，β=μ/h=40 Pa·s/m。相同外部驱动不强制相同实际外功或载荷面积。',
 '继承当前受支持四面体边界/裂面可达性/遮挡点探测代理，辅助四面体不引入材料刚度和质量。峰值假定外部速度25 m/s、Re_h上界625；这是高幅值机械剪切需求代理，不能换算为已验证微泡声流功率，也不是低Re流场求解。',
 '仅一次U和一次H，顺序执行，无材料/幅值/随机种子扫描。原有max Δd≤0.1事件细分和相对能量余额≤1%准入保持；没有新增内存或墙钟预算。',
 '|案例|状态|实际终时 ms|接受步|拒绝预测|目标块末态/历史最大|最大自由块|细屑/脱落质量|','|---|---|---:|---:|---:|---:|---:|---:|'])
 for m in CASES:
  r=q[m];s=r['native'];o=r['objective'];fraction=o['fine_debris_mass_fraction_detached'];lines.append(f"|{m}|{s['status']}|{s['time_s']*1e3:.8f}|{s['accepted_steps']}|{s['predictive_rejections']}|{o['n_fragments_20_200']}/{s['max_medium_components']}|{max(r['final_free_sizes'],default=0)}|{fraction if fraction is None else format(fraction,'.7%')}|")
  lines.append(f"{m}：CPU={s['cpu_s']:.3f}s，wall={s['wall_s']:.3f}s，peakRSS={s['peak_RSS_bytes']/2**20:.3f}MiB，状态文件={s['saved_state_count']}。实际驱动功={r['work']['drive_J']*1e9:.8g}nJ，拖曳={r['work']['drag_J']*1e9:.8g}nJ，净功={r['work']['actual_net_J']*1e9:.8g}nJ；最大相对余额={s['max_relative_energy_residual']:.7%}。停止原因={s['stop_reasons']}。")
  if m=='biphasic' and (O/'USER_RESUME_STATE.json').exists():
   resume=json.loads((O/'USER_RESUME_STATE.json').read_text());paused=resume['pause_record_to_resume_elapsed_s']
   lines.append(f"H遵用户指令SIGSTOP暂停后从同一进程SIGCONT恢复，未新增正式尝试。原始wall包含暂停；暂停记录至恢复间隔{paused:.3f}s，扣除此已记录间隔后约{s['wall_s']-paused:.3f}s。暂停记录略晚于实际SIGSTOP，差值不是精确纯计算墙钟；CPU不含暂停消耗，成本对比图采用CPU。暂停与恢复身份记录见USER_PAUSE_STATE.json及USER_RESUME_STATE.json。")
  ue=r['unload_energy']
  if ue is not None:lines.append(f"{m}已观测卸载窗口{ue['observed_start_s']*1e3:.6f}–{ue['observed_end_s']*1e3:.6f} ms：动能从{ue['start_kinetic_J']*1e9:.7g}到{ue['end_kinetic_J']*1e9:.7g} nJ，新增损伤耗能{ue['damage_increment_J']*1e9:.7g} nJ，驱动/拖曳增量{ue['drive_work_increment_J']*1e9:.7g}/{ue['drag_work_increment_J']*1e9:.7g} nJ，接触创建能量增量{ue['contact_creation_increment_J']*1e9:.7g} nJ。这是仍可有残余动能和继续断裂的有限动态观察窗口，不是静力平衡验证。")
  first_free=r['first_free_release'];lines.append(f"{m}首次真正游离："+(f"{first_free['birth_s']*1e3:.8f} ms，{first_free['nodes']}节点，component_id={first_free['component_id']}" if first_free else 'NOT_OBSERVED'))
  release=r['first_large_release'];release_text=f"{release['birth_s']*1e3:.8f} ms，{release['nodes']}节点，原始成员集合component_id={release['component_id']}" if release else 'NOT_OBSERVED'
  lines.append(f"{m}首次>200点真实游离分量：{release_text}。末态总游离质量/初始={o['total_detached_mass_fraction']:.7%}，仍附着比例={o['attached_body_fraction']:.7%}；最大分量/全部游离质量={o['largest_free_fragment_fraction_detached']}；末态目标块总质量={o['mass_target_kg']:.9g} kg。细屑比例以全部游离材料为分母；大块释放可能使该比例下降，不能把分母增加误读为细屑数量减少。")
  ancestor=large_lineage[m]
  if ancestor['status']!='NOT_OBSERVED':
   lines.append(f"{m}首次大块的原始{ancestor['parent_nodes']}个节点在末态分布于{len(ancestor['final_descendants'])}个真实后代，最大后代{ancestor['largest_descendant_nodes']}点，保留原始大块质量的{ancestor['largest_descendant_fraction_of_parent']:.5%}，其中目标区间后代{ancestor['target_descendant_count']}个。逐后代原始ID、内部存活键和相区组成见FIRST_LARGE_COMPONENT_DESCENDANTS.json；大于200点的凝聚体仍是真实游离块，只是不计入指定目标区间。")
 for filename,label in [('PEAK_LOAD_CASE_COMPARISON_READBACK.json','2.4 ms恒载结束'),('TIME_2P550_MS_NEAREST_CASE_COMPARISON_READBACK.json','2.55 ms附近的最近真实接受态'),('UPPER_DRIVE_OFF_CASE_COMPARISON_READBACK.json','2.8 ms驱动关闭瞬间')]:
  stage_path=O/filename
  if not stage_path.exists():continue
  stage=json.loads(stage_path.read_text());lines.append(label+'的阶段证据（各组实际时间分别列出，不插值；独立于末态结论，U相区标签仍只是空间协变量）：')
  for mode in CASES:
   r=stage[mode];c=r['components'][0];sm=r['metrics']
   lines.append(f"{mode}实际t={r['time_s']*1e3:.6f} ms，最大的七个游离分量={r['largest_free_sizes'][:7]}；最大块{c['nodes']}节点、{c['internal_alive_bonds']}条存活内部键，其中强区节点{c['strong_nodes']}，逐强域点数={c['strong_domain_node_counts']}，存活跨相键{c['live_internal_cross_phase_bonds']}、直接跨强域strong–strong键{c['live_internal_interdomain_strong_bonds']}，激活加载面{c['active_facets']}个，完整原始跨块边界无存活键。总游离质量占初始{sm['total_detached_mass_fraction']:.5%}，最大块及细屑分别占游离质量{sm['largest_free_fragment_fraction_detached']:.5%}/{sm['fine_debris_mass_fraction_detached']:.5%}。多个强域仍归于同一个材料分量时，不能按相区标签拆成多个碎块。源状态哈希见{filename}。")
 lines+=['能量账本 Wdrive+Wdrag=ΔU+ΔK+D+ΔUc−Jc+R；R始终单独记录，不并入D。每个接受步的标量账本、连接分量及损伤增量均保留；完整FP64状态按规定间隔和关键事件保存，数量在上表逐例列出。case/final.restart是终止时原生状态，若出现数值拒绝则不能将它当作last_accepted重启点，应同时保留failed_unaccepted与last_accepted证据且不得直接续算被拒绝状态。若未完成3.6ms，视频明确标注INCOMPLETE，最后接受态不冒充正式最终态。',
 '## 最终五项科学回答']
 sites=sp['biphasic']['early_spatially_separated_softening_sites'];first=sp['biphasic']['first_events']
 bulk_onsets={m:sp[m]['first_events']['bulk_softening'] for m in CASES}
 lines.append('区分最早全材料事件与最早bulk事件：'+ '；'.join(f"{m}首次bulk软化t={bulk_onsets[m]['time_s']}s，原始键ID={bulk_onsets[m]['original_bond_ids']}" for m in CASES)+f"。两组首次bulk原始键集合相同={bulk_onsets['uniform']['original_bond_ids']==bulk_onsets['biphasic']['original_bond_ids']}；全材料首次事件可因interface与bulk的先后关系而改变。多处软化存在于两组，不能表述成均匀材料没有多处起始位置。")
 first_count_ge5_text=(str(obj['first_count_ge5_time_s'])+' s') if obj['first_count_ge5_time_s'] is not None else 'NOT_OBSERVED'
 off_summary=(f"驱动关闭时目标块{obj['drive_off_target_count']}个，其中到最后接受态保持完全相同原始成员集合的{obj['drive_off_target_survivors_to_end']}个，仍处于目标区间的后代共{obj['drive_off_target_descendant_targets_at_end']}个。" if obj['drive_off_time_reached'] else '尚未到达2.8 ms驱动关闭时刻，卸载存活未观测，不能以零计数推断未发生。')
 lines+=[f"1. **是否出现多个自然裂纹起源？** H在首次bulk软化后0.1ms内，统一采用参考中点距离>2horizon=301µm的贪心分离规则，得到{len(sites)}个空间软化位置；U对应{len(sp['uniform']['early_spatially_separated_softening_sites'])}个。首次软化原始键集合是否相同：{sp['difference']['first_softening_same_original_bonds']}。这是多个空间损伤起始位置的证据，不能据此认证多条独立连续裂纹。H首次完整失效的材料对：{list(first['failure']['phase_pairs'])}。",
 '上述0.1 ms窗口分别从各自首次bulk软化开始，发生在不同绝对时刻和不同名义驱动水平；它们用于描述早期空间起始位置，不能作为同一牵引下独立连续裂纹数量的比较。',
 f"2. **是否贯通三维材料分离？** H末态共有{len(h['final_free_sizes'])}个与锚层不可达的真实材料分量。每个分量均通过全部原始存活bulk/interface图判定，故其与固定主体之间在该离散图上确已分离；连续三维裂面是否贯通尚不能仅由点图认证。强区内部连通与剩余材料桥逐区见STRONG_DOMAIN_RETENTION.csv，H剩余跨相存活键总数{sp['biphasic']['remaining_live_cross_phase_bonds']}，直接连接不同强域的strong–strong存活键为{sp['biphasic']['remaining_direct_strong_to_other_strong_domain_bonds']}。存活连接的当前几何伸长、历史κ/s0及尚未起裂的数量另见SPATIAL_DAMAGE_ANALYSIS.json；这些是几何和历史证据，并非逐键力测量。另对末态及保存的最大目标计数态中的全部目标块和最大游离块核验完整原始边界：跨边界存活键必须为零，原有缺口与自然失效分别统计，见GRAPH_SEPARATION_BOUNDARY_AUDIT.json。未做最小割优化，也未将failed-bond网络当成自由块。",
 f"3. **至少5个20–200点凝聚块是否实现？** H末态目标块{obj['n_fragments_20_200']}个、全程最大{maximum}个，目标质量占初始材料{obj['mass_fraction_target_initial']:.7%}；{off_summary}首次计数达到5块的时间={first_count_ge5_text}，已观测驱动关闭窗口内的最低目标计数={obj['minimum_target_count_during_observed_driveoff']}。UNLOAD_TARGET_LINEAGE.json分别记录成员不变、子块数和最大后代质量保留率，不把少量掉点自动等同于凝聚主体完全丧失。完成时程后的至少5块判定（卸载寿命另列，不附加隐藏计数门槛）：{obj['full_target']}。每个目标块的内部存活键、原始节点、出生/最后观测、卸载存活和点云形状见TARGET_FRAGMENT_EVIDENCE.csv；每个接受态核验其内部存活键数至少为节点数−1。同时比较参考形态与当前形态：弯曲细条的当前点云可能有三维秩，不能据此把原始细条/薄片解释为紧实块体；节点数准入不自动保证块状形貌。点云为链/片状时不得冒称已解析实体三维碎块。不能按相区或空间邻近合并节点。",
 f"4. **与同几何均匀材料相比？** U/H末态目标块={u['objective']['n_fragments_20_200']}/{obj['n_fragments_20_200']}，最大自由块={max(u['final_free_sizes'],default=0)}/{largest}节点，细屑质量比={u['objective']['fine_debris_mass_fraction_detached']}/{fines}。完整7档尺寸及质量分布见CSV与FIG_08；每个接受态记录14–19与20–33的独立区间，20和200均计入目标。U/H二次分裂事件数={u['secondary_count']}/{h['secondary_count']}，具体parent/children和birth/寿命由原始成员集合追踪。这里只能说明这一冻结空间实现的结果，不能宣称统计普适性。",
 f"5. **下一项最有证据支持的修改？** {nextstep} 这些是待检验机制，不是本轮已证明的混合模式模型。",
 '## 证据边界与复现',
 'U的分区统计仅借用H的参考区域标签作空间协变量，U的全部bulk仍为均匀卡。强区内部是否保持连接、弱区是否多处损伤、剩余跨相桥等分别量化，断键数量增加不等于可独立辨认的裂纹数量；交错连续裂纹表面不能由当前点图唯一恢复。界面参数、接触精度、整体强度、Gc均未实验标定。参考节点区域体积分数与组织材料真实相比例不可直接等同。',
 '目标尺寸块也逐块检查受支持加载面数、面积、表面节点覆盖与实际表面功率，见TARGET_SURFACE_SUPPORT.json；不能仅因节点数达到20–200就假设其表面载荷已解析。对没有受支持3D辅助四面体表面的细屑/细链，当前流体加载仍未解析；必须保留覆盖率证据，不能把未表示的表面力解释为物理零阻力。原始预切重新张开与自然新裂面加载分别统计，不能把二者混称新裂纹增益。',
 'PNG及矢量PDF由冻结输入、明确标记的离线验证和真实FP64接受态生成。第二视频左侧保留完整节点与细屑；右侧只显示整个真实≥20节点游离分量，在尚无该类分量时显示附着材料相区。该显示筛选不拆分或合并材料分量，统计始终来自完整原始存活图。视频位移×1，相机与全景世界坐标范围固定，U/H各标实际保存时间；以最近接受态对齐，不插值物理轨迹。视频中3节点及以上游离分量的颜色按原始成员谱系追踪：仅节点数严格超过父块一半的子块继承父色，其余子块使用新色；精确成员集合与分量ID始终独立保存，见VIDEO_COMPONENT_COLOR_LINEAGE.json。1–2节点细屑统一用较小的淡棕色点显示，全部原始身份仍被记录。颜色延续不表示节点集合完全不变，也不改变碎块计数。图中球点仅为显示符号，不是重建组织外壳。近邻断键只作示意，碎块连通用全部原始存活键，游离块内部显示全部存活材料键。',
 '所有路径/源码/库哈希见SOURCE_FREEZE、BUILD_AUDIT和最终DELIVERY_MANIFEST；历史production和已有正式结果保持。复现命令见REPRODUCE.md，单轮正式运行数已用完，任何后续物理试验必须建立新的运行标识。']
 lines+=mechanism_tables(q,sp)
 lines.append('末态目标块加载表示诊断：'+json.dumps(loading_representation,ensure_ascii=False)+'。零激活面不能单凭零力判为物理无阻力；候选面、开放面、遮挡和支持四面体数分别见TARGET_SURFACE_SUPPORT.json。原生冻结表面算法保持，本轮没有据结果修改载荷或追加运行。')
 report=markdown_blocks(lines);name='REVIEW_CLOT_BIPHASIC_HETEROGENEOUS_SHEAR_FRAGMENTATION_V2.md';(O/name).write_text(report)
 repro=f'''# Reproduction

Project: /workspace/y_clot/project. Activate /workspace/y_clot/env.sh.

Sources: analysis/{TASK}/. Runtime/source/input hashes: verification/{TASK}/{RUN}/SOURCE_FREEZE.json.
The final isolated library is under verification/{TASK}/{RUN}/build. The capacity-only prior comparison library is under baseline_capacity. Build commands and upstream object hashes are retained.

Completed order: generate_biphasic_material.py → preflight_figures.py (P0 displayed) → build_capacity_baseline.py → build.py → uniform_probe.py old_uniform/new_uniform (separate libraries) → checks.py → path_audit.py → freeze_sources.py → run_cases.py --name uniform → run_cases.py --name biphasic → analyze.py → figures.py → render_videos.py → report.py → final_delivery.py.

Do not rerun formal cases in this output directory. The sole two permitted formal attempts are preserved. Reproduction that changes physical trajectories requires a new RUN and new authorization/scope; read-only plotting does not.

Read-only plotting on the existing server environment:
```bash
source /workspace/y_clot/env.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /workspace/y_clot/project/analysis/{TASK}
python analyze.py
python stage_comparison.py --key peak_load
python stage_comparison.py --time-ms 2.55 --nearest
python stage_comparison.py --key upper_drive_off
python preflight_figures.py
python figures.py
python render_videos.py --smoke
python render_videos.py
python report.py
```

Y_CLOT_PROJECT_ROOT overrides the root for replay in an environment containing numpy/scipy/pandas/matplotlib/PyVista/VTK/Pillow/FFmpeg. Native execution additionally needs mpi4py, the LAMMPS Python binding and hashed native build dependencies; LD_LIBRARY_PATH must select the isolated library. No whole-history dependency directory was duplicated.

Frame indices: videos/*_FRAME_TIMES.csv. Full FP64 snapshots: states/<case>/keys/*.npz. Every accepted scalar/graph and exact damage increments remain in each case. All geometry uses SI internally. Original IDs start at 1.
''';(A/'REPRODUCE.md').write_text(repro)
 blocks=['<!doctype html><html lang="zh"><meta charset="utf-8"><title>双相血栓剪切碎裂 V2</title><style>body{font:17px system-ui;max-width:1300px;margin:auto;padding:24px;line-height:1.6;background:#fafbfc;color:#253b49}img,video{width:100%;background:white}section{margin:28px 0;border:1px solid #c5cdd3;padding:18px;background:white}a{color:#175b84}pre{white-space:pre-wrap}h1{font-size:28px}</style>',f'<h1>双相血栓 · 相同均匀剪切</h1><p>{html.escape(conclusion)}</p><p>10强区 · 33.329%强相 · 1000 Pa · 单一探索性空间实现 · 未标定机械剪切代理</p>']
 blocks.append('<style>table{border-collapse:collapse;width:100%;font-variant-numeric:tabular-nums}th,td{padding:10px;border:1px solid #c5cdd3;text-align:left}th{background:#edf2f6}.scroll{overflow-x:auto}aside{padding:18px;background:#fff2e5;border-left:4px solid #a47642}</style><section><h2>两组实际结果</h2><div class="scroll"><table><tr><th>案例</th><th>状态／终时</th><th>目标块：末态／历史最大</th><th>目标质量／初始</th><th>总游离质量／初始</th><th>最大游离块</th><th>细屑／游离质量</th></tr>')
 for mode in CASES:
  r=q[mode];o=r['objective'];fine=o['fine_debris_mass_fraction_detached'];fine_text=f'{fine:.2%}' if fine is not None else '无游离材料'
  blocks.append(f"<tr><td>{'U · 均匀' if mode=='uniform' else 'H · 双相'}</td><td>{r['native']['status']} / {r['native']['time_s']*1e3:.4f} ms</td><td>{o['n_fragments_20_200']} / {r['native']['max_medium_components']}</td><td>{o['mass_fraction_target_initial']:.2%}</td><td>{o['total_detached_mass_fraction']:.2%}</td><td>{max(r['final_free_sizes'],default=0)} 节点</td><td>{fine_text}</td></tr>")
 blocks.append('</table></div><p>目标区间包含20和200节点；游离块按全部原始存活材料键及锚点可达性判定。节点数与图连通性不自动认证紧实三维形貌或完整流体受力。</p><aside>'+html.escape(nextstep)+'</aside></section>')
 for file,label in [(O/name,'完整科研报告'),(O/'FINAL_STATUS.json','研究状态'),(O/'COMPUTE_COST.json','计算成本'),(A/'REPRODUCE.md','复现说明'),(O/'TRUE_FRAGMENT_HISTORY.csv','每接受态碎块指标'),(O/'ENERGY_LEDGER.csv','完整能量账本'),(O/'GRAPH_SEPARATION_BOUNDARY_AUDIT.json','全部原始边界键分离核验'),(O/'UNLOAD_TARGET_LINEAGE.json','卸载目标块后代与质量保留'),(O/'FIRST_LARGE_COMPONENT_DESCENDANTS.json','首次大块的后代去向'),(O/'COHESIVE_SECONDARY_SPLITS.json','凝聚子块之间的二次分裂'),(O/'TARGET_SURFACE_SUPPORT.json','目标块表面载荷表示')]:blocks.append(f'<a href="{os.path.relpath(file,F)}">{label}</a> · ')
 blocks.append('<details><summary>冻结输入、材料场与逐块证据</summary><ul>')
 for filename in ['SOURCE_FREEZE.json','GEOMETRY_FREEZE.json','BIPHASIC_MATERIAL_SPEC.json','BIPHASIC_FIELD.npz','STRONG_DOMAINS.csv','NODE_MATERIAL_PHASE.csv','PAIR_MATERIAL_PROPERTIES.csv','FROZEN_SHEAR_DRIVE.json','VARIATIONAL_ENERGY_AUDIT.json','PRESCRIBED_PATH_ENERGY.json','TARGET_FRAGMENT_EVIDENCE.csv','STRONG_DOMAIN_RETENTION.csv','SPATIAL_DAMAGE_ANALYSIS.json','FRAGMENT_SIZE_DISTRIBUTION.csv','FRAGMENT_MASS_DISTRIBUTION.csv','PEAK_LOAD_CASE_COMPARISON_READBACK.json','TIME_2P550_MS_NEAREST_CASE_COMPARISON_READBACK.json','UPPER_DRIVE_OFF_CASE_COMPARISON_READBACK.json','END_STATE_CASE_COMPARISON_READBACK.json','UNLOAD_ENERGY_DIAGNOSTIC.json']:
  target=O/filename
  if target.exists():blocks.append(f'<li><a href="{os.path.relpath(target,F)}">{filename}</a></li>')
 blocks.append('</ul></details>')
 for p in sorted((F/'videos').glob('*.mp4')):blocks.append(f'<section><h2>{p.stem}</h2><video controls preload="metadata" poster="videos/{p.stem}_COVER.png" src="videos/{p.stem}.mp4"></video><a href="videos/{p.stem}_FRAME_TIMES.csv">真实接受态时间索引</a></section>')
 for p in sorted((F/'figures').glob('*.png')):blocks.append(f'<section><h2>{p.stem}</h2><a href="figures/{p.stem}.pdf">矢量PDF</a><img loading="lazy" src="figures/{p.name}"></section>')
 blocks.append('<section><h2>科学解释与限制</h2><pre>'+html.escape(report)+'</pre></section></html>');(F/'OPEN_RESULTS.html').write_text('\n'.join(blocks));(O/'OPEN_RESULTS.html').write_text('<!doctype html><meta charset="utf-8"><a href="'+os.path.relpath(F/'OPEN_RESULTS.html',O)+'">本轮完整科研图、视频与证据入口</a>')
 print('REPORT_DONE',conclusion,flush=True)
if __name__=='__main__':main()
