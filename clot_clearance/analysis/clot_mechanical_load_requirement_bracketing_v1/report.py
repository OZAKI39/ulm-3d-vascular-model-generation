from common import *
from prescribed_traction import MODES,END
from analyze import collect,read
import re

def main():
 cases,rows,brackets=collect();resources=read(O/'RUN_RESOURCES.json');done=[z for z in cases if z['status']=='COMPLETED'];failed=[z for z in cases if z['status']!='COMPLETED'];valid_keys={(z['params']['mode'],z['params']['peak_pa']) for z in done};unresolved=[z for z in failed if (z['params']['mode'],z['params']['peak_pa']) not in valid_keys]
 def bracket(key):return {q['load_mode']:dict(no_event_lower_pa=q[key+'_negative_lower_pa'],observed_event_upper_pa=q[key+'_positive_upper_pa'],scope='Same 2.5 ms frozen protocol; finite sample bracket, not universal critical Pa') for q in brackets}
 status=dict(MECHANICAL_BRACKETING_STATUS='COMPLETED_LIMITED_ADAPTIVE_NATIVE_CAMPAIGN',SOFTENING_ONSET_RANGE=bracket('softening'),FIRST_FAILURE_RANGE=bracket('failure'),COHESIVE_CHUNK_RANGE=bracket('multinode'),FRAGMENTATION_REGIME={m:[dict(peak_pa=z['params']['peak_pa'],classification=z['classification'],history_max_nodes=z['history_max_free_nodes']) for z in cases if z['params']['mode']==m and (z['status']=='COMPLETED' or (z['params']['mode'],z['params']['peak_pa']) not in valid_keys)] for m in MODES},NUMERICAL_ADMISSION_STATUS=dict(completed=len(done),failed_attempts=len(failed),recovered_attempts=len(failed)-len(unresolved),unresolved_conditions=len(unresolved),energy_gate=.01,unresolved_cases=[z['tag'] for z in unresolved]),PHYSICAL_REALIZABILITY_STATUS='MECHANICAL_REQUIREMENT_ONLY; no validated acoustic streaming, jet or direct bubble feasibility',NEXT_RECOMMENDED_MECHANISM='TO_BE_SELECTED_FROM_ACTUAL_DATA')
 maxcase=max(cases,key=lambda z:z['history_max_free_nodes']);pull=[z for z in done if z['params']['mode']=='IDEAL_LOCAL_PULL'];imp=[z for z in done if z['params']['mode']=='LOCAL_IMPINGEMENT']
 if any(z['history_max_free_nodes']>=34 for z in pull):nextstep='LOCALIZED_OPENING_OR_DIRECT_BUBBLE_CLOT_TRACTION_WITH_PHYSICAL_FORCE_AND_WORK_BUDGET'
 elif any(z['history_max_free_nodes']>=34 for z in imp):nextstep='LOCAL_IMPINGEMENT_FORCE_AND_WORK_REALIZABILITY'
 elif any(z['history_max_free_nodes']>=34 and z['params']['mode']=='SURFACE_SHEAR' for z in done):nextstep='LOCALIZED_HIGH_INTENSITY_SHEAR_FROM_RESOLVED_JET_OR_DIRECT_BUBBLE_ACTION_WITH_FORCE_WORK_BUDGET'
 else:nextstep='LOCALIZED_TRACTION_GEOMETRY_AND_FRACTURE_PATH_DIAGNOSIS_BEFORE_FULL_V3_OR_BUBBLE_OPTIMIZATION'
 status['NEXT_RECOMMENDED_MECHANISM']='DEFORMATION_AWARE_TRACTION_AND_POST_DETACHMENT_UNLOADING_CONTROL_BEFORE_'+nextstep
 status['COHESIVE_CHUNK_RANGE']={'multinode_ge2':bracket('multinode'),'observed_amplitudes_ge14_pa':{m:sorted(set(r['peak_traction_pa'] for r in rows if r['load_mode']==m and r['max_detached_nodes']>=14)) for m in MODES},'observed_amplitudes_ge34_pa':{m:sorted(set(r['peak_traction_pa'] for r in rows if r['load_mode']==m and r['max_detached_nodes']>=34)) for m in MODES},'observed_amplitudes_ge114_pa':{m:sorted(set(r['peak_traction_pa'] for r in rows if r['load_mode']==m and r['max_detached_nodes']>=114)) for m in MODES},'scope':'2-node components are tracked but count as unresolved debris, not macroscopic cohesive chunks'}
 dump(O/'FINAL_STATUS.json',status)
 def range_text(key):
  names=dict(SURFACE_SHEAR='表面剪切',LOCAL_IMPINGEMENT='局部压入',IDEAL_LOCAL_PULL='理想拉脱');parts=[]
  for q in brackets:
   lo=q[key+'_negative_lower_pa'];hi=q[key+'_positive_upper_pa']
   value=f'({lo:.6g}, {hi:.6g}] Pa' if lo is not None and hi is not None else f'已在 {hi:.6g} Pa 观察到，未测得无事件下界' if hi is not None else f'完整观测窗至 {lo:.6g} Pa 未观察到事件' if lo is not None else '未识别'
   parts.append(names[q['load_mode']]+'：'+value)
  return '；'.join(parts)
 def mode_summary(m):
  z=[q for q in brackets if q['load_mode']==m][0]
  return f"{m}: 无软化完整窗 {z['no_damage_observed_pa']} Pa；观察到软化 {z['softening_observed_pa']} Pa；完整断裂 {z['failure_observed_pa']} Pa；多节点游离 {z['multinode_observed_pa']} Pa；二次裂解 {z['secondary_observed_pa']} Pa；数值未完成 {z['unresolved_pa']} Pa。"
 s=['# CLOT_MECHANICAL_LOAD_REQUIREMENT_BRACKETING_V1\n',f'运行标识：{RUN}。规定牵引机械需求试验，所有工况标记 `PRESCRIBED_TRACTION_MECHANICAL_REQUIREMENT` / `MECHANICAL_REQUIREMENT_ONLY`。\n','## 冻结基线与定义\n','同一1792材料点、336原始底部固定锚点、95666参考键、2736原预切键。材料与库SHA见 BASELINE_AND_SPECIMEN_FREEZE.json。bulk G=940 Pa、K=33333.333333 Pa、s0=0.15、sc=0.18；interface 0.140625/0.16875；固定锚连接1.0/1.2。dx=50 µm、δ=150.5 µm、ρ=1000 kg/m³。没有新预切、预设目标块或更换材料。\n','3080模型是双端规定移动、无独立锚点/界面卡的另一试样；其完整2.5 ms轨迹仅有主体3076点及4单点，并不是已证实能自然释放大块的锚定试样。此处选用最近的已有1792+336替代体，同时用于用途A与B。历史位移开裂试验曾形成最终1078点自由分量，但其顶夹具释放与本轮始终自由的受载顶点不同，不作为本轮牵引结果。\n','载荷时程明确适配为参考面上规定牵引：2 ms半余弦上升，加0.5 ms峰值保持。这个归一化时程/观察窗源于已验证三维预裂试样；不是声称保持历史240 µm位移边界。所有本轮幅值仅改变Pa。采用固定参考法向/切向与材料面片，面积也是参考面积，属于nominal dead traction；大转动后不再等同当前面上的纯切向/法向压力。已经脱离的原受载节点仍受载，这会影响二次裂解。\n','剪切：91个顶面四边形，2.275e−7 m²，+X。局部压入/拉脱：相同左侧42面片，1.05e−7 m²，−Z/+Z。四角各1/4积分权重；非负、无有符号外插。力矩和功率配对采用相同面心/节点插值。IDEAL_LOCAL_PULL另标 IDEAL_MECHANICAL_BOUNDARY。\n','积分仍为原生VV、κ/d不可逆、完整endpoint-star、全失效排除、原始ID及非键合接触；dt=0.25 µs，预测事件细分至0.03125 µs，max Δd≤0.1。候选在原生提交前因预测Δd拒绝；其他检查失败则停止该工况并保存未接受候选和最后接受态，不将其伪装为回滚后的成功续算。相对能量门槛1%，没有放松物理或数值门槛。\n','启动恢复：100 Pa原0.25 µs尝试在3.25 µs因绝对余额约1.35e−26 J、归一化1.3226%停止。独立保存失败后，只将最初0–2/2–4/4–8 µs的dt改为0.03125/0.0625/0.125 µs，随后恢复0.25 µs。原相对归一化和1%门槛未变；后续幅值统一采用此数值协议。原失败与恢复工况都保留，曲线优先显示同一物理工况的完成尝试。见 STARTUP_RECOVERY_FREEZE.json 和源码diff。\n','位置核验修正：1000 Pa在2.421375 ms候选处，展开式VV预测坐标比原生半步速度形式差1 ULP（3.4694e−18 m）；使用同一半步速度代数顺序后该候选差异精确为0。修正的是独立核验表达式，原生积分器/力学/绝对2e−18 m门槛未改。原失败接受前缀保留；续算先按原dt序列逐步重放并要求坐标、速度、力、κ/d、能量、计数器和时间逐位一致，再继续。见 PREDICTOR_IDENTITY_FIX.json 与续算 PREFIX_REPLAY_AUDIT.json。\n','## 实际结果\n','|模式|峰值Pa|数值状态|首软化ms|首断裂ms|历史/末态最大自由分量|末态类型|\n|---|---:|---|---:|---:|---:|---|']
 for z in cases:
  e=z['events'];fmt=lambda k:format(e[k]['time_s']*1e3,'.6g') if k in e else '未观察到'
  s.append(f"|{z['params']['mode']}|{z['params']['peak_pa']:.8g}|{z['status']}|{fmt('first_softening')}|{fmt('first_failure')}|{z['history_max_free_nodes']}/{z['last_row']['largest_detached_nodes']}|{z['classification']}|")
 s+=['\n首事件位置、原始atom/bond ID、参考与实际坐标、真实伸长及近锚/近表面标记均见 FRACTURE_EVENTS.json。初始预切完全排除。近表面是参考包围盒50 µm判据，不声称它能识别全部内部裂面。未发生事件用null，未完成数值窗不能建立无事件下界。\n','## 十三个科学问题\n','1. **软化量级？** '+range_text('softening')+'。仅对当前2.5 ms、固定区域和参考面积牵引协议成立。','2. **首完整断裂量级？** '+range_text('failure')+'。软化与完整失效的事件时间分别保存；同一粗夹区不意味着两个物理阈值相等。','3. **是否符合100–300 Pa估计？** 本剪切协议在56.2341 Pa已出现软化和完整失效，首事件峰值需求低于100–300 Pa粗估。100 Pa局部拉脱只软化、局部压入无损伤，说明方向、面积与时程会改变需求。事件发生时的瞬时Pa另列于主表；这些结果不构成普适血栓强度。','4. **三种模式区别？** '+ ' '.join(mode_summary(m) for m in MODES)+' 剪切作用面积较大；压入与拉脱面积相同而方向相反。逐档的实际损伤/连通性、冲量与功记录于主表。相同峰值Pa并非等力、等功或等应变率试验。','5. **真实多节点游离？** '+('观察到；原始存活bulk连接与任意存活interface锚连接共同判定。' if maxcase['history_max_free_nodes']>=2 else '在已接受轨迹中未观察到；仅限本试样、所测区域、载荷与时间窗。'),'6. **最大节点数及凝聚性？** '+f"历史最大 {maxcase['history_max_free_nodes']} 点，发生于 {maxcase['tag']}。每个拓扑事件记录自由成员原始ID、内部初始存活键数/当前存活键数及比例，见逐case的 DETACHED_COMPONENT_HISTORY.json；最大块接受态另存 max_free.npz。1–2节点另记未解析碎屑质量，不能因其存活一根键称宏观凝聚块。",'7. **过强载荷与二次碎裂？** '+ '; '.join(f"{z['tag']}: {len(z['secondary_fracture_events'])}次原自由成员分裂" for z in cases if z['secondary_fracture_events'])+('。没有观察到原自由分量二次分裂；不外推不存在。' if not any(z['secondary_fracture_events'] for z in cases) else '。完整ID父子对应保留，接触近邻从不当作连接。是否随Pa单调加剧需看有限采样，数值失败不是粉碎。'),'8. **响应单调？** 首软化/首失效只在实际有序样本支持时夹区；碎块尺寸并不预设单调，未对最大碎块做二分优化。下方单调性表列出已完成工况的观测比较。','9. **实际功与幅值？** 每步W增量为0.5(F_old+F_new)·(x_new−x_old)，不是末态力乘位移。下表给出真实功和冲量；不同模式功差异不能归因于方向的单独作用，因为加载面积也不同。','10. **能否转用于V3？** 不能直接搬运临界Pa。当前体尺寸、预裂尖端、底部平锚、顶面受载与真实1927点V3的曲率、异质连接、复杂载荷区不同；有限厚度和δ/dx也未作连续体收敛。没有新V3、CFD或Gc标定。','11. **普通Streaming差距？** 唯一直接可核实的上一轮机械剪切上限为0.079083543 Pa，是平面Couette模型输出，不是已测微泡Streaming实验应力。对本轮每个首软化有事件工况的峰值做数值比值如下；其时程也不同，不能当成声流不可实现性的严格倍数或临床结论。当前资料没有真实微泡到血栓表面牵引的已验证转换，因此实验差距尚未测定。声压幅值不能直接当作平均机械剪切Pa。','12. **更强空化/射流/直接泡—栓作用是否有理由？** 若普通平面剪切只进入弱弹性区，而规定牵引已经显示有限损伤需求，建立更强机制的实际力、面积、持续时间和能量预算有研究理由。此结果不证明某个空化或射流方案已经能实现这些牵引，更不证明安全性或疗效。','13. **下一轮优先？** 先做受控的小试样对照，建立变形后表面方向/面积与脱离后卸载规则，并检查等力或等功条件下的凝聚性。316.23 Pa末态原顶面法向已经旋转约94.8°，冻结参考切向几乎成为当前法向；高载荷原面片破裂后也不能再当成连续流固界面。完成这项核验后，再以本轮成功剪切工况的实际力、面积、脉冲和功为靶值，评估更强的微射流或直接泡—栓机械作用。现有证据不支持继续优化0.079 Pa附近的低Re平板流来追求本轮的块体目标，也不证明微射流已经可实现目标。不要启动高维随机载荷或完整V3扫描。\n','## 实际功、资源与有限采样单调性\n']
 for m in MODES:
  cc=sorted([z for z in done if z['params']['mode']==m],key=lambda z:z['params']['peak_pa']);metrics={}
  for key in ['max_tensile_kappa','new_softened','largest_detached_nodes','W_J']:
   vals=[z['last_row'][key] for z in cc];metrics[key]=dict(values=vals,nondecreasing=bool(np.all(np.diff(vals)>=-1e-18)))
  s.append(m+': '+json.dumps(metrics,ensure_ascii=False))
 for z,row in zip(cases,rows):
  s.append(f"- {z['tag']} 历史最大块首次达到时内部初始存活键保留比例={row['history_largest_internal_alive_fraction']}；≥14/34/114点历史最多分量数={row['history_max_N_ge14']}/{row['history_max_N_ge34']}/{row['history_max_N_ge114']}；≥14点父块二次分裂事件数={sum(len(e['parent_atom_ids'])>=14 for e in z['secondary_fracture_events'])}。")
  q=z['last_row'];ratio=z['params']['peak_pa']/.079083543 if 'first_softening' in z['events'] else None
  s.append(f"- {z['tag']}: W={q['W_J']:.9g} J，I={q['impulse_ns']:.9g} N·s，最大接受能量余额={z['max_relative_energy_residual']*100:.6g}%，相对上轮0.07908 Pa的有软化幅值比={ratio}。")
 chunk=read(O/'COHESIVE_CHUNK_AUDIT.json')[maxcase['tag']]
 if chunk['observed']:s.append(f"\n最大块逐ID追踪：{maxcase['tag']} 的 {chunk['history_largest_nodes']} 点分量在 {chunk['first_time_s']*1000:.6f} ms形成，形成时内部存活比例 {chunk['survival_at_formation']:.9g}；末态同一组原始成员仍连通={chunk['same_original_members_connected_at_final']}，可观察保持时长={chunk['observed_persistence_s']} s。同一工况发生二次裂解的最大父分量只有 {chunk['largest_secondary_parent_nodes']} 点，不能把该二次碎裂分类解释成最大块整体粉碎。")
 s.append(f"\n原生求解累计墙钟约 {resources['native_wall_s_including_failed_preflights']:.1f} s（含计时字段检查失败重放的LAMMPS日志墙钟{resources['additional_failed_preflight_native_wall_s']} s、精度1 s）、已测CPU合计 {resources['native_cpu_s']:.3f} s（该失败重放的CPU/RSS未被进程退出前的记录器捕获，未填造数值，因此CPU不是所有尝试的完整总计）、峰值RSS {resources['peak_RSS_bytes']/1024**3:.4f} GiB；{resources['accepted_steps']} 接受记录步（含续算重放前缀）、{resources['event_substeps']} 事件细分接受步。每个工况单rank单线程，无GPU；独立法向模式最多两组并行，必要重放可与其他工况重叠。没有新装依赖。资源是本轮各原生进程实测累计值，RSS是单进程峰值，不含人工分析等待，也不宣称与历史V3等成本可比。")
 if failed:
  s.append('\n保留的失败尝试（同条件后续已恢复者，不计为最终未解决物理工况）：')
  for z in failed:s.append(f"- {z['tag']}: last accepted {z['time_s']:.9g} s; {z['stop_reasons']}; 保留 keys/failed_unaccepted.npz 与 keys/last_accepted.npz，不能计为完整窗材料粉碎。")
 s.append('\n受载分量分配与游离动能见 LOAD_PARTITION_AUDIT.json；该表使用实际末态外力和速度，不把全部外功当成断裂耗能。')
 partition=read(O/'LOAD_PARTITION_AUDIT.json')[maxcase['tag']]
 s.append(f"最大块所属工况的末态总外力为 {partition['total_force_n']:.9g} N，其中 {partition['force_fraction_on_detached']*100:.6g}% 施加于自由分量；最大块本身的末态外力为 {partition['largest_component_force_n']:.9g} N。自由节点动能为 {partition['kinetic_energy_detached_j']:.9g} J。原受载点脱离后持续加速与大块凝聚性必须分别解释，不能把全部外功视作大块破坏所需能量。")
 s.append('\n当前面法向旋转与原面片面积变化诊断见 LOAD_DEFORMATION_AUDIT.json。破裂后的原四边形仅作为离散受载身份诊断，不能再当作连续真实流固界面。')
 s+=['\n## 交付与复现\n','MECHANICAL_LOAD_REQUIREMENT_TABLE.csv 是逐工况主表；LOAD_BRACKETING_SUMMARY.csv 是有限样本夹区。事件/分量/完整能量/资源与状态JSON，以及6组主图、1张全域投影补图的PNG/SVG和静态索引均在本轮新目录。证据索引引用逐case原始记录，避免重复存储全量FP64状态。执行路径见 REPRODUCE.md；源代码、库及数据SHA见清单。\n','数值门槛通过不等于已证明时间步长或空间离散收敛。材料宏观强度、Gc、界面尾部和真实泡—栓耦合仍未认证。本轮完成后停止。']
 document='\n\n'.join(part.strip() for part in s)+'\n'
 document=re.sub(r'(?<=\|)\n\n(?=\|)', '\n', document)
 (O/'REVIEW_CLOT_MECHANICAL_LOAD_REQUIREMENT_BRACKETING_V1.md').write_text(document)
 print(json.dumps(status,ensure_ascii=False))
if __name__=='__main__':main()
