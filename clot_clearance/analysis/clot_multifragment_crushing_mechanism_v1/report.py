from common import *
from analyze import read
def main():
 names=['uniform','heterogeneous'];c={k:read(O/'cases'/k/'ANALYSIS.json') for k in names};s={k:read(O/'cases'/k/'SUMMARY.json') for k in names};cr={k:read(O/'cases'/k/'CRACK_NETWORK_AUDIT.json') for k in names};r=read(O/'RUN_RESOURCES.json');u=c['uniform'];h=c['heterogeneous'];complete=all(q['status']=='COMPLETED' for q in c.values());anytarget=any(q['exploratory_target']=='EXPLORATORY_MULTIFRAGMENT_TARGET_MET' for q in c.values())
 status=dict(INDENTER_CONTACT_STATUS='PASS_NATIVE_CURRENT_COORDINATE_CONTACT',UNIFORM_CRUSHING_STATUS=u['classification'],HETEROGENEITY_IMPLEMENTATION_STATUS='PASS_ORIGINAL_ID_ONSET_FIELD_SORTING_RESTART_ENERGY_PARITY',MULTIPLE_CRACK_NUCLEATION_STATUS={k:dict(discrete_network_births=cr[k]['n_births'],max_simultaneous=cr[k]['max_simultaneous_networks'],discrete_network_merge_events=cr[k]['n_merge_events'],continuous_crack_surfaces='NOT_RECONSTRUCTED_OR_CERTIFIED') for k in names},MULTIFRAGMENT_FORMATION_STATUS={k:q['exploratory_target'] for k,q in c.items()},FRAGMENT_COHESION_STATUS={k:dict(history_max_medium_components=q['history_max_N20_200'],final_medium_components=q['final']['N20_200'],final_largest_nodes=q['final']['largest_nodes'],final_components_file=f'cases/{k}/ANALYSIS.json') for k,q in c.items()},POST_UNLOAD_STABILITY_STATUS={k:dict(medium_cohort=q['withdrawal_medium_cohort'],medium_cohort_remaining_connected=q['withdrawal_medium_cohort_surviving'],contact_after_withdrawal=q['contact_after_scheduled_withdrawal'],actual_contact_free_tail_s=q['actual_contact_free_tail_s'],claim='CONNECTED_THROUGH_OBSERVED_WINDOW only') for k,q in c.items()},NUMERICAL_ADMISSION_STATUS=dict(status='PASS_EXISTING_GATES' if complete else 'NUMERICALLY_UNRESOLVED',completed_cases=sum(q['status']=='COMPLETED' for q in c.values()),relative_energy_gate=.01,case_maxima={k:q['max_relative_energy_residual'] for k,q in c.items()},time_and_space_convergence='NOT_ESTABLISHED'),V3_TRANSFER_READINESS='NOT_READY_FOR_FULL_V3_PHYSICAL_TRANSFER',SCIENTIFIC_CONCLUSION='EXPLORATORY_TARGET_OBSERVED_IN_FROZEN_COMPARISON' if anytarget else 'EXPLORATORY_MULTIFRAGMENT_TARGET_NOT_OBSERVED_IN_FROZEN_COMPARISON')
 dump(O/'FINAL_STATUS.json',status)
 def stats(k):
  q=c[k];z=q['final'];return f"{k}：首软化 {q['first_softening_s']} s；首完整失效 {q['first_failure_s']} s；历史最多20–200点块 {q['history_max_N20_200']}；末态 {int(z['N20_200'])}；历史/末态最大块 {q['history_max_free_nodes']}/{int(z['largest_nodes'])} 点；末态最大块/初始材料 {100*z['largest_fraction']:.4g}%；1–2点碎屑/脱离质量 {('未脱离，比例未定义' if z['debris12_detached_fraction'] is None else format(100*z['debris12_detached_fraction'],'.5g')+'%')}。"
 def event_place(k):
  e=s[k]['events'].get('first_failure');return '未观察到' if not e else f"原始bond ID {e['original_bond_ids']}，参考中心m {e['reference_midpoints_m']}"
 def cohesion(k):
  entries=[q for q in c[k]['final_components'] if q['nodes']>=2]
  if not entries:return k+'末态没有至少两节点的内部连通自由分量'
  sizes=[q['nodes'] for q in entries];survival=[q['internal_survival'] for q in entries if q['internal_survival'] is not None]
  return f"{k}末态{len(entries)}个至少两节点自由分量，节点数范围{min(sizes)}–{max(sizes)}，内部原始连接存活率范围{min(survival):.6g}–{max(survival):.6g}"
 def criterion(k):
  q=c[k];z=q['final'];return f"{k}: {q['exploratory_target']}；峰值中等块数 {q['history_max_N20_200']}，末态中等块质量/初始质量 {100*z['mass_fraction20_200']:.5g}%"
 lines=['# CLOT_MULTIFRAGMENT_CRUSHING_MECHANISM_V1',f'运行：{RUN}（UTC标识）。研究记录日期：2026-10-17，Europe/Zurich。两组同一运动学加载对照，标签 `KINEMATIC_IMPACT_SURROGATE`、`EXPLORATORY_UNCALIBRATED`。','## 冻结模型和实际实现','1792材料点+336底部固定锚点、95666原始键和2736原始预切保持原字节。G=940 Pa，K=33333.333333 Pa，ρ=1000 kg/m³，dx=50 µm，δ=150.5 µm。bulk s0/sc=0.15/0.18；interface 0.140625/0.16875；固定锚连接1/1.2。没有增加初始切口、随机删键、修改参考weighted volume或质量。',
 '压头为半径182.818 µm的有限球，圆形投影面积等于上一轮LOCAL_IMPINGEMENT的42面片面积1.05e−7 m²，中心沿用这些面片的参考中心。初始间隙25 µm；0–0.2 ms平滑接近，0.2–1.2 ms压入240 µm，保持至1.4 ms，1.4–2.4 ms撤回初始位置，继续至3.6 ms。行程取已有240 µm位移尺度，约试样高度30%；1 ms压入约1.21个剪切波穿越时间。不是把位移换算成1000 Pa，也不是具有真实质量/动量的锤击。',
 '复用原生LAMMPS `fix indent sphere`：δc=max(R−|x−c|,0)，Ui=ΣKind δc³/3，Fi=Kind δc²(x−c)/|x−c|。Kind=57973.58577 N/m²，使0.1dx重叠处的切线刚度等于既有粒子接触刚度0.5797358577 N/m。原生协调器仅新增在物理SI力转换之后的单压头准入；没有Python力回调、原参考节点死牵引或人为排除脱离颗粒。接触与释放依赖当前距离；再次撞上球的脱离碎片仍可受力。',
 '参考：[LAMMPS fix indent官方说明](https://docs.lammps.org/fix_indent.html)。实际兼容性另由服务器锁定源码及新库检查确认；不以网页当前版本代替源码审计。势能、反力及单位均用本轮实际原生输出核对。',
 '唯一非均匀场采用seed=20261017、lc=250 µm、高斯平滑白噪声。bulk键中点场减均值并按最大绝对值归一化，s0_e=0.15(1+0.10q_e)，范围0.135–0.1647022265；α是10%最大有界偏差，不是10%标准差；本场bulk s0实际变异系数约4.09%。lc/dx=5、lc/δ≈1.661；有限域与非局部平均会影响效果。sc=0.18固定，interface不变。随机场在正式计算前冻结且未更换；它不是实测血栓缺陷场，也不保证Gc不变。',
 '原生参数按原始bond ID存储，输入逆序仍正确；全局参考表与atom sorting/ownership槽位分离，新增restart V3保存逐键s0，兼容读取旧V2。一次组合排序/restart短检查实际改变2126个atom槽位并保持参数、κ/d及耗能一致。完整弹性响应与均匀体逐位一致；固定损伤力与变分势能梯度检查通过。没有重新运行MPI性能矩阵。',
 '正式积分继续使用原生VV及预测事件子步：基步0.25 µs、预测事件细分下限0.03125 µs、单步Δd≤0.1、相对能量门槛1%。预测拒绝发生在原生run之前；已失败候选不得计入接受态。没有数值阻尼。门槛通过不等于时间/空间收敛。',
 '## 两组实际结果',stats('uniform'),stats('heterogeneous'),criterion('uniform'),criterion('heterogeneous'),
 '|工况|峰值反力 µN|最大压头重叠/dx|末态外功 nJ|末态动能 nJ|损伤耗能 nJ|最大相对能量余额 %|', '|---|---:|---:|---:|---:|---:|---:|']
 for k in names:
  q=c[k];z=q['final'];lines.append(f"|{k}|{q['peak_reaction_N']*1e6:.8g}|{q['max_overlap_dx']:.7g}|{z['W_J']*1e9:.8g}|{z['K_J']*1e9:.8g}|{z['D_J']*1e9:.8g}|{q['max_relative_energy_residual']*100:.7g}|")
 lines+=['## 十六个科学问题',
 '1. **压头是否通过当前接触传力？** 是。未接触/释放力为0，接触方向、作用反作用、势能梯度与移动压头功均通过新接口检查；实际反力和接触数逐步保存。',
 '2. **脱离材料是否避免永久加载？** 是。没有受载原始节点表。材料点只有进入球体罚接触区才受力；脱离后仍可能真实再接触，不能以拓扑状态人为关闭碰撞。',
 '3. **均匀体如何破坏？** '+stats('uniform')+'分类：'+', '.join(u['classification'])+'。',
 '4. **非均匀场改变首完整断裂位置吗？** '+('没有改变。' if s['uniform']['events'].get('first_failure',{}).get('original_bond_ids')==s['heterogeneous']['events'].get('first_failure',{}).get('original_bond_ids') and 'first_failure' in s['uniform']['events'] else '原始键身份对比如下。')+'均匀：'+event_place('uniform')+'；非均匀：'+event_place('heterogeneous')+'。首软化位置另见FRACTURE_EVENTS.json与图3。仅一组随机场对照，不推断总体统计规律。',
 '5. **多个独立裂纹萌生区域？** '+ '; '.join(f"{k}：自然完整失效bulk键网络共记录{cr[k]['n_births']}次离散网络出生，最多同时{cr[k]['max_simultaneous_networks']}个" for k in names)+'。以共享原始端点定义网络，初始预切排除。这是离散损伤区诊断，尚不等同经过分辨率验证的连续裂纹数量。',
 '6. **裂纹是否贯通？** '+ '; '.join(f"{k}发生{cr[k]['n_merge_events']}次原始失效键网络合并" for k in names)+'。网络合并与存活材料图的真实分离分别记录；不能只凭断键点靠近就声称形成连续贯通裂面。',
 '7. **至少5个20–200节点游离块？** '+criterion('uniform')+'；'+criterion('heterogeneous')+'。完整探索目标还要求最大块占比、碎屑比例及撤回后连接保持，不能只看数量。',
 '8. **最大游离块占比？** '+ '; '.join(f"{k}末态{100*c[k]['final']['largest_fraction']:.5g}%，历史最大{100*c[k]['history_max_free_nodes']/1792:.5g}%" for k in names)+'。',
 '9. **碎屑比例？** '+ '; '.join(stats(k) for k in names)+'分母是实际脱离质量；没有脱离时留空，不解释为零阈值。',
 '10. **内部连接保留？** '+ '; '.join(cohesion(k) for k in names)+'。每个接受态对每个真实自由分量保存内部原始存活键数、存活率、平均d、质心速度和内部相对动能。1–2节点属于未解析碎屑，即使一根键存活也不能称宏观凝聚块。具体分量见MULTIFRAGMENT_SIZE_DISTRIBUTION.csv和图6。',
 '11. **卸载后保持连通吗？** '+ '; '.join(f"{k}撤回结束时20–200点分量{len(c[k]['withdrawal_medium_cohort'])}个，其中{c[k]['withdrawal_medium_cohort_surviving']}个原始成员组持续连通至末态；撤回后再接触={c[k]['contact_after_scheduled_withdrawal']}；末尾连续无接触观察窗={c[k]['actual_contact_free_tail_s']} s" for k in names)+'。只使用CONNECTED_THROUGH_OBSERVED_WINDOW，不称永久稳定。',
 '12. **非均匀性是否改善多块碎裂？** '+f"历史最多中等块数由{u['history_max_N20_200']}变为{h['history_max_N20_200']}，末态由{int(u['final']['N20_200'])}变为{int(h['final']['N20_200'])}；末态中等块质量比例由{100*u['final']['mass_fraction20_200']:.5g}%变为{100*h['final']['mass_fraction20_200']:.5g}%。实际输入功同时由{u['final']['W_J']:.9g}变为{h['final']['W_J']:.9g} J。需同时考虑碎屑、内部存活率及观察时长，不能将同轨迹称为等功比较。",
 '13. **若不能多块压碎，原因更支持什么？** 已有接触反力、软化/完整失效、局部变形与碎片图共同限定判断。两组实际软化/失效时间和中等块数量见上表；存在局部材料损伤不等于已输入足以贯通多个块体的机械功。若目标未达到，本轮能支持的结论是该冻结加载、lc和拉伸历史损伤规则的组合没有组织出目标分块。单一行程和单一随机场不能独立识别加载不足、局部化不足或本构局限。LOCAL_DEFORMATION_AUDIT给出起裂附近保留邻域的Green应变、非仿射残差及材料键virial/参考体积诊断；该应力排除接触，不是经过标定的局部Cauchy应力。',
 '14. **值得研究剪切/压缩混合失效吗？** 当前κ只累计正键伸长；压入只可经侧向拉伸、几何剪切等间接触发该损伤律。结果可为混合失效机制提出可检验问题，但尚不足以授权替换本构或宣称新失效准则必然改善。下一轮优先核查实际起裂区拉伸/压缩/剪切路径和接触重叠误差，再决定是否需要经标定的混合准则。',
 '15. **能否迁移真实1927节点V3？** 尚无充分物理依据。平底固定锚、已有预裂纹、刚性球压入、单一未标定随机场与真实弯曲血管/微泡载荷均有差异；接触分辨率、时间/空间收敛、材料耗能和真实力传递链未完成标定。本轮不运行新V3。',
 '16. **资源消耗？** '+f"两组正式驱动（含原生步进、在线核验及记录）累计墙钟{r['native_wall_s']:.3f} s、CPU{r['native_cpu_s']:.3f} s、单进程峰值RSS{r['peak_RSS_bytes']/1024**3:.5f} GiB。每组1 rank×1线程，顺序执行。隔离构建另用{r['build_wall_s']:.3f} s，成功关键检查{r['key_check_wall_s']:.3f} s；早期变量语法失败预检单独保留，其未计时部分不伪造为0。没有安装新依赖、CFD或完整V3。",
 '## 能量、接触精度与解释边界',
 'W由接受态之间的压头中心路径积分∂Ui/∂c·dc得到，在接触进入/离开根处分段做Gauss积分。检查4点与8点积分差；不使用末态反力×位移，也不把材料点接触功再加一次。完整账本为 W=ΔU+ΔK+D+ΔUc−Jc+ΔUi+R，数值阻尼功为0。Jc仅处理新断键时原粒子接触排除解除产生的势能偏移。',
 '相对能量门槛使用历时最大能量尺度：max_history(max(|W|, |ΔU|+K+D+Uc+Ui))，不是末态外功作分母。末态|R/W|另为'+ '; '.join(f"{k} {abs(c[k]['final']['R_J']/c[k]['final']['W_J'])*100:.7g}%" for k in names)+'；该区分不改变原有接受门槛或原始余额。',
 '原生力/势能及余额通过说明实现一致，不自动证明罚接触的几何精度。图2与上表报告真实最大压头重叠，应与50 µm点距和接触半径共同审查；未做刚度或网格收敛扫描。',
 '原始日志FREE_EVOLUTION是规定运动时程的阶段标签，表示压头已回到初始位置，并不保证无接触；反弹或脱离材料仍可碰到静止球。真正无压头力的区段必须依据实际contact_nodes=0判定。碎屑比例同时参考脱离质量和初始质量；描述性分类仅在脱离量≥初始质量10%且1–2点占脱离量>20%时称EXCESSIVE_PULVERIZATION，否则将微小碎屑脱落称SURFACE_EROSION，并另标DEBRIS_DOMINATED_DETACHMENT。该分类约定不改变预先冻结的5块探索目标。',
 '所有碎片由真实存活材料连接和锚点可达性定义，没有接触合并、目标块分组、增加预切或随机删键。20点不等于充分解析的200 µm实体；≥34/114点及≥200/300 µm参考体积等效尺寸另行记录。这些不是医学安全标准。',
 '## 输出与下一步',
 'OPEN_RESULTS.html索引8组主图（图4分两组、另附完整末态投影），CSV/JSON保留实际分量历史、事件、能量、反力和资源。必要FP64接受态、逐步损伤增量、原生日志、失败预检及源码SHA均保留。REPRODUCE.md说明复现入口。',
 '本轮最有价值的结论来自冻结对照，无论是否达到5块都不换种子或提高行程。下一轮先审查损伤网络为何贯通或未贯通、罚接触重叠及卸载后轨迹，再设计一个具有明确可辨识性的局部机制对照。真实微泡、微射流或超声到压头轨迹/力/功的物理映射仍缺失。完成后停止，等待科学审查。']
 for k in names:
  z=c[k]['final'];thresholds=['N_ge34','N_ge114','N_ge200um','N_ge300um'];values=', '.join(f"{key}: {int(z[key])}/{c[k]['history_max_size_threshold_counts'][key]}" for key in thresholds)
  lines.insert(lines.index('## 输出与下一步'),f"{k}额外尺寸阈值（末态/历史最多）：{values}。末态脱离总量{int(z['detached_nodes'])}/1792点，即初始质量{z['detached_nodes']/1792*100:.6g}%；尺寸等效直径按原始体积计算，不表示已经解析的连续碎块。")
 for name in names:
  diag=read(O/'cases'/name/'LOCAL_DEFORMATION_AUDIT.json').get('first_failure')
  if diag:
   local=[q for q in diag['local_affine'] if q['resolved']]
   if local:
    eig=np.array([q['green_principal_strains'] for q in local]);mean=np.array([q['material_virial_mean_tension_positive_Pa'] for q in local]);nonaff=np.array([q['affine_fit_relative_residual'] for q in local])
    lines.insert(lines.index('## 输出与下一步'),f"{name}首失效附近的实际诊断：保留邻域Green主应变逐分量范围{eig.min(0).tolist()}至{eig.max(0).tolist()}，材料键virial平均正拉负压指标{mean.min():.7g}至{mean.max():.7g}Pa，非仿射拟合相对残差{nonaff.min():.7g}至{nonaff.max():.7g}。局部键伸长与邻域平均应变并非同一量；压缩背景中的局部拉伸断键不证明已有压缩损伤能力，也不能仅凭无中等块认定拉伸损伤律错误。")
 if u['first_softening_s'] is not None and h['first_softening_s'] is not None:
  delay=(h['first_softening_s']-u['first_softening_s'])*1e6
  lines.insert(lines.index('## 输出与下一步'),f"非均匀场相对均匀体的首次软化时间差为{delay:.6g}µs（正值表示更晚）。首次软化原始键ID：均匀{s['uniform']['events']['first_softening']['original_bond_ids']}，非均匀{s['heterogeneous']['events']['first_softening']['original_bond_ids']}。只对这组冻结场和轨迹成立，不解释为所有缺陷场都会延迟或提前起裂。")
 timing=read(O/'TIMESTEP_BOUNDARY_AUDIT.json')
 for k in names:
  t=timing['cases'][k]
  lines.insert(lines.index('## 输出与下一步'),f"{k}实际时间步审计：接受{t['accepted_steps']}步，非浮点残余步的最小dt={t['min_non_roundoff_dt_s']:.9g}s；另有{len(t['boundary_residual_steps'])}个规定边界对齐残余步，实际最小dt={t['min_actual_accepted_dt_s']:.9g}s。这些末端浮点残余步按原始日志保留，不归类为突破损伤预测细分下限；详见TIMESTEP_BOUNDARY_AUDIT.json。")
 context=read(O/'HISTORICAL_LOADING_CONTEXT.json')
 lines.insert(lines.index('## 输出与下一步'),f"与上一轮局部1000Pa死牵引相比：旧外功{context['old_W_J']*1e9:.7g}nJ，动能{context['old_K_J']*1e9:.7g}nJ，损伤耗能{context['old_D_J']*1e9:.7g}nJ；本轮均匀/非均匀外功分别只有旧值的{context['current']['uniform']['W_ratio_to_old']*100:.5g}%/{context['current']['heterogeneous']['W_ratio_to_old']*100:.5g}%。两轮时长、力与面积路径不同，不能把碎片变化全部归因于接触机制或非均匀性，也不能把旧死载产生的高动能当成必须消耗的断裂功。")
 text='\n\n'.join(lines)+'\n';text=text.replace('|\n\n|','|\n|');(O/'REVIEW_CLOT_MULTIFRAGMENT_CRUSHING_MECHANISM_V1.md').write_text(text);print('REPORT_DONE',status['SCIENTIFIC_CONCLUSION'])
if __name__=='__main__':main()
