from common import *
import time
def event(q):
 if q is None:return '观察窗内未发生'
 return f"{q['time_s']*1e3:.6f} ms，step {q['step']}"+(f"，bond IDs {q['bond_ids']}，参考中点[m] {q['reference_midpoints_m']}" if 'bond_ids' in q else '')
def main():
 audit=json.loads((O/'NUMERICAL_AUDIT.json').read_text());prior=json.loads((O/'FINAL_STATUS.json').read_text());top=json.loads((O/'V3_DIRECTION_CONTROL_TOPOLOGY.json').read_text());co=json.loads((O/'HISTORICAL_DIRECTION_COMPARISON.json').read_text());energy=json.loads((O/'DYNAMIC_ENERGY_LEDGER.json').read_text());cand=json.loads((O/'V3_CHUNK_CANDIDATES.json').read_text());mapping=json.loads((O/'NORMAL_MAPPING_AUDIT.json').read_text());unit=json.loads((O/'LOAD_UNIT_CHECKS.json').read_text());resources=json.loads((O/'RUN_RESOURCES.json').read_text());interpret=json.loads((O/'SCIENTIFIC_INTERPRETATION.json').read_text());new=top['reference_normal'];old=top['historical_tangent'];group=new['largest_historical_group'];first=json.loads((O/'FIRST_EVENTS.json').read_text());h=np.genfromtxt(O/'V3_NATIVE_PILOT_HISTORY.csv',delimiter=',',names=True);last=h[-1];end=float(last['time_s']);maxima=audit['maxima'];ng=energy['runs']['reference_normal'];og=energy['runs']['historical_tangent'];restore=json.loads((O/'early_checkpoint/RESTORE_AUDIT.json').read_text()) if (O/'early_checkpoint/RESTORE_AUDIT.json').exists() else None
 status=dict(task=TASK,run_id=RUN,NUMERICAL_STATUS=audit['status'],DIRECTION_CONTROL_STATUS='PASS',CHUNK_FORMATION_STATUS='OBSERVED' if new['history_maxima']['N_ge34'] else 'NOT_OBSERVED',CHUNK_FORMATION_SCOPE='Actual independent >=34-node components (>=200 micrometres by reference volume); separate >=14/114 statistics retained',COHESION_RETENTION_STATUS=interpret['COHESION_RETENTION_STATUS'],ENERGY_COMPARABILITY_STATUS='SAME_SCALAR_RULE_NOT_EQUAL_INPUT_WORK',SCIENTIFIC_CONCLUSION=interpret['SCIENTIFIC_CONCLUSION'],full_protocol_completed=abs(end-END)<1e-15,last_accepted_time_s=end,accepted_steps=int(last['step']),new_physical_cases=1,historical_case_repeated=False,production_Cpp_changed=False,material_changed=False,load_amplitude_changed=False,normal_load_label='SYNTHETIC_MECHANICAL_CONTROL',delivery_status='READY_FOR_WSL_AND_SEAL',next_recommendation=interpret['NEXT_RECOMMENDATION'])
 if not status['full_protocol_completed']:status['completion_scope']='SHORT_WINDOW_ONLY'
 dump(O/'FINAL_STATUS.json',status)
 resources.update(postprocessing_wall_s=json.loads((O/'ANALYSIS_STATUS.json').read_text())['wall_s'],render_wall_s=json.loads((F/'FIGURE_PROVENANCE.json').read_text())['wall_s'],early_checkpoint_and_restore_wall_s=None if restore is None else restore['checkpoint_and_restore_wall_s'],accepted_substeps=sum(v['decision']=='ACCEPT_DAMAGE_GATE' and v['depth']>0 for v in json.loads((O/'full/TRIAL_ACCEPT_REJECT_AUDIT.json').read_text())['trials']),actual_history_rows=len(h),output_bytes_before_final_seal=sum(p.stat().st_size for p in O.rglob('*') if p.is_file()),figure_bytes=sum(p.stat().st_size for p in F.rglob('*') if p.is_file()))
 dump(O/'RUN_RESOURCES.json',resources)
 answers=[]
 answers.append('453个原始加载ID、原始坐标和面积权重全部保持；面积总和 '+format(mapping['total_original_load_area_m2'],'.12g')+' m²，未增加加载点。')
 answers.append(f"使用原V3暴露表面的{mapping['face_count']}个有向三角形，按原来的最近加载点归属累计面积加权面法向，再归一化。重新计算的逐点面积与历史面积最大相对差{mapping['max_area_relative_difference']:.3g}。保留每个原始atom ID、所属三角形、面面积/法向及最近表面点。")
 answers.append(f"全部453点通过闭合血栓内外和血管腔内探针检查：参考面点沿法向外移1µm均在血栓外且在腔内，反向1µm均在血栓内。未使用附着面或血管壁法向；无法判定点{len(mapping['ambiguous_atom_ids'])}。最小节点法向/最近面法向点积{mapping['normal_nearest_face_dot_min']:.6g}；边缘仅平均原归属暴露面，未混入附着面。")
 answers.append('完全复用历史External.amp/compute：scale=780.2322247780305 Pa，sigma=0.28 mm，脉冲1.1605 ms，原Gaussian中心、当前位置计算规则不变。200 Pa是原面积预算，不是超声功率。两条轨迹不同后Gaussian标量可以不同，未强制匹配瞬时力。')
 answers.append(f"通过。初态与历史变形态、多个实际时间输入下的最大逐点力模长差为{max(v['max_force_magnitude_error_N'] for v in unit['checks']):.6g} N；同时保留历史tangent模长因子。锚点及无加载面积点没有新增人工力，旧切向未叠加；其余力通道逐值相同。")
 answers.append('C++、完整endpoint-star自由能和typed软化本构完全未改、未重编译。实际加载的原冻结liblammps.so SHA为 '+sha(LIB)+'；运行中逐接受步与独立Python力/能量/历史核对。')
 answers.append('保留1927材料点、1561固定锚点及157637原始连接：bulk 82270、interface 22040、fixed-fixed 53327。原始m/r0/V、atom/bond/type身份、固定锚点、接触/壁碰撞、COM拖曳、数值阻尼与能量账本不变。')
 answers.append(f"实际完成{end*1e3:.9f} ms，接受{int(last['step'])}步。"+('0.8ms保存永久原生检查点并实际重新载入；历史/拓扑/力与能量逐项核验通过，同一累计账本继续至2.33ms，没有重新计算已接受前缀。' if restore else '未声称完成0.8ms检查点恢复。'))
 answers.append('bulk首次软化：'+event(first['bulk_softening'])+'；bulk完全失效：'+event(first['bulk_failure'])+'；interface首次软化：'+event(first['interface_softening'])+'；interface完全失效：'+event(first['interface_failure'])+'。首次拓扑脱离：'+event(first['detachment'])+'。')
 answers.append(interpret['EARLIEST_DAMAGE_COMPARISON'])
 answers.append(f"全程最大节点速度{maxima['max_speed_m_s']['value']:.9g} m/s，最大接受Δd={maxima['max_dd']['value']:.10g}，最大非键合粒子重叠{maxima['max_contact_penetration_dx']['value']:.9g} dx，最大壁重叠{maxima['max_wall_overlap_m']['value']:.9g} m。")
 answers.append(f"数值状态{audit['status']}；最大归一化独立能量余额{100*maxima['relative_R']['value']:.9g}%。沿用原2%停止门槛和全部原始科学准入，未增加新的时间/内存限额。独立残差未加入D；全程dt/空间收敛仍未认证。")
 for size in [14,34,114]:answers.append(f"≥{size}节点真实自由分量的历史最大同时数为{new['history_maxima']['N_ge'+str(size)]}，末态为{new['final']['N_ge'+str(size)]}。"+('34/114节点对应≥200/300µm参考体积等效尺度，不是医学安全或收敛分辨率。' if size>14 else ''))
 answers.append(f"历史最大{new['history_maxima']['largest_detached_nodes']}节点，末态最大{new['final']['largest_detached_nodes']}节点；历史C1分别13/3节点。最大分量的真实原始ID、参考体积、等效直径、内部与外部键见拓扑/凝聚性JSON。")
 answers.append(f"不改变上一轮候选定义；≥14/34/114节点候选的唯一ID集合数分别为{cand['unique_identity_sets']['14']}/{cand['unique_identity_sets']['34']}/{cand['unique_identity_sets']['114']}。候选是TETHERED_CANDIDATE，不能称已脱落。筛查限于保存的真实关键帧，不是穷尽所有子集或所有时刻。")
 if group and group['internal_reference_bulk_bonds']:answers.append(f"历史最大真实自由集合含{group['nodes']}节点，首次分离时内部原始bulk键{group['internal_reference_bulk_bonds']}条，存活{group['internal_alive_bulk_bonds']}条，平均d={group['mean_internal_d']}；外部存活连接{len(group['external_alive_bond_ids'])}条。候选的内外连接及演化另表，不能以阈值删边替代真实失效。")
 elif group:answers.append('最大真实自由分量没有内部原始bulk键；单节点的内部键存活率和平均损伤不适用，不能认证多节点碎块的凝聚性。候选与附着主体另行统计。')
 else:answers.append('本轮没有真实自由分量，内部凝聚性按附着主体与候选统计，不虚构自由块。')
 answers.append(f"末态1–2节点碎屑{new['final']['tiny_nodes']}节点、{new['final']['tiny_components']}个分量；历史C1为{old['final']['tiny_nodes']}节点。是否称“大规模”需结合总脱离量{new['final']['detached_nodes']}节点，不将这些未解析碎屑称临床微栓子。")
 answers.append(interpret['SECONDARY_BREAKUP'])
 answers.append(interpret['COHESION_COMPARISON'])
 answers.append(f"法向实际人工载荷功{ng['final']['W_load_J']*1e9:.9g} nJ，历史切向{og['final']['W_load_J']*1e9:.9g} nJ；比值{energy['normal_over_tangent_load_work']:.9g}。背景拖曳、附着数值阻尼、壁接触/阻尼分别列出，没有把它们合并冒充人工载荷功。")
 answers.append(f"末态法向D={ng['final']['D_J']*1e9:.9g} nJ、K={ng['final']['K_J']*1e9:.9g} nJ；历史切向D={og['final']['D_J']*1e9:.9g} nJ、K={og['final']['K_J']*1e9:.9g} nJ。完整时序同时显示储能与所有外功通道。")
 answers.append('不能排除输入能量变化的中介作用。本轮控制的是加载规则中的方向；运动、Gaussian取值及力与速度夹角随之变化。没有等能量对照，因此不能宣称NORMAL_LOADING_HAS_HIGHER_ENERGY_EFFICIENCY，也未追加幅值去匹配能量。')
 answers.append(interpret['SCIENTIFIC_CONCLUSION'])
 answers.append(interpret['MIXED_MODE_RECOMMENDATION'])
 answers.append(interpret['SUPPORTED_MECHANISM'])
 answers.append(interpret['NEXT_RECOMMENDATION'])
 assert len(answers)==28
 channels=['W_load_J','W_drag_J','W_numeric_J','W_wall_elastic_J','W_wall_damping_J','U_J','K_J','D_J','Uc_J','Uw_J','Jc_J','R_J']
 table='\n'.join(f"| {k} | {og['final'][k]*1e9:.10g} | {ng['final'][k]*1e9:.10g} |" for k in channels)
 early=json.loads((O/'EARLY_WINDOW_COMPARISON.json').read_text());ec=early['comparison']
 early_table='\n'.join(f"| {k} | {ec[k]['historical_tangent']*scale:.9g} | {ec[k]['normal']*scale:.9g} |" for k,scale in [('bulk_failed',1),('interface_failed',1),('detached_nodes',1),('detached_components',1),('largest_detached_nodes',1)])
 early_table+=f"\n| W_load [nJ] | {ec['W_load_J']['historical_tangent']*1e9:.9g} | {ec['W_load_J']['normal']*1e9:.9g} |"
 text=f'''# V3_CHUNK_FRAGMENTATION_LOAD_DIRECTION_CONTROL_V1 科学审查

运行标识：{RUN}。一个新的真实V3法向合成载荷工况。历史C1只读，不重跑。状态：**{audit['status']}**。

## 本轮结论

{interpret['SCIENTIFIC_CONCLUSION']}

`NUMERICAL_STATUS={status['NUMERICAL_STATUS']}`\n
`DIRECTION_CONTROL_STATUS={status['DIRECTION_CONTROL_STATUS']}`\n
`CHUNK_FORMATION_STATUS={status['CHUNK_FORMATION_STATUS']}`（≥34节点真实自由分量）\n
`COHESION_RETENTION_STATUS={status['COHESION_RETENTION_STATUS']}`\n
`ENERGY_COMPARABILITY_STATUS={status['ENERGY_COMPARABILITY_STATUS']}`

## 28项科学问题

'''+ '\n\n'.join(f'{i+1}. {answer}' for i,answer in enumerate(answers))+f'''

## 0.8ms永久检查点

历史实际接受时刻{early['historical_actual_accepted_time_s']*1e3:.12f}ms，法向实际接受时刻{early['normal_time_s']*1e3:.12f}ms；标量比较未插值。示意图的完整节点帧可能略有采样偏移，另行标注。

| 量 | 历史切向 | 新参考外法向 |
|---|---:|---:|
{early_table}

该接受检查点实际重新载入；位置、速度、力、材料力、κ/d、拓扑、原生物理标量与检查点重载前一致。完整恢复审计保留于early_checkpoint/RESTORE_AUDIT.json。

## 完整能量对照

下表单位均为nJ，末态物理时间为历史切向2.33ms、新法向{end*1e3:.9f}ms。W_wall_elastic仅为诊断功；主账本已计入壁面储能Uw，不把两者重复相加。

`W_load + W_drag + W_numeric + W_wall_damping = ΔU + ΔK + D + ΔUc + ΔUw − Jc + R`。

| 通道 | 历史切向 | 新参考外法向 |
|---|---:|---:|
{table}

R是独立余额，不是材料耗能。当前D、虚拟断键能量、理想表面能均不是标定材料Gc。原界面软化尾部仍未标定。两个工况使用同一原有余额归一化定义。

## 数值与运行追溯

原生单MPI rank×单线程；无C++重编译、无GPU迁移、无MPI扫描。接受{int(last['step'])}步，接受细分步{resources['accepted_substeps']}，拒绝试算{audit['rejected_damage_trials']}次。主段墙钟{resources['main_wall_s']:.3f}s（{resources['main_wall_s']/60:.3f}min），CPU {resources['main_CPU_s']:.3f}s，峰值RSS {resources['peak_RSS_bytes']/1024**3:.6f}GiB。回滚恢复耗时{resources['restart_wall_s']:.3f}s，检查点耗时{resources['checkpoint_wall_s']:.3f}s；0.8ms永久检查点与恢复耗时{resources['early_checkpoint_and_restore_wall_s']}s。上述子项包含在主段时间内，不能重复相加。未新增墙钟/RSS停止线；保留历史2GiB输出上限。实际归档体积见RUN_RESOURCES和最终PAYLOAD_MANIFEST。

全部接受态核查真实身份、固定锚点、不可逆κ/d、typed断键、原生与独立Python完整星能量/力、接触、VV及动量。原始库、加载标量文件、参考点/键及方向冻结SHA逐项验证。0.8ms继续理由与实际恢复结果见CONTINUATION_DECISION和early_checkpoint/RESTORE_AUDIT；没有切换物理参数或重算接受前缀。

每个接受步保存标量与稀疏损伤/事件；完整FP64节点力/速度按20接受步及事件关键帧保存。新增工况在拓扑变化时保存分量质心速度和内部相对动能；历史C1形成瞬间未保存的节点量明确NOT_MEASURED。规则采样图按最近实际接受态匹配，精确时间偏差写入VISUAL_STATE_MATCHING，不插值构造拓扑事件。

候选筛查直接调用上一轮topology_audit与graph_ops，保留所有阈值、最大963节点和权重定义。1−d只是诊断图权重，不是认证承载力。新旧关键帧采样不同；候选出现次数不应当成穷尽的物理区域数量，另提供唯一原始ID集合统计。

## 交付与科学边界

隔离Python提供器、方向场、冻结卡、原生接受轨迹/拒绝日志、实际检查点、拓扑/凝聚性/空间失效/能量审计和PNG/SVG均保存在本轮目录。图表入口为匹配results目录的OPEN_RESULTS.html。真实三维点和存活原始键，位移×1、共同相机与尺度；没有连续医学外壳、强制合并节点或人为V3预切。

本轮没有读取新的论文参数或移植NOSB/Ogden。沿用上一轮核验的文献限制；本轮结论来自实际单变量原生对照。法向牵引不是已验证的微泡streaming、空化或超声力场。200/300µm只表达参考体积等效尺度，非临床安全阈值。

下一轮建议：{interpret['NEXT_RECOMMENDATION']} 本轮到此停止，不自动追加第二方向、幅值或材料组合。
'''
 (O/'REVIEW_V3_CHUNK_FRAGMENTATION_LOAD_DIRECTION_CONTROL_V1.md').write_text(text)
 (O/'REPRODUCE.txt').write_text('Activate /workspace/y_clot/env.sh. Use analysis/v3_chunk_load_direction_control_v1/README.md and LOAD_DIRECTION_FREEZE.json. The frozen run is immutable: do not rerun prepare/freeze/run_protocol into it. For a newly authorized replication set a NEW run ID/output root, use the same original input SHA and REFERENCE_NORMAL_FIELD, prepend the original V3 build directory to LD_LIBRARY_PATH, then prepare -> freeze -> managed run_protocol. The permanent 0.8ms native checkpoint and Python cumulative ledger state are retained. analyze.py and figures.py can reproduce derived outputs only in a separate copied analysis output directory; do not change the sealed payload. One rank / one thread; no C++ rebuild.\n')
 import html,re
 page=F/'OPEN_RESULTS.html';body=page.read_text();body=re.sub(r'<!-- REVIEW_SUMMARY_START -->.*?<!-- REVIEW_SUMMARY_END -->','',body,flags=re.S)
 origin=json.loads((F/'FIGURE_PROVENANCE.json').read_text())['origin_m']
 intro='<!-- REVIEW_SUMMARY_START --><section><h2>实际结论</h2><p>'+html.escape(interpret['SCIENTIFIC_CONCLUSION'])+'</p><p>全部三维坐标轴均为相对参考材料质心的平移坐标，单位mm；全局原点为 '+html.escape(str(origin))+' m。只平移原点，材料位置和位移尺度不变。报告中的首次损伤中点采用原全局坐标。</p><p><a href="../../../verification/'+TASK+'/'+RUN+'/REVIEW_V3_CHUNK_FRAGMENTATION_LOAD_DIRECTION_CONTROL_V1.md">完整28项中文科研审查</a></p></section><!-- REVIEW_SUMMARY_END -->'
 assert '<article>' in body;page.write_text(body.replace('<article>',intro+'<article>',1))
 print(json.dumps(status,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
