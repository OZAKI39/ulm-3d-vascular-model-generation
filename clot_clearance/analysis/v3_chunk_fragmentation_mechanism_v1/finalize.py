from common import *
import shutil,time,re
def main():
 small=json.loads((O/'SMALL_SPECIMEN_VALIDATION.json').read_text());assert small['status']=='COMPLETED' and len(small['cases'])==2
 cand=json.loads((O/'V3_CHUNK_CANDIDATES.json').read_text());comp=json.loads((O/'COHESIVE_ADHESIVE_COMPETITION.json').read_text());mixed=json.loads((O/'MIXED_MODE_DIAGNOSTIC.json').read_text());contact=json.loads((O/'CONTACT_0376MS_AUDIT.json').read_text());group=json.loads((O/'LOCAL_FRAGMENT_COHESION.json').read_text());base=json.loads((O/'BASELINE_AUDIT.json').read_text());h=np.genfromtxt(B/'V3_NATIVE_PILOT_HISTORY.csv',delimiter=',',names=True)
 c0,c1=small['cases'];g=next(s for s in group['states'] if s['step']==3282);peak=max((sum(t['below_history_onset_and_shear_ge_015'] for t in f['typed']),f['step']) for f in mixed['frames']);syntheticmax=max(c['maximum_free_nodes'] for c in small['cases']);syntheticchunks=syntheticmax>=34
 firstdesign=json.loads((O/'benchmark/DESIGN_FREEZE.json').read_text());finaldesign=json.loads((O/small['design_freeze']).read_text())
 for key in ['hypothesis','cases','only_changed_factor','material','geometry','loading','integration','library_sha256']:assert firstdesign[key]==finaldesign[key],key
 unchanged_inputs={}
 for name in ['opening_z','shear_x']:
  unchanged_inputs[name]={}
  for filename in ['reference.npz','reference.native','body.data']:
   h0=sha(O/'benchmark'/name/'provenance'/filename);h1=sha(O/small['benchmark_root']/name/'provenance'/filename);assert h0==h1;unchanged_inputs[name][filename]=h1
 dump(O/'IMPLEMENTATION_REPAIR_AND_FREEZE_IDENTITY.json',dict(status='PASS',physical_design_unchanged=True,exact_input_sha256=unchanged_inputs,first_failure='second NVE fix rejected at grip release',fix='unfix old integration and grip motion; use exactly one NVE fix covering all free material',second_failure='grouped native replay changed floating time accumulation; strict continuation gate rejected it',recovery='Replay original run 1 sequence; x/v/F/material_F/kappa/d/time/topology/K/Uc exactly match original accepted prefix',failed_sources='benchmark_attempt01_source and benchmark_attempt02_source; exact original frozen hashes checked',no_material_contact_or_amplitude_tuning=True))
 oldcase=json.loads((O/'benchmark/opening_z/SUMMARY.json').read_text());log=(O/'benchmark_corrected_v2/opening_z/native.log').read_text();hh,mm,ss=map(int,re.findall(r'Total wall time: (\d+):(\d+):(\d+)',log)[-1]);replay_failed=hh*3600+mm*60+ss;replay=json.loads((O/small['benchmark_root']/'opening_z/PREFIX_REPLAY_AUDIT.json').read_text())
 cost=dict(original_failed_prefix_wall_s=oldcase['wall_s'],rejected_grouped_replay_native_wall_s_rounded=replay_failed,corrected_opening_wall_s_including_exact_replay=c0['wall_s'],exact_replay_included_in_corrected_opening_wall_s=replay['replay_wall_s'],shear_case_wall_s=c1['wall_s'],total_benchmark_attempt_wall_s_approx=oldcase['wall_s']+replay_failed+c0['wall_s']+c1['wall_s'],benchmark_peak_RSS_bytes=max(oldcase['peak_RSS_bytes'],c0['peak_RSS_bytes'],c1['peak_RSS_bytes']),P0_P1_wall_s=base['wall_s'],topology_wall_s=cand['wall_s'],spatial_audit_wall_s=comp['wall_s'],mixed_mode_wall_s=mixed['wall_s'],MPI_ranks=1,threads=1,GPU_used=False,arbitrary_wall_or_RSS_stop=False,offline_peak_RSS='NOT MEASURED',note='Approximate benchmark total includes both rejected attempts; accepted-prefix replay is included once, within corrected opening wall. Code development and transfer time are not included.')
 dump(O/'RUN_RESOURCES.json',cost)
 # No result-driven constitutive modification or second V3 is inferred from a toy test.
 nextchange='仅改变真实V3外载方向：参考切向改为冻结的参考法向，保持同一坐标/时刻的逐点力幅值、Gaussian空间包络、脉冲、材料、界面和接触不变，先做短窗机制对照，检验载荷路径能否从表层侵蚀转向凝聚开裂；法向工况须明确标为合成对照，不称为已确认的微泡载荷。'
 status=dict(status='COMPLETED_P0_P4_NO_NEW_V3',ROOT_CAUSE_STATUS='SUPPORTED_HYPOTHESIS',CHUNK_FORMATION_STATUS='NOT_OBSERVED',CHUNK_FORMATION_SCOPE='REAL_V3: no >=200 um component/candidate found in archived states and exact accepted topology',SYNTHETIC_CHUNK_FORMATION_STATUS='OBSERVED' if syntheticchunks else 'NOT_OBSERVED',synthetic_maximum_free_nodes=syntheticmax,NEXT_RECOMMENDED_CHANGE=nextchange,real_V3_new_cases=0,material_law_changed=False,native_Cpp_changed=False,baseline_replayed=False,small_specimen_physical_cases=2,small_specimen_failed_implementation_prefix_retained=True,scientific_limits=['No unique causal identification of V3 fragmentation','No calibrated 3D mixed-mode criterion','No calibrated bulk/interface Gc','No horizon/size or full time-step convergence','No clinical embolus size inference','Synthetic loading is not physiological flow'],final_V3_detached_nodes=int(h['detached_nodes'][-1]),final_V3_detached_components=int(h['detached_components'][-1]),final_V3_largest_nodes=int(h['largest_detached_nodes'][-1]),historical_V3_largest_nodes=int(h['largest_detached_nodes'].max()),large_tethered_status=cand['large_chunk_status'],small_benchmark_max_relative_R={c['case']:c['max_relative_energy_residual'] for c in small['cases']})
 dump(O/'FINAL_STATUS.json',status)
 dump(O/'P5_ADMISSION_DECISION.json',dict(decision='NOT_RUN',new_real_V3_cases=0,reason='P2 supports local cohesion loss and secondary breakup, but has not isolated a unique cause. The direction test uses prescribed grips, whereas V3 uses spatially localized surface forces; it establishes model capability without validating a unique real-V3 mechanism. Objective shear diagnostics do not establish a calibrated failure blind spot. A real-V3 normal-traction control must first freeze direction/sign and equal point-force-magnitude semantics. Therefore the joint admission conditions are not all met in this round.',unchanged=['real geometry','material card','contact stiffness','load amplitude/pulse','accepted historical baseline']))
 table='\n'.join(f"| {c['case']} | {c['accepted_steps']} | {c['time_s']*1e3:.6g} | {c['last_row']['bulk_failed']} | {c['last_row']['interface_failed']} | {c['maximum_free_nodes']} | {100*c['max_relative_energy_residual']:.6g}% | {c['wall_s']:.2f} |" for c in small['cases'])
 details=[]
 for c in small['cases']:
  a=c['spatial_analysis'];q=c['quantitative_response'];frac=a['first_50us_notch_fraction'];details.append(f"- **{c['case']}**：首本体软化 {a['first_bulk_soft_time_s']*1e3:.6g} ms；随后50 µs内软化键落在裂尖δ邻域的比例 {frac*100:.4g}%；末态原裂纹面尚存 {a['final_remaining_cross_plane_bonds']}/{a['initial_remaining_cross_plane_bonds']} 条原始未预切跨面连接。峰值顶端反力 {q['peak_abs_top_reaction_N']*1e6:.6g} µN，末态 W/U/K/D/Uc/R = {q['final_W_J']*1e9:.6g}/{q['final_U_J']*1e9:.6g}/{q['final_K_J']*1e9:.6g}/{q['final_D_J']*1e9:.6g}/{q['final_contact_U_J']*1e9:.6g}/{q['final_R_J']*1e9:.6g} nJ，K/U={q['final_K_over_U']:.6g}。")
  chunks=[v for v in a['free_components'] if v['nodes']>=34]
  if chunks:
   z=max(chunks,key=lambda v:v['nodes']);details.append(f"  末态最大真实自由块为{z['nodes']}点，内部原始bulk键存活{z['internal_alive_bulk']}/{z['internal_reference_bulk']}（{100*z['internal_alive_fraction']:.5g}%），原始内部预切{z['internal_initial_cut']}条，外部存活连接{z['external_alive']}。失去锚点连接的首接受时刻与夹具释放后的自由统计时刻分别见chunk_time_semantics，不能混称同一断裂事件。")
 review=f'''# V3_CHUNK_FRAGMENTATION_MECHANISM_V1 科学审查

运行标识：{RUN}。历史唯一对照为 `lps_variational_softening_v3_integration_and_single_pilot_v1/20261014T232612Z`。本轮完成只读轨迹分析、客观混合模式模块及两组实际原生小试样；没有新真实V3运行，也没有更改生产材料律或接触刚度。

## 直接结论

`ROOT_CAUSE_STATUS=SUPPORTED_HYPOTHESIS`：证据支持自由表面受载区的局部凝聚性先受损、而界面与主体大部分连接仍保留，随后小分量带着内部相对动能脱离并继续裂解。尚未通过唯一物理因素对照确认V3的主导原因。

`CHUNK_FORMATION_STATUS=NOT_OBSERVED`：此状态明确指真实V3中≥200 µm离散体积尺度的游离块及本轮筛查的大块候选，不否认历史13节点小分量。独立小试样的≥34节点自由块状态为 **{status['SYNTHETIC_CHUNK_FORMATION_STATUS']}**，最大 {syntheticmax} 节点；两者不能混为一项真实V3成果。

## 十项必答问题

1. **真实V3有没有内部凝聚、尚未脱离的大块？** 在17个完整接受关键帧的完整度阈值连通分量与谱切割搜索中，≥34/114节点合格候选均为0，故 `NO_COHERENT_CHUNK_CANDIDATE`。这不是对所有时刻、所有子集的数学穷尽证明。实际连接图另在全部10503个接受状态（含t=0）逐项核对，整个历史最大仅13节点，末态3节点。

2. **有哪些真实桥接键？** 唯一满足冻结筛查条件的是1.000 ms的三节点集合 `[1625,1632,1649]`，等效体积直径约89.47 µm。它由 bulk bond **90571** 拴住，没有存活界面外连；该键 κ=0.1698109341、d=0.6999879333，已经进入拉伸软化，并于1.00425 ms自然失效。到末态这三个原始点又各自分开。候选的完整外连ID、类型、κ/d、伸长及演化均在JSON中；没有把诊断割边当成真实断键。

3. **是否存在内部凝聚性过早丧失的证据？** 对历史最大13节点集合，0.554125 ms脱离时42条原始内部连接仅13条存活，内部平均d=0.690476；0.5745 ms分成9/3/1节点。故对这个实际小分量，内部凝聚性在脱离前已明显下降，并出现脱离后继续碎裂。不能由此说整个血栓都已粉碎：末态仍有1747材料点附着，其内部原始bulk键存活率99.33%，还有21620条界面连接存活。全体bulk完全失效比例约11.63%，主要损伤集中在局部受载/脱离区域。

4. **本体/界面失效顺序支持什么？** bulk首次软化/失效为0.12025/0.12775 ms，interface为0.3225/0.449625 ms，首次拓扑脱离0.330125 ms。首本体失效两端距原自由表面约42.65/12.21 µm，原始邻居44/39（完整点阵122），m比值0.3363/0.3100。界面首次软化前的626个早期断键端点中95.85%处于自由表面一个δ内；其原始m均值比0.7301，全材料为0.7910。结合时空位置，支持受载表层本体损伤先行，而不是先整体界面剥离。未识别出可认证的连续内部裂面或单一裂纹前缘。表面载荷、几何集中与邻域缺失彼此混杂，不能直接判为数值表面错误。

5. **0.376 ms接触造成细碎化了吗？** 接受态已从无损归档恢复且SHA与原件完全一致。574/1429距离28.1995 µm，重叠0.4360095dx，排斥力12.6385 µN，均仍附着。固定参考邻域在首次非零接触前已断273条bulk键，在峰值前又已断1467条；峰值附近0.5 µs窗新增2条。这说明破坏早于峰值，不能将峰值直接认定为起因。574/1429局部仿射拟合detF约−0.1413/0.8813，前者提示局部折叠/强非仿射，不能解释成负物理密度；当前δ球点数比为1.000/1.022，也不构成连续体密度测量。强接触可能影响后续运动，尚缺接触单因素对照。没有改接触刚度。

6. **伸长历史有可量化Mode-II盲区吗？** 有可量化的“剪切较强但历史伸长未达起点”连接：已保存关键帧最多 {peak[0]} 条（step {peak[1]}，约0.330125 ms，含bulk/interface）。0.15剪切只是比较尺度，并非临界值；有限刚体转动诊断误差 {mixed['objectivity_max_abs_error']:.3g}。三节点候选的桥键局部星秩不足，剪切必须标 **NOT MEASURED**，且它已被现有拉伸判据识别为软化。因此没有证明大块被“漏检剪切桥”拴住，也未认证一个物理Mode-II失效盲区。完整ΔU_virtual经128次独立整模型重算核对，包含两端星耦合，不能等同Gc或独立单键耗能。

7. **Karmakar哪些可借鉴？** 已读取用户提供的18页全文并核验关键公式/图。可直接借鉴不可逆历史、预裂纹反力实验、本体与粘附分开校准、horizon和异质性限制的审查方法。NOSB/Ogden应力退化、界面弹簧刚度、长度比1.38/2.27及horizon缩放不能直接移植到当前完整LPS星能量；见 `LITERATURE_TO_CODE.md`。Zhang–Qiao参考连接方向在其算例中已知，本轮三维极分解指标不是其critical skew材料律。

8. **小试样形成了真实多节点自由块吗？** 两组实际结果见下表，均以存活bulk/interface、锚点可达性和夹具状态判定。两组最大自由分量为 {c0['maximum_free_nodes']} / {c1['maximum_free_nodes']} 节点；≥34节点状态为 {status['SYNTHETIC_CHUNK_FORMATION_STATUS']}。初始预裂纹是几何定义；计算期间只由原生历史自然断键，未指定最终碎块或合并粒子。夹具在2 ms释放时不删除材料键，释放前含夹具点的分量不计自由碎块。首50 µs裂尖比例、原裂纹面剩余桥接及空间分区见下文，不能仅把损伤云称连续裂纹。

9. **是否运行新真实V3？** 没有。小试样已检验现有模型的能力，但位移夹持与真实V3的空间局部表面力并不等价，反力/外功也不同；不能把方向差异直接当作真实V3的唯一因果证明。下一步法向对照还需冻结方向、符号及逐点力幅值一致的定义。数值一致性虽通过，但P5联合准入尚不满足。保留原52分钟基线，不重复运行。

10. **下一轮唯一主要改动？** {nextchange} 这是待证伪的加载路径假设，不是本轮已确认的真实微泡作用；不同时修改力幅值/空间包络、接触、界面或失效判据。现有模型已经在独立试样形成凝聚块体，因而本轮没有证据要求优先扩大软化参数或更换失效律。

## 实际小试样与计算证据

14×8×16材料点（1792）+三层固定锚点（336），dx=50 µm、δ=150.5 µm，原始星m保留未预切邻域；2736条初始预切不计自然损伤。两组唯一差别是顶面位移方向z/x；相同240 µm净位移、2 ms半余弦、2 ms释放夹具后再积分0.5 ms。所有材料参数、界面、接触、几何及积分准入相同。没有流场，也没有额外阻尼。是有限尺寸动态试样，不是准静态韧度实验。

| 案例 | 接受步 | 完成ms | 本体自然断键 | 界面断键 | 历史最大自由点数 | 最大相对R | 本次墙钟s |
|---|---:|---:|---:|---:|---:|---:|---:|
{table}

{chr(10).join(details)}

反力、外功和惯性随方向改变是实际响应；两组不是等能量输入的对照。`SMALL_SPECIMEN_VALIDATION.json` 保存每个真实自由分量的内部/外部连接及原始ID。能量账本为 W=ΔU+ΔK+D+ΔUc−Jc+R，R独立保留。每步比对Python完整星力/能量、κ/d、接触、VV自由DOF方程、原始身份与固定锚点，末态原生gather拓扑一致。全程时间步收敛和材料Gc仍未认证。

两组都首先在预裂纹尖端附近软化，但法向首50µs的损伤并非全部局限于裂尖，原z=0面也仍有跨面连接；自然分离路径不能等同于平直裂面已经全部贯通。法向块体保留高比例内部连接，足以反驳“现有材料律必然只能生成1–2点碎屑”。切向试样同样发生裂尖损伤和自然断键，故不能说当前伸长历史完全不响应剪切。定量模式韧度优劣仍缺独立标定基准。

小试样的残差分母是历史最大max(|W|, |ΔU|+K+D+Uc, 1e−15J)，两个方向一致；真实V3因存在多个正/负外功通道使用其原有累计尺度，不能仅横比两个百分比认证谁更精确。接触对计数包括r≈dx的近零力几何候选，不等于所有候选都发生有意义的碰撞。

原始history中的reaction_N是顶端点集内力/接触合力的负值，仅夹持期间是机械反力项；实际驱动器力还包含规定运动的惯性项，释放后驱动器力为0。修正语义后的力列另存derived_boundary_forces.csv，图中使用实际驱动器力，不将释放后的点集内力称为外力。外功账本一直包含规定DOF动能变化，未事后修改能量定义。

首次法向运行在2 ms释放夹具时因重复NVE积分fix被协调器拒绝，前缀保留于 `benchmark/`。修正只将释放后的自由材料点纳入同一NVE组；零步释放检查通过后，批量重放因浮点时间累计差异未通过严格续算门槛，保留于 `benchmark_corrected_v2/`。改为原逐条run 1顺序，在 `benchmark_corrected_v3/` 重放8001个接受步，x/v/总力/材料力/κ/d/时间/键数/K/Uc逐值一致后继续。scalar数组中的计时字段允许不同，物理累计D单独核对。两个失败版本的驱动源码按原冻结SHA完整保留。没有放宽续算门槛或把失败候选作为接受态。小试样只有两种冻结物理配置，接口失败/重放不是第三种物理对照。

实际两种对照及失败/重放的合计墙钟约{cost['total_benchmark_attempt_wall_s_approx']:.2f}s（{cost['total_benchmark_attempt_wall_s_approx']/60:.3f}min），其中原失败前缀{cost['original_failed_prefix_wall_s']:.3f}s、被拒绝的批量重放原生日志约{replay_failed}s；表内法向墙钟已含{replay['replay_wall_s']:.3f}s逐条重放，不能再重复相加。原生单MPI rank×单线程、无GPU，实测峰值RSS约{cost['benchmark_peak_RSS_bytes']/1024**2:.3f}MiB。未设置4GiB或1800s停止条件。离线图分析/空间分析/混合模式分别耗时{cand['wall_s']:.3f}/{comp['wall_s']:.3f}/{mixed['wall_s']:.3f}s；离线峰值RSS未测量。详见RUN_RESOURCES.json。

## 二次裂解补充诊断

13节点集合在0.55425 ms已脱离：内部相对动能 {g['K_internal_J']*1e9:.6g} nJ、星储能 {g['star_U_J']*1e9:.6g} nJ，K_internal/U≈{g['K_internal_J']/g['star_U_J']:.4g}。同一时刻基于力计算半步速度的载荷相对功率 {g['load_relative_power_W']:.6g} W，COM拖曳相对功率约0，附着数值阻尼已关闭。它解释了为什么“发生拖曳”并不等于“内部相对运动已被耗散”，但这些是采样状态的瞬时诊断，不是全部时程外功，更不证明修改拖曳即可得到大块。

## 数据完整性、文献、图形与限制

P1核对1414个原基线文件及完整原始ID；527个常规节点帧、17个全损伤关键帧、每接受步稀疏d事件与接触记录存在。普通逐节点力并非每接受步都保存；完整κ也不是每一步全量保存。额外5个已接受before-state从原无损拒绝试步归档恢复，绝不使用rejected trial充作真实轨迹。13节点形成的精确拓扑可从事件恢复，其形成瞬间完整逐节点力未声称存在。

图表位于对应 `results/{TASK}/{RUN}/index.html`，全部三维点/原始存活键、位移×1、固定相机和共同尺度；实际脱离、候选、离线诊断、合成试样分别标识。大图中固定支架与界面线省略仅为清晰度，不参与错误拓扑归类。无虚构连续裂面。图形审阅记录保留于 `FIGURE_QA.json`。

P0未删除仍用于冻结复现/构建的旧源码；新分析共享 `common.py`/`graph_ops.py`，原生桥复用冻结类。无C++或材料律修改。`PROTECTED_INPUTS.json`、原始清单重核及本轮SHA清单形成追溯链。原有G0理想瞬断标称值不作为本轮Gc；虚拟能量、D、R均未改名为材料断裂韧度。无临床微栓子、真实微泡耦合或AI训练资格声明。
'''
 (O/'REVIEW_V3_CHUNK_FRAGMENTATION_V1.md').write_text(review)
 for name in ['P0_CLEANUP_AUDIT.md']:
  assert (O/name).exists()
 for name in ['LITERATURE_TO_CODE.md','LITERATURE_INPUTS.json']:shutil.copy2(A/name,O/name)
 protected=json.loads((O/'PROTECTED_INPUTS.json').read_text());dump(O/'FINAL_SOURCE_IDENTITY.json',dict(native_library_sha256=sha(LIB),frozen_reference_model_sha256=sha(R/'analysis/softening_feasibility_v1/model.py'),current_analysis_files={p.name:sha(p) for p in A.iterdir() if p.is_file()},no_Cpp_change=True))
 print(json.dumps(status,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
