from common import *
import html

def pct(x):return '未定义（无脱离质量）' if x is None else f'{100*x:.5g}%'

def main():
    data={m:json.loads((O/p).read_text()) for m,p in zip(CASES,['SINGLE_VORTEX_NATIVE_RESULT.json','CHECKERBOARD_NATIVE_RESULT.json'])}
    status=json.loads((O/'FINAL_STATUS.json').read_text())
    workdiag=json.loads((O/'WORK_MATCH_DIAGNOSTICS.json').read_text())
    status['formal_case_names']={'counter_rotating':'CASE_B_WORK_MATCHED_COUNTER_ROTATING','checkerboard':'CASE_C_WORK_MATCHED_CHECKERBOARD_VORTICES'}
    for m,p in zip(CASES,['SINGLE_VORTEX_NATIVE_RESULT.json','CHECKERBOARD_NATIVE_RESULT.json']):
        data[m]['formal_case_name']=status['formal_case_names'][m]
        ob=data[m]['objective'];ob['total_detached_mass_kg']=sum(data[m]['final_free_sizes'])*RHO*V
        ob['largest_free_component']=max(data[m]['final_free_sizes'],default=0)
        ob['largest_free_fraction_of_detached']=ob['largest_free_fragment_fraction_detached'];ob['fines_1_2_mass_fraction']=ob['fine_debris_mass_fraction_detached']
        dump(O/p,data[m])
    topology=json.loads((O/'FRAGMENT_TOPOLOGY_HISTORY.json').read_text())
    for m in CASES:topology[m]['final_metrics']=data[m]['objective']
    dump(O/'FRAGMENT_TOPOLOGY_HISTORY.json',topology)
    status.update(scientific_target='EXPLORATORY_TARGET_NOT_MET',work_comparison='NOT_AN_EQUAL_WORK_COMPARISON',surface_repair='PARTIAL_GEOMETRIC_PROXY; true full wetted surface and fine-fragment hydrodynamics remain unresolved')
    dump(O/'FINAL_STATUS.json',status)
    def each(fn):return '\n\n'.join(('B 单对反向涡旋' if m==CASES[0] else 'C 棋盘涡旋')+'：'+fn(q) for m,q in data.items())
    rows=[]
    for m,q in data.items():
        w=q['work'];ob=q['objective'];s=q['native']
        rows.append(f"|{'B 单对' if m==CASES[0] else 'C 棋盘'}|{s['status']}|{w['actual_net_J']*1e9:.6f}|{w['relative_error']*100:+.3f}%|{ob['n_fragments_20_200']}|{max(q['final_free_sizes'],default=0)}|{pct(ob['fine_debris_mass_fraction_detached'])}|{100*s['max_relative_energy_residual']:.5g}%|")
    answers=[]
    def add(title,body):answers.append(f'### {len(answers)+1}. {title}\n\n{body}\n')
    add('是否实现新开放裂面的流体可达性近似？','实现了 QUASISTATIC_FLUID_ACCESSIBILITY_PROXY。实际点阵分成9828个辅助四面体、20558个唯一三角形。六条原始键存活、当前体积方向正常的四面体提供材料支撑；失去支撑本身不代表流体空间。局部连接已断开且体积增量对应开口>1 µm的单元，才可作为间隙；经断开且面积张开的三角通道连到外部，并检查当前实体遮挡。新面必须有存活四面体支撑、对侧实际节点法向距离增加>1 µm。1–5 µm采用平滑开口权重。辅助单元不参与材料求解，不改变节点质量或任何材料键。')
    add('多少新裂面获得受力表示？',each(lambda q:f"累计首次受载的新内部三角面 {q['surface']['unique_ever_loaded_new_faces']} 个，其中自然新断裂形成的面 {q['surface']['unique_ever_loaded_natural_new_fracture_faces']} 个、原预切口重新张开暴露的面 {q['surface']['unique_ever_loaded_reopened_precut_faces']} 个；瞬时最多 {q['surface']['max_accessible_new_faces']} 个，最终 {q['surface']['final']['accessible_new_faces']} 个；最大新增可达面积 {q['surface']['max_new_area_m2']:.8g} m²。这是有节点与键支撑的离散代理面数，不是独立连续裂纹条数。"))
    add('原始面损失后是否仍有表面承担载荷？',each(lambda q:f"末态原始可达面 {q['surface']['final']['original_faces']}，新可达面 {q['surface']['final']['accessible_new_faces']}。本轮包含有支撑的侧面和顶面，初始面积 {q['work']['initial_area_m2']:.8g} m²，与旧顶面182片/2.275e−7 m²不同；不能用旧末态87/182=47.802%作为新模型的覆盖率分母。"))
    add('真实游离块是否继续受力？',each(lambda q:f"最终 {q['surface']['final']['true_free_components']} 个真实游离分量中，{q['surface']['final']['free_components_with_surface']} 个有已表示表面，合计 {q['surface']['final']['free_surface_faces']} 面。按当前位置、法向及相对速度计算牵引；脱离不是关闭受力的条件。离开驱动窗口或2.8 ms后停止主动驱动，已表示面仍有被动阻力。若该例有表面的游离分量数为0，说明正式轨迹未产生可加载的三维自由块，不能把896点已知游离夹具的通过冒充主试样中的观测。"))
    add('还有多少未表示的流体可达表面？',each(lambda q:f"最终闭合／未判定候选面 {q['surface']['final']['closed_unresolved_faces']}，候选面积中未判定比例 {pct(q['surface']['final']['unresolved_candidate_area_fraction'])}；没有已表示面片的游离分量 {q['surface']['final']['free_components_without_surface']} 个、共 {q['surface']['final']['free_nodes_in_components_without_surface']} 节点。真实连续表面缺失比例为 UNDETERMINED，不能拿候选面积比例冒充。1–2点细屑没有假闭合外壳，也未加点粒子阻力，其实际水动力仍未解决。"))
    add('是否引入重复计力？','每个排序后三节点集合只有一个几何面ID；原始面与新内部面互斥。游离标签是同一面的拓扑属性，不另加一份力。三角形重心力按正1/3分配到原始节点，关键态合力、当前力矩及伴随功率逐项验证。新暴露面的工作量含原预切口张开与自然新断裂，不能全部归于新生成裂纹；它已包含在总外功内。')
    add('单对涡旋是否匹配历史工作预算？',f"{data[CASES[0]]['work']['status']}。\n\n"+each(lambda q:f"实际净外功 {q['work']['actual_net_J']*1e9:.9g} nJ，驱动 {q['work']['drive_J']*1e9:.9g} nJ，被动阻力 {q['work']['drag_J']*1e9:.9g} nJ；相对25.698197 nJ偏差 {q['work']['relative_error']*100:+.6g}%。"))
    add('棋盘是否匹配相同外功？',f"{data[CASES[1]]['work']['status']}，实际偏差 {data[CASES[1]]['work']['relative_error']*100:+.6g}%。匹配同时要求完整3.6 ms协议和绝对相对误差≤5%。两个新模式采用相同表面算法和控制规则；历史A只有工作预算意义。")
    add('若未匹配，原因及限制是什么？',each(lambda q:f"状态 {q['work']['status']}；幅值上下界触发 {q['work']['cap_saturated_steps']} 个接受步，变化率限制触发 {q['work']['slew_limited_steps']} 步。冻结上限6000 Pa、变化率8e6 Pa/s、固定空间窗口和2.8 ms停止驱动共同限制补功。控制器对后段阻力使用瞬时线性衰减代理；实际破坏、暴露面积、速度及非线性弹性恢复可偏离预测。当驱动功率/当前幅值≤1e−14时，冻结控制器退回3200Pa而非反解非正增益；这种规则也可能阻碍补功，不把控制器看成最优。各限制是否主导须结合 ACTUAL_LOADING_HISTORY.csv 和 WORK_MATCH_DIAGNOSTICS.json，不从偏差直接推断微泡能力。"))
    add('实际峰值剪切和冲量？',each(lambda q:f"峰值实际面剪切 {q['work']['actual_peak_shear_Pa']:.8g} Pa；∫Σ|面片切向力|dt={q['work']['facet_shear_impulse_Ns']:.8g} N·s；节点力模长积分 {q['work']['node_force_norm_impulse_Ns']:.8g} N·s；峰值节点力 {q['work']['peak_node_force_N']:.8g} N；最大表面速度 {q['work']['max_surface_speed_m_s']:.8g} m/s。模长积分不等于净矢量冲量。"))
    add('棋盘实际有多少单元？','实际顶面650×350 µm，冻结方格尺度250 µm。9个连续流函数单元与顶面相交，仅3个涡心位于顶面内，只有中间1格完整，另8格被边界裁切；不能称为9个完整涡旋。中心、尺寸、裁切面积见 SHEAR_PARAMETERIZATION.json。该文件的 rotation_sign 是流函数叶瓣符号：+Z法向下，+号对应顺时针，−号对应逆时针；相邻格交替。有限区域平滑窗口和当前法向影响实际速度。')
    add('裂纹空间分布是否不同？',each(lambda q:f"最终自然本体失效 {q['native']['last_row']['bulk_failed']:.0f} 键、界面 {q['native']['last_row']['interface_failed']:.0f} 键；自然失效键最终区域 {q['damage_regions']['final_regions']} 个。Figures 2/4提供实际节点损伤和开口代理面，未将所有失效键画成连续裂面。"))
    add('是否有多个独立损伤区域？',each(lambda q:f"自然完全失效键共享端点的区域诞生 {q['damage_regions']['birth_count']} 次，最多同时 {q['damage_regions']['max_concurrent_regions']} 个，合并 {q['damage_regions']['merge_count']} 次。这个指标不等同于贯通裂面或游离材料分量。"))
    add('是否观察到真实二次碎裂？',each(lambda q:f"真实游离父分量分裂事件 {q['secondary_count']} 次。父子component ID、原始成员及发生时间见 SECONDARY_FRAGMENTATION.json；仅按实际存活材料图进一步分裂统计；小细屑的拆分不等同于凝聚大块受流体剪切后的二次分裂。"))
    add('是否形成20–200节点凝聚块？',each(lambda q:f"历史最多 {q['native']['max_medium_components']} 个，最终 {q['objective']['n_fragments_20_200']} 个；34–200节点 {q['objective']['n_fragments_34_200']} 个，≥114节点 {q['objective']['n_fragments_ge_114']} 个。总脱离质量={q['objective']['total_detached_mass_kg']:.8g} kg；目标质量/初始质量={pct(q['objective']['mass_fraction_target_initial'])}，目标质量/脱离质量={pct(q['objective']['mass_fraction_target_detached'])}。内部存活键、原始成员、三维形态及寿命逐分量保存。"))
    add('是否达到至少5块及完整标准？',each(lambda q:f"最终目标块 {q['objective']['n_fragments_20_200']} 个，3.4–3.6 ms相同成员保持连通的目标块 {q['objective']['cohort_target_survivors_to_end']} 个；完整标准={q['objective']['full_target']}。还要求最大块≤脱离质量50%、细屑≤20%。本轮不额外加入旧版10%脱离质量门槛。"))
    add('细屑比例是否下降？',each(lambda q:f"细屑/脱离质量={pct(q['objective']['fine_debris_mass_fraction_detached'])}，细屑/初始质量={pct(q['objective']['fine_debris_mass_fraction_initial'])}。历史A为4.8309%、旧B为100%，但表面模型不同，跨轮差异不能单独归因于涡旋数量。"))
    add('最大块是否仍占主导？',each(lambda q:f"最终最大块 {max(q['final_free_sizes'],default=0)} 节点，占脱离质量 {pct(q['objective']['largest_free_fragment_fraction_detached'])}；历史最大 {q['native']['history_max_free_nodes']} 节点。全部最终尺寸：{q['final_free_sizes']}。仍存活的非局部材料连接保留大块整体性，没有进行后处理拆块。"))
    add('新裂面受力有何可解释贡献？',each(lambda q:f"新增暴露面（含原预切口重新张开与自然新断裂）的合计实际净外功 {q['surface']['contribution']['actual_net_work_on_new_accessible_faces_J']*1e9:.8g} nJ，面片状态提供开放及加载直接证据。没有运行同一模型关闭新面受力的额外对照，因此不能定量宣称碎裂变化中有多少由新面单独造成。"))
    add('是否支持少量BoTorch优化？',f"判定：{status['optimization_admission']}。非细屑主导的目标带块可提供后续有限研究信号，仍需审查等功偏差及未表示表面。本轮优化trial=0，不自动启动下一轮或增加新物理工况。")
    add('真实微泡声流可实现性有哪些限制？',each(lambda q:f"实际最大外缘速度 {q['work']['max_outer_speed_m_s']:.7g} m/s，Re_h最大 {q['work']['max_Re_h']:.7g}。")+'\n\nμ=.001 Pa·s、h_eff=25 µm的薄层滑移律是机械代理；高Re不是低Re Stokes认证。每个有向面使用n×∇ψ，不保证对应一个全域不可压缩三维流场。没有泡一阶/二阶流、渗润时程或双向FSI。支撑单元、局部面积/体积开口、有限遮挡探针和参考邻接连通都有分辨率局限；未重建表面不能称真实零水动力。')
    add('两次正式协议成本？',each(lambda q:f"墙钟 {q['native']['wall_s']:.3f} s，CPU {q['native']['cpu_s']:.3f} s，峰值RSS {q['native']['peak_RSS_bytes']/2**20:.3f} MiB；接受步 {q['native']['accepted_steps']}，实际时间 {q['native']['time_s']*1e3:.6g} ms。")+'\n\n每例1 rank×1线程，两例并行，墙钟之和不是项目经过时间。新路径微步、初始化故障及后处理另列。没有4 GiB或1800 s限制。首次初始化因统计未排除锚点，各留下一个未接受微步且没有损伤变化；另一次哈希防护在动力学开始前发现生成diff的冻结顺序问题。失败目录和日志完整保留；恢复的是同一两组冻结协议，没有利用结果重新标定。')
    text=f'# CLOT_SHEAR_MULTIVORTEX_FRAGMENTATION_V2\n\nRun {RUN}。科学标签 **{LABEL}**。原生协议状态 **{status["status"]}**。\n\n'
    text+='本轮未实现等外功对照，也未达到多块碎裂目标。已实现的表面路径是部分、保守的几何代理；实际游离细屑的水动力仍未解决。当前不支持启动BoTorch。\n\n'
    text+='|工况|原生状态|净外功 nJ|目标偏差|最终20–200块|最大块|细屑/脱离质量|最大能量余额|\n|---|---|---:|---:|---:|---:|---:|---:|\n'+'\n'.join(rows)+'\n\n'+'\n'.join(answers)
    text+='\n## 工作预算诊断\n\n'+'\n\n'.join(f'{m}：实际保持6000Pa上限 {d["upper_cap_duration_s"]*1e3:.6g} ms；驱动功率为负 {d["negative_drive_power_duration_s"]*1e3:.6g} ms；瞬时单位幅值增益≤1e−14约 {d["nonpositive_or_tiny_unit_power_duration_s"]*1e3:.6g} ms。2.8–3.6ms实际净功变化 {d["off_to_end_actual_work_J"]*1e9:.7g} nJ，停驱时预测拖曳损失 {d["off_state_drag_forecast_J"]*1e9:.7g} nJ。' for m,d in workdiag.items())+'\n'
    text+='\n\n两组在2.8ms停驱前已显著欠功。卸载段实际阻力损失与预测仅相差约0.008nJ，远小于11.76–16.35nJ的最终预算缺口，所以该预测误差不是主要原因。幅值上限、非正瞬时驱动增益、固定空间作用域/包络与破坏后的表面变化限制了本次控制；没有额外对照可分摊各因素的因果比例。不能据此比较等功条件下哪种涡旋更优。下一步应先审查表面覆盖及可达的共同工作预算，不盲目增大牵引。\n'
    text+='\n## 数值与数据边界\n\n1792材料点+336锚点，95666原始键/2736预切，dx50 µm、horizon150.5 µm、G940 Pa、K100000/3 Pa、ρ1000，bulk .15/.18、interface .140625/.16875均冻结。native库SHA '+sha(LIB)+'。没有改材料、随机弱区、人为新增主试样裂纹或接触连通语义。\n\n'
    text+='SUMMARY内继承的field_parameters.nominal_shear_Pa=1000只是未用于实际幅值的基函数默认值，其静态applicability上界不能代表本轮载荷；实际幅值由冻结控制器给出，逐接受态记录，以上实际速度和Re优先。\n\n'
    text+='原生VV与第二半步隐式一致速度保持。预测Δd>0.1时在提交前缩步；几何重算不改历史，累计曝光仅在接受后提交。若提交后检查失败则输出failed_unaccepted和last_accepted并停止，不声称原生状态已回滚。每接受步保存标量、损伤变更和真实分量记录；关键时刻及约50 µs间隔保存FP64 x/v/κ/d/节点力/面片状态，并非每个接受步都保存完整坐标数组。\n\n'
    text+='外功以接受端点真实力平均点乘实际native位移，驱动/阻力分开。W=ΔU+ΔK+D+ΔUc−Jc+R，R独立保留，不混入损伤耗能。源码哈希、稳定ID、精确锚点、不可逆损伤、能量恒等式、关键态力重算均独立读回。参考态几何代理未认证为唯一连续裂面或真实润湿面积。\n\n'
    text+='[科研图表](../../../results/'+TASK+'/'+RUN+'/OPEN_RESULTS.html) · [数值验证](DELIVERY_VALIDATION.json) · [重现说明](REPRODUCE.md)。\n'
    (O/'REVIEW_CLOT_SHEAR_MULTIVORTEX_FRAGMENTATION_V2.md').write_text(text)
    reproduce=f'''# Reproduce {TASK}

Accepted run {RUN} is immutable. Use a NEW task/run directory; never overwrite its results.

Activate /workspace/y_clot/env.sh. Set single-thread BLAS/OpenMP and LD_LIBRARY_PATH to the frozen library directory. Library SHA: {sha(LIB)}.

Chain: shear_field_library.py → current_surface.py → native.py → run_cases.py → analyze.py.
work_controller.py changes one global amplitude only.

Change isolated RUN before reproducing. Copy exact original inputs; run prepare.py, checks.py, native_checks.py. Freeze and read back every SHA. Geometry fixtures are not fragmentation evidence; native checks take one microstep per mode.

Run the TWO specified modes once each under supervisor, from the identical original state. No additional uniform case or outcome-based calibration trials.

Read-only postprocessing: python analyze.py; python figures.py; python report.py.

Initialization failures and source records are retained in startup_failure_00/ and preflight_failure_01/. Do not hide process attempts.

Full accepted scalar/damage/component records and FP64 key states are retained. Work is actual endpoint force dot displacement; no nominal-area substitution or curve scaling. Unknown true surface area remains unknown.
'''
    (O/'REPRODUCE.md').write_text(reproduce);(A/'REPRODUCE.md').write_text(reproduce)
    images=sorted(F.glob('Fig*.png'));body='\n'.join(f'<section><h2>{html.escape(p.stem)}</h2><a href="{p.stem}.svg"><img src="{p.name}" alt="{html.escape(p.stem)}"></a></section>' for p in images)
    index='<!doctype html><html lang="zh"><meta charset="utf-8"><meta name="viewport" content="width=device-width"><title>Exposed surface / multivortex shear</title><style>body{font:17px system-ui;max-width:1500px;margin:35px auto;padding:0 24px;color:#21333f;background:#f9fbfc}img{width:100%;background:white}section{margin:35px 0}a{color:#17658c}</style><h1>当前可达裂面与等外功多涡旋剪切</h1><p>SHEAR_MECHANICAL_REQUIREMENT_ONLY · 真实原生状态 · 位移×1 · 连续表面完整性未认证</p><p><a href="../../../verification/'+TASK+'/'+RUN+'/REVIEW_CLOT_SHEAR_MULTIVORTEX_FRAGMENTATION_V2.md">22问报告</a> · <a href="../../../verification/'+TASK+'/'+RUN+'/WORK_MATCH_RESULT.csv">实际功匹配</a></p>'+body+'</html>'
    (F/'OPEN_RESULTS.html').write_text(index)
    (O/'OPEN_RESULTS.html').write_text('<!doctype html><meta charset="utf-8"><meta http-equiv="refresh" content="0;url=../../../results/'+TASK+'/'+RUN+'/OPEN_RESULTS.html"><a href="../../../results/'+TASK+'/'+RUN+'/OPEN_RESULTS.html">Open figures</a>')
    print('REPORT_COMPLETE')

if __name__=='__main__':main()
