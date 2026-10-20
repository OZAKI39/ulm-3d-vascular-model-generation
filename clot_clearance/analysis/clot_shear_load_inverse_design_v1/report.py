"""Chinese research record generated from accepted evidence, no invented success."""
from analyze import *
import html
def fmt(v,scale=1):return '未观察到／无定义' if v is None else f'{v*scale:.8g}'
def main():
    q=read(O/'ANALYSIS.json');summ=read(O/'NATIVE_SHEAR_CASE_RESULTS.json');status=read(O/'FINAL_STATUS.json');bridge=read(O/'CRACK_BRIDGE_AUDIT.json')['cases'];near=read(O/'NEAR_SEPARATION_AUDIT.json');res=read(O/'RUN_RESOURCES.json');surface=read(O/'CURRENT_SURFACE_TRAJECTORY_AUDIT.json');valid=read(O/'MODEL_VALIDITY_WINDOWS.json');checks=read(O/'KEY_CHECKS.json')
    coverage=read(O/'FRAGMENT_SURFACE_COVERAGE_AUDIT.json')['cases'];finalcoverage={k:next(z for z in v if z['key']=='final') for k,v in coverage.items()}
    names={'uniform':'A 单向随动剪切','counter_rotating':'B 双涡旋随动剪切'};answers=[]
    def add(question,answer):answers.append(f'### {len(answers)+1}. {question}\n\n{answer}\n')
    def each(fn):return '；'.join(names[k]+'：'+fn(k,v) for k,v in q.items())
    add('是否实现两种实际可运行的切向剪切场？','是。`ShearField` 提供单向场与相反符号双 Gaussian 流函数场，均实际进入 native LAMMPS 两组计算。7个参数为名义应力上限、平面方向、中心x/y、宽度、涡间距、驱动时长；当前冻结1000 Pa、0 rad、(0,0)、125 µm、300 µm、2.8 ms。单向场不使用涡间距/宽度，固定共同作用区域；本轮没有进行优化搜索。')
    add('是否根据当前坐标和法向更新？','是。182个原始顶面三角形、112个原始材料ID，每次外力回调重新计算质心、面积、单位法向和顶点平均速度。τ=(µ/h)Pₜ(u_outer−v_surface)，µ=.001 Pa·s、h=25 µm，取自既有独立Stokes模板，非泡模型校准。u_outer是在假设剪切层外缘的速度，不是把无滑移界面速度差误当速度梯度。位置使用实际x(t)，位移显示倍率1。')
    add('哪些新裂面具备流体可达性？','本轮没有证明任何新裂面的流体可达性，统一 `NEW_CRACK_FACE_WETTING_UNRESOLVED`。仅跟随原始顶面；三条原始边键仍存活、面积/质量非退化的三角形保留。失去支撑的面片退出载荷离散，不改变材料键。保留规则也不构成流体可达性的证明。其他原外表面和新裂面未构建；没有保留面片的碎片，其水动力是未解决，而非真实零力结论。'+each(lambda k,v:f'最终保留原参考受载面积{v["final"]["active_reference_area_fraction"]*100:.5g}%，有保留面片的游离分量{finalcoverage[k]["free_components_with_retained_faces"]}个、完全无保留面片的分量{finalcoverage[k]["free_components_without_retained_faces"]}个/共{finalcoverage[k]["free_nodes_in_components_without_retained_faces"]}节点'))
    add('是否修复局部映射方向问题？','新小试样路径使用面片原始顶点和正1/3权重，没有V3有符号外插；离散总力、当前力矩、伴随虚功精确匹配。实际关键态最大守恒误差：'+each(lambda k,v:f'力{max(z["force_error_N"] for z in surface[k]):.4g} N、力矩{max(z["moment_error_Nm"] for z in surface[k]):.4g} N·m、功率{max(z["power_error_W"] for z in surface[k]):.4g} W、局部力/正权重预算≤{max(z["local_force_vs_positive_budget_max"] for z in surface[k]):.6g}')+'。初始平面双涡旋节点合力相对解析节点方向最大差4.6266°（只统计解析幅值>峰值5%）；相邻面片平均产生离散误差，非网格收敛证明。节点分摊面积仍有2/3–4/3差异。旧1927点V3仿射路径未修改、未通过准入，不能称其历史问题已经普遍修复。')
    add('两种模式实际输入多少功？',each(lambda k,v:f'W={v["final"]["W_J"]*1e9:.8g} nJ，其中外缘驱动{v["final"]["W_drive_J"]*1e9:.8g} nJ、被动拖曳{v["final"]["W_drag_J"]*1e9:.8g} nJ；实际相对滑移面剪切峰值{v["peak_shear_Pa"]:.7g} Pa')+'。1000 Pa约束的是静止面的驱动上限，不是两例相同的实际相对剪切峰值或外功。2.8 ms后驱动严格为0，保留的原面片仍受速度相关被动拖曳。')
    eventlines=[]
    for k in CASES:
        parts=[]
        for ev in ['first_softening','first_failure']:
            e=summ[k]['events'].get(ev)
            parts.append(ev+'未观察到' if e is None else f'{ev}={e["time_s"]*1e3:.8g} ms，原始bond IDs={e["original_bond_ids"]}，参考中点µm={np.round(np.array(e["reference_midpoints_m"])*1e6,5).tolist()}，当前中点µm={np.round(np.array(e["actual_midpoints_m"])*1e6,5).tolist()}')
        eventlines.append(names[k]+'：'+'；'.join(parts))
    add('首次软化和首次断裂位置是否不同？','\n\n'.join(eventlines)+'。真实原始ID及局部仿射伸长诊断另见 FRACTURE_EVENTS.json / LOCAL_DEFORMATION_AUDIT.json。')
    add('是否产生不同的空间裂纹网络？',each(lambda k,v:f'自然bulk/interface失效网络诞生{v["inclusive_failed_material_networks"]["n_births"]}次、最多同时{v["inclusive_failed_material_networks"]["max_simultaneous_networks"]}个；bulk诞生中心最大间距{v["crack_networks"]["max_birth_centroid_separation_m"]*1e6:.7g} µm；最终本体/界面失效{int(v["final"]["bulk_failed"])}/{int(v["final"]["interface_failed"])}条')+'。网络由自然失效原始边共享端点定义，初始2736预切不计；锚端点保留但不冒充材料点。')
    add('双涡旋是否促进多区域裂纹贯通？',each(lambda k,v:f'失效网络合并{v["inclusive_failed_material_networks"]["n_merge_events"]}次，真正游离{int(v["final"]["detached_nodes"])}节点/{int(v["final"]["components"])}分量')+'。合并数量不是三维闭合断面或凝聚块的证据；只能结合下述实际存活图判断。当前两个固定场不能证明所有双涡旋参数都优于单向场。')
    add('是否产生20–200节点真实游离块？',each(lambda k,v:f'历史最大目标带数量{v["history_max_N20_200"]}、最终{int(v["final"]["N20_200"])}；历史最大游离分量{v["history_max_free_nodes"]}节点')+'。必须通过原始存活材料键不能到达任何固定锚点，contact或几何接近不增加材料连通。尺寸、参考体积等效直径、三轴跨度、内部连接密度和寿命均保留。20–200点对应约168.4–362.8 µm参考体积等效直径，仅为离散工程探索尺度，非医学阈值。')
    add('是否达到至少5个目标碎块？',each(lambda k,v:f'全部冻结探索标准通过={v["exploratory_target_met"]}；3.4 ms目标成员集合至3.6 ms仍连通的数量={v["common_cohort_medium_surviving"]}')+'。还要求脱离质量≥10%、最大块≤脱离质量50%、1–2点细屑≤20%。驱动卸载后的观察保留拖曳，寿命为有限窗口右删失，不代表永久稳定。')
    add('1–2节点细屑占比是否下降？',each(lambda k,v:f'最终细屑/脱离质量={fmt(v["final_metrics"]["fine_debris_mass_fraction_detached"],100)}%，细屑/初始质量={v["final_metrics"]["fine_debris_mass_fraction_initial"]*100:.7g}%')+'。无脱离时比例为null，不能写成0%改善。此前死载反向剪切为91.8919%细屑，但空间场、速度闭合和包络均不同；不能将跨轮差异单独归因于表面随动。')
    add('是否仍有跨区非局部连接？','\n\n'.join(names[k]+'：'+ '；'.join(f'{r["region"]} {r["candidate_nodes"]}节点，存活跨界{r["surviving_crossing_bonds"]}键，其中历史未达起点{r["historically_below_onset"]}、当前压缩{r["currently_compressive"]}、长度≥2dx {r["long_bonds_ge_2dx"]}，锚可达{r["anchor_reachable_nodes"]}节点' for r in bridge[k]['final']['regions']) for k in CASES)+'。A的诊断子区域虽不再锚可达，仍由存活跨界键连到更大的游离块，不能把它们另计为独立碎块；B仍有大量锚可达节点。三个历史诊断区域不是载荷掩码或预切碎片。CSV逐条保留原始bridge ID和状态，额外保存真实三维锚路径与固定源端最小割代理；不执行删键。')
    add('是否有三维凝聚的近分离候选？','\n\n'.join(names[k]+'：'+ '；'.join(f'{stage}审查{len(v["candidates"])}个去重诊断分区，满足冻结近分离代理门槛{v["passing"]}个，其中当前三维点云{sum(c["near_separation_proxy_pass"] and c["current_shape"]["shape_class"]=="THREE_DIMENSIONAL_POINT_CLOUD" for c in v["candidates"])}个，进一步满足附着、20–200点、内部连通和三维形态的候选{v['attached_target_band_3D_candidate_count']}个' for stage,v in near[k].items()) for k in CASES)+'。代理门槛：至少3点且不超过半体、参考bbox占比≥.2、内部存活≥80%、边界完整性损失≥50%且超过内部损伤、加权边界/度≤.15。完整原始成员、厚度/形状、内部存活和锚可达均保存。它不是目标碎块证据，也不是唯一可能的候选搜索。')
    add('结果是否支持继续优化方向、位置及时序？',f'当前证据分级：`{status["OPTIMIZATION_READINESS"]}`。优先比较实际目标带数量和真实近分离结构；仅κ、断键数或单表面细屑不作为有价值的多块目标信号。本轮保留7维有界参数接口和既有优化器环境，但不因接口存在就追加搜索。')
    a_peak=next(z for z in surface['uniform'] if z['key']=='peak_load');b_peak=next(z for z in surface['counter_rotating'] if z['key']=='peak_load')
    add('未达目标的主要原因更可能是什么？','证据须分层：单向与涡旋加载改变不同区域的正伸长，存活三维桥接限制真正自由；表面边键失效还会降低载荷覆盖。'+each(lambda k,v:f'表面支撑首次损失={fmt(valid[k]["surface_support_loss"]["time_s"] if valid[k]["surface_support_loss"] else None,1e3)} ms、末态受载原面积{v["final"]["active_reference_area_fraction"]*100:.5g}%')+f'。在{HOLD*1e3:g}ms，A保留{a_peak["active_surface_faces"]}面片中{a_peak["retained_faces_outside_drive_region"]}片在驱动区外，驱动节点力模长和{a_peak["sum_drive_node_force_N"]:.7g}N、拖曳{a_peak["sum_drag_node_force_N"]:.7g}N。B保留{b_peak["active_surface_faces"]}面片、区外{b_peak["retained_faces_outside_drive_region"]}片。这是按当前空间重算，不是脱离开关。内点无需直接成为受载顶点：A history中{int(q["uniform"]["final"]["free_nodes_without_resolved_surface"])}个游离非表面顶点包含受力大块内部节点，只有{finalcoverage["uniform"]["free_components_without_retained_faces"]}个无面片分量/共{finalcoverage["uniform"]["free_nodes_in_components_without_retained_faces"]}节点完全缺少已表示作用面。相对速度拖曳、空间离开加载区和当前法向都会降低/改变实际剪切，不再按死牵引无限加速。仅两例无法唯一确定“伸长型本构不能表征剪切断裂”；混合模式本构需要独立实验/校准证据，不能为达标直接替换。')
    add('是否启动下一阶段BoTorch？',f'本轮新增trial=0，正式原生工况恰好2组。就目前保存证据，准入状态为 `{status["OPTIMIZATION_READINESS"]}`。未识别成功需求或可靠近分离信号时停止扩大优化；若出现信号也应先检查失去表面支撑导致的模型局限，再定义少量有界后续实验。未运行20–40次扫描。')
    add('载荷属于机械需求还是真实微泡声流？','所有结果标记 `SHEAR_MECHANICAL_REQUIREMENT_ONLY`。1000 Pa、µ=.001 Pa·s、h=25 µm给出25 m/s，Reₕ上限625、宽度尺度Re上限3125，不能解释为原4 mm/s低速Stokes流。双Gaussian平面流函数核心区在初始平面无散；固定区域平滑截断及变形曲面上只是局部切向运动学近似，没有求解三维不可压缩流体边值问题。新加载虽解决跟随与相对速度，未认证薄层近似在高Re和大变形下成立。')
    add('对接声学可实现性还需哪些真实微泡数据？','需要泡与血栓/管壁的统一坐标和动态间隙、泡数量/尺寸/壳层及相位、声频和复声压/速度、径向及非球形振荡、平移速度、介质µ/ρ、壁面边界及空间分辨率；再求/验证近壁一阶速度梯度与二阶平均流/边界驱动。输出当前表面外缘速度或应力张量、压力、时空插值与有效域，并提供动量/功率和来源校验。声压、R(t)或泡质心速度不能直接替代这些量。')
    add('何时可迁移到真实1927节点V3？','至少先在小样形成可解释的目标块或近分离结构，解决实际局部表面映射方向及精度、当前原面/新裂面和碎片流体覆盖、真实载荷输入与适用范围，并确认时间步/材料/空间分辨率足以支持目标。V3原有有符号仿射映射未通过本轮准入；本轮不运行完整V3、不重建CFD、不改native材料内核。')
    add('实际计算时间和内存？',each(lambda k,v:f'墙钟{summ[k]["wall_s"]:.3f} s、CPU {summ[k]["cpu_s"]:.3f} s、峰值RSS {summ[k]["peak_RSS_bytes"]/2**20:.3f} MiB、接受步{summ[k]["accepted_steps"]}')+f'。两组并行，每例1 rank×1线程；累计算例墙钟{res["native_wall_s_sum"]:.3f} s不是项目经过时间，CPU合计{res["native_cpu_s"]:.3f} s；单例最大RSS {res["max_single_case_peak_RSS_bytes"]/2**20:.3f} MiB，两例各自峰值和上界{res["sum_individual_peak_RSS_upper_bound_bytes"]/2**20:.3f} MiB，非同步测得峰值。必要新接口检查另计1.607 s；全部已记录的后处理与重绘成本、CPU及RSS见 AGGREGATE_COMPUTE_COST.json。该汇总不包含未计时的环境读取、网络传输或人工审图，不伪称项目端到端CPU总量。没有4 GiB或1800 s硬限制。')
    intro=f'''# CLOT_SHEAR_LOAD_INVERSE_DESIGN_V1

Run `{RUN}`。两种当前表面随动剪切的原生小试样对照。科学标签：**{LABEL}**。目标状态：**{status['MULTIFRAGMENT_FORMATION_STATUS']}**。

{status['SCIENTIFIC_CONCLUSION']}

| 工况 | 首次软化(ms) | 首次断裂(ms) | 最终游离节点 | 历史最大块 | 历史最大20–200块数 | 外功(nJ) | 最大能量余额(%) |
|---|---:|---:|---:|---:|---:|---:|---:|
'''
    for k,v in q.items():intro+=f'|{names[k]}|{fmt(v["first_softening_s"],1e3)}|{fmt(v["first_failure_s"],1e3)}|{int(v["final"]["detached_nodes"])}|{v["history_max_free_nodes"]}|{v["history_max_N20_200"]}|{v["final"]["W_J"]*1e9:.8g}|{v["max_relative_energy_residual"]*100:.6g}|\n'
    tail=f'''
## 冻结与数值边界

1792材料点、336锚点、95666原始键、2736原预切。G=940 Pa、K=100000/3 Pa、ρ=1000 kg/m³、dx=50 µm、δ=150.5 µm；bulk s0/sc=.15/.18，interface .140625/.16875，均匀材料及接触参数未修改。原生库SHA `{sha(LIB)}`，实际版本20260930。材料尚未针对本任务的剪切断裂实验标定。

dt=0.25 µs、事件最小0.03125 µs，预测Δd>0.1的步长在native提交前缩减。沿用原生Velocity-Verlet，新增速度相关外力在第二半步中以收敛固定点解隐式一致速度；每接受步验证预测x、κ/d、VV速度、力端口和稳定ID。新路径两种模式各一次微小原生步与同态restart力重算通过，没有重复旧MPI/material/full-restart矩阵。

W=W_drive+W_drag=ΔU+ΔK+D+ΔUc−Jc+R。外功以实际native drift位移积分，面片边失效根分段，Gauss4/8和端点功率积分差均保留。没有将数值R归入损伤耗能，没有额外数值阻尼。当前x×F力矩和角冲量逐步更新。1%能量准入只针对这一数值模型，不是高Re薄层物理认证。

几何诊断分别记录位移>0.1h、法向转角>15°、面积倍率超[.5,2]、支撑损失及Re>1首次时刻，见 MODEL_VALIDITY_WINDOWS.json。它们是保守近似诊断而非已校准统一有效窗；触发后不隐藏数据，也不把运行完成视为物理有效。原外表面支持丢失后，不能从未表示裂面的零外力推出真实水动力为零。

## 文献与用户讨论材料

- 已完整读取用户 `/home/lzy/projects/FEM网格claude.pdf`（3页）。采纳先研究剪切需求的思路；本轮明确要求优先，未采用其弱格子、12组筛选、20–40轮扫描或完整V3建议。
- [Mobadersany & Sarkar公开论文](https://arxiv.org/abs/1802.07212)研究脉动泡附近墙面声流和切应力，支持涡流形态的研究动机。其轴对称环状涡并非本轮两个平面Gaussian涡，不能移植数值或声称逐式复现。
- [Doinikov等双泡研究摘要](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/abs/acoustic-microstreaming-produced-by-two-interacting-gas-bubbles-undergoing-axisymmetric-shape-oscillations/9DDC94B3F63A9BF4C1162C474A402456)描述给定形变模态振幅下的一阶/二阶相互作用；本轮只核对摘要，未复现其流场。
- [LAMMPS fix external官方文档](https://docs.lammps.org/fix_external.html)：当前节点ID/坐标外力接口；本轮同时核对本地定制native力密度到N的转换顺序。

完整每步拓扑在 `CONNECTED_COMPONENT_HISTORY.jsonl.gz`，成员注册、自然断裂原始ID、FP64接受态、能量/表面覆盖CSV均保留。查看 [全部科研图表](../../../results/{TASK}/{RUN}/OPEN_RESULTS.html) 与 `DELIVERY_VALIDATION.json`。新旧源码diff和历史依赖SHA保留；没有覆盖旧正式结果。
'''
    text=intro+'\n'+'\n'.join(answers)+tail;path=O/'REVIEW_CLOT_SHEAR_LOAD_INVERSE_DESIGN_V1.md';path.write_text(text)
    figures=read(O/'FIGURE_MANIFEST.json')['figures'];items=''.join(f'<section><h2>{html.escape(f["name"])}</h2><p>{html.escape(f["caption"])}</p><a href="{f["name"]}.svg"><img src="{f["name"]}.png"></a></section>' for f in figures)
    page='<!doctype html><html lang="zh"><meta charset="utf-8"><title>Current-surface shear comparison</title><style>body{font:18px system-ui;max-width:1600px;margin:25px auto;padding:20px;background:#f3f5f8;color:#263544}section{padding:20px;background:white;margin:25px 0}img{width:100%}p{line-height:1.6}</style><h1>当前表面随动剪切：单向与双涡旋</h1><p>SHEAR_MECHANICAL_REQUIREMENT_ONLY。原生真实坐标与存活拓扑；新裂面润湿未解决。</p>'+f'<p>{html.escape(status["MULTIFRAGMENT_FORMATION_STATUS"])} · <a href="../../../verification/{TASK}/{RUN}/{path.name}">20问中文报告</a></p>'+items+'</html>'
    (F/'OPEN_RESULTS.html').write_text(page);(O/'OPEN_RESULTS.html').write_text(f'<!doctype html><meta charset="utf-8"><a href="../../../results/{TASK}/{RUN}/OPEN_RESULTS.html">打开全部科研图表</a>');print('REPORT_DONE',len(answers))
if __name__=='__main__':main()
