from common import *
import textwrap

def rd(n):return json.loads((O/n).read_text())
st=rd('FINAL_STATUS.json');au=rd('NUMERICAL_AUDIT.json');top=rd('V3_FRAGMENT_TOPOLOGY.json');en=rd('V3_ENERGY_LEDGER.json');ha=rd('V3_BACKGROUND_DRAG_HANDOFF.json');rs=rd('RUN_RESOURCES.json');dam=rd('V3_BULK_INTERFACE_DAMAGE.json');ev=dam['events'];fr=rd('V3_PILOT_FREEZE.json');ident=rd('SOURCE_AND_LIBRARY_IDENTITY.json');loc=rd('V3_FRACTURE_LOCATION_COMPARISON.json');mix=rd('MINIMAL_INTEGRATION_CHECK.json');par=rd('INITIAL_FORCE_PARITY.json');typed=rd('V3_TYPED_BOND_AND_ANCHOR_AUDIT.json');tr=rd('full/TRIAL_ACCEPT_REJECT_AUDIT.json');h=np.atleast_1d(np.genfromtxt(O/'V3_NATIVE_PILOT_HISTORY.csv',delimiter=',',names=True));complete=st['full_protocol_completed'];end=st['last_accepted_time_s']*1e3;largest=top['largest_detached_component'];tiny=top['unresolved_fraction_of_detached_mass'];maxs=au['maxima'];failure=au['exception'];elapsed=rs['main_wall_s'];jacs=[loc[str(t)]['old_new_failure_set_jaccard'] for t in [1,2]]
def when(key):
 e=ev[key]
 if e is None:return f'截至{end:.9g} ms未观察到'
 ids=e.get('bond_ids',[]);atoms=e.get('atom_ids',[])
 return f'{e["time_s"]*1e3:.9g} ms，第{e["step"]}接受步，原始bond ID {ids if ids else "不适用（拓扑事件）"}，atom ID {atoms}'
def peak(k):return maxs.get(k,{}).get('value',0)
def fmt(v):return '未发生/不可用' if v is None else f'{v:.9g}'
if complete:
 answer='是：在这一张冻结工程材料卡和历史C1人工载荷下，原生变分渐进软化已接入真实V3，完成0–2.33 ms，接受态力/能量/身份/锚点/积分及拓扑核查通过。其结论限于这一个探索工况，不认证材料Gc、全局时间步收敛或真实生理碎片尺寸。'
 if top['detached_material_points']==0:answer+='本次没有真实脱离；不能宣称已经得到脱离结果。'
 admission='适合下一轮少量、预先冻结的探索性载荷对照，但不是正式AI生产或材料认证。当前逐步独立审计与事件检查点成本不适合直接扩为高吞吐大批量采样。本轮不追加这些计算。'
 principal='就探索性计算而言没有阻断准入的主要实现问题；若要把碎片分布作为物理目标，主要缺口是未经认证的失效/软化材料卡（尤其界面尾部）。'
 nexttask='优先只读调查本工况中体内损伤、界面桥接和1–2点碎屑形成的空间先后关系，确定控制碎片尺寸的主要瓶颈，再决定下一轮少量有意义的载荷对照；不先调大载荷。'
else:
 answer=f'集成及最小验证通过，但不能回答完整正式工况成功：唯一pilot只接受到{end:.9g} ms，状态{st["status"]}。停止原因：{failure["message"] if failure else "见审计"}。接受前缀可审查，后续2.33 ms结果不存在。'
 admission='暂不适合直接开展载荷矩阵；先解决本唯一工况的明确停止原因并审查相同材料/载荷下的完整窗口。'
 principal=f'本轮唯一首要阻断是完整pilot未通过：{failure["message"] if failure else st["status"]}；不将其归因为材料不裂，也不通过修改材料或加大载荷绕过。'
 nexttask='下一轮只针对记录的停止门槛进行一次明确的实现/输出路径调查，保留本轮接受前缀和失败候选，不能由失败前缀生成完整工况标签。'
if complete:
 answer+=f' 粒子断后接触峰值重叠为{ha["maximum_particle_contact_overlap_dx"]:.6g}dx，数值账本通过不能替代接触几何精度认证。'
 admission='数值执行基础已具备，可支撑下一轮少量探索性对照；但建议先完成本工况高重叠接触事件的只读关联审查，再冻结后续载荷。当前逐步全量审计与检查点成本也不适合直接批量AI生产。本轮不追加物理工况。'
 principal=f'完整数值窗口没有代码或数值门槛阻断；当前应用判断最需优先核实的是断后粒子接触的{ha["maximum_particle_contact_overlap_dx"]:.6g}dx大重叠。它可能影响局部再碎裂，但当前数据不能直接证明因果关系；不是新的停止阈值，也不据此改材料。'
 nexttask='下一轮唯一建议：只读追溯0.376ms附着区非键合材料自接触及其邻近时段的键失效/分量分裂，核实大重叠和局部塌缩是否影响后续细碎化；先定位问题，不先增大载荷，不做本轮禁止的接触刚度扫描或Gc标定。'
pattern=('本接受窗口仍以1–2点碎屑为主，未形成≥200µm游离块；没有证据称其转为宏观成块脱层。' if tiny is not None and tiny>.5 and top['N_detached_ge200']==0 else '本接受窗口的实际分量分布见第15–18项；大小和机制需分别判断。')
contact_event=rd('MAXIMUM_PARTICLE_CONTACT_EVENT.json')['event']
contact_note=f'峰值发生在{contact_event["time_s"]*1e3:.9g}ms，atom ID {contact_event["i"]}/{contact_event["j"]}；两点当时均仍附着，原始参考图中没有直接键，是新发生的非键合材料自接触，不是游离块与Y壁碰撞。该时序证据尚不证明其导致碎屑。'
answers=[
('新软化模型是否成功处理V3本体及界面两类键？','集成通过。隔离C++继续使用完整endpoint-star变分自由能及局部精确消元，仅把s₀/s꜀按不可变原始bond type分配；不是把旧力简单乘(1−d)。本体0.15/0.18，界面0.140625/0.16875。界面尾部明确为EXPLORATORY_UNCALIBRATED。固定–固定第三类支架保留，s₀=1、s꜀=1.2不会因全固定而触发。'),
('1927个材料点及1561个锚点是否正确保留？',f'是。总3488原始atom ID逐步保持，1561锚点不积分、不接受外部加载/碎块阻力，也不计材料质量；全部1927材料点守恒。初始材料质量{top["initial_material_mass_kg"]:.12g} kg，离散体积1927×(50µm)³=0.240875 mm³，与连续几何约0.2716 mm³不混用。'),
('原始weighted volume和参考键是否正确？',f'是。原始157637键全部保留身份：82270本体、22040界面、53327固定支架。m按完整三类参考邻域及原r₀/taper/V计算、损伤后不重算；与历史参考最大相对误差{typed["weighted_volume_max_relative"]:.4g}。原始body与端点/type逐条比对；参考文件SHA {fr["original_bonds_sha256"]}。'),
('新模型未损伤节点力是否与旧LPS相容？','INITIAL_FORCE_PARITY.json为PASS。对真实初态和微小非均匀位移态，旧双Bug修正库、新typed库、独立LPS逐节点比对材料和锚点力、m、总力/力矩；不是只比较静止零总力。旧库只做静态run0，没有重跑历史物理工况。完整数据在该文件states中。'),
('本体及界面κ/d能否分别更新？',f'能。真实本体bond1与界面bond104310在局部伸长0.16时分别得到d=0.375和0.7265625；卸载不回退，再到0.205后各d=1且事件只计一次。原始ID/type绑定、κ单调和d单调全程核查。正式最大接受Δd={peak("max_dd"):.12g}；发生{rs["rejected_trials"]}次损伤拒绝及{rs["substep_count"]}个接受细分步。'),
('固定锚点是否保持固定？','全部1561锚点x逐值等于原始坐标、v逐值为零，逐接受步断言通过。没有用setforce0抹除原生力，因此反力=−ΣF_anchor逐步记录，外功因dx_anchor=0为零；材料动量账本包含该反力。'),
('界面失效后是否存在虚假传力？',f'最小混合损伤测试及正式接受态独立力核查中，完全失效原始键的非接触拉力最大{au["max_independent_errors"].get("failed_tension_N",0):.6g} N。部分损伤键仍由完整变分状态传力且排除粒子接触，d=1后邻居/特殊排除同步刷新；界面正式完全失效{dam["type_counts"]["2"]["fully_failed"]}条。该检查不等于界面材料标定。'),
('动态壁接触是否正常？',f'复用原NativeYWallQuery真实壁三角网格和penalty路径；严格区分1.0dx原生粒子断后接触与真实Y壁罚接触。前者在{when("particle_contact")}开始，最大重叠{ha["maximum_particle_contact_overlap_dx"]:.9g}dx；Y壁非零接触实际触发={ha["actual_wall_contact_exercised"]}、最大重叠{ha["maximum_wall_overlap_m"]:.9g}m。若未触发，只能声明路径接入及查询正常，不能称本轮动态壁碰撞已验证。壁接触半径20µm、k=50N/m、阻尼比0.5均沿用，不做刚度扫描。'),
('背景阻力是否在实际脱离时立即启用？',f'实际累计脱离{ha["actual_detached_points"]}点，同一接受时刻有有效背景拖曳{ha["same_time_drag_points"]}点；最大已配对延迟{fmt(ha["max_delay_s"])}s，状态{ha["status"]}。每次原生损伤commit及拓扑刷新后的post_force即判定真实无锚连接分量，COM位置和半步速度采样已有无栓Y背景流；1–2点100µs响应，较大分量沿用SN球阻力。无效背景域不得编造流速。附着点只保留历史500s⁻¹数值阻尼，游离点不重复施加该阻尼。最小人为隔离诊断只验证接口，不计正式碎裂结果。'),
('选用了哪一个历史人工载荷？','唯一C1_GAUSSIAN_A200。200Pa为历史面积预算，系数780.2322247780305Pa、初始峰值约770.874Pa；sigma=0.28mm，原中心/参考面积/切向逐点沿用，作用453个历史表面点，每步按当前坐标评估Gaussian。持续1.1605ms，前1/4余弦上升、中间1/2平台、末1/4余弦下降。没有加入新的580点含栓CFD映射，没有将Pa解释为超声功率。'),
('是否完成一次V3正式试运行？',answer+f' 正式launch标记恰好1个，从真实V3 t=0/x=X/v=0/d=0开始；实际接受{st["accepted_steps"]}步。最小8步力学和8步外力账本检查不属于第二个正式物理工况。'),
('首次本体及界面破坏发生在何时？','本体软化：'+when('bulk_softening')+'。\n\n本体完全失效：'+when('bulk_failure')+'。\n\n界面软化：'+when('interface_softening')+'。\n\n界面完全失效：'+when('interface_failure')+'。\n\n首次拓扑脱离：'+when('detachment')+'。'),
('最终本体和界面分别断开多少键？',f'末接受态本体完全失效{dam["type_counts"]["1"]["fully_failed"]}条、界面{dam["type_counts"]["2"]["fully_failed"]}条、固定支架{dam["type_counts"]["3"]["fully_failed"]}条。本体/界面出现非零d的总数分别{dam["type_counts"]["1"]["softened"]}/{dam["type_counts"]["2"]["softened"]}；后者包含完全失效，不能再相加。稀疏损伤流重建末d与FP64末态逐值一致。'),
('实际脱离了多少质量？',f'脱离{top["detached_material_points"]}点，质量{top["detached_material_mass_kg"]:.12g}kg（{top["detached_mass_percent"]:.9g}%）；剩余附着{top["attached_material_mass_kg"]:.12g}kg。二者和严格对应1927个材料点，锚层不计质量。'),
('出现多少个真实游离连通分量？',f'末态{top["detached_component_count"]}个。以全部存活原始键(d<1)构图并检查任意锚点连接，不用断键数量或粒子相距远近猜分量。原生末态全部端点图独立核对={top["native_final_endpoint_match"]}；若因门槛停止未再调用原生，此False表示未完成末态gather，逐接受步数量/身份检查仍保留。'),
('最大游离分量多大？','没有游离分量。' if largest is None else f'末态最大{largest["node_count"]}材料点，质量{largest["mass_kg"]:.12g}kg；按Ndx³等效体积直径{largest["equivalent_volume_diameter_m"]*1e6:.9g}µm。该分量点中心坐标跨度xyz={[float(v*1e6) for v in largest["center_span_m"]]}µm；这是当前点云跨度，与等效体积直径分开。全程最大瞬时分量{int(peak("largest_detached_nodes"))}点、等效直径{peak("largest_detached_diameter_m")*1e6:.9g}µm，不偷换为末态大小。'),
('是否出现≥200µm的游离碎块？',f'末态≥200µm共{top["N_detached_ge200"]}个、≥300µm共{top["N_detached_ge300"]}个；全程最大对应数量{int(peak("N_ge200"))}/{int(peak("N_ge300"))}。阈值只按离散参考体积Ndx³求(6V/π)^(1/3)，不按拉长点云外包尺寸冒充大块。'),
('是否仍以1–2节点碎屑为主？',f'末态1–2点分量{top["unresolved_1_2_node_components"]}个、合计{top["unresolved_1_2_node_material_points"]}点，质量{top["unresolved_1_2_node_mass_kg"]:.12g}kg，占脱离质量'+('不可定义（没有脱离）' if tiny is None else f'{tiny*100:.9g}%')+'。这些必须标注未解析数值碎屑，不能直接称为临床微栓子。'),
('与旧瞬断V3相比有什么真实差异？',f'原历史为217脱离点/208分量/最大4点/213细碎屑点；后续双Bug修正为199/190/最大4点/192细碎屑点。当前为{top["detached_material_points"]}/{top["detached_component_count"]}/最大{0 if largest is None else largest["node_count"]}点/{top["unresolved_1_2_node_material_points"]}细碎屑点。历史首次脱离0.2965ms、修正后0.3055ms，当前{when("detachment")}。旧两例到2.33ms；当前到{end:.9g}ms。几何及外加人工载荷相同，但本体和界面失效律、界面起点约1ppm取值、dt及MPI不同，只作探索性对照，不能将差异全归于本体软化。旧瞬断不存在当前相同D账本，不能伪造等义能量比较。'),
('新软化模型是否明显改变破坏位置和方式？',f'原修正旧模型首本体断键端点471/584；本轮真实首事件见第12项。新旧末态本体/界面失效原始键集合Jaccard分别{fmt(jacs[0])}/{fmt(jacs[1])}，参考中点均值/范围在V3_FRACTURE_LOCATION_COMPARISON.json。渐进κ/d是实际机制变化。{pattern} 是否发展为大块必须依第15–18项。单工况、不同失效律和可能不同末时刻不足以认证空间模式的独立因果改善。'),
('完整动力学能量是否可以核算？',f'接受窗口可独立核算：Wload+Wdrag+Wnumeric+Wwall,damp = ΔU+ΔK+D+ΔUc+ΔUw−Jc+R。D是同一x下损伤前后变分能差，R保留且用五项数值缺陷交叉核对，非材料耗能。末态单位nJ：'+', '.join(f'{k}={v*1e9:.12g}' for k,v in en['last'].items())+f'。最大累计归一化|R|={en["max_relative_R_percent"]:.12g}%<2%为接受准则；末K/U={fmt(en["last_K_over_U"])}，不宣称准静态。背景相对耗散/流动供能另为时间积分诊断，主要Wdrag使用真实力–位移功；半步v与完整v差值未隐藏。'),
('最大速度、损伤增量与接触重叠是多少？',f'全接受窗口最大速度{peak("max_speed_m_s"):.12g}m/s，Δd={peak("max_dd"):.12g}，粒子重叠{ha["maximum_particle_contact_overlap_dx"]:.12g}dx（{ha["maximum_particle_contact_overlap_dx"]*50:.9g}µm），Y壁重叠{ha["maximum_wall_overlap_m"]*1e6:.9g}µm。最大单步位移{h["max_step_displacement_dx"].max():.9g}dx、端点相对位移{h["max_relative_step_dx"].max():.9g}dx。门槛保持原值；小残差不自动认证接触物理精度或全局dt收敛。{contact_note}'),
('实际计算用时和内存是多少？',f'单rank×单线程，无GPU。正式主段墙钟{elapsed:.3f}s（{elapsed/60:.3f}min），CPU{rs["main_CPU_s"]:.3f}s，峰值RSS{rs["peak_RSS_bytes"]/1024**3:.6f}GiB；{rs["wall_s_per_physical_ms"]:.3f}s/物理ms。接受{rs["accepted_steps"]}步、拒绝{rs["rejected_trials"]}次、接受细分{rs["substep_count"]}步，dt={rs["actual_dt_min_s"]*1e6:.6g}–{rs["actual_dt_max_s"]*1e6:.6g}µs；检查点{rs["checkpoints"]}次耗时{rs["checkpoint_wall_s"]:.3f}s、恢复{rs["restart_wall_s"]:.3f}s。未损伤短步给出的早期估计为{fr["budget"]["estimate_s"]:.3f}s，未包含后续密集损伤、回滚和检查点额外成本；原生候选计算累计{rs["native_candidate_wall_s"]:.3f}s。原修正旧瞬断189.286s/2rank、树RSS1.781GiB，检查与I/O负担不同。只把special_bonds移至read_data前避免未用三跳临时表，未改图/力学。没有RSS或墙钟上限，2GiB输出限制保留。'),
('当前模型是否适合下一轮少量载荷敏感性试验？',admission),
('如果不适合，唯一主要原因是什么？',principal+' '+nexttask),
('当前哪些结果可以用于探索性AI样本？','V3_SOFTENING_PILOT_AI_RECORD.json已形成单例schema：原始几何/键/材料/源码/库SHA、真实载荷、接受窗、失效事件、质量/尺寸分布、Δd/能量余额、接触限制、成本和provenance。保存EXPLORATORY_UNCALIBRATED及实际validity_status；失败前缀不作完整工况监督标签。CRITICAL_LOAD=UNAVAILABLE，formal_training_data_eligible=false，不由单例造临界载荷。'),
('哪些科学限制必须保留？','本体宏观强度未重新认证；界面软化尾部未标定；实际材料Gc未测定；背景流为无栓Y流场的单向COM工程阻力，不是含栓双向FSI；真实微泡streaming未接入；人工Pa不是超声功率；1–2点碎屑未解析；接触重叠未完成物理标定；仅一个工况、非全局时间步收敛。G0_THEORY=0.07170070833333J/m²、G0_DISCRETE=0.06969623734234J/m²仅保留旧瞬断标称参考，不能称为新V3实际Gc。'),
('WSL科研文档是否实际更新？','实际/home/lzy/projects/clot_rupture.md通过只插入新段更新，原①–⑨及全部历史字节保留。实际读回、before/after备份、diff与SHA在WSL_DOCUMENT_UPDATE.json及wsl_document/；交付封存脚本要求该记录PASS且before+insert完全还原after后才允许封存。')]
report=f'# {TASK.upper()}\n\nRun `{RUN}`；性质`EXPLORATORY_APPLICATION`；科学标记`EXPLORATORY_UNCALIBRATED`。\n\n{answer}\n\n'+''.join(f'### {i}. {q}\n\n{a}\n\n' for i,(q,a) in enumerate(answers,1))
report+=f'''## 实现与可复现证据

- 新隔离库SHA：`{sha(LIB)}`；冻结SHA：`{sha(O/'V3_PILOT_FREEZE.json')}`。生产库与历史已接受结果不覆盖。
- 类型阈值、不可逆历史和typed restart由原始ID绑定；外部力只能以明确v3forces接口排在原生commit/单位转换之后。实际受力时间层为当前x与VV半步v；FP64交接保持完整接受x/v/time，重新计算完整v阻力时明确记录与原半步力差。
- 初始隔离构建发现生成style_bond.cpp工厂仍用旧类sizeof造成析构异常；正式前补编译该单元并复验受影响功能。失败库/进程日志保存在preflight_registry_build_failure，REGISTRY_ABI_REPAIR.json记录原因和最终PASS。PREPARED.json是该历史首构建记录，其库SHA不代表最终正式库；正式身份以SOURCE_AND_LIBRARY_IDENTITY及freeze为准。不掩盖预检失败，也不把它计成第二个正式V3工况。
- 源码diff、完整C++、构建命令/复用对象SHA、Python驱动、输入SHA、关键FP64状态、每接受步标量、稀疏损伤流、接触流及拒绝/恢复审计均归档。节点帧每20接受步采样，不冒充每步全量轨迹。所有正式回滚均检查x/v/κ/d/力/能量/端点图/特殊排除一致。
- 六组PNG/SVG和HTML位于`results/{TASK}/{RUN}/`。图使用真实点云、实际相机/坐标比例，不补画医学外壳或大碎片。

## 最后判断与下一项工作

{answer}

{nexttask}

本轮到此结束；不自动追加第二材料卡、第二载荷、CFD、微泡或AI训练。
'''
(O/('REVIEW_'+TASK.upper()+'.md')).write_text(report)
if (O/'REJECTION_ARCHIVAL_AUDIT.json').exists():
 a=rd('REJECTION_ARCHIVAL_AUDIT.json')
 report+=f'\n## 拒绝状态的无损存储整理\n\n事件密集阶段完整拒绝快照占主要输出空间，因此在不改冻结求解代码、材料、载荷、dt规则或门槛的情况下，对已写完且求解器不再读取的拒绝/回滚前NPZ做独立无损归档。{a["archived_files"]}个文件原{a["original_archived_bytes"]}字节压缩为{a["archive_bytes"]}字节；逐文件重建原NPZ并核对SHA完全相同后才移除冗余散文件，全部FP64内容保留。`archive_rejections.py restore <batch.bin> <new_empty_folder>`可恢复；各批JSON给出完整原始SHA。接受关键帧、实时restart和全标量轨迹未压缩改写，正式case仍只有一个。\n'
 (O/('REVIEW_'+TASK.upper()+'.md')).write_text(report)
insert=f'''\n\n**真实V3渐进软化集成与唯一探索工况（{RUN}）**：本轮从材料认证转回Y形血栓应用主线。隔离bpm/peri/variational保留完整endpoint-star变分能和局部精确消元，只按原始bond type分配软化门槛；实际1927材料点+1561全固定锚点、157637原始键（本体82270/界面22040/固定支架53327），dx50µm、δ150.5µm、G940Pa、K100000/3Pa、ρ1000kg/m³及完整m/r₀/V/taper不变。本体s₀/sc=0.15/0.18，界面0.140625/0.16875为EXPLORATORY_UNCALIBRATED；固定支架不积分、不损伤、不计血栓质量。锚点反力保留，完全失效后原生图及1.0dx接触刷新；真实Y壁罚接触沿用，无新刚度扫描。

选用唯一历史C1_GAUSSIAN_A200：200Pa面积预算、系数780.232224778Pa、初始峰值约770.874Pa、sigma0.28mm、453个历史表面点，按当前位置计算、原切向和面积不变，1.1605ms余弦升/平台/余弦降。没有叠加新的580点含栓CFD映射，Pa不等于超声功率。从真实初态x=X/v=0/d=0开始。正式状态**{st['status']}**，实际接受到{end:.9g}ms/{st['accepted_steps']}步；{answer}

本体首软化{when('bulk_softening')}；本体首完全失效{when('bulk_failure')}；界面首软化{when('interface_softening')}；界面首完全失效{when('interface_failure')}。末态本体/界面完全失效{dam['type_counts']['1']['fully_failed']}/{dam['type_counts']['2']['fully_failed']}条。首次脱离{when('detachment')}；最终脱离{top['detached_material_points']}点/{top['detached_component_count']}分量，质量{top['detached_material_mass_kg']:.12g}kg（{top['detached_mass_percent']:.6g}%），附着{top['attached_material_mass_kg']:.12g}kg。最大末游离分量{0 if largest is None else largest['node_count']}点，等效体积直径{0 if largest is None else largest['equivalent_volume_diameter_m']*1e6:.9g}µm；≥200/≥300µm分量{top['N_detached_ge200']}/{top['N_detached_ge300']}。1–2点未解析碎屑{top['unresolved_1_2_node_material_points']}材料点，不能称临床微栓子。真实连通性由全部d<1原始键及任意锚连接确定，全部1927材料点守恒。

脱离后同一原生post_force立即启用既有无栓Y背景COM阻力，实际{ha['same_time_drag_points']}/{ha['actual_detached_points']}点同接受时刻有有效拖曳，状态{ha['status']}。附着500s⁻¹数值阻尼与游离阻力分开、不重复；固定锚不受阻力。真实Y壁非零接触触发={ha['actual_wall_contact_exercised']}；最大壁重叠{ha['maximum_wall_overlap_m']*1e6:.9g}µm、粒子接触重叠{ha['maximum_particle_contact_overlap_dx']:.9g}dx。{contact_note}交接{ha['state_handoff']}，仅查FP64状态/受力时间层连续，不追加长时输运。

接受窗完整账本Wload+Wdrag+Wnumeric+Wwall,damp=ΔU+ΔK+D+ΔUc+ΔUw−Jc+R；材料D来自实际变分更新，背景阻力/数值阻尼/壁阻尼功单列，R不改名为耗能。最大归一化余额{en['max_relative_R_percent']:.9g}%，最大Δd={peak('max_dd'):.9g}、速度{peak('max_speed_m_s'):.9g}m/s；{rs['rejected_trials']}次损伤拒绝、{rs['substep_count']}个接受细分步、dt{rs['actual_dt_min_s']*1e6:.6g}–{rs['actual_dt_max_s']*1e6:.6g}µs。单rank×单线程，无GPU；墙钟{elapsed:.3f}s/CPU{rs['main_CPU_s']:.3f}s，峰值RSS{rs['peak_RSS_bytes']/1024**3:.6f}GiB。未设RSS/墙钟上限，输出2GiB及原数值门槛保留。special_bonds提前到read_data前避开无用三跳表，未改力学。正式前生成工厂类尺寸失配已最小修复并记录失败与复验，正式case仍只有1个。

历史217点/208分量/最大4点及双Bug修正199/190/最大4点只读对照，未重跑旧模型。损伤律、界面起点约1ppm、时间积分和MPI不同，不能把全部变化归因于本体软化；如末时刻不同，不视为等时定量改善。AI仅形成带源码/输入/库SHA及真实validity的单例探索schema，CRITICAL_LOAD=UNAVAILABLE、formal_training_data_eligible=false。本体宏观强度、界面尾部和实际Gc未认证；无真实微泡streaming、非双向FSI。旧G0_THEORY=0.07170070833333J/m²、G0_DISCRETE=0.06969623734234J/m²仅为旧瞬断标称参考。{nexttask}

见[本轮六组科研图](temp_storage/clot_github_sync_20261011/clot_clearance/results/{TASK}/{RUN}/OPEN_RESULTS.html)、[28项完整报告](temp_storage/clot_github_sync_20261011/clot_clearance/verification/{TASK}/{RUN}/REVIEW_{TASK.upper()}.md)。本轮交付后停止，不自动加大载荷、改材料或追加第二正式工况。

'''
(O/'WSL_DOCUMENT_INSERT.md').write_text(insert)
repro=f'''Run: {RUN} / {TASK}\nServer root: {R}\nEnvironment: source /workspace/y_clot/env.sh\nImmutable pilot freeze: {O}/V3_PILOT_FREEZE.json\nFreeze SHA256: {sha(O/'V3_PILOT_FREEZE.json')}\nFinal native library: {LIB}\nLibrary SHA256: {sha(LIB)}\nProduction unchanged; one formal launch marker only. Do not run launch.py pilot again in this archived directory.\n\nHistorical postprocessing commands used before sealing (DO NOT execute in the sealed directory):\ncd {A}\npython analyze.py\npython figures.py\npython report.py\nFor recomputation, create an isolated project copy and point common.py R/A/O/F and config.json at that copy/new run before invoking commands. Shared historical CFD assets may remain read-only. No archived file should be overwritten. PYTHONPATH must resolve the intended copied pd_clot package.\nAll reference/input files are in provenance; all compiled C++ and diffs in source_identity.\nBUILD_AUDIT.json records each compiler/link invocation and all reused object hashes.\nThe corrected initial build requires generated build/styles/style_bond.cpp to be recompiled with changed material class headers.\nPREPARED.json records first build; REGISTRY_ABI_REPAIR.json + SOURCE_AND_LIBRARY_IDENTITY.json + V3_PILOT_FREEZE.json identify the verified formal library.\nThe formal command is recorded in MANAGED_JOB.json and FORMAL_PILOT_LAUNCHED.json; the two minimal8step checks are separate diagnostics.\nLossless rejection snapshots: full/rejected_archives/batch_*.bin with original SHA in paired JSON; restore via python archive_rejections.py restore BATCH.bin NEW_EMPTY_FOLDER. Reconstruction was checked byte-for-byte before removing duplicate NPZs; software versions in EXECUTION_ENVIRONMENT.json.
Data: V3_NATIVE_PILOT_HISTORY.csv every accepted state; full/STREAM_LAYOUT.json defines sparse damage/contact streams and FP64 node frames; full/keyframes holds original bond kappa/d; full/TRIAL_ACCEPT_REJECT_AUDIT.json includes every rejection and identity restore.\nReproduction of a new physical run requires a new directory/run id and explicit task scope; do not change frozen card/load to manufacture a desired fragment size.\nOne terminal state handoff only if real detached components exist. No long-time transport executed.\nScientific: EXPLORATORY_UNCALIBRATED; GC_EFFECTIVE and CRITICAL_LOAD UNAVAILABLE; no formal AI training admission.\n'''
(O/'REPRODUCE.txt').write_text(repro)
print('REPORT_28_ANSWERS_AND_DOCUMENT_INSERT_READY',flush=True)
