"""Generate the Chinese review from sealed protocols and computed evidence."""
from fragment_topology_analysis import *
def pct(value):return '未定义（无脱离质量）' if value is None else f'{100*value:.6g}%'
def main():
    cases=read(O/'ANALYSIS.json');resources=read(O/'RUN_RESOURCES.json');validation=read(O/'DELIVERY_VALIDATION.json');assert validation['status']=='PASS'
    resources['validation_wall_s']=sum(q['validation_wall_s'] for q in validation['cases'].values())
    early=O/'SIMULTANEOUS_EARLY_REDUCTION.json'
    if early.exists():resources['early_read_only_reduction_wall_s']=read(early)['wall_s']
    dump(O/'RUN_RESOURCES.json',resources)
    old=read(O/'CRACK_BRIDGE_AUDIT.json')['cases'];after=read(O/'CRACK_BRIDGE_AFTER_LOADING.json')['cases'];baseline=read(O/'BASELINE_MANIFEST.json');checks=read(O/'KEY_CHECKS.json');build=read(O/'BUILD_AUDIT.json')
    labels={'simultaneous':'同步','staggered':'错时'};summaries={k:read(O/'cases'/k/'SUMMARY.json') for k in cases};fragmenthistory=read(O/'FRAGMENT_HISTORY.json')
    success=[labels[k] for k,v in cases.items() if v['exploratory_target_met']]
    observed=[labels[k] for k,v in cases.items() if v['history_max_N20_200']>=2]
    conclusion=('、'.join(success)+'满足冻结的探索目标；仍仅适用于本离散试样与机械替代加载。') if success else ('、'.join(observed)+'观察到多个目标带真实分量，但没有工况满足全部冻结约束。' if observed else '两组均未观察到两个及以上20–200节点真实独立块体；多块碎裂科学目标未达到。')
    status=dict(
        CRACK_BRIDGE_AUDIT_STATUS='PASS_ORIGINAL_ID_READ_ONLY_DIAGNOSTIC',
        MULTISITE_LOADING_STATUS='PASS_NATIVE_TWO_SOURCE_MULTISITE_MECHANICAL_SURROGATE',
        MULTIPLE_CRACK_NUCLEATION_STATUS={k:dict(bulk_only=v['crack_networks'],bulk_plus_interface=v['inclusive_failed_material_networks']) for k,v in cases.items()},
        CRACK_PERCOLATION_STATUS={k:dict(final_target_band_components=v['final_metrics']['fragment_count_target_band'],remaining_upper_region_bridges=after[k]['final']['regions'][0]['surviving_crossing_bonds'],scope='Free original-material graph is decisive; failed-edge mergers alone do not prove separation') for k,v in cases.items()},
        MULTIFRAGMENT_FORMATION_STATUS={k:dict(target_met=v['exploratory_target_met'],classification=v['classification'],history_max_N20_200=v['history_max_N20_200']) for k,v in cases.items()},
        TARGET_SIZE_DISTRIBUTION_STATUS={k:v['target_criteria'] for k,v in cases.items()},
        FINE_DEBRIS_FRACTION={k:dict(of_detached=v['final_metrics']['fine_debris_mass_fraction_detached'],of_initial=v['final_metrics']['fine_debris_mass_fraction_initial']) for k,v in cases.items()},
        FRAGMENT_COHESION_STATUS={k:dict(final_ge14=v['final_metrics']['N_ge14'],all_observed_ge14=len(fragmenthistory[k]['candidates_ge14']),evidence=f'cases/{k}/FRAGMENT_HISTORY.json',unresolved_two_node_pairs_are_not_continuum_blocks=True) for k,v in cases.items()},
        POST_UNLOAD_STABILITY_STATUS={k:dict(status='NOT_APPLICABLE_NO_TARGET_BAND_COHORT' if not v['common_cohort_medium_ids'] else 'FINITE_CONNECTIVITY_OBSERVATION_ONLY',common_cohort_surviving=v['common_cohort_medium_surviving'],actual_contact_free_tail_s=v['actual_contact_free_tail_s'],secondary_after_common_cohort=v['secondary_after_common_cohort'],scope='Scheduled return and connectivity observation are not a guaranteed force-free window; stationary actuators can be recontacted') for k,v in cases.items()},
        NUMERICAL_ADMISSION_STATUS=dict(status='PASS_TWO_COMPLETED_ACCEPTED_TRAJECTORIES',validation='DELIVERY_VALIDATION.json',relative_energy_gate=.01,maxima={k:v['max_relative_energy_residual'] for k,v in cases.items()}),
        PHYSICAL_FEASIBILITY_STATUS='UNVALIDATED_BUBBLE_STREAMING_OR_MICROJET; frictionless kinematic mechanical surrogate only',
        SCIENTIFIC_CONCLUSION=conclusion)
    dump(O/'FINAL_STATUS.json',status)
    lines=['# CLOT_MULTIFRAGMENT_BREAKUP_V1 中文科研审查',f'正式运行 `{RUN}`；全部科学判断基于实际原生接受态。{conclusion}',
        '本轮完成软件与数值交付，不以多断键、单一大块或失效簇数量代替多块碎裂。只运行冻结的两组均匀材料小试样；未追加第三工况、完整V3、CFD或参数扫描。',
        '## 基线与冻结协议',
        f"沿用1792材料点、336固定锚点、95666原始键、2736原预切。dx=50 µm，horizon=150.5 µm，G=940 Pa，K=33333.333 Pa，ρ=1000 kg/m³；bulk s0/sc=0.15/0.18，interface=0.140625/0.16875。无新增缺口或随机缺陷，无人为删除桥接键。旧库SHA `{baseline['baseline_library_sha256']}`；新隔离库SHA `{build['library_sha256']}`。旧库未替换。",
        '两个外部表面作用点为(−325,−75,175) µm、(325,75,−175) µm，分别朝+X、−X运动；Y偏置和跨切口上下错位产生不同方向的法向接触响应。半径182.8183198 µm，罚系数57973.58577 N/m²；初始间隙25 µm、压入240 µm。每源0–0.2 ms接近、0.2–1.2 ms压入、保持至1.4 ms、2.4 ms撤回。同步同时启动；错时第二源延后0.4 ms。延迟按两表面点间剪切波传播时间的一半取整，在正式计算前冻结。两组均观察至3.6 ms，采用共同2.8 ms群组追踪至末态。',
        '接触依当前坐标决定，释放后力消失，重新接触可再受力。两个源分别保存中心、速度、接触节点、反力、势能、冲量、接触时长和外功；没有恢复原参考节点dead traction。无摩擦球体仅施加径向法向力；偏置运动不等于已验证的切向摩擦牵引。该定义与[LAMMPS fix indent](https://docs.lammps.org/fix_indent.html)的球面罚接触形式一致。',
        'A/B唯一改变冻结的启动时间；实际反力、接触面积代理和外功可以不同。历史单压头位置与方向也不同，因此它只是历史上下文，不能作为仅改变源数量的单因素对照。',
        '此前1000 Pa规定参考表面剪切曾产生882节点游离块、内部连接存活约96.43%，它证明单个凝聚性大块可以脱离；这不是本轮要求的多块目标。该历史参考牵引仍有大变形跟随及脱离后持续受力局限，本轮没有重跑或恢复那种dead traction。',
        '## 固定多块目标与形态定义',
        '真实碎块由所有存活原始bulk/interface键及锚点可达性定义；非键合接触不计材料连接。成功须共同2.8 ms状态及末态同时满足：≥5个20–200节点块、脱离质量≥初始10%、最大块≤脱离质量50%、1–2节点细屑≤脱离质量20%，且至少5个共同群组的原始成员保持连通至末态。10%与50%是事先声明的工程门槛。各项独立保存，不使用加权总分掩盖失败。',
        '分组1–2、3–13、14–33、34–113、114–200、≥201节点仅用于当前离散分辨率。体积=N dx³，质量=ρN dx³，等效直径=(6N dx³/π)^(1/3)。每个接受态保留实际跨度、回转半径、PCA主轴/RMS长宽厚比、点云秩与内部存活键密度。退化厚度或比值记未定义；不补造平滑壳层。两个节点即使键100%存活也仍是未解析碎屑。',
        '## 实际结果',
        '| 工况 | 首软化/首完整失效 ms | 末态脱离点数 | 历史最大块/末态最大块 | 历史最多/末态20–200块 | ≥34/≥114块（末态） | 细屑/脱离质量 | 总外功 nJ |',
        '|---|---:|---:|---:|---:|---:|---:|---:|']
    for k,v in cases.items():
        e=v['final'];m=v['final_metrics'];lines.append(f"| {labels[k]} | {v['first_softening_s']*1e3:.6g}/{v['first_failure_s']*1e3:.6g} | {int(e['detached_nodes'])} | {v['history_max_free_nodes']}/{int(e['largest_nodes'])} | {v['history_max_N20_200']}/{m['fragment_count_target_band']} | {m['N_ge34']}/{m['N_ge114']} | {pct(m['fine_debris_mass_fraction_detached'])} | {e['W_J']*1e9:.8g} |")
    lines+=['','完整CSV、分量原始成员、逐接受态历史、损伤增量和关键FP64状态均保存。图的位移倍率为1；主相机固定、裁剪点数明示，并另给完整末态投影。桥接和起裂位置比较图明确使用参考坐标承载真实末态图信息。','## 18项科学问题直接回答']
    def both(fn):return '；'.join(labels[k]+'：'+fn(k,v) for k,v in cases.items())+'。'
    oldflap=old['uniform']['final']['regions'][0]
    answers=[
      ('之前的多处失效网络为何没有形成独立材料块？',f"历史均匀工况上半部896节点诊断区仍有{oldflap['surviving_crossing_bonds']}条真实跨区存活键，其中{oldflap['right_of_notch_tip']}条位于切口尖端右侧、{oldflap['long_bonds_ge_2dx']}条长≥2dx；存在完全避开|Y|≤50 µm中心观察层的锚点存活路径。37次失效簇出生和31次合并并没有切断材料图。"),
      ('有多少真实存活连接阻止主要候选裂缝贯通？','历史均匀末态三个诊断区（上半部、加载表面帽、左切口唇）跨区存活键分别为3461、2961、3439；非均匀为3467、2959、3384。'+both(lambda k,v:'新末态为'+str([q['surviving_crossing_bonds'] for q in after[k]['final']['regions']]))+'这些区域不是预设碎块边界；最小割是有界诊断代理，不能解释成应删除键或实际Gc。'),
      ('多区域加载是否改变首次起裂位置？',both(lambda k,v:'首软化原始键ID '+str(summaries[k]['events']['first_softening']['original_bond_ids'])+'，参考中点µm '+str((np.asarray(summaries[k]['events']['first_softening']['reference_midpoints_m'])*1e6).tolist()))+'历史首软化在(−175,0,375) µm；LOCAL_DEFORMATION_AUDIT保留局部应变与非仿射诊断，不能将键virial认作已标定Cauchy应力。'),
      ('是否出现多个空间独立裂纹起始区域？',both(lambda k,v:f"离散失效网络出生{v['crack_networks']['n_births']}次，最大并存{v['crack_networks']['max_simultaneous']}个，出生质心最大间距{v['crack_networks']['max_birth_centroid_separation_m']*1e6:.6g} µm")+'这支持离散失效区的空间分散；不等价于已解析多条连续裂面。'),
      ('多裂纹是否形成完整分离路径？',both(lambda k,v:f"失效簇合并{v['crack_networks']['n_merges']}次，末态目标带真实分量{v['final_metrics']['fragment_count_target_band']}个，上半部仍有{after[k]['final']['regions'][0]['surviving_crossing_bonds']}条桥接")+'只有零外部存活材料边界和失去锚点路径的分量才算完成分离；不能凭红色断键区相交判断。'),
      ('是否形成两个以上真实多节点游离块？',both(lambda k,v:f"历史最大同时≥2节点分量{v['history_max_size_threshold_counts']['N_ge2']}个，≥14节点{v['history_max_size_threshold_counts']['N_ge14']}个；末态分别{v['final_metrics']['N_ge2']}和{v['final_metrics']['N_ge14']}个")+'必须将2节点对与具有厚度的中等块体分开。'),
      ('是否形成至少5个20–200节点碎块？',both(lambda k,v:f"历史最多{v['history_max_N20_200']}个，末态{v['final_metrics']['fragment_count_target_band']}个，完整冻结目标={v['exploratory_target_met']}")),
      ('有多少块达到≥34节点？',both(lambda k,v:f"历史同时最大{v['history_max_size_threshold_counts']['N_ge34']}，末态{v['final_metrics']['N_ge34']}")),
      ('有多少块达到≥114节点？',both(lambda k,v:f"历史同时最大{v['history_max_size_threshold_counts']['N_ge114']}，末态{v['final_metrics']['N_ge114']}")),
      ('最大游离块占总脱离质量多少？',both(lambda k,v:pct(v['final_metrics']['largest_free_fragment_fraction_detached'])+'；占初始材料'+pct(v['final_metrics']['largest_free_fragment_fraction_initial']))),
      ('1–2节点碎屑占总脱离质量多少？',both(lambda k,v:pct(v['final_metrics']['fine_debris_mass_fraction_detached'])+'；占初始材料'+pct(v['final_metrics']['fine_debris_mass_fraction_initial']))),
      ('中等块体卸载后是否保持内部连接？',both(lambda k,v:f"共同2.8 ms的20–200节点群组{len(v['common_cohort_medium_ids'])}个，原成员连通至末态{v['common_cohort_medium_surviving']}个；共同窗口后二次裂解{v['secondary_after_common_cohort']}次；最后连续无源接触窗口{v['actual_contact_free_tail_s']*1e6:.6g} µs")+'若没有中等块，稳定性应写“不适用”，不能因为没有分裂事件便称稳定。'),
      ('同步与错时哪种更容易形成多条贯通裂纹？',both(lambda k,v:f"历史最多目标带块{v['history_max_N20_200']}，最终脱离质量{pct(v['final_metrics']['total_detached_mass_fraction'])}，实际W={v['final']['W_J']*1e9:.7g} nJ")+'仅两条冻结轨迹，不能按失效簇出生数推断普遍优劣；应按真实中等分量、桥接剩余和独立质量门槛比较。'+conclusion),
      ('限制来自空间加载、非局部桥接还是损伤准则？',both(lambda k,v:f"上半部末态{after[k]['final']['regions'][0]['surviving_crossing_bonds']}条桥接中，当前压缩{after[k]['final']['regions'][0]['currently_compressive']}条、历史κ未越过s0为{after[k]['final']['regions'][0]['historically_below_onset']}条、长≥2dx为{after[k]['final']['regions'][0]['long_bonds_ge_2dx']}条")+'这些重叠统计直接支持残余空间覆盖和正伸长历史不足的解释，不能唯一识别本构错误。当前正键伸长损伤下，压缩连接不会仅因高接触压力自动失效；亦不能仅凭此否定整个材料模型。'),
      ('是否值得引入文献混合模式剪切/拉伸损伤？','值得作为后续独立、可标定且受能量约束的假设比较，前提是先把尚存桥接的局部剪切/张开响应与数据对应，再验证强度、混合模式耗能及不可逆变分一致性。本轮没有实现或验证新混合模式律，不能把更易碎当作模型正确。'),
      ('是否足以支持真实V3迁移？','不足以作为真实V3多块碎裂机制已经成立的证据。当前只有固定网格/接触罚参数的小试样、两种时间组织；仍需机制有辨识力的验证、分辨率/接触敏感性与真实几何的独立准入。没有运行完整V3。'),
      ('距离微泡Streaming/微射流物理可实现性还缺什么？','缺少微泡或微射流到两作用点的牵引方向、面积、峰值、时间谱、持续时长和总冲量/功的独立物理映射，以及流体耦合、载荷可达性和实验标定。这里的无限惯性规定轨迹球体不是微泡；本轮没有证明流动能实现这些反力、行程和输入功。'),
      ('实际耗时和内存？',f"两组累计墙钟{resources['native_wall_s']:.3f} s、过程CPU {resources['native_cpu_s']:.3f} s，单进程峰值RSS {resources['peak_RSS_bytes']/1024**3:.6f} GiB。"+both(lambda k,v:f"{resources['cases'][k]['accepted_steps']}接受步，wall={resources['cases'][k]['wall_s']:.3f} s，CPU={resources['cases'][k]['cpu_s']:.3f} s，RSS={resources['cases'][k]['peak_RSS_bytes']/1024**2:.3f} MiB")+'1 rank×1线程顺序运行；正式driver时间包括逐步核验和输出。离线桥接、编译与后处理单列，不重复计算历史仿真成本。')]
    for i,(question,answer) in enumerate(answers,1):lines.append(f'{i}. **{question}**\n\n   {answer}')
    lines+=['## 能量、接触与数值边界',
        '账本W0+W1=ΔU_LPS+ΔK+D+ΔU_pair−J_pair_activation+ΔU_tools+R。损伤耗能按同一当前构形损伤更新前后势能差核算；无数值阻尼。R保持独立，不伪装为材料耗散。相对残差的分母是累计能量尺度的历史最大值，不是末态净外功。',
        both(lambda k,v:f"两源末态功={[v['final'][f'site{j}_W_J']*1e9 for j in range(2)]} nJ；两源峰值反力={np.asarray(v['peak_site_reactions_N'])*1e6} µN；最大重叠/dx={v['max_site_overlap_dx']}；最大相对能量残差={v['max_relative_energy_residual']:.8g}"),
        '罚接触允许有限重叠；实际重叠已记录，不能称完全刚性不可穿透。未进行接触刚度或网格收敛扫描，因此不能把一次离散轨迹解释为物理极限或真实血栓强度认证。',
        '原生修改限于物理单位转换后接纳两个独立命名接触源；没有改变变分材料核、原始ID历史、MPI或restart机制，因此未重复完整旧回归矩阵。新接口检查包括无接触/释放、径向力、势能梯度、独立及叠加反力/势能/外功。两组每步检查原生VV预测、不可逆历史、损伤增量≤0.1、原始键计数、接触加和和1%能量门槛。交付时以稀疏原始ID损伤记录重放全部接受态分量成员及关键状态。',
        '一次预检fixture曾把固定材料点留在球心穿越路径上，跨越r=0势能尖点，长区间单段功积分未通过。失败原值及日志保留在PRECHECK_FIXTURE_RECOVERY和checks_attempt01目录；在正式运行前改为接触进入附近的短路径梯度/功检查，未修改物理协议或放松容差。正式积分使用≤0.25 µs的接受步和接触进入/退出分段，同时监测4/8点积分差和能量残差；不据此声称跨球心长步积分已验证。',
        '本轮仅把终点循环浮点余量停止门槛从1e−15 s改为1e−14 s，避免上一轮末尾约1e−15 s的多余接受步。记录实际原生时钟，不将其重写为理论时刻；不改变VV或损伤事件规则。预测损伤子步下限0.03125 µs，非结果驱动调整。',
        '历史非均匀残余桥接统计首次离线版本误用了均匀s0比较κ阈值，已在加载冻结前改用原始ID冻结场并保留初版及修正记录。存活图、桥接键数和均匀基线不受影响。',
        '## 图论代理与真正材料分离',
        '失效网络主对照沿用历史纯bulk网络定义，以便37/31等计数可直接对照；历史两组自然interface失效均为0，因此P0的bulk图已包含所有自然失效材料边。新工况另保存完整bulk+interface失效图，保留固定锚端点与每条原始键的类型、出生/合并和全部失效ID。两种视图不混用计数。完整图末态与独立稀疏连通计算、原始接受态失效ID逐项核对。'+both(lambda k,v:f"完整图新增bulk/interface失效键{v['inclusive_failed_material_networks']['new_failed_bulk_bonds']}/{v['inclusive_failed_material_networks']['new_failed_interface_bonds']}，出生/合并{v['inclusive_failed_material_networks']['n_births']}/{v['inclusive_failed_material_networks']['n_merge_events']}"),
        '初始历史6个状态×3个诊断区域，另对两组新工况各取最高冻结目标状态及末态×3区域；共30次有界最小割。源核≥34节点，固定锚点/指定远端材料为sink，容量max(1,round(1024(1−d)))仅作用于存活原始键。残余可达性与割容量独立核对；从不向物理模型写回删键。整数容量计算采用[SciPy maximum_flow](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csgraph.maximum_flow.html)。三个几何诊断区不是已形成块体。',
        '要从连接主体形成多个中等块，必须在多个有实际体积的区域周围形成闭合的材料分离边界，消除所有bulk/interface跨界路径（包括非局部长键和三维绕行），同时保留区域内部连通。仅增加表面断键、失效网络出生或接触压力不满足这一图条件。本轮加载是对此条件的一次有界检验，不预设分割面，也不保证目标形态。',
        '## 交付与复现',
        '源码：analysis/clot_multifragment_breakup_v1；证据：本目录；图像：results/clot_multifragment_breakup_v1/'+RUN+'。打开OPEN_RESULTS.html查看八组主题的13幅PNG/SVG。REPRODUCE.md列出冻结步骤、隔离库、正式运行与只读后处理。SOURCE_DIFF、SHA清单、逐步日志、实际关键状态和WSL文档前后diff随交付保存；可重建编译对象与动态库留服务器，仅同步其SHA及构建命令。',
        '本轮结束后停止；不自动继续V3、混合损伤开发或微泡验证。']
    blocks=[]
    for line in lines:
        if line.startswith('|') and blocks and blocks[-1].startswith('|'):blocks[-1]+='\n'+line
        else:blocks.append(line)
    (O/'REVIEW_CLOT_MULTIFRAGMENT_BREAKUP_V1.md').write_text('\n\n'.join(blocks)+'\n')
    print('REPORT_DONE',conclusion)
if __name__=='__main__':main()
