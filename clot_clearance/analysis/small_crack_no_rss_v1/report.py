from common import *
st=json.loads((O/'FINAL_STATUS.json').read_text());e=json.loads((O/'DYNAMIC_ENERGY_LEDGER.json').read_text());n=json.loads((O/'NUMERICAL_VALIDITY.json').read_text());r=json.loads((O/'RUN_RESOURCES.json').read_text());ev=json.loads((O/'FIRST_SOFTENING_AND_FAILURE.json').read_text());loc=json.loads((O/'CRACK_LOCALIZATION_AUDIT.json').read_text());top=json.loads((O/'NATIVE_COMPONENT_TOPOLOGY.json').read_text());bridge=json.loads((O/'CROSS_PLANE_BRIDGE_AUDIT.json').read_text());area=json.loads((O/'CRACK_AREA_IDENTIFIABILITY.json').read_text());contact=json.loads((O/'CONTACT_STREAM_AUDIT.json').read_text())
regions={'NOTCH_TIP_INTERIOR':'裂尖厚度中央','NOTCH_TIP_SIDE_SURFACE':'裂尖侧表面','CRACK_FLANK':'裂纹侧面','GRIP':'加载端','OTHER_FREE_SURFACE':'其他自由表面','REMOTE_INTERIOR':'远处内部'}
def eventtext(event):
 if event is None:return '在已接受窗口内未观察到；不能补造事件时间或位置。'
 return f'{event["time_s"]*1e3:.9g} ms，第{event["step"]}接受步；区域'+ '、'.join(regions[z] for z in sorted(set(event['regions'])))+'；原始bond ID '+','.join(map(str,event['bond_ids'][:20]))+('等' if len(event['bond_ids'])>20 else '')+f'，共{len(event["bond_ids"])}条同一步事件。完整参考/实际中点、κ/d及实际伸长见FIRST_SOFTENING_AND_FAILURE.json。'
def yn(z):return '是' if z else '否'
endms=st['actual_end_s']*1e3;last=e['last'];energytext='，'.join(f'{k[:-2]}={v*1e9:.12g} nJ' for k,v in last.items())
lead=loc['notch_softening_lead_over_grip_s'];onsettext=('是，裂尖先于加载端进入软化。' if loc['notch_softening_before_grip'] else '没有证实裂尖先于加载端进入软化。')+(f'领先{lead*1e6:.9g}µs。' if lead is not None else '加载端未出现事件或缺少可比较事件。')
tip_failed=sum(z['naturally_failed'] for z in loc['regions'][:2]);tip_D=sum(z['star_center_D_J'] for z in loc['regions'][:2]);geometric=json.loads((O/'GEOMETRIC_BRIDGING_REVIEW.json').read_text())
finding=f'裂尖首先软化并断键，但后期出现广泛的其他区域损伤：尖端两区完全失效{tip_failed}/{st["natural_failed"]}条，端点star中心归属的尖端D占全局{100*tip_D/last["D_J"]:.5g}%。末态候选开口带仅50µm宽，未通过100µm持续前缘推进准则；仍有{st["live_cross_bonds"]}条跨面桥接、一个3080点连通体，horizon无桥接保守子集为空。该50µm/单体结论限定本次实际结果。'
questions=[
('最终采用什么三维几何？','沿用14×10×22胞元中心三维矩形试样，700×500×1100 µm胞元体积尺寸，参考点中心跨度650×450×1050 µm；不是完整Y形血栓。'),
('实际节点和原始键数？','3080原始点、144448原始参考键，dx=50 µm、δ=150.5 µm；未重新生成或粗化。'),
('初始裂纹、前缘和预切？','原始z=0，x从−350至−50 µm，贯穿y±250 µm；原始线段精确交点决定3524预切。初始投影0.15 mm²，不算自然新增失效、D或新面积；m仍包含全部原始参考邻域。'),
('相对历史试样的设计依据？','剩余韧带2.66δ、厚度3.32δ、裂尖到受控面3.49δ，预注册区分裂尖中央/侧面与加载端。M5B/M5C曾显示端部竞争与大量桥接；本轮没有通过改几何消除这种竞争，必须按实际结果判断，不能由设计就宣称裂尖优先。'),
('加载和时程？',f'仍为上下各侧120 µm半余弦2 ms张开、0.5 ms保持；计划2.5 ms。实际完成到{endms:.9g} ms。最大每侧规定速度0.0942478 m/s；从原始预切t=0新运行，没有把旧0–50 µs接到新轨迹。'),
('自由和受控DOF？','上下各140点仅规定z，全部点x/y自由；2800内部点z也自由，8960自由DOF/280规定DOF，原生VV与fix move NULL分量积分不重叠。边界CPP与前轮已验证库完全相同。'),
('第一次软化发生在哪里？',eventtext(ev['first_softening'])),
('第一次完整断键发生在哪里？',eventtext(ev['first_complete_failure'])),
('裂尖是否先于加载端进入损伤？',onsettext+'当前分类为'+loc['classification']+'，仅描述首次起裂与已观测的集中损伤云，不表示后续全部损伤都位于裂尖；详细区域竞争见CRACK_LOCALIZATION_AUDIT.json。'),
('是否观察局部损伤集中？',yn(loc['localized_cloud_observed'])+'。采用冻结的≥20自然损伤键、≥60%在尖端径向δ、≥6相连50 µm体素并跨≥2dx及接近原前缘的判据。集中损伤云不是唯一连续裂面。'),
('是否发生连续裂纹前缘推进？',yn(st['crack_propagation_candidate'])+'（按冻结候选平面的无桥接+实际开口+连续缺口前缀，且中央厚度条带在两个已保存时刻推进≥2dx）。本次末态仅前进一格50µm，未达到100µm；首个开口无桥接单元出现在1.375ms，不能称为已建立连续裂纹传播。'),
('是否仍有存活跨面键？',f'原始8100跨面记录，预切3524后4576存活；末接受态存活{bridge["last_live_cross"]}，新增完全失效跨面{bridge["natural_failed_cross"]}，其中存活长键r0>dx为{bridge["last_long_live_cross"]}。'),
('实际连通分量？','节点数为'+str(top['component_sizes'])+'。图使用全部原始d<1键，部分软化仍保留连接；每个分量质量、实际坐标包络、质心、最大点间距和atom ID见NATIVE_COMPONENT_TOPOLOGY.json。原生末态端点集合核对：'+str(top['native_final_endpoint_match'])+'；null代表停止后未采集，不表示失配。'),
('加载端、自由表面或接触影响？',f'首次失效在裂尖侧表面，后期加载端及其他自由表面大量损伤，不能解释为纯裂尖扩展。实际最大接触重叠{contact["max_overlap_dx"]:.9g}dx，约{contact["max_overlap_dx"]*DX*1e6:.6g}µm，达到点距的{100*contact["max_overlap_dx"]:.6g}%，不能排除对局部动力学的影响。有限厚度/自由表面与horizon限制保留，不另调接触刚度。'),
('预裂纹以外其他损伤？',f'不跨候选原始面的新增完全失效{loc["off_plane_natural_failed"]}条；逐区域自然软化/失效数见区域表。这是相互作用统计，不能直接当作另一连续裂面面积。'),
('新增分离区域可追溯吗？','交点、存活长键、对应上下节点开口及从缺口延伸的每个50 µm网格单元均可追溯；实际网格代理见PLANE_HISTORY.npz与CRACK_PLANE_HISTORY.csv。物理新生分离区域尚未独立认证。'),
('是否能测量新投影分离面积？',f'当前{area["status"]}，PROJECTED_SEPARATION_AREA=null。末态无桥接/开口/连续前缀的网格代理为{area["notch_connected_grid_proxy_m2"]:.12g}m²，保守horizon无桥接子集为{area["conservative_grid_proxy_m2"]:.12g}m²；它们不是断键数×dx²，也不能自动称唯一真实裂面。'),
('有限horizon不确定性？','保留δ=150.5 µm、50 µm单元、全部长程跨面记录及前后真实开口。严格子集要求整个单元的δ邻域无桥接。冻结FP64 floor交点在精确网格边界的归属存在相邻单元不确定性；GEOMETRIC_BRIDGING_REVIEW.json另将边界交点同时计入所有相邻单元检查虚假空单元，未改变原规则或原生轨迹。'),
('完整能量账本？',f'末接受{endms:.9g}ms：{energytext}。初态U/K/Uc=0，W=ΔU+ΔK+D+ΔUc−Jc+R。W仅受控z分量反力功加受控质量动能变化；自由横向不算规定运动功。末态K/U={e["last_K_over_U"]}，保持段受控速度已为零而自由节点动能仍显著，因此归为动态损伤响应，不是准静态试验。'),
('裂纹附近与其他区域耗能？','按相同endpoint-star中心分区精确求和；见下表。分区只是声明的耦合能量归属约定，不是唯一逐键裂面能；全局D与分区和差为'+f'{e["regional_sum_error_J"]:.12g}J。'),
('协议相关分离耗能？','G_SEPARATION_PROTOCOL=UNAVAILABLE。尚未同时认证连续新生物理面积、裂纹推进和可用的裂纹相关耗能归属；不得将全局D除以任意截面。'),
('为什么仍不能称材料Gc？','GC_EFFECTIVE=UNAVAILABLE。数值积分/拓扑通过、出现失效或出现多个分量均不等于实验材料Gc；还受动态惯性、非局部厚度和边界效应及分离面积定义限制。'),
('宏观强度是否不变？','没有宏观峰值强度独立认证；G/K/s0/sc未改仅保证输入数学参数保持，不能替代宏观试验。'),
('材料是否保持？','G940 Pa、K100000/3 Pa、ρ1000 kg/m³、s0.15、sc.18、α0、dx50 µm、δ150.5 µm、完整m/r0/V/taper及1.0dx接触完全不变。'),
('时间步和事件积分是否通过？',f'已接受{n["accepted_steps"]}步，损伤拒绝{n["rejected_steps"]}次，实际dt范围{n["dt_min_s"]*1e6:g}–{n["dt_max_s"]*1e6:g}µs，最大接受Δd={n["max_dd"]:.12g}≤0.10。完整协议状态{n["status"]}，全局时间收敛仍未认证；并未取消任何数值门槛。'),
('原始身份/状态错位？','所有接受步核对3080 atom ID及144448原始bond ID，κ/d与独立解比较；完整d更新流重建末态及所有采样时刻分量。前50 µs的201行全部标量及11帧x/v/力与上一运行逐值一致；不声称比较了未保存的逐步完整κ。'),
('最大能量余额？',f'最大接受归一化余额{e["max_relative_R_percent"]:.12g}%，最大绝对R={e["max_absolute_R_J"]:.12g}J；原累计归一化和2%门槛保持。若存在未接受候选，其余额不混入合法历史。'),
('1.0dx接触是否干扰解释？',f'记录几何候选{contact["geometric_candidate_records"]}条，其中承载r<dx记录{contact["force_bearing_records"]}条、恰好r=dx零力记录{contact["zero_force_touching_records"]}条；最大重叠{contact["max_overlap_dx"]:.12g}dx。保留接触限制，不另调刚度。'),
('实际计算成本？',f'本次主段墙钟{r["main_wall_s"]:.6f}s、CPU{r["main_CPU_s"]:.6f}s、峰值RSS{r["peak_RSS_bytes"]/1024**3:.6f}GiB；1rank×1线程。原生候选累计{r["native_candidate_wall_s"]:.6f}s，回滚{r["rollback_wall_s"]:.6f}s。按用户新指令取消4GiB门槛，无替代RSS上限；原1800s墙钟及2GiB输出预算保持。未改高阶邻接初始化顺序、未新编译或重跑P1。'),
('科研图是否交付？','六组PNG/SVG及OPEN_RESULTS.html；同几何图1复用并标明，其余按本次真实结果生成。全部点阵、事件和曲线来自实际数据，无人工医学外壳，无位移放大。'),
('生产库与历史是否不变？','复用前轮隔离库SHA730e11651aee6533b4fb7d38c8341e06ad853bf3397c4908fc6ef905d4e3ac90，无C++重编译；生产库、env.sh、旧研究库及上一157文件清单保持。只在新目录改路径与RSS条件，活跃入口推进时备份。'),
('WSL记录是否实际更新？','本轮交付时实际更新/home/lzy/projects/clot_rupture.md，记录取消内存限制的用户授权和本次结果，保留上次资源失败事实及全部历史；before/after/diff/SHA、并发写入检查和实际读回见WSL_DOCUMENT_UPDATE.json。'),
('主要科学发现或失败原因？',finding+f'完整协议完成={st["full_protocol_completed"]}；停止分类={st["stop_category"]}，实际墙钟{r["main_wall_s"]:.6f}s超原1800s。预算每200接受步检查，故超过约26s才停止；输出仍低于2GiB，未触发内存或数值门槛。尚缺0.1ms保持段，不推断其结果。'),
('下一轮优先工作？','优先审查裂尖损伤为何在50µm局部开口后停滞，而端部/其他自由表面继续耗能，以及15.44%点距接触重叠的影响。本次没有数值准入失败，不应优先放宽数值门槛。若需要完整2.5ms协议，应另行明确从已保存2.4ms原生检查点补足0.1ms的恢复方案，不能将当前状态冒充完成；本轮不自动追加求解或改几何/材料。')]
assert len(questions)==34
text=f'''# 取消4 GiB上限后的同工况继续运行审查

任务：{TASK}；Run `{RUN}`。用户明确指令：“不设置 4 GiB 预算，继续进行运行”。

**实际结果：{st['status']}；{loc['classification']}。** 已接受到{endms:.9g} ms；自然软化{st['natural_softened']}条，新增完全失效{st['natural_failed']}条，末态存活跨面{st['live_cross_bonds']}条。当前未认证唯一裂面面积，分离能与材料Gc不可用。

{finding} 几何边界归属审查后仍为10个开口单元、0.025mm²网格代理，不能当作可信新生裂面。末态K/U={e['last_K_over_U']:.6g}，属于明显动态响应。数值检查通过仅限已接受的0–2.4ms，未完成最后0.1ms保持；本次停止原因是原墙钟预算，不是RSS或物理数值门槛。

## 本次授权与实施范围

取消应用层4 GiB RSS停止条件，不设替代内存门槛；原几何、材料、载荷、接触、库、事件细分和全部数值门槛不变。原1800s墙钟及2GiB输出门槛保持。只新增一次本工况原生运行；不扫描参数、几何、MPI或GPU。

上一轮50µs FP64不是完整原生restart，原生检查点只有t=0，因此本次从原始预切t=0重新推进，没有拼接前缀。新路径前201行标量及11个完整节点帧与旧0–50µs逐值一致。P1数值检查复用其相同配置和库的已验证结果，无新增三步或额外物理路径。

前一资源失败记录完整保留，其预检资源遗漏不被追溯改写；本次不再把同样峰值视为4GiB违例。当前资源成本、数值停止与科学结论分开报告。

## 实测汇总

| 项目 | 结果 |
|---|---|
| 实际物理时间 | {endms:.9g} / 2.5 ms |
| 接受 / 损伤拒绝 | {n['accepted_steps']} / {n['rejected_steps']} |
| 最大接受Δd | {n['max_dd']:.12g} |
| 最大归一化余额 | {e['max_relative_R_percent']:.12g}% |
| 自然软化 / 新完全失效 | {st['natural_softened']} / {st['natural_failed']} |
| 初始切断 | 3524（不计自然耗能） |
| 存活跨面键 | {st['live_cross_bonds']} |
| 分量节点数 | {top['component_sizes']} |
| 主段墙钟 / CPU | {r['main_wall_s']:.6f} / {r['main_CPU_s']:.6f} s |
| 峰值RSS | {r['peak_RSS_bytes']/1024**3:.6f} GiB，无RSS停止上限 |
| 实际停止原因 | {n['stop_reasons']} |
| 新裂面 / 分离能 / Gc | UNRESOLVED / UNAVAILABLE / UNAVAILABLE |

## 34项科研回答

'''
text+='\n\n'.join(f'### {i}. {q}\n\n{a}' for i,(q,a) in enumerate(questions,1))
text+='\n\n## 真实区域统计\n\n| 预注册区域 | 自然软化（含完全失效） | 自然完全失效 | star中心D / nJ |\n|---|---:|---:|---:|\n'
for z in loc['regions']:text+=f'| {regions[z["region"]]} | {z["naturally_softened"]} | {z["naturally_failed"]} | {z["star_center_D_J"]*1e9:.12g} |\n'
text+='''

## 最终回答与归档

本次可以依据实际事件回答损伤发生在哪里、是否扩展及是否形成图分量，但**当前仍不能宣布已建立空间可辨识的唯一自然裂面—耗能对应，也不能报告已标定材料Gc**。初始预切面积和断键数不能替代新生裂面。旧`G0_THEORY=0.07170070833333 J/m²`与`G0_DISCRETE=0.06969623734234 J/m²`仍是旧瞬断LPS标称值。

[六组真实科研图](OPEN_RESULTS.html)。可复核数据包括完整CSV、FP64节点帧/关键态、原始d变更与事件、原生restart、跨面/区域/分量/能量审计；缺失的候选或末态原生检查均明确标出。交付后停止，不自动开始另一工况。
'''
(O/('REVIEW_'+TASK.upper()+'.md')).write_text(text)
insert=f'''\n**取消4GiB内存停止条件后的同工况继续运行（{RUN}）**：用户明确指令“不设置4GiB预算，继续进行运行”。本次只取消应用层RSS上限，保留原1800s墙钟/2GiB输出及全部损伤、能量、运动学门槛；3080点/144448原始键/3524预切、700×500×1100µm几何、完整m/r0/V/taper、G940Pa、K100000/3Pa、ρ1000kg/m³、s₀=0.15、s_c=0.18、α=0、1.0dx接触和120µm每侧2ms半余弦+0.5ms保持均不变。复用前轮混合DOF隔离库，无新C++编译或P1物理重跑。旧50µs记录不是完整原生restart，因此从原始预切t=0新运行；前201行标量及11帧x/v/力与旧0–50µs逐值一致，不拼接轨迹。\n\n本次实际到{endms:.9g}ms，接受{n['accepted_steps']}步/损伤拒绝{n['rejected_steps']}次，dt范围{n['dt_min_s']*1e6:g}–{n['dt_max_s']*1e6:g}µs，最大接受Δd={n['max_dd']:.9g}；完整协议状态 **{st['status']}**，停止原因{n['stop_reasons']}。首次软化：{eventtext(ev['first_softening'])}首次完整失效：{eventtext(ev['first_complete_failure'])}区域分类 **{loc['classification']}**；自然软化{st['natural_softened']}、新增完全失效{st['natural_failed']}，末态跨面存活{st['live_cross_bonds']}，分量点数{top['component_sizes']}，全部3080原始ID保留。\n\n末接受能量：{energytext}，最大归一化余额{e['max_relative_R_percent']:.9g}%（原门槛2%不变）。K/U={e['last_K_over_U']}，未认证准静态；区域D按已声明endpoint-star中心分配，不能冒充唯一逐键裂面功。实际承载接触记录{contact['force_bearing_records']}，零力r=dx记录{contact['zero_force_touching_records']}，最大重叠{contact['max_overlap_dx']:.9g}dx。主段墙钟{r['main_wall_s']:.3f}s/CPU{r['main_CPU_s']:.3f}s，峰值RSS{r['peak_RSS_bytes']/1024**3:.5f}GiB；本次已按用户授权取消4GiB停止条件，不追溯覆盖上次资源失败与漏检记录。\n\n裂纹候选平面前缘推进判据触发={st['crack_propagation_candidate']}，新生唯一物理裂面尚未认证，**CRACK_AREA_UNRESOLVED；G_SEPARATION_PROTOCOL=UNAVAILABLE、GC_EFFECTIVE=UNAVAILABLE**。宏观强度、全局时间收敛和V3生产仍未认证。旧341点700µs成功基线、旧瞬断 **G0_THEORY=0.07170070833333J/m²、G0_DISCRETE=0.06969623734234J/m²**及所有历史保持。见[本次六组真实图](temp_storage/clot_github_sync_20261011/clot_clearance/results/{TASK}/{RUN}/OPEN_RESULTS.html)、[34项科研审查](temp_storage/clot_github_sync_20261011/clot_clearance/verification/{TASK}/{RUN}/REVIEW_{TASK.upper()}.md)。本次只完成用户授权的同工况追加运行和结果交付；不自动放宽任何数值门槛或再改几何/材料/载荷。\n\n'''
insert+=f'本次科学区分：{finding} 动能/储能约{e["last_K_over_U"]:.6g}，应归为动态损伤响应；接触重叠峰值为{100*contact["max_overlap_dx"]:.6g}%点距，不能排除其局部影响。边界归属审查保留10个开口网格单元、0.025mm²代理，严格horizon无桥接子集仍为零，不能用该代理或总D计算Gc。原墙钟预算每200接受步检查，本次在1826.010s附近归档停止，数值已接受到2.4ms，尚缺0.1ms保持段；此处不声明完整协议通过。\n\n'
(O/'WSL_DOCUMENT_INSERT.md').write_text(insert)
print('REPORT_AND_WSL_INSERT_WRITTEN')
