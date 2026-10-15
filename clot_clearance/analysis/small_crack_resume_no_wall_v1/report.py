from common import *
def read(n):return json.loads((O/n).read_text())
st=read('FINAL_STATUS.json');n=read('NUMERICAL_VALIDITY.json');e=read('DYNAMIC_ENERGY_LEDGER.json');r=read('RUN_RESOURCES.json');ev=read('FIRST_SOFTENING_AND_FAILURE.json');loc=read('CRACK_LOCALIZATION_AUDIT.json');top=read('NATIVE_COMPONENT_TOPOLOGY.json');bridge=read('CROSS_PLANE_BRIDGE_AUDIT.json');area=read('CRACK_AREA_IDENTIFIABILITY.json');contact=read('CONTACT_STREAM_AUDIT.json');chain=read('CONTINUATION_CHAIN_AUDIT.json')
parent=json.loads((PARENT/'FINAL_STATUS.json').read_text());last=e['last'];plane=bridge['last_plane'];names={'NOTCH_TIP_INTERIOR':'裂尖厚度中央','NOTCH_TIP_SIDE_SURFACE':'裂尖侧表面','CRACK_FLANK':'裂纹侧面','GRIP':'加载端','OTHER_FREE_SURFACE':'其他自由表面','REMOTE_INTERIOR':'远处内部'}
energytext='，'.join(f'{k[:-2]}={v*1e9:.12g}nJ' for k,v in last.items())
isolated=[z['atom_ids'][0] for z in top['components'] if z['node_count']==1];tipD=sum(z['star_center_D_J'] for z in loc['regions'][:2]);lead=loc['notch_softening_lead_over_grip_s']
def et(q):
 if q is None:return '未观察到。'
 return f'{q["time_s"]*1e6:.9g}µs，原始bond ID '+','.join(map(str,q['bond_ids']))+'；'+','.join(sorted(set(names[z] for z in q['regions'])))+'。参考与实际中点、真实伸长及κ/d均保存在事件JSON。'
finding=f'完整2.5ms原生非均匀动力学数值准入通过；裂尖优先软化并出现局部损伤，但持续前缘推进判据未通过。末态候选开口带仅{plane["max_front_advance_m"]*1e6:g}µm，跨面桥接仍有{bridge["last_live_cross"]}条，其中{bridge["last_long_live_cross"]}条为r₀>dx长键。{len(top["components"])}个连通分量的节点数为{top["component_sizes"]}，小分量均为单节点，不是沿候选裂面形成的两个宏观分离体。新裂面与分离能仍不可认证。'
qa=[
('最终采用什么三维几何？','同一14×10×22胞元中心矩形试样，胞元体积尺寸700×500×1100µm；点中心跨度650×450×1050µm。没有新几何或完整Y形血栓求解。'),
('节点数和原始键数？','3080原始节点、144448原始参考键，全部ID保留，dx=50µm、δ=150.5µm。'),
('初始预裂纹、前缘及预切数？','原始z=0，x∈[−350,−50]µm、贯穿y±250µm。按原始线段交点预切3524键；初始投影0.15mm²，既不计入自然失效，也不计入后续D和新增面积。m保留全部原始参考邻域。'),
('相对历史小试样的设计依据？','剩余韧带2.66δ、厚度3.32δ、裂尖到受控面3.49δ，预注册区分裂尖中央、侧表面与端部；吸取M5B/M5C加载端竞争和桥接教训，但仍需承认有限厚度和horizon效应，不能由设计本身认证裂面。'),
('真实加载及本次时程？',f'上下每侧120µm半余弦2ms张开，再保持0.5ms。每侧最大规定速度0.0942478m/s；本次从2.4ms原生restart继续到{st["actual_end_s"]*1e3:.9g}ms，仅新增{chain["additional_steps"]}步，没有重跑t=0或延长冻结物理窗口。'),
('自由和受控DOF？','上下各140点只规定z，所有点x/y自由；2800内部点z也自由。8960自由DOF、280规定DOF，原生VV与fix move NULL分量没有重复积分。'),
('首次软化在哪里？',et(ev['first_softening'])),
('首次完整断键在哪里？',et(ev['first_complete_failure'])),
('裂尖是否先于夹持区损伤？',f'是，裂尖首次软化比加载端早{lead*1e6:.9g}µs。首次事件从父运行无改写继承，继续按预注册区域判断；后期大量其他区域损伤不能归入裂尖传播。'),
('是否观察局部损伤集中？',f'是，{ev["first_localized_cloud"]["time_s"]*1e6:g}µs触发冻结判据；当时{ev["first_localized_cloud"]["damaged_bonds"]}个自然损伤键，裂尖占比{ev["first_localized_cloud"]["notch_fraction"]:.0%}，{ev["first_localized_cloud"]["connected_voxels"]}相连50µm体素。其含义是局部损伤云，不是唯一连续裂面。'),
('是否连续推进裂纹前缘？',f'未通过冻结的两个采样时刻、中央条带推进至少2dx=100µm判据。末态无桥接/实际开口/缺口连续前缀只延伸{plane["max_front_advance_m"]*1e6:g}µm；不能把1.375ms首次无桥接开口单元称为连续裂纹传播。'),
('是否仍有跨面桥接？',f'8100条原始跨面记录，预切后4576存活；末态仍有{bridge["last_live_cross"]}条，新增完全失效跨面{bridge["natural_failed_cross"]}条。全部参考长度均纳入，其中{bridge["last_long_live_cross"]}条存活长程桥接。'),
('真实连通分量是什么？',f'{top["component_sizes"]}，单节点原始atom ID为{isolated}。末态原生端点集合与独立全部d<1原始键图完全一致；部分软化键仍保持连接。主分量点中心包络为'+ '×'.join(f'{v*1e3:.6g}' for v in top['components'][0]['center_span_m'])+'mm。单节点包络为零不代表物理碎片体积为零，其指定体积仍为dx³；逐分量质量、位置和ID已保存。'),
('端部、表面或接触是否有影响？',f'首次损伤在裂尖侧表面，后期加载端及其他自由表面大量损伤。最大接触重叠{contact["max_overlap_dx"]:.9g}dx，即{100*contact["max_overlap_dx"]:.6g}%点距，不能排除其局部影响；本轮没有修改接触刚度或增加诊断工况。'),
('预裂纹之外是否出现损伤？',f'不跨原候选面的自然完全失效有{loc["off_plane_natural_failed"]}条。区域数量和D分配见下表；这些是相互作用统计，不能直接称为另一连续裂纹或面积。'),
('新增分离区域可追溯吗？','所有原始交点、跨面桥接、对应上下节点真实开口和缺口前缀单元均可追溯。当前只能确认局部网格代理和4个单节点分量，不能据此认证可用于裂面耗能的连续新分离面。'),
('能否测量投影分离面积？',f'CRACK_AREA_UNRESOLVED；PROJECTED_SEPARATION_AREA=null。末态网格代理{area["notch_connected_grid_proxy_m2"]*1e6:g}mm²，严格horizon无桥接子集{area["conservative_grid_proxy_m2"]*1e6:g}mm²。0.025mm²来自10个满足几何规则的50µm单元，绝不是断键数×dx²，也不被冒充唯一裂面。'),
('有限horizon不确定性？','δ=150.5µm大于仅50µm的候选推进带。保守子集要求整个单元的δ邻域无存活桥接；精确网格边界交点另按所有相邻单元计数，仍为10个代理单元、零保守单元。冻结floor规则没有被结果后更改。'),
('完整能量账本是多少？',f'末态{energytext}。W=ΔU+ΔK+D+ΔUc−Jc+R，以原始t=0为能量起点；W只包含受控z反力功和受控质量动能变化。末态K/U={e["last_K_over_U"]:.9g}，保持段自由动能仍显著，属于动态损伤响应，不是准静态试验。'),
('损伤耗能来自哪些区域？',f'按一致endpoint-star中心分区，尖端两区D占全局{100*tipD/last["D_J"]:.6g}%。分区和与原生全局D之差{e["regional_sum_error_J"]:.9g}J；这只是声明的状态型耦合能量归属，不是唯一逐键裂面功。'),
('协议相关分离耗能可用吗？','G_SEPARATION_PROTOCOL=UNAVAILABLE。连续推进、可信新增面积和裂纹相关耗能归属未同时满足；不能用总D除以任意截面。'),
('为什么不是材料Gc？','GC_EFFECTIVE=UNAVAILABLE。完整数值准入、多个图分量和局部损伤云都不能替代唯一物理裂面与实验材料认证；动态惯性、接触和非局部边界效应仍在。'),
('宏观强度是否保持不变？','未独立认证。保持输入G/K/s0/sc不等于证明宏观峰值强度不变。'),
('材料及模型是否保持？','G=940Pa、K=100000/3Pa、ρ=1000kg/m³、s0=0.15、sc=0.18、α=0、dx=50µm、δ=150.5µm，以及完整m/r0/V/taper、1.0dx接触和变分软化函数全部不变；未新编译C++。'),
('时间步和事件积分是否通过？',f'完整0–2.5ms已接受{n["accepted_steps"]}步，损伤拒绝{n["rejected_steps"]}次；dt={n["dt_min_s"]*1e6:g}–{n["dt_max_s"]*1e6:g}µs，最大Δd={n["max_dd"]:.12g}≤0.10。保持原事件细分机制和全部数值门槛；未认证全局dt收敛。'),
('状态或原始身份是否错位？','恢复时x/v/κ/d/步数/时刻与2.4ms完整检查点一致，原生端点图、全部原始bond ID和m复核通过。逐步独立力/能量/κd检查继续通过；旧9601行CSV、节点流和损伤流前缀逐字节保持，未重复或跳过接受步，恢复时R跳变量为零。'),
('最大能量余额？',f'全协议最大归一化余额{e["max_relative_R_percent"]:.12g}%，最大绝对R={e["max_absolute_R_J"]:.12g}J，低于原2%门槛。R独立保留，不转为物理D。'),
('1.0dx接触是否影响解释？',f'全程几何候选{contact["geometric_candidate_records"]}条时刻记录，承载r<dx记录{contact["force_bearing_records"]}条，r=dx零力记录{contact["zero_force_touching_records"]}条；最大并发{contact["max_force_bearing_pairs"]}对，最大重叠{contact["max_overlap_dx"]:.9g}dx。数值账本通过不等于接触物理精度已标定。'),
('实际计算成本？',f'本次400步续算墙钟{r["segment_wall_s"]:.6f}s、CPU{r["segment_CPU_s"]:.6f}s、峰值RSS{r["segment_peak_RSS_bytes"]/1024**2:.6f}MiB；1rank×1线程。与父运行合计墙钟{r["protocol_total_wall_s"]:.6f}s、CPU{r["protocol_total_CPU_s"]:.6f}s，全协议峰值{r["protocol_peak_RSS_bytes"]/1024**3:.6f}GiB。取消墙钟和RSS上限，无替代上限；2GiB输出门槛保留。'),
('科研图是否生成？','六组PNG/SVG及OPEN_RESULTS.html：图1复用同一冻结几何；事件图、演化图、桥接图、反力/能量图和末态图使用完整连续轨迹及实际末态，不放大位移、不造医学外壳。'),
('生产库和历史是否保持？','复用SHA730e11651aee6533b4fb7d38c8341e06ad853bf3397c4908fc6ef905d4e3ac90隔离库；生产库和旧物理归档不变。父运行125文件在恢复前、链核验和交付前核对；仅当前开发入口在备份后推进。此前内存和墙钟停止历史不被追溯抹去。'),
('WSL文档是否实际更新？','本轮交付时实际插入/home/lzy/projects/clot_rupture.md，保留此前2.4ms未完成事实及所有历史章节。before/after/diff/SHA、并发检查及实际读回见WSL_DOCUMENT_UPDATE.json；不把旧瞬断G0改成新软化Gc。'),
('最重要的发现？',finding+f'本次额外100µs新增{st["natural_failed"]-parent["natural_failed"]}条自然完全失效，并形成4个单节点分量；跨面桥接数量没有进一步减少。完成时程只提升完整数值验收范围，不改变裂面/耗能不可辨识结论。'),
('下一轮最值得做什么？','优先审查为何裂尖局部开口止于50µm，而端部/其他表面持续耗能，并评估接触重叠对小分量的影响；据此决定独立的几何/载荷或接触验证。当前不自动延长物理时间、不调材料、不运行更大血栓或CFD/微泡/AI工况。')]
assert len(qa)==34
text=f'''# 取消墙钟预算后的原生续算与完整时程审查

任务：{TASK}；Run `{RUN}`。用户授权：“不设置1800秒墙钟预算，继续进行运行”。

**{st['status']}：完整2.5ms时程已完成。** 本次仅从父运行的完整2.4ms原生检查点追加400步，原物理设计、材料与全部数值门槛不变。墙钟与RSS上限均为null，没有替代上限；2GiB输出门槛保留。没有重跑前2.4ms、没有第二种参数工况。

{finding}

## 连续性与数值验收

- 检查点恢复时未推进物理步；节点位置/速度、κ/d、原始ID、m及原生存活端点图核对通过。
- 累计W/D/Jc、接触功、能量归一化尺度、四项数值缺陷、区域D和原有最大误差恢复，原t=0能量基准保持；R跳变量为0。
- 父运行0–2.4ms原始CSV/FP64/损伤/接触流保持精确前缀；新接受步9601–10000连续追加，原首次事件不重计。
- 最大Δd={n['max_dd']:.12g}，最大归一化余额{e['max_relative_R_percent']:.12g}%，原生末态端点图与独立原始图核对通过。没有数值门槛放宽。

| 项目 | 完整协议结果 |
|---|---|
| 实际物理时间 | 2.5ms（原计划全部完成） |
| 接受 / 损伤拒绝 | {n['accepted_steps']} / {n['rejected_steps']} |
| 自然软化 / 完全失效 | {st['natural_softened']} / {st['natural_failed']} |
| 初始预切 | 3524，单独排除 |
| 分量节点数 | {top['component_sizes']} |
| 存活跨面键 / 长键 | {bridge['last_live_cross']} / {bridge['last_long_live_cross']} |
| 末态K/U | {e['last_K_over_U']:.9g}，动态响应 |
| 本次续算墙钟 / 峰值RSS | {r['segment_wall_s']:.6f}s / {r['segment_peak_RSS_bytes']/1024**2:.3f}MiB |
| 连续协议累计墙钟 / 峰值RSS | {r['protocol_total_wall_s']:.6f}s / {r['protocol_peak_RSS_bytes']/1024**3:.6f}GiB |
| 裂面 / 分离能 / Gc | UNRESOLVED / UNAVAILABLE / UNAVAILABLE |

## 34项科研回答

'''
text+='\n\n'.join(f'### {i}. {q}\n\n{a}' for i,(q,a) in enumerate(qa,1))
text+='\n\n## 实际区域统计\n\n| 预注册区域 | 自然软化（含失效） | 自然完全失效 | star中心D/nJ |\n|---|---:|---:|---:|\n'
for z in loc['regions']:text+=f'| {names[z["region"]]} | {z["naturally_softened"]} | {z["naturally_failed"]} | {z["star_center_D_J"]*1e9:.12g} |\n'
text+='\n完整数值协议通过与连续裂纹面积—耗能认证分开验收。现有模型出现可追溯的裂尖局部损伤和4个单节点分量，但仍未建立满足冻结准则的连续前缘推进、可信新增分离面积及对应耗能。旧G0_THEORY=0.07170070833333J/m²、G0_DISCRETE=0.06969623734234J/m²仍仅为旧瞬断LPS标称值。\n\n[六组科研图](OPEN_RESULTS.html)。本次到原定2.5ms结束，不继续扩大规模或延长时间。\n'
(O/('REVIEW_'+TASK.upper()+'.md')).write_text(text)
insert=f'''\n**取消1800秒墙钟预算后的原生续算与完整时程（{RUN}）**：根据用户“不设置1800秒墙钟预算，继续进行运行”的指令，墙钟和RSS上限均取消，不设替代上限；2GiB输出及原损伤、能量、运动学门槛保留。直接恢复上轮2.4ms完整原生检查点与累计W/D/Jc、接触功、归一化尺度、四项数值缺陷、区域D和事件记录，R跳变量为0；原0–2.4ms轨迹与事件保持精确前缀，不重跑、不重计。原3080点/144448键/3524预切、700×500×1100µm几何、G940Pa、K100000/3Pa、ρ1000kg/m³、s₀=0.15、s_c=0.18、δ=150.5µm、dx=50µm、1.0dx接触及120µm每侧2ms半余弦+0.5ms保持均不变，无C++重编译。\n\n本次追加400步，完成**2.5ms完整协议，FULL_PROTOCOL_NUMERICAL_PASS**；累计10000接受步/0损伤拒绝，dt=0.25µs，最大Δd={n['max_dd']:.12g}，最大归一化能量余额{e['max_relative_R_percent']:.12g}%（原2%门槛）。全部原始节点和键身份保持，原生末态端点图与独立全部存活原始键图核对通过。首软化747.5µs、首完全失效799.75µs，均为裂尖侧表面原始键68923/77580；裂尖软化领先加载端91.25µs。完整末态自然软化{st['natural_softened']}、完全失效{st['natural_failed']}，分量点数{top['component_sizes']}，四个单节点ID{isolated}；新增100µs内完全失效增加{st['natural_failed']-parent['natural_failed']}条。\n\n{finding} 候选网格代理仍为0.025mm²、严格horizon无桥接子集为0；**CRACK_AREA_UNRESOLVED，G_SEPARATION_PROTOCOL=UNAVAILABLE，GC_EFFECTIVE=UNAVAILABLE**。末态{energytext}，K/U={e['last_K_over_U']:.9g}，属于动态损伤响应；最大接触重叠{contact['max_overlap_dx']:.9g}dx，不能排除其局部影响。区域D只按endpoint-star中心分配，不是唯一逐键裂面功。\n\n本次续算墙钟{r['segment_wall_s']:.3f}s/CPU{r['segment_CPU_s']:.3f}s、峰值RSS{r['segment_peak_RSS_bytes']/1024**2:.3f}MiB；与原2.4ms段合计墙钟{r['protocol_total_wall_s']:.3f}s、峰值RSS{r['protocol_peak_RSS_bytes']/1024**3:.6f}GiB。保留此前50µs内存停止及2.4ms墙钟停止的历史事实，不追溯改写。完整数值通过不认证宏观强度、全局dt收敛、实验材料Gc或V3生产适用性。旧341点700µs成功基线、旧瞬断G0_THEORY=0.07170070833333J/m²、G0_DISCRETE=0.06969623734234J/m²均保留。见[完整时程科研图](temp_storage/clot_github_sync_20261011/clot_clearance/results/{TASK}/{RUN}/OPEN_RESULTS.html)、[34项完整审查](temp_storage/clot_github_sync_20261011/clot_clearance/verification/{TASK}/{RUN}/REVIEW_{TASK.upper()}.md)。到原定2.5ms后停止，不自动延长加载或扩大试样。\n\n'''
(O/'WSL_DOCUMENT_INSERT.md').write_text(insert)
(O/'REPRODUCE.txt').write_text(f'''Native continuation of the same frozen specimen from 2.4ms to END=2.5ms. No wall-clock or RSS ceiling; output and numerical gates unchanged.
Parent: {PARENT}
Parent manifest SHA: {sha(PARENT/'PAYLOAD_MANIFEST.json')}
Library: {LIB}, SHA {sha(LIB)}. No new C++ build.
Activate source /workspace/y_clot/env.sh. Export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1.
Set LD_LIBRARY_PATH to include {LIB.parent}; set PYTHONPATH to include {NATIVE_PREV/'source/python'}.
Recorded physical command: python -u {A/'launch.py'} full. Archive guards prohibit overwriting this completed run.
prepare_resume.py verifies the parent payload and checkpoint; it copies the old accepted trajectory without hardlinks and freezes new driver identities before the native resume. The restored native state and complete cumulative mechanics ledger are checked before the first new step; original t0 remains the energy origin.
Postprocessing order: chain_audit.py, analyze.py, geometry_review.py, figures.py, report.py. Run only in an isolated redirected copy if reproducing; these scripts write output.
NATIVE_RESTART_RESTORE_AUDIT.json, LEDGER_RESTORE_AUDIT.json and CONTINUATION_CHAIN_AUDIT.json prove the checkpoint handoff and exact old raw prefixes. Gzip damage/contact files have concatenated members; standard Python gzip.open reads the complete stream. KEYFRAME_INDEX contains historical stop keys and new final keys, with no repeated accepted node frame.
Full raw FP64 layout is in full/STREAM_LAYOUT.json. Every accepted scalar and d change is retained; subthreshold kappa is only in complete saved states. No unsaved field is reconstructed by interpolation.
Latest native checkpoint is at 2.5ms; original physical protocol is complete. Further physical time would be a changed protocol and was not executed.
WSL update is performed once by update_wsl_document.py with backup, SHA guard and actual readback. seal.py verifies provenance and protects the immutable parent evidence before updating the active entry and writing the payload manifest.
''')
print('REPORT_34_AND_WSL_INSERT_READY')
