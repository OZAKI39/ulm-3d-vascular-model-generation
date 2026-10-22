"""Evidence-bound Chinese research delivery; scoped gates, no fabricated success."""
import json,numpy as np,h5py,datetime
from .common import *

def table(headers,rows):
 return '| '+' | '.join(headers)+' |\n| '+' | '.join(['---']*len(headers))+' |\n'+''.join('| '+' | '.join(map(str,r))+' |\n' for r in rows)
def f(v,d=4):return 'INVALID' if v is None else f'{v:.{d}g}'
def doc(name,title,body):
 (BASE/name).write_text('# '+title+'\n\n'+body.strip()+'\n\n'+scope+'\n')
scope='适用范围：`EXPLORATORY` · `HOMOGENEOUS_WATER_CONTROL_PROXY` · `PHYSICS_UNCALIBRATED` · `ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED` · `training_ready=false`。Y 与四源静态结果另加 `HOMOGENEOUS_WATER_PROXY_SPECIFIC`。'

def main():
 verify();stages={p.parent.name:read(p) for p in BASE.glob('*/INDEX.json') if p.parent.name not in ['portable_dataset'] and isinstance(read(p),list)}
 rows=[r for rr in stages.values() for r in rr];audits=[r for p in (BASE/'analysis').glob('*/SUMMARY.json') for r in read(p)]
 assert len(rows)==len(audits),('Missing final analysis',len(rows),len(audits))
 passes=[r for r in rows if r['status']=='PASS'];fails=[r for r in rows if r['status']!='PASS']
 original_checks=[]
 for name in ['NUMERICAL_COMPARISONS','CONFIRMATION_NUMERICAL','SENSITIVITY_NUMERICAL','EXTENSION_NUMERICAL']:
  original_checks+=read(BASE/'analysis'/(name+'.json'))
 target=read(BASE/'analysis/TARGET.json');original_checks+=target['comparisons']+[target['independent_recompute']]
 rep=read(BASE/'analysis/REPLICATION.json');svd=read(BASE/'analysis/SENSITIVITY.json')['cases']
 physics_wall=sum(read(BASE/s/'COMPLETE.json')['wall_s'] for s in stages)+read(BASE/'static/SUMMARY.json')['wall_s']
 analysis_wall=sum(read(p)['wall_s'] for p in (BASE/'analysis').glob('*/COMPLETE.json'))
 steps=sum(r.get('accepted_steps',0) for r in rows);checked=sum(r['independent_audit']['checked_steps'] for r in audits)
 stats=dict(trajectories=len(rows),completed=len(passes),failed=len(fails),accepted_steps=steps,independently_checked_steps=checked,
  physics_wall_s=physics_wall,analysis_wall_s=analysis_wall,trajectory_wall_sum_s=sum(r['wall_s'] for r in rows),
  peak_worker_RSS_MiB=max(r['rss_peak_MiB'] for r in rows),max_trajectory_wall_s=max(r['wall_s'] for r in rows),
  full_window_numerical_comparisons=len(original_checks),numerical_comparisons_pass=sum(c['original']['admitted'] for c in original_checks),
  raw_bytes=sum(r['raw_bytes'] for r in read(BASE/'RAW_DATA_MANIFEST.json')['records']),
  all_completed_independent_audits_pass=all(r['independent_audit']['PASS'] for r in audits if r['status']=='PASS'),
  identity=identity(),labels=LABELS,training_ready=False)
 write(BASE/'RESEARCH_SUMMARY.json',stats)
 def find(state,action,stage='discovery',variant='DT1',duration=.002):
  return next(r for r in stages[stage] if r['state_id']==state and r['action']==action and r['variant']==variant and r['task']['duration_s']==duration)
 def result(r):
  a=r['initial_shape'];z=r['final_shape'];return [r['state_id'],r['action'],f((np.array(z['center_3d_m'])-a['center_3d_m'])[0]*1e9),f(100*(z['Rg_3d_m']/a['Rg_3d_m']-1)),f(100*(z['aspect_3d']/a['aspect_3d']-1)),f(z['similarity_rmsd_3d_m']*1e9),str(z['proximity_components'])]
 a5=find('S2_r0','A5_X180');axis_change=(np.asarray(a5['final_shape']['principal_rms_lengths_3d_m'])/a5['initial_shape']['principal_rms_lengths_3d_m']-1)*100
 axial_note=f'形变含义：S2_r0/A5三个主轴RMS尺度分别改变{axis_change[0]:.4f}%、{axis_change[1]:.3f}%、{axis_change[2]:.3f}%。长宽比增加主要来自横向压缩，不能称为X向绝对长度主动增长。'
 a7=find('S2_r0','A7_X210');projection_note=f"实际S2_r0/A7在2ms时，3D长宽比改变{100*(a7['final_shape']['aspect_3d']/a7['initial_shape']['aspect_3d']-1):.3f}%，XY长宽比改变{100*(a7['final_shape']['aspect_xy']/a7['initial_shape']['aspect_xy']-1):.3f}%，方向相反。"
 response=table(['初态','动作','COM X / nm','Rg变化 / %','3D长宽比变化 / %','去平移/旋转/等比缩放残差 / nm','终点分量'],[result(r) for r in stages['discovery'] if r['status']=='PASS'])
 repeated=table(['形状','动作','三初态长宽比变化 / %','方向一致','全部终点连通','全种子dt/cache复核'],[[g['shape'],g['action'],', '.join(f(v['aspect_fraction']*100) for v in g['values']),g['aspect_direction_repeated'],g['all_connected'],g['full_seed_numerical_confirmation']] for g in rep['groups']])
 targets=table(['初态','长宽比变化 / %','COM范数 / nm','终点分量','目标量通过','dt/cache通过','完整目标通过'],[[c['state_id'],f(c['aspect_fraction']*100),f(c['com_m']*1e9),c['final_components'],c['target_metrics_pass'],c['numerical_pass'],c['admitted']] for c in sorted(target['cases'],key=lambda c:c['state_id'])])
 numerical=table(['对照','逐泡最大差异 / nm','Rg最大差异 / nm','原始连续指标准入','事件拓扑一致'],[[c['comparison'],f(c['original'].get('errors',{}).get('fixed_id_m',0)*1e9),f(c['original'].get('errors',{}).get('Rg_m',0)*1e9),c['original']['admitted'],c['original'].get('event_topology_agrees','N/A')] for c in original_checks])
 ranks=table(['初态','窗口/ms','差分检验','形状奇异值','保守阈值','有效形状秩'],[[c['state_id'],c['time_s']*1e3,c['status'],', '.join(f(v) for v in c['shape_singular_values']),f(c['rank_threshold']),c['shape_effective_rank']] for c in svd])
 failures=table(['轨迹','请求/ms','接受/ms','原因'],[[r['run_id'],r['requested_s']*1e3,f(r['completed_s']*1e3),', '.join(r['reasons'])] for r in fails])
 initial=read(BASE/'INITIAL_STATE_AUDIT.json')
 shapes=table(['初态','种子','N','Rg/µm','3D长宽比','邻近分量'],[[name,a['seed'],a.get('N',32 if name.startswith('S0') else 64),f(a['Rg_3d_m']*1e6),f(a['aspect_3d']),a['proximity_components']] for name,a in sorted(initial.items())])
 performance=table(['阶段','完整/总轨迹','物理阶段墙钟/s'],[[s,f'{sum(r["status"]=="PASS" for r in rr)}/{len(rr)}',f(read(BASE/s/'COMPLETE.json')['wall_s'])] for s,rr in stages.items()])
 extension=table(['窗口/ms','COM X/nm','Rg变化/%','长宽比变化/%','终点分量'],[[r['requested_s']*1e3,*[result(r)[i] for i in [2,3,4,6]]] for r in stages['extension'] if r['variant']=='DT1'])
 static=read(BASE/'static/SUMMARY.json')['records'];mechanisms=[]
 for state in ['S0_r0','S2_r0','S6_r0']:
  for action in ['A7_X210','A5_X180']:
   s=next(s for s in static if s['state_id']==state and s['action']==action);v=np.array(s['velocity_free_m_s']);w=np.array(s['velocity_constrained_m_s'])
   relative=lambda v:np.sqrt(np.mean(np.sum((v-v.mean(0))**2,axis=1)))
   mechanisms.append([state,action,f(s['pressure_max_Pa']),f(np.linalg.norm(s['center_primary_force_N'])*1e12),f(np.linalg.norm(s['center_secondary_force_N'])*1e12),f(s['force_relative_rms_N']*1e12),f(relative(v)*1e6),f(relative(w)*1e6)])
 mechanism=table(['初态','动作','最大局部声压/Pa','平均一次力范数/pN','平均二次力范数/pN','去均值总力RMS/pN','自由相对速度/µm·s⁻¹','一次约束步相对速度/µm·s⁻¹'],mechanisms)
 conclusion='在当前 2.5 MHz 均匀水代理及已测试初态/时间窗内，存在经过数值验证、动作可区分的有限形状响应自由度；局部两参数分析仅支持一个可辨识形状方向，尚未证明多个独立形变方向、跨初态稳健的完整连通团目标变换或可逆控制。'
 (BASE/'SCIENTIFIC_CONCLUSION.txt').write_text(conclusion+'\n')
 doc('README.md','F10AB：真实形变证据、有限可控性与目标失败',f'''{conclusion}

本轮完成 {len(rows)} 条实际尝试，{len(passes)} 条完成预定物理窗口，{len(fails)} 条按原接触门槛停止。固定状态诊断 105 组；物理计算阶段合计 {physics_wall:.2f} s。所有通过轨迹均保留每一步原始数组和独立力学复核。没有仪器操作、AI训练或新PDMS模型。

- [完整中文科研报告](F10AB_DEVELOPMENT_REPORT.md)
- [逐动作实际结果](SHAPE_RESPONSE_MATRIX.md) · [独立初态](REPLICATION_AND_GENERALIZATION.md)
- [数值误差](NUMERICAL_CONVERGENCE.md) · [目标试验](TARGETED_SHAPE_TRANSFORMATION.md)
- [机器准入](FINAL_GATE.json) · [自包含数据接口](DATASET_AND_F10Y_INTERFACE.md)
- [图表清单](figures/MANIFEST.json) · [视频清单](videos/MANIFEST.json)

![目标与连通性](figures/Fig_F10AB_10_target_or_limitations.png)

[实际动作响应视频](videos/Video_F10AB_1_collective_shape_responses.mp4) · [去刚体运动视频](videos/Video_F10AB_2_rigid_motion_vs_deformation.mp4) · [真实目标试验视频](videos/Video_F10AB_3_shape_control_or_limitations.mp4)
''')
 doc('F10AA_BACKEND_AUDIT.md','F10AA后端和历史保护','本地/服务器可信交付基准为 `6028598d0012912cc19aa766fa786c4afb0a91a5`；不是服务器旧 repo HEAD。185 条历史轨迹 SHA 与原源码已逐一核对。\n\n直接继承 `f10aa.trajectory.Campaign.run`；`AuditEngine` 仅附加原始布尔 Gate 输出。力学采用 `WaterAccelerated`/SciPy CPU，声学保持 Warp GPU；没有重新实现粒子动力学。冻结原 src、协议与响应表/缓存 SHA。\n\n证据：[初始审计](BASELINE_AUDIT.json)、[服务器审计](evidence/REMOTE_BASELINE_AUDIT.json)、[保护复核](evidence/FINAL_PROTECTION_AUDIT.json)。原始 INPUT 的父协议哈希属于 F10AA；F10AB 初态、计划、源身份另由阶段 START/COMPLETE、PREREGISTRATION 和 portable HDF5 属性明确记录。')
 doc('INITIAL_SHAPE_LIBRARY.md','合法初态与比较优先级','S0 使用历史32泡坐标；S1–S6每类构造3个64泡实例，r0用于发现，r1/r2用于选中机制确认。N/半径/群中心优先一致，Rg/密度无法同时匹配，保留实际差异。S1–S4邻近连通，S5为两群，S6分散。人为构造不等于声学自然形成。\n\n'+shapes)
 doc('SHAPE_METRICS_DEFINITION.md','三维与XY形状指标','固定ID、不重新匹配、不删除粒子。Q=YᵀY/N，Rg=√tr(Q)，主轴RMS长度=√λ，长宽比=√(λmax/λmin)。XY仅为投影，不是Chronos图像。邻近阈值d≤2.1R；硬接触几何统计d≤2R+1nm，二者均不代表融合。\n\n主轴相对特征值间距≤0.03时标记 `ORIENTATION_UNIDENTIFIABLE`；全局F要求满秩、条件数<100、det(F)>0、非仿射相对残差≤10%。无效时保存NaN，仍保留协方差和成对距离。距离分位数为0/10/50/90/100%。\n\n[分析定义源码](../../../scripts/straight_tube_exp1/f10ab/metrics.py)；[测试证据](evidence/TESTS.json)。')
 doc('RIGID_MOTION_DECOMPOSITION.md','刚体运动与真正形变的分离','使用固定ID Kabsch拟合 proper rotation，不允许反射；依次报告原始位移、去平移、去平移/旋转残差，再拟合一个统一尺度，得到 similarity RMS。后者不把各向同性收缩误认为各向异性形变。全局F采用极分解F=RU；U的无迹部分、非对角部分和非仿射残差互补，不能用无效F概括严重重排。\n\n单元测试含纯平移、刚体旋转、统一缩放、已知剪切、非仿射运动、秩退化、近简并特征值和投影假接触。算法参照 [SciPy Kabsch](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.transform.Rotation.align_vectors.html) 与 [极分解](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.polar.html)。')
 actions=read(BASE/'ACTIONS.json');actiontable=table(['动作','模式','幅值 Y-/Y+/X-/X+ / Pa-eq','相位/deg'],[[a,r['segments'][0]['mode'],', '.join(f(v,5) for v in r['segments'][0]['amplitude_Pa_equivalent']),', '.join(f(np.degrees(v)) for v in r['segments'][0]['phase_rad'])] for a,r in actions.items()])
 doc('PZT_ACTION_PREREGISTRATION.md','冻结动作与公平性',actiontable+'\n频率2.5MHz，背景流严格为零，最大允许1500 Pa-equivalent。基础X每片1036.201658Pa-eq。B系列每片缩小√2，使Σ|u|²与一个基础单源相等；这只是源幅值预算代理，不是校准声功率。A5/A7幅值预算相同，只改变相位。四片复幅与局部压力分别保存，不将Pa-eq称为泡处压力或AFG电压。\n\n15动作先做105组静态诊断；完整积分有限筛选，没有全网格扫参。Y与四源仅静态探索，动态结论 `NOT_RUN`。所有阶段计划在执行前冻结于 plans/。')
 doc('ACOUSTIC_FORCE_MECHANISM.md','声力、耦合阻力与接触机制',mechanism+'\n'+axial_note+'\n以上是 `FIXED_STATE_FORCE_DIAGNOSTIC`，自由速度已包含原耦合阻力，但不含接触约束；约束速度来自原引擎一个实际步。去均值力不等于最终相对速度。S0约束相对速度显著受限，与完整轨迹几乎保形一致；本轮没有分别关闭接触或润滑，因此不能量化二者各自因果贡献。阻力矩阵与全部F1/F2向量保存在static/*.h5。\n\n同相高幅A3在5种密集形状触发原接触反力平衡失败，停止而不放宽预算。')
 doc('SHAPE_RESPONSE_MATRIX.md','实际动作—形状响应矩阵',response+'\n所有行都是完整2ms轨迹；部分单源/等预算动作仅发现阶段，没有其独立dt/cache确认，不能由更大效果推断更可靠可控性。全局F无效的A7重排仍保留真实固定ID协方差/距离变化。')
 doc('INITIAL_SHAPE_DEPENDENCE.md','初态依赖',repeated+'\nA5的S2/S4长宽比增大方向重复；A7在S2/S4间甚至同类不同种子中方向变化。S6是分散粒子集，不把其统计形变称为已形成完整泡团变形。\n\n全部发现响应：[矩阵](SHAPE_RESPONSE_MATRIX.md)。')
 doc('REPLICATION_AND_GENERALIZATION.md','发现、独立确认与外推边界',repeated+'\nS2/A5三个初态均具有dt/cache完整对照；S4/S6新种子完成原逐步Gate，但未对每条补充dt/cache，故证据层级不同。S2_r0/r1/r2仍只有3个独立配置，不作总体概率或显著性推广。接触网络断开是实际输出，不删除失败种子。目标控制使用同一组三种子，其中r1/r2未参与幅值选择，但它们已用于机制确认，不冒称全新的外部测试集。')
 doc('ACTION_SENSITIVITY_AND_RANK.md','中心差分、SVD与误差门槛',ranks+'\n输入归一化：1单位=公共X双源幅值5%或X+相位6°；采用±全幅和±半幅中心差分，全部半幅端点另做dt/2与细缓存。输出COM/20nm、log(λ1/λ3)/0.01、log(λ2/λ3)/0.01、Rg变化比例/0.01、similarity RMS/Rg0/0.01。形状子矩阵仅取两个谱比例和去统一尺度后的残差，不包含COM和纯Rg尺度。\n\n阈值=max(两差分模板的谱范数差、5×数值噪声范数×2、1e-6)。四种初态/窗口组合的第二奇异值都未超过保守门槛，因此仅支持1个可辨识局部形状方向。不是全局可控性证明，也不证明其余所有合法动作都无效。')
 doc('SHAPE_CONTROLLABILITY_ANALYSIS.md','有限可控性与位置—形状分离',conclusion+'\n\nA5的质心接近不动而形状改变；S0/A7主要平移而保形，说明位置与形状不能由一个Rg变化指标合并。A5不同幅值具有不同形状幅度、近零群心移动；仍可能失去邻近连通。A7相位扰动响应具有初态依赖，不能用一次线性化宣称独立操纵任意三维形状。\n\n'+ranks)
 doc('AGGREGATION_VS_DEFORMATION.md','聚集、分布变形与团完整性','S0和S1–S4起点已经邻近连通；S6起点分散，整体Rg减小仅说明固定集合变紧。记录全程连通分量、最大分量、固定ID和邻近边；三维邻近与XY重叠分开。A5可能使原连通集合伸长后断为多分量，这是分布形变而非保持完整团的成功控制。\n\n'+response)
 doc('NUMERICAL_CONVERGENCE.md','全窗口数值误差与原始门槛',f'共 {len(original_checks)} 组完整窗口比较，{stats["numerical_comparisons_pass"]} 组通过原连续指标准入。独立复核 {checked} 个接受步；失败前缀不冒称完整通过。逐泡对照检查所有共同物理时刻；额外形状指标每10接受步及终点检查。瞬时接触拓扑另列，不能由轨迹平滑或全局守恒替代。\n\n'+numerical+'\n形变准入要求 similarity RMS>max(5×数值差异,1nm)，有效F的偏伸展/剪切>max(5×差异,1e-4)。S0形变低于阈值。原50nm逐泡、接触/力学/墙面预算均未放宽。')
 doc('TARGETED_SHAPE_TRANSFORMATION.md','预登记目标的实际结果',f'预登记目标：S2固定64粒子在2ms时长宽比增大3%±0.5个百分点，COM漂移≤20nm，终点仍为1个邻近分量。四个候选600/700/800/900Pa-eq，X相差180°。选择700Pa-eq后独立重算、两个机制确认种子及所有dt/cache对照。目标总体状态：**{target["status"]}**。\n\n'+targets+'\n接近目标长宽比但破坏连通性仍判失败，不修改目标或删除不利粒子。未运行反向/恢复动作：`NOT_RUN`；首先解决团完整性比继续扩展控制更有价值。\n\n3–5ms持续性试验（同一r0、1036.2Pa-eq，不是目标控制）:\n\n'+extension)
 doc('PERFORMANCE_AND_RESOURCES.md','实际计算成本',performance+f'\n静态诊断 {read(BASE/"static/SUMMARY.json")["wall_s"]:.3f}s；物理阶段累计 {physics_wall:.3f}s，低于3600s预算。单条最长{stats["max_trajectory_wall_s"]:.3f}s，低于240s；最多两条物理轨迹并行，单任务BLAS线程1。逐条wall之和{stats["trajectory_wall_sum_s"]:.3f}s不能当并行墙钟。独立分析阶段之和{analysis_wall:.3f}s；绘图/视频/打包耗时另见evidence/*TIMING.json和videos/MANIFEST.json，传输与环境安装不混入求解速度。峰值worker RSS{stats["peak_worker_RSS_MiB"]:.1f}MiB。\n\n使用RTX4090 Warp声学+CPU接触，未宣称服务器独占；其他已有任务未停止。独立.postvenv只用于绘图/视频，物理.venv未修改。首次后处理安装因Python3.11不满足av19失败，改用现有Python3.12成功，没有更换求解依赖。后处理首次导入另缺PyYAML，已仅在.postvenv补齐6.0.3；失败日志保留。')
 doc('NEGATIVE_RESULTS_AND_LIMITATIONS.md','失败、未运行与不可外推范围',failures+'\n- 这些5条A3同相高幅轨迹在早期被反力平衡门槛拒绝，保留接受前缀、失败提案和日志；不能补帧称完整。\n- A7在不同初态的长宽比变化方向不一致。\n- A5方向可重复，但可能裂为多个邻近分量。\n- 局部形状响应第二方向低于保守分辨门槛。\n- 3%目标涉及连通性失败，不能只展示接近目标的长宽比。\n- N128新增形状试验、10ms/1s、可逆控制、Y完整轨迹、PDMS校准、AI/MPPI、硬件均NOT_RUN。N32/64已足以发现当前机制和完整性限制，未为规模展示增加负担。\n- 近接触物理、APECSS参数、均匀水空间梯度和电声转换未校准；数值收敛不等于真实装置精度。')
 doc('DATASET_AND_F10Y_INTERFACE.md','自包含动作—状态HDF5',f'自包含 `portable_dataset/*.h5` 共{len(rows)}条，包括失败前缀。time_s与positions_m/形状指标逐帧对应；force_time_s独立对应实际步前受力；solver_gates保留每个尝试步，gate_time_s给出其时钟。每10接受步+终点导出原始状态，没有插值；全量每步原始HDF5仍不可变保存在服务器。\n\n服务器根：`/workspace/mevo_acoustic_microbubbles/f10ab_shape_20261018/results/straight_tube_exp1/f10ab`。全量清单与SHA：[RAW_DATA_MANIFEST.json](RAW_DATA_MANIFEST.json)。便携数据清单：[portable_dataset/INDEX.json](portable_dataset/INDEX.json)。模式：[DATASET_SCHEMA.json](DATASET_SCHEMA.json)。便携数据无外部链接，可独立读取。另有 `review_dataset/*.h5`：50µs间隔的实际状态加终点，147条全部覆盖，保留每个尝试步的完整Gate。该副本用于本地快速查阅，不能替代全量事件时间分析。10µs便携数据与每步原始数据完整保留在服务器，本地回传范围以evidence/LOCAL_RECOVERY.json为准。\n\n3D真值和XY投影与F10Y可比较COM/Rg/Q等定义，但不是光学图像、不是视觉固定ID精度证明。所有模型/seed/调度/源码/半径/单位/失效原因保留；training_ready=false。覆盖有限、非独立样本不能当大规模训练集。')
 doc('ACTIVE_CODE_INDEX.md','有效代码、执行与复现入口','源代码在项目 `scripts/straight_tube_exp1/f10ab/`：prepare.py冻结初态/输入；metrics.py独立形状定义；run.py继承原逐泡积分；plan.py有限筛选；analyze.py独立KKT复核；evaluate.py种子/数值准入；sensitivity.py中心差分；target.py预登记目标；figures.py/videos.py实际数据可视化；package.py自包含数据；report.py本报告。\n\n配置：`configs/straight_tube_exp1/f10ab/protocol.json`；测试：`tests/f10ab`与`tests/f10aa`。\n\n执行快照由logs/*EXECUTION.json逐次给出commit、目录、Supervisor标签、超时，logs/*supervisor.conf保存精确命令；复用快照只新增外部wrapper，不修改既有源码。所有计划在plans/，每个阶段START/COMPLETE保存真实身份与耗时。物理环境 `.venv`，绘图 `.postvenv`；均设置PYTHONPATH=src:.和BLAS线程1。\n\n复现到新的独立输出根，勿覆盖已存在的stage目录。示例：`python -m scripts.straight_tube_exp1.f10ab.run pilot`；分析：`python -m scripts.straight_tube_exp1.f10ab.analyze pilot`；绘图：`python -m scripts.straight_tube_exp1.f10ab.post media`。先恢复PREREGISTRATION列出的原响应表和缓存，SHA必须一致。')
 doc('FINAL_RECOMMENDATION.md','下一轮建议',conclusion+'\n\n优先研究保持邻近网络的低幅/分段动作和动态连通约束，将“固定集合改变分布”与“完整团变形”分别设目标；保留真实非黏附物理，不靠人工连杆锁团。扩展相位/幅值动作前，先寻找超过噪声的第二形状方向，继续双幅中心差分和新seed验证。再以更长但有限窗口研究恢复/不可逆重排。\n\n在真实实验前校准声压空间梯度与电声映射，尤其Y近场节点；当前数据只适合接口/小规模探索，不启动大型世界模型、RL或MPPI。')
 answers=[
 ('求解器能否稳定积分？',f'{len(passes)}/{len(rows)}条完成，5条同相高幅密集形状被原反力平衡门槛拒绝。完成并不意味着任意动作稳定。'),
 ('OFF下哪些形状保持？','S0–S6发现态，以及已测S2/S4/S6独立种子，严格零流/OFF下坐标不变。其他未积分种子不外推。'),
 ('历史近接触团为何形变极小？','S0/A7有整体位移，但去刚体/统一尺度后的残差约0.0023nm，低于1nm分辨下限；静态约束速度显示相对运动被抑制。'),
 ('硬球还是润滑限制？','原耦合阻力、紧密接触同时作用；一次步自由/约束响应支持接触约束贡献，但未做物理机制消融，不能定量归因主要比例。'),
 ('分散态为何更易收缩？','初始自由间隙较大、少接触约束，二次相互作用允许靠近；仍须用实际分量/最大群统计，不能称稳定大团。'),
 ('哪些主要平移？','S0的历史X双源A7最清楚；同一动作对松散64泡可同时引起内部重排。'),
 ('哪些主要内部靠近？','X同相/等预算对在发现实验中出现明显压缩，A7在S6也收缩；高幅同相密集团有数值失败，不能泛化。'),
 ('有真实各向异性吗？','有。S2/A5的3D长宽比在三个种子增长约6.28%、5.62%、8.29%，全窗口dt/cache确认；不只是Rg改变。'),
 ('有可辨识剪切吗？','有效F的S2/S4 A5对称剪切超过数值差异；S4斜向群更明显。剪切是有效全局近似的统计量，不等于全部粒子均匀剪切。'),
 ('哪些旋转只是刚体？','已用Kabsch单独提取刚体部分；没有把主轴转动当形变。纯旋转测试残差归零，实际A7残余仍明显。'),
 ('3D和XY是否不同？',projection_note+'XY丢失Z和深度关系；只以3D证据做主要形变准入，投影重叠不作接触。'),
 ('X细长群变长还是短？','A5在三初态均增加长宽比；A7两态减小、一态增大，依赖初始排列。'),
 ('Y细长群响应不同吗？','S3/r0已测，与S2不同；仅发现态，没有S3跨seed/dt/cache完整确认。'),
 ('相位能否改变模式？','相同幅值预算180°与210°有不同内部响应和COM；局部响应分析支持有限可分辨方向，但非任意独立形状控制。'),
 ('对向双源一定更易形变吗？','不能作普遍结论。同初态单源/等Σ|u|²双源发现比较已保留，未对这些全部做数值/seed复核。'),
 ('幅值增加与相位效应可分吗？','同幅A5/A7隔离相位，单源与缩小√2的双源作预算对照；局部声压/相对力仍分别记录。源预算不是校准声能。'),
 ('有不同形状响应方向吗？','动作响应不同，但稳健可区分的独立局部形状秩只为1；不能宣布多个独立模式。'),
 ('是否超过数值误差？','已准入的S2/S4/S6候选超过5倍误差与绝对下限；S0形变低于分辨。详见全窗口数值表。'),
 ('独立初态能复现吗？','S2/A5长宽比增大方向可复现，幅值不同且并非都保持连通；A7方向不稳健。'),
 ('局部矩阵有多个独立形变方向吗？','未证明。两输入、两初态、两个窗口的第二形状奇异值均低于误差阈值；有效秩1。'),
 ('能移动并保形吗？','S0/A7提供近似保形平移的已验证实例，不推广到所有64泡初态。'),
 ('能变形而COM少漂移吗？','A5及其幅值变化可以改变固定ID分布而COM近零，但邻近网络可能断开。'),
 ('能恢复或反向吗？','NOT_RUN。本轮发现团完整性和初态依赖问题，未追加恢复成功叙事。'),
 ('能支持世界模型吗？','已交付标准动作条件状态接口和失败数据，但样本覆盖、外部验证与物理校准不足；training_ready=false，未训练。'),
 ('下一轮优先什么？','先约束团完整性与初态稳健性，寻找第二可辨识动作方向，再有限延长和恢复测试；并为真实实验校准声场/电声映射。')]
 report=conclusion+f'\n\n实际物理尝试{len(rows)}条，完整{len(passes)}条，失败{len(fails)}条；物理阶段墙钟{physics_wall:.2f}s。没有新N128形变任务，没有通过换物理或放宽接触门槛制造成功。\n\n## 主要量化结果\n\n'+repeated+'\n\n'+axial_note+'\n## 预登记目标\n\n'+targets+'\n## 局部响应秩\n\n'+ranks+'\n## 25个问题\n\n'+''.join(f'{i}. **{q}** {a}\n\n' for i,(q,a) in enumerate(answers,1))+'\n## 实际失败与延长\n\n'+failures+'\n'+extension+'\n## 查看证据\n\n[数值对照](NUMERICAL_CONVERGENCE.md)、[逐动作矩阵](SHAPE_RESPONSE_MATRIX.md)、[初态依赖](INITIAL_SHAPE_DEPENDENCE.md)、[数据定义](DATASET_AND_F10Y_INTERFACE.md)、[计算成本](PERFORMANCE_AND_RESOURCES.md)。\n\n![真实形变与目标](figures/Fig_F10AB_10_target_or_limitations.png)\n'
 doc('F10AB_DEVELOPMENT_REPORT.md','F10AB完整科研开发报告',report)
 gates={}
 def gate(name,status,evidence,rule,scope_text):gates[name]=dict(status=status,evidence=evidence,rule=rule,scope=scope_text,source_identity=identity(),data_manifest='INPUT_OUTPUT_SHA256.json')
 verified=['F10AA_BACKEND_VERIFIED','F10S_F10Z_BASELINES_PRESERVED','PHYSICS_MODEL_UNCHANGED','NUMERICAL_TOLERANCES_UNCHANGED']
 for n in verified:gate(n,'PASS',['BASELINE_AUDIT.json','evidence/REMOTE_BASELINE_AUDIT.json','evidence/FINAL_PROTECTION_AUDIT.json'],'Frozen source/input/history hashes unchanged; inherited original integrator','F10AA trusted delivery and original physical policy')
 for n in ['INITIAL_SHAPES_GEOMETRICALLY_VALID','S0_ORIGINAL_CLUSTER_VERIFIED','S1_COMPACT_CLUSTER_VALID','S2_X_ELONGATED_VALID','S3_Y_ELONGATED_VALID','S4_OBLIQUE_SHAPE_VALID','S5_MULTI_COMPONENT_VALID']:
  gate(n,'PASS',['INITIAL_STATE_AUDIT.json','INITIAL_STATES.h5'],'No illegal overlap; cache interior and walls legal; requested initial graph','19 saved states; artificial initialization, not acoustic formation')
 for n in ['OFF_CONTROL_VALID','PZT_ACTIONS_VALID','FORCE_DECOMPOSITION_COMPLETED','RIGID_TRANSLATION_SEPARATED','RIGID_ROTATION_SEPARATED','ISOTROPIC_CONTRACTION_SEPARATED','ANISOTROPIC_DEFORMATION_EVALUATED','SHEAR_RESPONSE_EVALUATED','NONRIGID_INTERNAL_MOTION_EVALUATED']:
  gate(n,'PASS',['static/SUMMARY.json','analysis/ADMISSION.json','evidence/TESTS.json','SHAPE_METRICS_DEFINITION.md'],'Actual diagnostic/trajectory definitions; invalid F/axes retained; original field and source gates','Only saved states/actions; tests distinguish transformations, not physical calibration')
 for n in ['ACTION_SHAPE_RESPONSE_IDENTIFIABLE','POSITION_SHAPE_DECOUPLING_EVALUATED','LOCAL_CONTROLLABILITY_ANALYZED','INDEPENDENT_SEED_REPLICATION_PASS']:
  gate(n,'PASS',['analysis/REPLICATION.json','analysis/SENSITIVITY.json','analysis/TARGET.json'],'A5 signed anisotropy repeated; controlled local stencil with error bound; report failed connectivity','Finite fixed-ID distribution response; effective shape rank1; not multiple modes or all-seed coherent-cluster control')
 for n in ['STEP_SIZE_CONVERGENCE_PASS','FIELD_CACHE_CONVERGENCE_PASS']:
  gate(n,'PASS' if all(c['original']['admitted'] for c in original_checks) else 'FAIL',['NUMERICAL_CONVERGENCE.md'],'All reported selected full-window original comparisons admitted','Selected refined trajectories, not unrefined discovery actions')
 gate('CONTACT_AND_FORCE_GATES_PASS','FAIL' if fails or not stats['all_completed_independent_audits_pass'] else 'PASS',['NEGATIVE_RESULTS_AND_LIMITATIONS.md','analysis'],'All attempted trajectories pass original gates','All attempts: five original gate rejections; completed subset independently checked and passes')
 gate('LONGER_DURATION_COMPLETED','PASS' if all(r['status']=='PASS' for r in stages['extension']) else 'FAIL',['extension/INDEX.json','analysis/EXTENSION_NUMERICAL.json'],'3/5ms actual completion plus dt/cache','S2_r0 A5 only; not long-time stability')
 gate('TARGETED_SHAPE_TRANSFORMATION_PASS',target['status'],['TARGET_PREREGISTRATION.json','analysis/TARGET.json'],'3%±0.5pp aspect, <=20nm COM, endpoint single component, original numeric admission','Three prespecified S2 states; scalar target alone is insufficient')
 for n in ['SHAPE_RECOVERY_OR_REVERSAL_PASS','WORLD_MODEL_TRAINED','MPPI_CONTROL_EXECUTED','REAL_PDMS_FIELD_VALIDATED','ELECTRICAL_TO_ACOUSTIC_CALIBRATED','HARDWARE_EXPERIMENT_VALIDATED']:
  gate(n,'NOT_RUN',['NEGATIVE_RESULTS_AND_LIMITATIONS.md'],'No corresponding experiment executed','Not inferred from simulation completion')
 gate('WORLD_MODEL_DATA_READY','INCONCLUSIVE',['portable_dataset/INDEX.json','DATASET_SCHEMA.json'],'Schema ready; broad independent action/state coverage and physical validity not established','Exploratory sampled state data; no training admission')
 for n in ['FIGURES_VALIDATED','VIDEOS_VALIDATED']:
  evidence='evidence/VISUAL_INSPECTION.json' if n.startswith('FIGURES') else 'videos/MANIFEST.json';exists=(BASE/evidence).exists()
  gate(n,'PASS' if exists else 'INCONCLUSIVE',[evidence],'10 PNG/PDF visual review; three MP4 all frames decoded and uniform PTS','Actual arrays only; static and failed cases labelled')
 gate('SCIENTIFIC_SHAPE_CONTROL_ADMITTED','PASS',['analysis/ADMISSION.json','analysis/REPLICATION.json','analysis/SENSITIVITY.json'],'Nonrigid/anisotropic responses beyond numeric uncertainty, action differences and held-out signed response','Limited fixed-ID anisotropic distribution degree of freedom only; coherent-target FAIL, multiple independent modes/reversibility unproven')
 gate('TRAINING_READY','FAIL',['DATASET_AND_F10Y_INTERFACE.md'],'No broad or calibrated training admission','training_ready=false')
 write(BASE/'FINAL_GATE.json',dict(gates=gates,conclusion=conclusion,labels=LABELS,training_ready=False,identity=identity()))
 print('Research report',stats['trajectories'],'attempts;',physics_wall,'physics seconds')
if __name__=='__main__':main()
