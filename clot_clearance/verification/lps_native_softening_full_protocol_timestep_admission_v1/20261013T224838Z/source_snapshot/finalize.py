"""Final reporting and scoped protection. Does not start native simulation."""
from common import *
import subprocess,shutil,datetime
a=json.loads((O/'FULL_PROTOCOL_NATIVE_TIMESTEP_AUDIT.json').read_text());en=json.loads((O/'FULL_PROTOCOL_DYNAMIC_ENERGY_LEDGER.json').read_text());co=json.loads((O/'HISTORICAL_DT_COMPARISON.json').read_text());ct=json.loads((O/'DYNAMIC_CONTACT_CLOSURE_AUDIT.json').read_text());top=json.loads((O/'FINAL_TOPOLOGY_AND_COMPONENTS.json').read_text());doc=json.loads((O/'WSL_DOCUMENT_UPDATE.json').read_text());res=json.loads((O/'RUN_RESOURCES.json').read_text());diag=json.loads((O/'FAILURE_EVENT_DIAGNOSIS.json').read_text());identity=json.loads((O/'SOURCE_AND_LIBRARY_IDENTITY.json').read_text());ev=json.loads((O/'FIRST_SOFTENING_AND_FAILURE_EVENTS.json').read_text());p0=json.loads((O/'ACTIVE_ENTRYPOINT_AUDIT.json').read_text());v=a['final_ledger'];maxima=a['maxima']
assert doc['status']=='PASS_ACTUAL_WSL_FILE_UPDATED' and sha(O/'clot_rupture.after.md')==doc['after_sha256']
protected=json.loads((O/'PROTECTED_FILES_BEFORE.json').read_text())['files'];rootindex=R/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.md';assert all(sha(p)==h for p,h in protected.items()),'Protected asset changed before intentional active-index update'
assert subprocess.check_output(['git','diff','--binary'],cwd=R)==(O/'tracked_diff.before.patch').read_bytes()
index=f'''# 当前clot原生软化开发入口

本轮：{TASK} / {RUN}（2026-10-14 Europe/Zurich）。

- **结论：从无损伤t=0开始，dt=0.25µs在第2562步/640.5µs因Δd=0.11363469>0.10停止。完整700µs FAILED_ADMISSION。**
- 最后通过全部冻结门槛的状态是第2561步/640.25µs；越界现场已保存，不继续积分。
- 已覆盖完整升载与峰值保持，卸载350–640.5µs仅部分覆盖；末保持NOT_RUN。
- 原生历史/力/能量/动量在实际窗口一致，能量相对余额最大0.372426823%；不认证全局时间收敛。
- 232.5µs首次真实闭合接触，1576个步有活动接触；最大软罚力重叠24.6649%dx，物理穿透精度未认证。
- 停止时4834条完全失效、7213条存活；仍是一个341点连通体。新裂面/Gc不可用。
- 先读verification/{TASK}/{RUN}/FINAL_STATUS.json和36项REVIEW报告。
- 新驱动/审计：analysis/full_protocol_timestep_v1/；run_protocol.py有单工况启动锁，禁止对封存目录重跑。
- 冻结工况：同目录verification下FULL_PROTOCOL_TIMESTEP_DESIGN_FREEZE.json及其SHA。
- 核心数据：trajectory.npz、damage_history_sparse.npz、keyframes/、actual_contact_pairs.npz、final_executed.restart。
- 关键失败：NUMERICAL_SAFETY_STOP_RECORD.json、FAILURE_EVENT_DIAGNOSIS.json；键1813，自由节点24/245，独立损伤公式完全吻合。
- 生产环境仍为source /workspace/y_clot/env.sh；本轮只在进程内优先隔离库。
- 复用库：{LIB}
- 库SHA：866065394c48ebc8a6d32f147bde992a325c2a503a49d50392243eee10fea9bd；本轮没有C++修改/编译。
- 实际C++：analysis/dynamic_softening_v1/bond_bpm_peri_variational.h/.cpp及fix_lps_variational_dynamic.h/.cpp。
- 内部style LPS/VARIATIONAL/DYNAMIC，POST_INTEGRATE一次提交，POST_FORCE将N/m³乘V为N。
- 冻结离线数学：analysis/softening_feasibility_v1/model.py，SHA 3eb1396b376cf30e40586c26ba649909ab9382d9eac9607666359258273acf74。
- 341点/12047键，左右49/47规定点、245自由点；rhoV=1.25e-10kg；G940、K100000/3Pa，s0=.15、sc=.18。
- 每侧60µm，半余弦升载300µs、保持50µs、卸载300µs、末保持50µs；冻结物理加载不变。
- 上轮依据：verification/lps_native_softening_dynamic_commit_and_small_crack_v1/20261013T214800Z/；旧0.5µs到122µs与105–122µs短窗口均保持。
- 本轮不重扫31001历史文件；按上轮195文件payload及直接依赖做限定SHA核查，删除0。
- 根活跃入口是本轮明确更新的可变文件；其上轮SHA副本保存在本轮ACTIVE_DEVELOPMENT_ENTRYPOINTS.before.md；旧manifest不改。
- 生产库SHA d1b2f53d6b942d1ca9277a8b8cd9c195ce75ea6f616681ec21c311e4faf3b41a，env.sh和已有tracked diff均保持。
- 六组PNG/SVG：results/{TASK}/{RUN}/OPEN_RESULTS.html。
- 实际WSL文档：/home/lzy/projects/clot_rupture.md，before/after/diff/SHA在本轮verification。
- 下一轮唯一优先：同一工况的局部损伤事件分步与力—功—能量一致积分，先设计/预注册；不是盲目全局dt=0.125µs。
- P4缺口试样NOT_RUN；不自动进行V3、CFD、微泡、AI、材料Gc标定或生产部署。
- 后续只读上述相关入口，无需展开全部历史报告。
'''
assert sha(rootindex)==protected[str(rootindex)];rootindex.write_text(index);(O/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.after.md').write_text(index)
changes=[]
for p,h in protected.items():
 actual=sha(p)
 if actual!=h:
  assert p==str(rootindex);changes.append(dict(path=p,before=h,after=actual,reason='Authorized mutable current entrypoint; prior snapshot preserved'))
assert len(changes)==1
dump(O/'HISTORICAL_PROTECTION_AFTER.json',dict(status='PASS_SCOPED_PROTECTION',verified_file_count=len(protected),unchanged_file_count=len(protected)-1,intentional_mutable_changes=changes,prior_manifest_sha256=sha(OLD/'PAYLOAD_MANIFEST.json'),historical_recursive_31001_scan=False,scope='Prior195-file sealed payload plus direct dependencies; reuse earlier whole-history certification. No claim of a fresh whole-repository rescan.',production_unchanged=True,environment_unchanged=True,tracked_diff_unchanged=True,deleted_files=[]))
p0['current_entrypoint_after_sha256']=sha(rootindex);p0['current_entrypoint_before_sha256']=protected[str(rootindex)];dump(O/'ACTIVE_ENTRYPOINT_AUDIT.json',p0)
items={
 'SOURCE_LIBRARY_IDENTITY':'PASS_REUSED_IDENTICAL_LIBRARY_NO_RECOMPILE',
 'INITIAL_UNDAMAGED_STATE':'PASS_ORIGINAL_T0_NO_RESTART',
 'PHYSICAL_LOADING_IDENTITY':'PASS_IDENTICAL_ANALYTIC700US_PROTOCOL',
 'FULL_DT_025_PROTOCOL':'FAILED_ADMISSION_STOP_AT640_5US',
 'ACTUAL_FREE_NODE_DYNAMICS':'PASS_OBSERVED_NATIVE245_FREE_NODES',
 'ACCEPTED_STEP_DAMAGE_COMMIT':'PASS_ONE_PER_EXECUTED_PHYSICAL_STEP',
 'RUN0_AND_TRIAL_READ_ONLY':'PASS_SHA_VERIFIED_PRIOR_EVIDENCE_PLUS_INITIAL_RUN0',
 'DAMAGE_HISTORY_IRREVERSIBLE':'PASS_OBSERVED_WINDOW_INCLUDING_PARTIAL_UNLOAD',
 'FIRST_SOFTENING':'OBSERVED_85US_BOND3393',
 'FIRST_COMPLETE_FAILURE':'OBSERVED_92_5US_SEVEN_BONDS',
 'MAX_DAMAGE_INCREMENT_GATE':'FAILED_ADMISSION_0_11363469187762688_GT_0_10',
 'NODE_STEP_DISPLACEMENT_GATE':'PASS_OBSERVED_WINDOW',
 'RELATIVE_ENDPOINT_DISPLACEMENT_GATE':'PASS_OBSERVED_WINDOW',
 'MOMENTUM_BALANCE':'PASS_OBSERVED_WINDOW',
 'COMPLETE_VARIATIONAL_ENERGY':'PASS_INDEPENDENT_EVERY_STEP',
 'DAMAGE_DISSIPATION_VALIDITY':'PASS_FIXED_NEW_COORDINATE_DROP_NO_DOUBLE_COUNT',
 'DYNAMIC_ENERGY_BALANCE':'PASS_OBSERVED_MAX0_372426823_PERCENT_NOT_FULL700US',
 'FORCE_SWITCH_TIME_LAYER_ERROR':'QUANTIFIED_DOMINANT_NUMERICAL_ERROR_NOT_MATERIAL_D',
 'UNLOADING_HISTORY_VALIDITY':'PASS_PARTIAL_UNLOAD_ONLY',
 'POSTBREAK_DYNAMIC_CONTACT':'OBSERVED_NUMERICAL_CHECKS_PASS_PHYSICAL_OVERLAP_NOT_CERTIFIED',
 'FULL_700US_COMPLETION':'FAILED_ADMISSION_2562_OF2800_STEPS',
 'GLOBAL_TIME_CONVERGENCE_NOT_CERTIFIED':'NOT_CERTIFIED',
 'NATURAL_CRACK_AREA_UNAVAILABLE':'UNAVAILABLE',
 'MATERIAL_GC_NOT_CERTIFIED':'NOT_CERTIFIED_GC_UNAVAILABLE',
 'MACROSCOPIC_STRENGTH_NOT_CERTIFIED':'NOT_CERTIFIED',
 'P4_NOT_RUN':'NOT_RUN',
 'FULL_V3_NOT_RUN':'NOT_RUN',
 'PRODUCTION_LIBRARY_UNCHANGED':'PASS',
 'HISTORICAL_RESULTS_PROTECTED':'PASS_SCOPED_SHA_CHECK_NO_NEW_FULL_RESCAN',
 'WSL_DOCUMENT_UPDATE':'PASS_ACTUAL_FILE_UPDATED'
};assert len(items)==30
status=dict(task=TASK,run_id=RUN,overall='DELIVERED_WITH_FAILED_FULL_PROTOCOL_ADMISSION',scientific_answer='否。0.25µs从原始t=0开始仍在640.5µs触发Δd门槛，不能安全完成700µs。已观察窗口内不可逆历史、原生力、动量及完整能量账本一致；动态接触真正发生，但较大罚力重叠的物理精度尚未认证。',program_finished_as_registered=True,full_protocol_safe=False,items=items,main_summary={q:a[q] for q in ['executed_steps','last_admitted_step','final_time_s','dt_s','stop_reasons','contact_status']},next_single_priority=diag['next_single_priority'],material_Gc='UNAVAILABLE',new_crack_area='UNAVAILABLE',macroscopic_strength='NOT_CERTIFIED',P4_future_admission='NOT_GRANTED',additional_physical_runs=0,CPP_modified=False,library_recompiled=False,production_deployed=False,WSL_document=doc)
dump(O/'FINAL_STATUS.json',status)
res['postprocess']=json.loads((O/'IMPLEMENTATION_CHECKS.json').read_text());res['final_render']=json.loads((O/'FIGURE_MANIFEST.json').read_text());res['render_attempts']=2
if (O/'FIGURE_RENDER_FIRST_ATTEMPT.json').exists():res['initial_render_before_layout_fix']=json.loads((O/'FIGURE_RENDER_FIRST_ATTEMPT.json').read_text())
res['output_bytes_before_final_report_seal']=sum(p.stat().st_size for p in O.rglob('*') if p.is_file())+sum(p.stat().st_size for p in F.rglob('*') if p.is_file());res['final_payload_size_source']='PAYLOAD_MANIFEST.json total_bytes (does not include itself/receipt)';res['orchestration_total_wall_s']='NOT_MEASURED';dump(O/'RUN_RESOURCES.json',res)
answers=[
('使用的是哪个原生软化共享库？',f'实际加载 `{LIB}`（运行映射为liblammps.so.0），SHA256 **{identity["library_sha256"]}**。源与库现场复核，复用上一轮已验收隔离库；本轮未重新编译。'),
('原始341节点、12047键和材料参数是否完全保持？','保持。左右规定节点49/47、自由节点245；ρV=1.2500000000000003×10⁻¹⁰kg，G=940Pa、K=100000/3Pa、ρ=1000kg/m³、dx=50µm、δ=150.5µm、s₀=.15、s_c=.18、α=0。初始X、m、原始键ID、r、taper、体积均以原文件及SHA核对。'),
('是否确实从原始无损伤状态开始？','是。正式进程只有一个Native实例，以原body.data从step0/time0开始；初始速度、κ、d、D、U、K及力为零，没有read_restart或损伤导入。初始重复run0逐值只读。'),
('完整加载协议是否仍为700µs？','是。完全相同的左右±半余弦解析曲线，每侧幅值60µm，升载300µs、保持50µs、卸载300µs、末保持50µs，y/z仅夹持点固定。所有实际边界坐标与解析式最大差1.0842×10⁻¹⁹m；最大速度差9.9920×10⁻¹⁶m/s。'),
('实际运行多少物理时间步？','2562个原生物理步，dt=0.25µs；共2563个保存状态含初态。第2562步是原生已提交的越界现场，不属于满足全部门槛的接受集；最后合格是第2561步。'),
('是否完成2800步？','没有。完整协议FAILED_ADMISSION，剩余238步未运行；没有事后继续或改dt拼接。'),
('若未完成，首次安全停止发生在什么时间？','第2562步、640.5µs。最后合格状态640.25µs。唯一触发理由MAX_DAMAGE_INCREMENT；NUMERICAL_SAFETY_STOP_RECORD.json保存两端速度、材料/总力、κ/d前后、实际伸长、完整能量以及停止restart。'),
('最大Δd是多少？是否始终小于等于0.10？',f'最大**0.11363469187762688**，并非始终≤0.10；原始键1813（节点24、245）越界。前一步d=0.0002808428166676996，后一步d=0.11391553469429458；κ=0.15000702139906827→0.15290300458732148。门槛从冻结到停止均未改变。'),
('首次自然软化的时间和bond ID是什么？','85µs，第340步，原始键3393；与旧dt=.5µs的实际物理事件时间一致。'),
('首次完全失效的时间和bond ID是什么？','92.5µs，第370步，共7条：2180、2194、3368、3393、3904、9197、10264。这是实际观察到的相同集合，没有强制匹配旧事件。'),
('是否发生任何损伤历史回退？','没有。全实际窗口（包括350–640.5µs部分卸载）的κ/d最大回退为0，与独立模型误差为0。稀疏d变化、逐步FP64坐标重建与完整关键帧逐值一致。'),
('是否重复提交损伤或失效事件？','没有。2562个实际步各提交一次，原生step、commit和token连续一致；4834条唯一原始键各发生一次完全失效。重复force/run0、trial和restart语义复用SHA校验过的旧证据，不虚称本轮重新跑了全面矩阵。'),
('最大实际速度和单步位移是多少？',f'最大速度{maxima["max_speed_m_s"]:.12g}m/s；最大步位移/dx={maxima["step_disp_dx"]:.12g}，即{maxima["step_disp_dx"]*DX*1e6:.9f}µm；相对端点步位移/dx={maxima["relative_endpoint_step_dx"]:.12g}。均通过各自冻结门槛5m/s、.05dx、.10dx。'),
('与前轮0.5µs重叠时程的差别有多大？',f'只对比已存0–122µs的245个相同物理时刻，无插值或旧重跑。最大坐标分量差{co["maximums"]["x_max_component_m"]*1e9:.9f}nm、最大节点向量差{co["maximums"]["x_max_vector_m"]*1e9:.9f}nm、最大速度分量差{co["maximums"]["v_max_component_m_s"]:.12g}m/s、最大d差{co["maximums"]["d_max_abs"]:.12g}。共同时间点完全失效ID集合全部相同；122µs均为520条。新/旧末态D差−3.548721799pJ；新R=−3.601696779pJ、旧R=−7.261265927pJ。新重叠窗口最大Δd=.0516474106784；旧.102554248321。旧105–122µs短restart窗口另行保留，不作为本轮完整新路径。'),
('是否可以认为0.25µs已经达到全程时间收敛？为什么？','不可以。本轮自身未完成700µs且Δd越界；没有更细完整路径参考，更没有连续时间极限研究。重叠窗口改善只说明该窗口时间分辨率敏感性。TIMESTEP ADMISSION ≠ GLOBAL TIME CONVERGENCE。'),
('外功、完整储能、动能和损伤耗能分别是多少？',f'640.5µs停止现场W={v["W_J"]*1e9:.12f}nJ，U={v["U_J"]*1e9:.12f}nJ，K={v["K_J"]*1e9:.12f}nJ，D={v["D_J"]*1e9:.12f}nJ；Uc={v["Uc_J"]*1e9:.12f}nJ，Jc=0。Kfree={v["Kfree_J"]*1e9:.12f}nJ。没有数值阻尼。'),
('全程最大能量余额是多少？',f'实际0–640.5µs最大绝对R为{en["maximum_absolute_residual_J"]*1e12:.12f}pJ，最大相对余额**{en["relative_residual_max"]*100:.9f}%**（第{en["max_relative_residual_step"]}步），低于2%。停止R={v["R_J"]*1e12:.12f}pJ。归一化完全继承：运行至本步的max[10⁻¹⁵J,|W|,|ΔU|+|ΔK|+D+|ΔUc|]，不重置；终尺度6.653193900330409nJ。不能把已观察窗口称完整700µs。'),
('损伤力切换项是否仍是主要误差来源？','是。累计VV缺陷+0.020985434641pJ、固定d空间求积−0.009172867557pJ、力切换+22.775050094316pJ、接触求积−0.014058011909pJ；负和重建R，最大未解释余量2.23245×10⁻²²J。力切换约占末态余额绝对值的100.01%，其余小项有符号抵消；不将该误差计入材料D，不宣称断裂整体二阶。'),
('动量平衡是否通过？',f'已观察窗口通过，最大逐步外部冲量余额{maxima["momentum_residual_Ns"]:.12g}N·s，低于10⁻¹⁸。外冲量包含规定边界ΔP与梯形材料反力。自由积分位置/速度方程最大差分别{maxima["x_equation_m"]:.6g}m、{maxima["v_equation_m_s"]:.6g}m/s。'),
('是否进入升载保持及卸载阶段？','进入。升载至300µs完成、峰值保持至350µs完成，卸载仅到640.5µs。300µs已有4245条失效，350µs为4296条；保持阶段新增51条，部分卸载新增538条。保持结束自由最大速度仍0.471727m/s；实际惯性和能量振荡不等于静态平衡。'),
('卸载是否发生非物理损伤愈合？','在已执行350–640.5µs范围没有κ/d回退，D单调非负。尚未执行末段卸载与末保持，不能认证完整循环；本协议没有再加载段，RELOADING=NOT_IN_FROZEN_PROTOCOL_NOT_RUN。'),
('已失效键是否出现不合理的重新传力？','没有发现。全时程完全失效原键两端标量拉力严格为0，最终原生存活端点对与原始ID推导拓扑一致，m不变。软化后状态重分配和接触压缩力不是失效键重新粘合。'),
('是否真正发生动态闭合接触？',f'发生。第930步、232.5µs，节点247/278、原始已失效键11121首次r<dx。{ct["active_steps"]}个实际步有活动接触，同时最多{ct["max_simultaneous_pairs"]}对，停止时{ct["last_pairs"]}对。独立接触力误差{ct["max_force_error_N"]:.6g}N、接触储能差{ct["max_contact_U_error_J"]:.6g}J；接触功−98.874275315pJ，Uc末态98.860217303pJ。'),
('接触与bond力是否重复？','没有存活原始键与接触pair重叠，已失效端点拉力为0；独立枚举全部实际接触对，逐对压缩方向及原生力/能量一致。最大重叠0.246649195614dx=12.33246µm，不能解释为硬壁非穿透或真实材料接触精度已通过。本轮提前冻结的新增几何保护是r/dx<0.10近重合或接触相对步位移>.10dx即停，未触发；它不是“允许90%穿透”的物理标定。接触罚刚度/重叠精度未标定，单列NOT_CERTIFIED。'),
('最终完全失效键数是多少？','停止时4834条完全失效，7213条存活；不是700µs合格末态数。'),
('试样是否保持连通？','保持。以全部存活原始键构图，停止态为一个341节点连通分量；接触不会把原始键重新连上。完整节点分量ID已保存。'),
('是否形成可辨识新裂面？','未建立可辨识新生连续裂面。点云损伤颜色、4834条失效和个别空间间隙不等于裂面面积；新裂面面积UNAVAILABLE，不用断键数×dx²替代。'),
('目前是否可以报告材料Gc？','不可以。Gc为UNAVAILABLE。G0_THEORY=0.07170070833333J/m²、G0_DISCRETE=0.06969623734234J/m²仍为旧瞬断LPS名义参考；旧运动学173.156623/230.776817nJ也不是本轮自然断裂材料Gc。未扫描参数或拟合面积。'),
('是否已证明宏观峰值强度不变？','没有。未损伤刚度、s₀和材料参数不变，并不意味着宏观峰值强度已测量或与材料实验相符。'),
('本轮是否修改任何原生C++或材料参数？','没有。四个研究C++、锁定Verlet/NVE/move/pair源码以及共享库与旧SHA一致，零编译。新增Python只负责冻结、逐步观测/安全停止、独立核查、保存、可视化和报告，没有用Python推进自由节点。'),
('计算耗时、CPU/RAM和输出规模是多少？',f'唯一主计算加逐步在线核查{res["main_wall_s"]:.6f}s，CPU {res["main_CPU_s"]:.6f}s；含最终压缩{res["wall_with_final_compression_s"]:.6f}s；峰值{res["peak_RSS_bytes"]/1024**2:.3f}MiB。1 MPI/1线程、不用GPU。实际CPU配额15.36核、RAM上限42.078GiB；冻结主预算300s/4GiB/600MiB，未越界。独立后处理和两次图形排版费用另计RUN_RESOURCES.json；精确最终文件数/字节见PAYLOAD_MANIFEST.json，不虚报未计量的总编排时间。'),
('是否生成必要真实科研可视化？','生成六组PNG/SVG与OPEN_RESULTS.html：实际加载/运动、损伤及门槛、六个真实关键点云、完整能量及缺陷、旧新共同时间对照、独立准入状态。图3固定坐标/相机，节点颜色来自全部原始邻键；为可读性仅画r₀≤1.01dx的边，失效虚线明确仅为参考身份痕迹。没有医学外壳、假末态或4K视频。'),
('生产库和全部历史结果是否保持不变？',f'生产库、env.sh、已有tracked diff、前轮正式科研数据和直接依赖保持。按已封存195项payload及直接资产共{len(protected)}个路径限定复核；不重新扫描31001历史文件，不能声称完成了新的全历史扫描。唯一已授权变动为根ACTIVE_DEVELOPMENT_ENTRYPOINTS.md，前版本已备份且旧manifest原样保留。P0只审查{len(p0["candidates"])}个上轮动态入口文件，安全删除0。'),
('WSL clot_rupture.md是否实际更新？',f'实际更新 `{doc["actual_path"]}`；旧全文在移除本次单一插入块后逐字节相同，①–⑨与全部历史数值保留。修改前检查最新SHA和并发状态，写后读回一致。'),
('备份、diff和SHA在哪里？',f'本轮verification：clot_rupture.before.md、clot_rupture.after.md、clot_rupture.diff、WSL_DOCUMENT_UPDATE.json。旧SHA `{doc["before_sha256"]}`，新SHA `{doc["after_sha256"]}`。冻结SHA、源码/库身份、原始输入、停止restart、原生事件和payload manifest均在同目录。'),
('下一轮唯一值得开展的开发任务是什么？','同一冻结341点路径的**局部损伤事件分步与力—功—能量一致时间积分研究**。键1813的原生与独立损伤公式完全相同，局部κ增量.002895983和近s₀处较陡d(κ)直接解释Δd=.11363469；两个端点为自由节点且越界时无直接接触力。因此优先研究事件分辨与积分耦合，尚无已证实本构代码bug；不能凭这次数据直接改材料或换模型。冻结局部平均伸长率得到约.2195µs的单事件分辨率估计只是诊断，并非新dt工况或保证。下一算法必须独立预注册，保留当前方程/损伤耗能定义并继续审查接触物理精度。本轮不实现/试跑下一算法，不自动全局dt=.125µs、不进入P4。')
]
assert len(answers)==36
report=f'''# LPS_NATIVE_SOFTENING_FULL_PROTOCOL_TIMESTEP_ADMISSION_V1

运行：`{RUN}`；报告日期2026-10-14（Europe/Zurich）。服务器：`{O}`。

**结论：否。当前隔离变分软化模型使用从t=0开始的0.25µs固定步长，仍不能在冻结门槛下安全完成700µs完整加载—卸载。第2562步、640.5µs的Δd=0.1136346919>0.10触发立即停止。最后通过全部门槛是640.25µs。**

程序按预注册停止，不代表完整协议通过。已观察区间内原生自由积分、不可逆历史、完整力/储能/耗能及动量一致，能量最大相对余额0.372426823%<2%。动态闭合接触真正发生，但软罚力最大重叠24.6649%dx的物理精度未认证。仍为一个341点连通体，新裂面面积/Gc不可用，不能准入P4缺口试样。

## 设计与实际执行

唯一物理改变为dt由0.5µs减至0.25µs；本构、参考数据、原始键、质量和解析加载不变。运行前冻结文件SHA为 `{sha(O/'FULL_PROTOCOL_TIMESTEP_DESIGN_FREEZE.json')}`。只运行一个原生实例、2562个真实时间步，没有dt/4、第二窗口、restart拼接、材料调参或阻尼。沿用已有源/库和只读证据，零C++重编译。P0限定核查，不做全历史清理。

κ/d在漂移后POST_INTEGRATE提交，先算同新坐标Uold再算Unew，随后更新全失效owned拓扑及接触屏蔽；POST_FORCE乘V，第二次NVE半步使用新力。`run 0`不提交。每步原生计数、token和物理时间均记录。能量始终按 `W=ΔU+ΔK+D+ΔUc−Jc+R`；R保留，−D拓扑项只显示不重复计入，能量尺度不在阶段转换时重置。

## 各阶段实际覆盖

| 阶段 | 实际范围 | 覆盖和现象 |
|---|---|---|
| 升载 | 0–300µs | 完成；85µs软化、92.5µs完全失效；300µs累计4245失效。 |
| 最大位移保持 | 300–350µs | 完成；新增51条失效，保持末自由速度约0.471727m/s，非静态平衡。 |
| 卸载 | 350–640.5µs | 部分执行，新增538条；终步Δd越界，已保存前后完整态。 |
| 卸载末9.5µs及最终保持 | 640.5µs之后 | NOT_RUN。650µs末保持入口未到达。 |
| 再加载 | 无此预注册阶段 | NOT_IN_FROZEN_PROTOCOL / NOT_RUN。 |

## 独立验收状态（30项）

| 项目 | 状态 |
|---|---|
'''+''.join(f'| `{k}` | `{val}` |\n' for k,val in items.items())+'\n## 36项逐项回答\n\n'+''.join(f'### {i}. {q}\n\n{answer}\n\n' for i,(q,answer) in enumerate(answers,1))
report+='''## 数据、复核与科学限制

`trajectory.npz`每步FP64只保存一份x/v/总力/材料力/标量；`damage_history_sparse.npz`按原始ID保存每次d变化时的κ/d，全部κ可由逐步坐标、初态和冻结最大历史公式精确重建，已与关键原生全态逐值核对。原始文件没有舍入/降采样。首次软化、首次失效、接触激活、阶段转换及停止前后关键帧另存完整状态；停止restart不被用作下一物理步的续算。

`actual_contact_pairs.npz`记录全部已观察接触对、真实距离、径向相对速度和逐对力；`FINAL_TOPOLOGY_AND_COMPONENTS.json`及NPZ记录真实存活图。逐步损伤/动量/能量检查在run1之间立即执行；没有先跑到700µs再事后标红。第一越界点在图上完整保留，后续未运行区间只作背景标记。

当前证据限定单MPI的小图架构；没有新MPI迁移、生产规模、几何收敛或接触物理参数认证。没有未解释原生—独立模型分歧，也没有足够证据把较大接触重叠直接断言为材料真行为。不能由首事件一致或残差较小推出严格时间收敛、宏观强度不变或材料Gc。

保护与交付见 `HISTORICAL_PROTECTION_AFTER.json`、`SOURCE_AND_LIBRARY_IDENTITY.json`、`PAYLOAD_MANIFEST.json`、`REPRODUCE.txt`、`WSL_DOCUMENT_UPDATE.json`。六组科研图见 `OPEN_RESULTS.html`。本轮交付后停止，等待审查，不自动部署或运行新的物理工况。
'''
(O/'REVIEW_LPS_NATIVE_SOFTENING_FULL_PROTOCOL_TIMESTEP_ADMISSION_V1.md').write_text(report)
(O/'REPRODUCE.txt').write_text(f'''READ-ONLY EVIDENCE FIRST. This sealed case must not be rerun in place.
1. Read FINAL_STATUS.json and the36-question Chinese report. Protocol admission FAILED at640.5us; only2562 steps executed.
2. Server environment: source /workspace/y_clot/env.sh
   Research-only LD_LIBRARY_PATH must prepend {OLD}/build
   Research-only PYTHONPATH must prepend {OLD}/source/python
   Do not edit env.sh or production library.
3. Exact existing library: {LIB}
   SHA256 {sha(LIB)}. The .so.0 runtime mapping has the same hash.
   Reused complete source: {OLD}/source; current four files plus frozen model in source_snapshot.
   No new compiler invocation. Local prior binary is retained under verification/lps_native_softening_dynamic_commit_and_small_crack_v1/20261013T214800Z/binary/liblammps.so.
4. Executed once: prepare.py -> run_protocol.py (exclusive FORMAL_CASE_STARTED marker) -> postprocess.py -> figures.py.
   prepare froze exact reference hashes, analytic curve, library/driver identity, gates and budget BEFORE case start.
   run_protocol used original body.data (no restart or history import). No time step beyond first violation.
   postprocess is saved-data-only. figures ran twice solely to improve label spacing; no second simulation.
5. Actual output: trajectory.npz (each physical step x,v,Ftotal,Fmaterial,scalars); sparse original-ID damage updates;
   keyframes/ (native FP64 at scientific events); final_executed.restart; step_scalars.csv; actual_contact_pairs.npz.
   Reconstruct kappa_n = max(kappa_(n-1),max((|x_j-x_i|-r0)/r0,0)) and frozen d(kappa) using source_snapshot/model.py.
   Sparse d and all keyframes were verified exactly. Do not replace native raw data with reconstructed rounded tables.
6. Read-only delivery check on either server or local clone:
   python analysis/full_protocol_timestep_v1/verify_delivery.py
   This checks the scoped payload manifest and hyperlinks; it does not run LAMMPS.
7. Future explicitly authorized replay must choose a NEW run_id/output sandbox and preserve this evidence.
   Do not bypass the exclusive started marker or overwrite design freeze. Current authorization was for one case only.
8. WSL document updated by update_wsl_doc.py on the actual local file; before/after/diff/hash retained here.
   Root current-entry index intentionally changed; old version retained. All other scoped protected SHA unchanged.
9. No dt0125, no P4 specimen, no V3/CFD/bubbles/AI, no production deployment, no Git commit/push.
''')
for p in A.glob('*.py'):shutil.copy2(p,O/'source_snapshot'/p.name)
print('FINAL_REPORT_AND30_STATUSES_WRITTEN',O,flush=True)
