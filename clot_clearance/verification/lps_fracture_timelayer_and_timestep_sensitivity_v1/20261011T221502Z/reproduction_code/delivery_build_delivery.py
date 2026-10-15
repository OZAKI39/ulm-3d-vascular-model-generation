import json,hashlib,shutil,html,re
from pathlib import Path
C=json.loads(Path(__file__).with_name('config.json').read_text());R=Path(C['local']);O=R/'verification'/C['task']/C['run'];F=R/'results'/C['task']/C['run'];CODE=R/'analysis/timelayer_v1'
def load(n):return json.loads((O/n).read_text())
def dump(n,d): (O/n).write_text(json.dumps(d,ensure_ascii=False,indent=2,allow_nan=False))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert load('HISTORICAL_PROTECTION_AFTER.json')['status']=='PASS'
assert load('SOURCE_AND_LIBRARY_IDENTITY.json')['status']=='PASS_ON_DISK_AND_ACTUAL_PROCESS_MAPS'
assert sha('/home/lzy/projects/clot_rupture.md')==load('DOCUMENT_UPDATE_AUDIT.json')['after_sha256']
items=[
('REFERENCE_LIBRARY_IDENTITY','PASS_ACTUAL_PROCESS_MAPS_SHA','两组实际加载最终双修复库，盘上与运行期身份均确认','SOURCE_AND_LIBRARY_IDENTITY.json'),
('SMALL_SPECIMEN_REPRODUCIBILITY','PASS_REFERENCE_IDENTITY_AUTHORIZED_DYNAMIC_ADAPTATION','原节点/键/材料一致；连续动力时程为已授权适配','SMALL_SPECIMEN_IDENTITY.json'),
('PHYSICAL_LOADING_EQUIVALENCE','PASS_SAME_T0_AND_PHYSICAL_TIME_FUNCTION','两组相同15 µm/10 µs半余弦加载，观察50 µs','TIMESTEP_DESIGN_FREEZE.json'),
('PREBREAK_TIMESTEP_RESPONSE','PASS_FROZEN_ENGINEERING_SCREEN','共同未断键时刻匹配，使用运行前冻结尺度','PREBREAK_TIMESTEP_RESPONSE.json'),
('FRACTURE_EVENT_ALIGNMENT','SAME_FIRST_AND_FINAL_BONDS_ONSET_SAMPLING_DIFFERS','首断202相同；5.00/4.75 µs，事件分组3/4','FRACTURE_EVENT_ALIGNMENT.json'),
('LAGGED_THETA_FORCE_DIFFERENCE','IDENTIFIED_NONZERO_NATIVE_MIXED_VERIFIED','拓扑力变和lag力差分离；瞬时lag力差仍存在','THETA_TIMELAYER_FORCE_AUDIT.json'),
('TIMELAYER_IMPULSE_PROXY','REDUCES_WITH_HALF_DT_ON_ACTUAL_PATH','首断自由节点代理比值0.463657；非反事实轨迹差','IMPULSE_AND_MOMENTUM_AUDIT.json'),
('ACTUAL_MOMENTUM_RESPONSE','PASS_NATIVE_VV_BALANCE_NOT_COUNTERFACTUAL','直接半步/整步速度及实际自由DOF动量重建通过','ACTUAL_MOMENTUM_EVENT_WINDOWS.json'),
('TOPOLOGY_ENERGY_RELEASE','SIGNED_IDENTIFIED_NOT_GC','同坐标有符号储能变化，未标定固有Gc','TOPOLOGY_ENERGY_RELEASE_AUDIT.json'),
('ENERGY_LEDGER_CONSISTENCY','PASS_FROZEN_SCREEN_WITH_SIGNED_DECOMPOSITION','余额小于冻结筛查线，同时披露各带符号残差项','TOPOLOGY_ENERGY_RELEASE_AUDIT.json'),
('POSTBREAK_TIMESTEP_SENSITIVITY','SAME_FINAL_TOPOLOGY_EVENT_GROUPING_DIFFERS','末态同五键、单341节点组件；仍有局部瞬时响应差异','POSTBREAK_COMPARISON.json'),
('TIMESTEP_CONVERGENCE_SCOPE','TIMESTEP_SENSITIVITY_NOT_FULLY_RESOLVED','仅两个步长和小试样，不能证明全局断裂收敛','POSTBREAK_COMPARISON.json'),
('FRACTURE_ALGORITHM_CHANGE_REQUIRED','NOT_SUPPORTED_BY_CURRENT_SMALL_SPECIMEN_EVIDENCE','保留EXPLICIT_LAGGED_THETA_SCHEME；未开发新算法','IMPLEMENTATION_CHECKS.json'),
('PRODUCTION_LIBRARY_UNCHANGED','PASS','生产库、指定原生源码及env.sh身份未变','SOURCE_AND_LIBRARY_IDENTITY.json'),
('HISTORICAL_RESULTS_PROTECTED','PASS_9477_REGISTERED_FILES','9477文件、792806399字节SHA零变化，范围见保护清单','HISTORICAL_PROTECTION_AFTER.json'),
('WSL_DOCUMENT_UPDATE','PASS_ACTUAL_WSL_UPDATED','实际WSL文件已更新、备份、diff和SHA已保留','DOCUMENT_UPDATE_AUDIT.json')]
status={'task':C['task'],'run_id':C['run'],'completion':'BOUNDED_AUDIT_COMPLETED_AND_STOPPED','delivery_date':'2026-10-12','timezone':'Europe/Zurich','scheme':'EXPLICIT_LAGGED_THETA_SCHEME','statuses':{k:{'status':v,'scope':note,'evidence':ev} for k,v,note,ev in items},'formal_completed_dynamic_cases':2,'native_instance_attempts':3,'zero_dynamic_step_setup_failures':1,'additional_dt_levels':0,'full_V3_runs':0,'native_CXX_changes':False,'production_changes':False,'Gc_calibrated':False,'next_priority_not_executed':'已有含栓CFD保守节点牵引的禁止断键原生完整态受载接口验证：单位、总力/力矩、约束反力、小变形响应'}
dump('FINAL_STATUS.json',status)
review=Path(__file__).with_name('review.md').read_text();table='| 检查项 | 状态 | 适用范围 |\n|---|---|---|\n'+'\n'.join(f'| `{k}` | `{v}` | {note} |' for k,v,note,ev in items);review=review.replace('{{STATUS_TABLE}}',table);name='REVIEW_LPS_FRACTURE_TIMELAYER_AND_TIMESTEP_SENSITIVITY_V1.md';(O/name).write_text(review)
assert len(re.findall(r'^### \d+\.',review,re.M))==28
perf=load('CAMPAIGN_PERFORMANCE.json');loops={}
for case in ['dt_0p5us','dt_0p25us']:
 logs=(O/case/'log.lammps').read_text();loops[case]=re.findall(r'Loop time of ([\d.eE+\-]+) on .* for (\d+) steps',logs)
resource={'formal_cases':{n:load(n+'/PERFORMANCE.json') for n in ['dt_0p5us','dt_0p25us']},'formal_outer_wall_s':sum(q['wall_s'] for q in perf['records']),'formal_process_tree_peak_RSS_bytes_sampled':max(q['process_tree_peak_RSS_bytes_sampled'] for q in perf['records']),'formal_CPU_user_plus_system_s':sum(load(n+'/PERFORMANCE.json')[k] for n in ['dt_0p5us','dt_0p25us'] for k in ['CPU_user_s','CPU_system_s']),'native_loop_time_and_steps_from_logs':loops,'setup_only_failed_attempt':load('setup_failure_1/CAMPAIGN_PERFORMANCE.json'),'analysis_wall_s':load('ANALYSIS_SUMMARY.json')['analysis_wall_s'],'last_figure_generation':load('FIGURE_PERFORMANCE.json'),'GPU_used':False,'MPI_ranks_per_case':1,'cases_sequential':True,'OMP_BLAS_threads':1,'memory_measurement':'Observed sampled process-tree RSS; no claim of a hard maximum','output_scope':'See FILE_MANIFEST.json: exact per-file bytes/count for this verification and results run, excluding the manifest, checksums and transfer-verification metadata','formal_stop_conditions_triggered':[],'failed_setup_condition':'Locked LAMMPS scalar min parser failure, accepted positive dynamic steps=0; preserved, equivalent expression re-frozen before formal runs','resources_reference':'PREFLIGHT.json','managed_completion':'MANAGED_COMPLETION.json'}
dump('RESOURCE_AND_COST.json',resource)
repro=O/'reproduction_code';repro.mkdir(exist_ok=True)
for p in CODE.iterdir():
 if p.is_file() and p.suffix in ('.py','.json','.md'):shutil.copy2(p,repro/p.name)
for p in ['remote.py','update_document.py','build_delivery.py','review.md']:
 shutil.copy2(Path(__file__).with_name(p),repro/('delivery_'+p))
(O/'REPRODUCE.txt').write_text('''本轮是已完成、已封存的有限审计。以下说明便于复核，不授权新增仿真。

1. 查看入口
  OPEN_RESULTS.html：四组实际科研图、历史静态补充图、中文28问及证据链接。
  FINAL_STATUS.json：各项独立结论和未解决范围。
  FILE_MANIFEST.json / SHA256SUMS：交付文件身份。

2. 原生身份与依赖
  source /workspace/y_clot/env.sh
  cd /workspace/y_clot/project
  Python /workspace/y_clot/envs/clot/bin/python，3.13.11；OpenMPI /workspace/y_clot/deps/openmpi，4.1.6。
  诊断库 verification/lps_native_state_identity_and_fracture_bugfix_v1/20261011T203015Z/build_owner_fix/liblammps.so
  库SHA 00b246ebe4224e648c27b74afb64fc22b9cea4ba6ffb23f7c40d96c37ca22ae9。
  复用上轮observer.so；本轮副本及源码在observer.so、provenance/，SHA见IMPLEMENTATION_CHECKS。
  独立势能依赖pd_clot/lammps_fragmentation/macro_lab_v1.py及analysis/material_macro_audit_v1/energy.so。
  能量源码/调用模块快照在provenance/；二进制身份见SOURCE_AND_LIBRARY_IDENTITY。
  完整求解器树和历史结果保留在原路径；本交付未重复复制巨大构建树。

3. 实际已执行的命令与环境
  dt_0p5us/RUN_COMMAND.json、dt_0p25us/RUN_COMMAND.json 保存真实argv、cwd、LD_LIBRARY_PATH和冻结SHA。
  实际driver.py参数分别为dt_0p5us与dt_0p25us。
  各case保存reference_initialization.in、dynamics.in、body.data、native_states.npz、事件和日志。
  native_states：FP64 x、force_N、v_half、v_full、theta_native、m；每步拓扑为packed mask。
  原始atom ID对应数组下标+1；原始bond ID对应reference.npz的键表行下标+1，未用当前局部列表序号代替。
  已接受状态包括初始步，分别101/201条。初始run 0用于参考初始化；动力积分中不插入run 0。
  不能重跑launch.py/start_managed.py/driver.py覆盖本轮；它们具有重复运行保护。
  setup_failure_1为零动力步失败历史，不是第三个成功算例，不删除或重命名为成功。

4. 只读式数据复核与离线后处理
  在项目根目录运行：sha256sum -c verification/lps_fracture_timelayer_and_timestep_sensitivity_v1/20261011T221502Z/SHA256SUMS
  原分析程序：analysis/timelayer_v1/static_analysis.py、analyze.py、additional_checks.py。
  原绘图程序：analysis/timelayer_v1/figure_static.py、figures.py。
  这些程序不启动LAMMPS动态实例，但会写分析输出；若需要再生成，应复制到新验证目录并把common.py/config.json的输出目录指向新目录。
  不要重新运行preflight.py或freeze.py来替换原设计及身份清单。
  reproduction_code/保留本轮执行和交付代码；原生C++没有改动，所以不存在新原生修复patch。
  导出图为PNG/SVG，不依赖网络资源，不需要视频播放器。

5. 文档与边界
  document_provenance保存实际WSL文件更新前后及diff；SHA见DOCUMENT_UPDATE_AUDIT.json。
  对历史静态试样增加连续时程的授权及数学冻结分别在DYNAMIC_ADAPTATION_AUTHORIZATION和TIMESTEP_DESIGN_FREEZE。
  明确区分正常拓扑释放、lag力差、原生实际动量、夹具功及有符号离散能量项；本轮不测量材料Gc。
  两级dt不足以证明全局收敛；不自动追加dt/4、完整V3或新theta算法。
''')
figs=[('Fig1_native_timelayer_force','图1：真实首断点阵与同坐标力分解','dt=0.50 µs，t=5.00 µs；参考原子4–90、原键202。当前x*固定，theta色标和力箭头尺度共同。'),('Fig2_timestep_comparison','图2：两种时间步的真实事件与响应','相同15 µm/10 µs加载及50 µs观察窗。每步原生状态，断键事件不插值。'),('Fig3_impulse_and_dynamics','图3：冲量代理与实际动力学','Force-difference impulse proxy ≠ measured velocity jump；规定节点90排除于自由DOF动量积分。'),('Fig4_energy_and_topology','图4：独立能量账本与末态拓扑','同五条断键、一个341节点连通体。带符号拓扑储能变化不是材料Gc。'),('Fig1_static_force_decomposition','补充：历史静态五键事件','复用历史坐标操作结果，不赋予不存在的动力学首断时间。')]
links=['FINAL_STATUS.json','TIMESTEP_DESIGN_FREEZE.json','SOURCE_AND_LIBRARY_IDENTITY.json','SMALL_SPECIMEN_IDENTITY.json','FRACTURE_EVENT_ALIGNMENT.json','THETA_TIMELAYER_FORCE_AUDIT.json','IMPULSE_AND_MOMENTUM_AUDIT.json','TOPOLOGY_ENERGY_RELEASE_AUDIT.json','TIMESTEP_COMPARISON.csv','IMPLEMENTATION_CHECKS.json','RESOURCE_AND_COST.json','DOCUMENT_UPDATE_AUDIT.json','FILE_MANIFEST.json','SHA256SUMS','DELIVERY_VERIFICATION.json','REPRODUCE.txt']
for dest in [O,F]:
 figprefix='' if dest==F else '../../../results/'+C['task']+'/'+C['run']+'/'
 evidenceprefix='' if dest==O else '../../../verification/'+C['task']+'/'+C['run']+'/'
 body='<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>LPS 断键时间层与时间步审计</title><style>body{font:17px/1.6 system-ui,sans-serif;max-width:1500px;margin:auto;padding:28px;color:#202c36;background:#fafbfc}h1{font-size:28px}h2{font-size:23px}figure{margin:28px 0;background:white;padding:14px;border:1px solid #d5dfe4}img{width:100%;height:auto}a{color:#14658b}table{border-collapse:collapse;width:100%;font-size:16px}td,th{padding:9px;border:1px solid #d5dfe4;text-align:left}.note{border-left:4px solid #14658b;background:#edf5f9;padding:15px}code{overflow-wrap:anywhere}.links{columns:2}</style><body>'
 body+='<h1>LPS 断键时间层与时间步敏感性审计</h1><p>Run ID：'+C['run']+' · 2026-10-12 · 真实小试样数据</p><p class="note">两次正式动力学适配已完成。dt 减半后首断冲量代理降为 0.46366 倍，最大能量余额残差从 0.75693 pJ 降至 0.30077 pJ；末态同五条键断裂、仍一个341节点组件。保留当前显式方案。两个步长不足以证明全局收敛，不能解释完整V3细碎屑成因。</p>'
 body+='<p><a href="'+evidenceprefix+name+'">完整中文28问报告</a> · <a href="'+evidenceprefix+'FINAL_STATUS.json">16项独立状态</a> · <a href="'+evidenceprefix+'REPRODUCE.txt">复现与依赖说明</a></p>'
 body+='<table><tr><th>指标</th><th>dt=0.50 µs</th><th>dt=0.25 µs</th></tr><tr><td>首断时刻</td><td>5.00 µs</td><td>4.75 µs</td></tr><tr><td>全窗自由节点累计冲量代理最大模</td><td>4.16178×10⁻¹³ N·s</td><td>2.10894×10⁻¹³ N·s</td></tr><tr><td>自由节点峰值动能</td><td>86.78255 pJ</td><td>86.74921 pJ</td></tr><tr><td>最终断键数 / 组件数</td><td>5 / 1</td><td>5 / 1</td></tr></table>'
 for stem,title,caption in figs:
  body+='<figure><h2>'+title+'</h2><a href="'+figprefix+stem+'.png"><img src="'+figprefix+stem+'.png" alt="'+html.escape(title)+'"></a><figcaption>'+caption+' <a href="'+figprefix+stem+'.png">PNG原图</a> · <a href="'+figprefix+stem+'.svg">SVG矢量图</a></figcaption></figure>'
 body+='<h2>审计与交付文件</h2><ul class="links">'+''.join('<li><a href="'+evidenceprefix+v+'">'+html.escape(v)+'</a></li>' for v in links)+'</ul><p>保留一次零动力步初始化失败记录；成功正式算例2次，原生实例尝试共3次。WSL科研文档已按SHA并发保护更新。计算已停止，本页为静态审查入口。</p></body></html>'
 (dest/'OPEN_RESULTS.html').write_text(body)
dump('FIGURE_REVIEW.json',{'review':'VISUALLY_INSPECTED_FOUR_GROUPS_AND_STATIC_SUPPLEMENT','latest_layout_adjustments':['Force-arrow scale annotation clear of coordinate labels','Five original bond IDs use separate categorical rows','Local native velocity near events uses matching y range','Required impulse-proxy warning uses ≠'],'data_modified_for_appearance':False,'plots':'Matplotlib real native points, true surviving/deleted bonds, native accepted states and independent same-coordinate diagnostics','no_AI_raster_generation':True})
print(json.dumps({'report_questions':28,'statuses':len(items),'report_bytes':(O/name).stat().st_size,'resource_wall_s':resource['formal_outer_wall_s'],'CPU_s':resource['formal_CPU_user_plus_system_s'],'native_loops':loops},ensure_ascii=False))
