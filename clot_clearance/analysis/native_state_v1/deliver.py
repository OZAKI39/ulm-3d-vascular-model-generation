"""Build Chinese review, static index, versioned statuses, and guarded actual WSL doc update."""
from pathlib import Path
import json,hashlib,difflib,html,re
C=json.loads(Path(__file__).with_name('config.json').read_text());R=Path.cwd();O=R/'verification'/C['task']/C['run'];F=R/'results'/C['task']/C['run'];P=O/'provenance';load=lambda name:json.loads((O/name).read_text());sha=lambda b:hashlib.sha256(b).hexdigest()
checks=load('IMPLEMENTATION_CHECKS.json');assert not checks['protected_changed'];gates=load('SMALL_TEST_GATES.json');assert gates['identity_pass'] and gates['independent_force_gradient_pass']
status={'WVOLUME_IDENTITY':'CONFIRMED_BUG_FIXED_MINIMAL_TEST_PASS','ATOM_SORT_SAFETY':'PASS_FINAL_REAL_NONIDENTITY_SORT','MPI_MIGRATION_SAFETY':'PASS_FINAL_49_ACTUAL_MIGRATIONS','RESTART_SAFETY':'PASS_FINAL_NEW_FORMAT_1_AND_2_RANK','LPS_FORCE_ENERGY_CONSISTENCY':'PASS_INDEPENDENT_GRADIENT_NATIVE_PE_TALLY_LIMITATION','FRACTURE_TIME_LAYER':'EXPLICIT_SCHEME_LIMITATION','POSTPROCESS_FRAGMENT_SEMANTICS':'PASS_ATTACHMENT_AWARE_V2_LEGACY_RETAINED','NATIVE_BUGFIX_VALIDATION':'PASS_SMALL_TESTS_FULL_V3_NOT_VALIDATED','V3_REPLAY_STATUS':'FAILED_PARTIAL_SINGLE_ATTEMPT_LAST_SAVED_STEP_320','MATERIAL_GC':'NOT_IDENTIFIED','PRODUCTION_LIBRARY_UNCHANGED':'PASS','HISTORICAL_DATA_PROTECTED':'PASS_2283_FILES_SHA256','MPI_FRACTURE_OWNERSHIP':'CONFIRMED_ORIGINAL_CRASH_FIXED_MINIMAL_TEST_PASS','DOCUMENT_UPDATE':'PENDING_GUARDED_WSL_UPDATE','NEXT_PRIORITY':'Next authorized round: exactly one unchanged C1 replay with final m+owned-endpoint guard build, full streamed ID/order/MPI/force logging; no material tuning','AUTOMATIC_NEXT_TASK':'STOPPED'}
# Read ALL latest actual WSL text and back it up before any modification.
doc=Path('/home/lzy/projects/clot_rupture.md');before=doc.read_bytes();before_sha=sha(before);text=before.decode();(P/'clot_rupture.before.md').write_bytes(before)
assert 'LPS_NATIVE_STATE_IDENTITY_AND_FRACTURE_BUGFIX_V1' not in text
addition='''**本轮原生求解器身份与断键审查（20261011T203015Z）**：已在实际生产库的341点小试样中复现参考加权体积 `wvolume` 的身份错配。真实排序后153个点错误；两MPI中49个点换rank，85个点的m错误。相应最大节点力误差约为独立参考力尺度的60.27%和63.40%。断后扩容、重新run及restart还会按剩余键重算m。这些是状态管理问题，不能通过改材料阈值补偿。

独立诊断版本把初始m及初始化标记交给LAMMPS正式的每原子属性存储，随ID排序、迁移、ghost通信、扩容和restart保存；断键后仍保留初始参考m。最终小试样m/theta/力约10⁻¹⁵相对误差，独立C++储能梯度与节点力通过冻结门槛，缺陷触发前的弹性状态与原版逐位一致。**生产源码和共享库未覆盖。** 原生PE×V的既有偏能计量限制仍存在，不用它计算Gc。

受控5键事件中，断键当步与同坐标断后重算的力差最大4.411 µN（31.85%参考力尺度）。原生力精确符合“断前theta、断后存活键”的显式顺序；下一次力计算恢复固定拓扑一致性。暂记为显式时间层限制，没有改变theta顺序或失效判据。估算ΔFΔt约2.21×10⁻¹² N·s，不是实测动力学速度跳变；正常拓扑储能释放另算。

唯一一次V3对照在后续跨MPI断键路径段错误停止：已保存到第320步的坐标、速度、拓扑与旧基线相同，43个完整力审计时刻到第327步也相同，首次bulk断键同为第241步。未到首次脱离，不能报告末态碎片改善。补充的小试样在**原生产库**也复现同一崩溃：`process_broken`把ghost索引当作本地键数组索引。最终独立版本增加“只修改本进程拥有端点”的两处边界检查，4条跨MPI场景断键、独立力及两rankrestart通过。**最终版本未再跑V3，完整动力学认证尚未完成。** G、K、s00、接触、200 Pa载荷和阻力参数都未改；材料Gc仍未识别。

'''
marker='# ⑥ 外加剪切：';assert marker in text;text=text.replace(marker,addition+marker,1)
stat='''**本轮统计语义修正**：上述历史217点、208个游离分量及真实轨迹完全不变。旧 `largest_free_fraction=0.88739` 实际是“全部分量中的最大值”，对应1710点附着主体；真正最大游离块是4点，质量占比0.2076%。旧≥3点分量数2包含附着主体，真正游离≥3点分量数为1。213个游离点属于1–2点未解析碎屑，不能称为已解析微栓。新增ATTACHMENT_AWARE_V2字段明确区分附着/游离，234个历史保存帧逐拓扑核对；旧字段数值、CSV/NPZ和AI标签保留。≥200/300 µm游离块仍为零。

'''
text=text.replace('# 碎块统计和本次结果\n\n','# 碎块统计和本次结果\n\n'+stat,1)
review=O/'REVIEW_LPS_NATIVE_STATE_IDENTITY_AND_FRACTURE_BUGFIX_V1.md'
text+='\n\n本轮 [原生身份与断键缺陷审查（含25项回答）]('+str(review)+') 和 [四组实际图及数据]('+str(F/'OPEN_RESULTS.html')+') 保存源码补丁、失败运行、独立验证与完整精度状态。此前静态加载接口建议暂缓；**下一轮唯一优先项**是使用最终“参考m身份＋owned端点检查”诊断库，对原C1工况做一次带逐步ID/迁移/力日志的完整V3验证。本轮不再追加运行，不自动调材料或启动CFD/微泡/AI任务。\n'
after=text.encode();assert doc.read_bytes()==before,'Concurrent WSL document edit detected';doc.write_bytes(after);assert doc.read_bytes()==after
(P/'clot_rupture.after.md').write_bytes(after);(P/'clot_rupture.diff').write_text(''.join(difflib.unified_diff(before.decode().splitlines(True),text.splitlines(True),str(doc)+' (before)',str(doc)+' (after)')))
update={'status':'PASS_ACTUAL_WSL_UPDATED','actual_path':str(doc),'before_sha256':before_sha,'after_sha256':sha(after),'concurrency_check':'SHA immediately before write matched latest full read','readback_equal':True,'preserved_numbered_titles':[x for x in re.findall(r'^# [①-⑨].*$',text,re.M)],'backup':str(P/'clot_rupture.before.md'),'diff':str(P/'clot_rupture.diff')};(O/'DOCUMENT_UPDATE_AUDIT.json').write_text(json.dumps(update,indent=2,ensure_ascii=False));status['DOCUMENT_UPDATE']=update['status'];(O/'FINAL_STATUS.json').write_text(json.dumps(status,indent=2,ensure_ascii=False))
resources=load('SUPERVISOR_JOB_RECORD.json');v3=resources['clot_native_state_replay_v1'];smallwall=sum(x['wall_s'] for x in checks['small_test_performance'].values());compilewall=sum(x['wall_s'] for x in checks['compile_resources']);fd=max(x['relative_error'] for r in load('NATIVE_FORCE_ENERGY_AUDIT.json')['records'] for x in r['FD']);fin=load('FRAGMENT_STATISTICS_SEMANTIC_AUDIT.json')
report=f'''# LPS_NATIVE_STATE_IDENTITY_AND_FRACTURE_BUGFIX_V1 — 中文科研审查

本轮确认并在**独立诊断构建**中修复两项实际原生缺陷：参考m身份/持久化，以及跨MPI断键访问ghost的本地键数组。最终小试样通过。**唯一V3对照以中间m修复版本运行，发生原生段错误；最终双修复版未追加V3。不能宣称完整生产验证通过，也没有大碎块改善证据。**

运行标识：{C['run']}。WSL工作树：`{R}`；服务器项目平铺目录：`/workspace/y_clot/project`。基于GitHub32e95bf快照；原生真实迁移提交519a5a，bond源码逐字节吻合锁定上游8de817dd，pair既有补丁一致。原生产库SHA仍为d1b2f53d6b942d1ca9277a8b8cd9c195ce75ea6f616681ec21c311e4faf3b41a。

## 25项明确回答

1. **wvolume是否错配？** 是。341点、12047初始键的同一三维非均匀邻域试样，参考m独立按初始ID/r0/taper/V求和。原生T0误差约1e−15；真实排序后153节点m错误，两rank迁移后85节点m错误，远超预注册2e−12门槛。
2. **历史V3默认排序是否发生？** 输入只设map array，没有关闭默认1000步排序。Verlet setup会调用Atom::sort；历史第1045步以后仍有明确ownership变化，结合Verlet源码支持后续排序调用。历史未保存本地顺序，**无法直接证明某一次具体非恒等排列或精确首次排序步**；本轮小试样保存了实际非恒等排列。
3. **历史是否有MPI迁移？** 是。壁回调按真实owned ID计数，rank0/rank1从1120/807变为1105/822，首次计数变化在第455步。小试样明确记录49个原子换rank及ID清单。
4. **m排序/迁移/restart是否正确？** 原版失败，断后扩容/重新run/restart亦有6节点m改变。最终版本真实排序、49点迁移、ghost、扩容、重新run及1/2rank新格式restart均保持初始m，最大相对误差约1.10e−15。旧损伤restart缺m信息不能自动补全。
5. **修了哪些C++？** bond_bpm_peri.cpp/.h：专用property/atom参考m和ready，跟随copy/grow/exchange/border/restart，刷新别名，不再销毁初始m。bond_bpm.cpp：process_broken两处端点访问限制为0≤index<nlocal。只在新source_owner_fix/build_owner_fix中实施。补丁没有材料、阈值、contact或theta公式修改。
6. **m/theta/F是否符合独立参考？** 是，最终无当步断键的静态、排序、迁移与断后固定拓扑状态均约1e−15。复用既有energy.cpp独立评估，所有方向差分最大相对误差{fd:.3g}，低于2e−6。
7. **断前弹性是否被改变？** 缺陷尚未触发时，三个配对T0的坐标、m、theta、力、存活键全部逐位一致。修复应改变身份已错乱后的结果，这不属于材料参数变更。
8. **断键同一步theta不一致？** 有。5条实际断键，最大Δtheta=0.06208、ΔF=4.41074 µN，约31.8526%独立节点力尺度。同坐标下一次力计算消除该差异。
9. **是明确bug还是时间层限制？** 测得力与源码明确的“断前theta＋断后存活键”预测一致至2.21e−16，判为EXPLICIT_SCHEME_LIMITATION。差异不小，不能宣称时间步收敛或忽略冲量；也不能把明确显式顺序直接改成另一算法。
10. **是否修theta顺序？** 没有。保留原失效判断和force pass。拓扑势能下降约5.9133e−11 J单独记录；ΔFΔt最大2.2054e−12 N·s，仅为估算，不是无积分测试中的实测速度跳变。
11. **参数和判据是否不变？** 是。G=940 Pa、K=33333.333333 Pa、bulk s00=.15、alpha0、dx50 µm、delta150.5 µm不变；V3原interface阈值、200 Pa Gaussian预算、0.5 µs、1.0dx接触、壁参数和拖曳完全保留。
12. **哪些字段语义有误？** largest_free_fraction/FI_free为全部分量最大值/其补量，名称不能解释成游离。debris_fraction及resolved_fragment_count是全部≤2/≥3节点分量，可能包含附着点；保留其原定义。新增8个明确字段并标记ATTACHMENT_AWARE_V2。旧classifier函数AST不变，新增failure_mode_detached_v2独立保存。
13. **真实游离质量和数量变化？** 没变。末态217点、2.7125e−8 kg、208个游离分量；201单点、6双点、1个4点。最大真实游离质量占比0.00207577；旧0.887390指1710点附着主体。旧≥3计数2→真正游离≥3计数1。213游离点为未解析碎屑。≥200/300 µm始终为零。
14. **统计修复是否改拓扑？** 没有。234保存帧重建与原component_id/attached逐项吻合，旧4项统计逐帧数值一致。5个人工图测试＋200/300µm边界与非等质量公式检查通过。独立新分类与旧泛化分类在{fin['new_vs_legacy_label_different_frames']}帧标签不同，显式标记POSTPROCESS_CLASSIFICATION_SEMANTICS_CHANGED；不是重写历史AI标签。
15. **跑了V3吗？** 跑了唯一一次正式尝试，中间m修复版，原C1完整输入字节一致、2MPI；发生SIGSEGV退出139，未完成4660步。最终owned guard加入后不再追加V3。
16. **V3比较结果？** 17完整帧至第320步（0.160 ms）坐标/速度/拓扑完全相同；第320步均30断键、vmax0.665175 m/s、KE约1.21625e−9 J、脱离0点。43完整力审计时刻至327步也完全相同。首次bulk均241；新轨迹未到旧首次脱离593或interface621。崩溃确切步未落盘，不能伪造。没有末态尺寸/锚反力改善结论。
17. **能说旧碎屑主要由bug造成吗？** 不能。P0确有风险，但完整修复版的V3因本轮预算未运行；部分对照发生在历史首次可证实迁移455步前且相同。旧217点/208分量负面结果保留。
18. **仍有哪些物理/数值限制？** 材料Gc未识别，瞬断不是耗散标定；显式theta lag未作时间步收敛；原生PE×V偏能计量不符合独立储能，已有文档限制保留；人工Gaussian非泡声场反演；背景无栓CFD、单向COM阻力、Gaussian与拖曳可能重叠为独立载荷问题。没有修改新表面、黏弹/疲劳/孔弹或微泡机制。
19. **必要原生与独立验证通过吗？** 最终小试样通过，包含额外证据驱动的原版两MPI断键崩溃复现和最终4键跨rank断裂/restart。完整V3最终版本未验证，不能总PASS。粒子pair源码字节未改；壁罚力离线公式检查通过。小试样pair zero不直接验证新断键接触，部分V3未有接触或拖曳事件，应标NOT_EXERCISED。
20. **计算数量/资源？** 12个小诊断场景执行（10成功，首次排序harness配置错误和原版跨MPI断裂2次失败），包含4次必要restart，共16次小型LAMMPS实例初始化；另1次V3，0次旧基线重跑。10成功小场景墙钟合计{smallwall:.3f}s；两次独立编译合计{compilewall:.3f}s，各8核，最大单子进程RSS约541144 KiB。V3 Supervisor记录墙钟{v3['stop']-v3['start']}s，崩溃导致CPU/峰值RSS未留存，标NOT_RECORDED，非0。CPU cgroup15.36核，内存cgroup45176848384字节；无需GPU。输出分项规模见IMPLEMENTATION_CHECKS，原生两套src/build留服务器，约260 MB，失败V3约40 MB。
21. **原始数据保护？** 2283个既有结果、生产源码与库/可执行文件SHA复查0变化。保留失败日志、原版诊断及中间构建。只按授权改后处理代码，旧分类函数数值保持。未push、未删历史数据、未停其他科研任务。
22. **科研图？** 四组PNG/SVG：实际ID m审查、断键时间层DIAGNOSTIC_ONLY、1927点历史真实附着/游离图、唯一V3部分失败对照。第四组明确失败边界，不画假定末态或改善。
23. **WSL文档是否真实更新？** 是，实际`/home/lzy/projects/clot_rupture.md`。完整读取、before备份、写前SHA并发校验、after备份、diff和读回均完成；保留①–⑨、旧表格、参数和Gc未识别结论。审计见DOCUMENT_UPDATE_AUDIT.json。
24. **证据在哪里？** 本目录是完整审查入口；BUGFIX_SOURCE_DIFF.patch为最终3文件补丁，BUGFIX_WVOLUME_ONLY.patch为中间修复；provenance/native_source、python_before及clot_rupture.before.md为备份；所有小状态reference/T0/事件/断后/增长/restart均全精度NPZ。SOURCE_AND_BINARY_IDENTITY、FINAL_PATCH_IDENTITY、复现说明和SHA清单随附。最终服务器构建在本run的build_owner_fix，生产入口env.sh未切换。
25. **下一轮唯一优先项？** 下一次授权后，用最终m＋owned guard诊断库，对原C1做**一次**完整V3回归，逐步落盘ID/排序/迁移/断键/力证据，首先确认跨rank稳定与首次有意义差异；不调材料追求大块。本轮停止。

## 分项状态

| 项目 | 状态 |
|---|---|
'''
report+=''.join(f'| `{k}` | `{v}` |\n' for k,v in status.items() if k not in ['NEXT_PRIORITY','AUTOMATIC_NEXT_TASK'])
report+='''
## 误差、使用和可移植性边界

m单位m⁴（本核omega=1/r），原生力乘V成为N。最大绝对误差、显式参考尺度、出错ID及邻居数、所有rank顺序/归属、ghost值在WVOLUME_IDENTITY_TEST/MPI_SORT_MIGRATION_AUDIT及NPZ。初始零力不使用虚构分母。参考几何/X/r0从独立初始模型计算，绝不以错误原生数组自身作真值。

最终修复只在当前等体积LPS、固定原子集合、锁定MPI库和源版本验证；未声称覆盖所有BPM派生模型、动态增删原子和跨MPI数restart。旧损伤restart丢失参考m的迁移须另有初始键资料。原生PE的单独差异早于本轮，已有clot_rupture说明；本轮force-energy PASS指独立变分势能，不把原生PE表述为PASS。

本轮两次额外定位是由真实错误触发的有界回归，未扫描材料或dt。所有停止/失败保留。原始生产代码和历史轨迹仍是原版本；最终诊断库尚未部署。

[实际源码路径审查](WVOLUME_SOURCE_AUDIT.md) · [各项状态](FINAL_STATUS.json) · [原生力与独立势能](NATIVE_FORCE_ENERGY_AUDIT.json) · [单次V3失败对照](V3_REPLAY_COMPARISON.json) · [统计语义](FRAGMENT_STATISTICS_SEMANTIC_AUDIT.json) · [复现](REPRODUCE.txt)
'''
review.write_text(report)
(O/'REPRODUCE.txt').write_text('''Use the locked server environment: source /workspace/y_clot/env.sh; cd /workspace/y_clot/project
Never run make_patch/build scripts against accepted output paths. For a fresh reproduction clone analysis/native_state_v1 into a new task directory and replace config.run plus the explicit frozen run paths in build/replay shell wrappers. Retain material/input constants and source hashes. Original production source must match SOURCE_AND_BINARY_IDENTITY.
1. capture.py records environment/SHA and protected baselines; DESIGN_FREEZE copied before tests.
2. Compile probe.cpp with mpicxx -std=c++17 -O2 -shared -fPIC -DLAMMPS_SMALLBIG -DOMPI_SKIP_MPICXX, exact source/src + third_party includes, and -L<production build> -llammps. Output probe_original.so in the new run only.
3. python -m analysis.native_state_v1.diagnose sort original
   mpiexec -np 2 python -m analysis.native_state_v1.diagnose migration original
   python -m analysis.native_state_v1.diagnose fracture original
4. Original cross-rank fracture reproduction: mpiexec -np 2 python -m analysis.native_state_v1.diagnose mpi_fracture original (known SIGSEGV; isolated process, writes only new run).
5. Final source changes: BUGFIX_SOURCE_DIFF.patch against an independent copy of the exact recorded production source. cmake Release BUILD_MPI=ON BUILD_SHARED_LIBS=ON PKG_BPM=ON PKG_MOLECULE=ON DOWNLOAD_POTENTIALS=OFF MPI_CXX_COMPILER=/workspace/y_clot/deps/openmpi/bin/mpicxx. Build --parallel 8 only after rechecking quota/jobs. Do not overwrite production or intermediate lib. Compile probe_final.so using final headers.
6. Set LD_LIBRARY_PATH=<new final build>:$LD_LIBRARY_PATH. Run sort final, 2-rank migration final, fracture final, 2-rank mpi_fracture final. audits.py compares full precision states against the locked existing independent energy.cpp/energy.so.
7. postprocess.py only reads historical accepted V3 and writes a new POSTPROCESS_OLD_NEW_COMPARISON. figures.py uses saved NPZ/JSON and generates static plots. final_checks.py rehashes protected records.
8. DO NOT rerun replay.py this round. One attempt already consumed and failed; final build has no full V3 trajectory. Next round needs a fresh independent case and streamed observer logs. The archived replay.py preserves exactly the attempted input and buffered observer which lost buffered state at native crash.
9. deliver.py is a one-time guarded actual WSL document update, not a generic rerunnable scientific test. The existing backup/diff/SHA are the documentary evidence; do not reapply duplicate text.
Dependencies: server private Python NumPy/SciPy/mpi4py/LAMMPS/VTK and existing C++ Energy library; plotting matplotlib. No AI-generated scientific imagery or coordinate editing.
''')
links=['FINAL_STATUS.json','WVOLUME_IDENTITY_TEST.json','MPI_SORT_MIGRATION_AUDIT.json','RESTART_IDENTITY_AUDIT.json','FRACTURE_TIME_LAYER_AUDIT.json','NATIVE_FORCE_ENERGY_AUDIT.json','FRAGMENT_STATISTICS_SEMANTIC_AUDIT.json','POSTPROCESS_OLD_NEW_COMPARISON.csv','V3_REPLAY_COMPARISON.json','IMPLEMENTATION_CHECKS.json','BUGFIX_SOURCE_DIFF.patch','SOURCE_AND_BINARY_IDENTITY.json','REPRODUCE.txt','DOCUMENT_UPDATE_AUDIT.json']
rel='../../../verification/'+C['task']+'/'+C['run']+'/'
page='<!doctype html><html lang="zh"><meta charset="utf-8"><title>LPS 原生缺陷审查</title><style>body{font:17px system-ui;max-width:1500px;margin:30px auto;padding:20px;color:#263747}img{width:100%;border:1px solid #bbb;margin-bottom:28px}a{color:#176293}li{margin:9px 0}.warning{background:#fff0da;padding:18px}</style><h1>LPS 原生身份与断键审查</h1><p class="warning">两项真实缺陷已在独立构建中修复，小试样通过；唯一 V3 运行失败。最终版本尚未完成 V3 回归，不能宣称大碎块改善。生产库及 2283 项历史文件未改变。</p><p><a href="'+rel+review.name+'">中文25项完整报告</a> · <a href="'+rel+'WVOLUME_SOURCE_AUDIT.md">源码行号审查</a></p><ul>'+''.join('<li><a href="'+rel+n+'">'+n+'</a></li>' for n in links)+'</ul>'
for name in ['Fig1_wvolume_identity','Fig2_fracture_time_layer','Fig3_fragment_semantics','Fig4_V3_partial_failure']:page+='<p><a href="'+name+'.svg">'+name+'（SVG）</a></p><img src="'+name+'.png" alt="'+name+'">'
(F/'OPEN_RESULTS.html').write_text(page+'</html>');(O/'OPEN_RESULTS.html').write_text('<!doctype html><html lang="zh"><meta charset="utf-8"><a href="../../../results/'+C['task']+'/'+C['run']+'/OPEN_RESULTS.html">打开四组图和中文报告</a></html>')
print('Delivered review and actual WSL update',update['before_sha256'],update['after_sha256'])
