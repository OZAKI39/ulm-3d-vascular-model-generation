"""Seal one typed V3 exploratory pilot, preserving the measured outcome."""
from common import *
import ast,difflib,re
from datetime import datetime,timezone
assert not (O/'PAYLOAD_MANIFEST.json').exists(),'Archive already sealed'
def rd(n):return json.loads((O/n).read_text())
st=rd('FINAL_STATUS.json');freeze=rd('V3_PILOT_FREEZE.json');au=rd('NUMERICAL_AUDIT.json');r=rd('RUN_RESOURCES.json');top=rd('V3_FRAGMENT_TOPOLOGY.json');en=rd('V3_ENERGY_LEDGER.json');ai=rd('V3_SOFTENING_PILOT_AI_RECORD.json')
assert sha(O/'V3_PILOT_FREEZE.json')==(O/'V3_PILOT_FREEZE.sha256').read_text().split()[0]=='49a534fc47ebe6827cc647f21855f0a55bd523549ef98ad5db52022a1108b65c'
for name,h in freeze['driver_sha256'].items():assert sha(A/name)==h,name
for name,h in freeze['input_sha256'].items():assert sha(O/'provenance'/name)==h,name
for p,h in rd('PROTECTED_BASELINES.json').items():assert sha(p)==h,p
assert sha(LIB)==freeze['library_sha256']=='da705cdd2358eab259d9f99ab402e743f6a64c829959eb094c148f8492552ae7'
identity=rd('SOURCE_AND_LIBRARY_IDENTITY.json')
for name,h in identity['cpp_after_sha256'].items():assert sha(A/name)==sha(O/'source_identity'/name)==h,name
assert identity['new_library_sha256']==sha(LIB)
prior=json.loads((BASELINE/'PAYLOAD_MANIFEST.json').read_text());assert sha(BASELINE/'PAYLOAD_MANIFEST.json')==(BASELINE/'PAYLOAD_MANIFEST.sha256').read_text().split()[0]
for name,z in prior['files'].items():assert (R/name).stat().st_size==z['bytes'] and sha(R/name)==z['sha256'],name
for p in A.glob('*.py'):ast.parse(p.read_text())
assert rd('PILOT_PROCESS.json')['returncode']==0,'Native launcher did not finish cleanly'
(O/'source_identity/PYTHON_NEW_FILES.diff').write_text(''.join(''.join(difflib.unified_diff([],p.read_text().splitlines(True),fromfile='/dev/null',tofile=str(p.relative_to(R)))) for p in sorted(A.glob('*.py'))))
for name in ['INITIAL_FORCE_PARITY.json','MINIMAL_INTEGRATION_CHECK.json','EXTERNAL_LEDGER_CHECK.json']:assert rd(name)['status']=='PASS',name
for name in ['HISTORICAL_LOAD_INTERFACE_CHECK.json','CONTACT_STREAM_AUDIT.json','V3_BULK_INTERFACE_DAMAGE.json','V3_FRAGMENT_TOPOLOGY.json']:assert rd(name)['status']=='PASS',name
assert rd('REGISTRY_ABI_REPAIR.json')['status']=='PASS_AFFECTED_INTEGRATION_RECHECK'
assert st['formal_cases']==1 and rd('FORMAL_PILOT_LAUNCHED.json')['allowed_cases']==1
assert freeze['budget']['RSS_limit_GiB'] is None and freeze['budget']['wall_limit_s'] is None
assert freeze['load']['case']=='C1_GAUSSIAN_A200' and freeze['load']['budget_Pa']==200
assert top['all_accepted_mass_conserved'] and top['all_material_points_preserved']
assert en['max_relative_R_percent']<=2 and au['maxima']['max_dd']['value']<=.10
assert all(q['status']=='PASS' and q['native_topology_exact'] and q['special_exclusions_exact'] for q in rd('full/TRIAL_ACCEPT_REJECT_AUDIT.json')['rollbacks'])
assert ai['formal_training_data_eligible'] is False and ai['CRITICAL_LOAD']=='UNAVAILABLE'
if st['full_protocol_completed']:
 assert st['status']=='FULL_PROTOCOL_NUMERICAL_PASS' and abs(st['last_accepted_time_s']-END)<1e-15
 assert top['native_final_endpoint_match'] and au['final_native_topology_match']
else:assert au['exception'] is not None and (O/'full/FAILURE_RECORD.json').is_file()
cp=rd('full/LATEST_NATIVE_CHECKPOINT.json')
for field in ['restart','state']:assert sha(cp[field])==cp[field+'_sha256']
doc=rd('WSL_DOCUMENT_UPDATE.json');assert doc['status']=='PASS' and doc['readback_exact'];assert sha(O/'wsl_document/clot_rupture.after.md')==doc['after_sha256']
assert (O/'wsl_document/clot_rupture.after.md').read_bytes().replace((O/'WSL_DOCUMENT_INSERT.md').read_bytes(),b'',1)==(O/'wsl_document/clot_rupture.before.md').read_bytes()
required=['REVIEW_'+TASK.upper()+'.md','FINAL_STATUS.json','V3_SOFTENING_INTEGRATION_PLAN.md','SOURCE_AND_LIBRARY_IDENTITY.json','V3_TYPED_BOND_AND_ANCHOR_AUDIT.json','V3_MATERIAL_CARD.json','INITIAL_FORCE_PARITY.json','V3_NATIVE_PILOT_HISTORY.csv','V3_BULK_INTERFACE_DAMAGE.json','V3_FRAGMENT_TOPOLOGY.json','V3_FRAGMENT_SIZE_DISTRIBUTION.csv','V3_ENERGY_LEDGER.json','V3_BACKGROUND_DRAG_HANDOFF.json','V3_SOFTENING_PILOT_AI_RECORD.json','RUN_RESOURCES.json','REPRODUCE.txt','OPEN_RESULTS.html']
for name in required:assert (O/name).is_file(),name
assert [int(i) for i in re.findall(r'^### (\d+)\.',(O/required[0]).read_text(),re.M)]==list(range(1,29))
assert len(list(F.glob('*.png')))==len(list(F.glob('*.svg')))==6
visual=rd('VISUAL_REVIEW.json');assert visual['status']=='PASS'
for name,digest in visual['png_sha256'].items():assert sha(F/name)==digest,name
assert rd('MANAGED_FINAL_STATUS.json')['statename']=='EXITED'
if (O/'full/rejected_archives').exists():
 archival=rd('REJECTION_ARCHIVAL_AUDIT.json');assert archival['status']=='PASS' and archival['byte_exact_NPZ_reconstruction_verified_before_removal']
 for p in (O/'full/rejected_archives').glob('batch_*.json'):
  q=json.loads(p.read_text());assert sha(p.parent/q['archive'])==q['archive_sha256']
# Explicitly retain reused source dependencies in addition to build-object hashes.
dependencies={}
for name in ['native_y_wall_query.py','native_y_wall_contact.py','background_flow_sampler.py','y_bifurcation_sources.py','y_mural_clot_geometry.py','early_flow_coupling.py','detached_fragment_metrics.py']:
 p=R/'pd_clot/lammps_fragmentation'/name;dependencies[str(p)]=sha(p)
for p in [BASELINE/'PAYLOAD_MANIFEST.json',SMALL/'BUILD_AUDIT.json',OLDV3/'FRAGMENT_ATTACHMENT_AWARE_V2.json',OLDV3/'FRACTURE_LOCATION_AND_BURST_AUDIT.json']:dependencies[str(p)]=sha(p)
dump(O/'REUSED_DEPENDENCIES.json',dict(files_sha256=dependencies,build_objects_in='BUILD_AUDIT.json',CFD_wall_assets_in='REUSED_CFD_ASSETS.json',source_not_modified=True))
assert (O/'REUSED_CFD_ASSETS.json').is_file()
# Preserve old entry bytes and only advance this pointer after all records exist.
entry=R/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.md';before=O/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.before.md';assert sha(entry)==sha(before)==prior['files']['ACTIVE_DEVELOPMENT_ENTRYPOINTS.md']['sha256']
new=f'''# 当前clot原生软化开发入口

当前：{TASK} / {RUN}。

- 正式性质EXPLORATORY_APPLICATION；状态**{st['status']}**，实际到{st['last_accepted_time_s']*1e3:.9g}ms，计划2.33ms。一次正式case，从真实V3初态开始；没有第二材料卡/载荷，没有重跑旧模型。
- 1927材料点+1561固定锚点、157637参考键，本体82270/界面22040/固定支架53327。保留完整m/r0/V/taper及原始atom/bond/type身份；材料质量守恒，锚点不计血栓质量。
- 本体G940Pa/K100000/3Pa/s0.15/sc.18；界面s0.140625/sc.16875，EXPLORATORY_UNCALIBRATED。完整endpoint-star变分力，局部精确消元，typed restart；固定支架无损伤。
- 新隔离库SHA {sha(LIB)}；冻结SHA {sha(O/'V3_PILOT_FREEZE.json')}。生产库/旧库不改。正式前生成工厂sizeof失配已补编译并复验，失败证据保留。
- C1_GAUSSIAN_A200：453历史表面点，200Pa面积预算、系数780.232224778Pa、sigma0.28mm、1.1605ms脉冲，当前位置/原面积切向。未叠加580点含栓CFD，人工Pa不是超声功率。
- 事件分辨VV：{st['accepted_steps']}接受步、{r['rejected_trials']}拒绝、{r['substep_count']}接受细分步；最大Δd{au['maxima']['max_dd']['value']:.9g}，最大归一化R{en['max_relative_R_percent']:.9g}%。没有RSS或墙钟上限；原数值门槛/输出2GiB保留。
- 末态本体/界面失效{st['actual_bulk_failed']}/{st['actual_interface_failed']}；真实脱离{top['detached_material_points']}点/{top['detached_component_count']}分量，质量{top['detached_material_mass_kg']:.12g}kg；最大{st['largest_detached_nodes']}点，≥200µm分量{top['N_detached_ge200']}。1–2点碎屑不能称临床微栓子。
- 同步真实存活图与锚连接；游离后即时复用无栓Y背景COM拖曳，附着仅原数值阻尼，非双向FSI。1.0dx粒子断后接触与实际Y壁罚接触分开；交接只查真实FP64状态，不追加长时输运。
- 单rank×单线程、无GPU；主段墙钟{r['main_wall_s']:.3f}s，峰值RSS{r['peak_RSS_bytes']/1024**3:.6f}GiB。special_bonds初始化顺序优化避免无用三跳表，图/力学不改。
- 原历史217/208与修正199/190只读探索对照；新旧损伤律、界面取值1ppm、dt及MPI不同，不作单一因果改善归因。
- 单案例AI schema已经输出；EXPLORATORY_UNCALIBRATED，CRITICAL_LOAD/GC_EFFECTIVE=UNAVAILABLE，宏观强度未认证，formal_training_data_eligible=false；真实微泡streaming未接入。
- 先读verification/{TASK}/{RUN}/FINAL_STATUS.json、28项REVIEW、V3_ENERGY_LEDGER.json、V3_FRAGMENT_TOPOLOGY.json与V3_PILOT_FREEZE.json。代码analysis/v3_softening_pilot_v1/；仅在新授权目录重现，禁止重跑/覆盖本归档。
- 六组PNG/SVG：results/{TASK}/{RUN}/OPEN_RESULTS.html，真实点云/事件/曲线，没有连续医学外壳。WSL实际clot_rupture.md已更新、备份、diff、SHA并读回。
- 上轮3080点2.5ms已接受基线20261014T223001Z及全部历史保留。旧G0标称值不是新V3实际Gc。
- 本轮计算和交付完成后停止；下一项调查见28项报告第25项，不自动加大载荷、改材料、运行CFD/微泡/正式AI或追加正式case。
'''
assert sha(entry)==sha(before);tmp=entry.with_suffix('.md.tmp');tmp.write_text(new);os.replace(tmp,entry);(O/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.after.md').write_text(new);(O/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.diff').write_text(''.join(difflib.unified_diff(before.read_text().splitlines(True),new.splitlines(True),fromfile='ACTIVE.before',tofile='ACTIVE.after')))
st.update(delivery_status='COMPLETE_ARCHIVED',WSL_document_updated_and_read_back=True,six_png_six_svg=True)
dump(O/'FINAL_STATUS.json',st)
# Keep AI validity as the scientific outcome, independent of delivery completion.
dump(O/'DELIVERY_VERIFICATION.json',dict(status='DELIVERY_VERIFIED',utc=datetime.now(timezone.utc).isoformat(),formal_cases=1,full_protocol_completed=st['full_protocol_completed'],scientific_status=st['status'],previous_manifest_files_verified_before_entry_update=len(prior['files']),only_previous_manifest_path_advanced='ACTIVE_DEVELOPMENT_ENTRYPOINTS.md; previous bytes preserved',frozen_runtime_input_card_load_unchanged=True,production_and_previous_libraries_unchanged=True,all_python_parse=True,six_PNG_six_SVG=True,WSL_actual_readback_verified=True,AI_single_exploratory_only=True,visual_review=visual))
files={}
for folder in [A,O,F]:
 for p in sorted(folder.rglob('*')):
  if p.is_file() and '__pycache__' not in p.parts and p.name not in ['PAYLOAD_MANIFEST.json','PAYLOAD_MANIFEST.sha256']:files[str(p.relative_to(R))]=dict(bytes=p.stat().st_size,sha256=sha(p))
files[str(entry.relative_to(R))]=dict(bytes=entry.stat().st_size,sha256=sha(entry))
assert sum(x['bytes'] for x in files.values())<2*1024**3,'Final retained payload exceeds the unchanged output budget'
dump(O/'PAYLOAD_MANIFEST.json',dict(task=TASK,run_id=RUN,files=files,file_count=len(files),total_bytes=sum(x['bytes'] for x in files.values()),self_and_sidecar_excluded=True,reused_dependencies=dependencies,scientific_status=st['status'],formal_cases=1))
(O/'PAYLOAD_MANIFEST.sha256').write_text(sha(O/'PAYLOAD_MANIFEST.json')+'  PAYLOAD_MANIFEST.json\n');print('SEALED',len(files),sha(O/'PAYLOAD_MANIFEST.json'),flush=True)
