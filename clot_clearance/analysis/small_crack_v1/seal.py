"""Seal negative scientific result without turning delivery success into physics success."""
from common import *
import difflib,ast,gzip
from datetime import datetime,timezone
assert not (O/'PAYLOAD_MANIFEST.json').exists(),'Already sealed; do not alter archived output'
status=json.loads((O/'FINAL_STATUS.json').read_text());assert status['resource_admission']=='FAIL' and not status['natural_crack_propagation_observed']
freeze=json.loads((O/'SMALL_CRACK_DESIGN_FREEZE.json').read_text())
assert sha(O/'SMALL_CRACK_DESIGN_FREEZE.json')=='f5c352f3466226c5a14855c88fa6d8b2fdead788960b1bf36c03a18125f03a75'
for n,h in freeze['driver_sha256'].items():assert sha(A/n)==h,n
for n,h in freeze['reference_sha256'].items():assert sha(O/'provenance'/n)==h,n
for p,h in json.loads((O/'PROTECTED_BASELINES.json').read_text()).items():assert sha(p)==h,p
prior=json.loads((PREV/'PAYLOAD_MANIFEST.json').read_text())
for n,m in prior['files'].items():assert sha(R/n)==m['sha256'],('previous artifact modified',n)
before=O/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.before.md';entry=R/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.md';assert sha(entry)==sha(before)
doc=json.loads((O/'WSL_DOCUMENT_UPDATE.json').read_text());assert doc['readback_exact'] and doc['status']=='PASS'
assert sha(O/'wsl_document/clot_rupture.after.md')==doc['after_sha256']
assert (O/'wsl_document/clot_rupture.after.md').read_bytes().replace((O/'WSL_DOCUMENT_INSERT.md').read_bytes(),b'',1)==(O/'wsl_document/clot_rupture.before.md').read_bytes()
for p in A.glob('*.py'):ast.parse(p.read_text())
required=['FINAL_STATUS.json','SOURCE_AND_LIBRARY_IDENTITY.json','SMALL_CRACK_DESIGN_FREEZE.json','GEOMETRY_AND_INITIAL_CUT_AUDIT.json','NATIVE_DYNAMIC_HISTORY.csv','FIRST_SOFTENING_AND_FAILURE.json','CRACK_LOCALIZATION_AUDIT.json','CROSS_PLANE_BRIDGE_AUDIT.json','CRACK_AREA_IDENTIFIABILITY.json','NATIVE_COMPONENT_TOPOLOGY.json','DYNAMIC_ENERGY_LEDGER.json','NUMERICAL_VALIDITY.json','RUN_RESOURCES.json','REPRODUCE.txt','OPEN_RESULTS.html','REVIEW_'+TASK.upper()+'.md','RESOURCE_STOP_DIAGNOSIS.json','FAILURE_RECORD_INTERPRETATION.json','EVENT_STREAM_FORMAT_AUDIT.json']
for n in required:assert (O/n).is_file(),n
assert len(list(F.glob('*.png')))==len(list(F.glob('*.svg')))==6
assert (O/'full/FAILURE_RECORD.json').is_file() and not (O/'full/final_native_topology.npz').exists()
assert len(gzip.open(O/'full/damage_updates.bin.gz','rb').read())==0
contact_raw=gzip.open(O/'full/contact_pairs.bin.gz','rb').read();contacts=np.frombuffer(contact_raw,dtype=[('step','<i4'),('i','<i4'),('j','<i4'),('distance_m','<f8')]);assert len(contacts)==296 and np.all(contacts['distance_m']==DX)
dump(O/'CONTACT_STREAM_AUDIT.json',dict(geometric_candidate_records=len(contacts),candidate_steps=[int(contacts['step'].min()),int(contacts['step'].max())],all_distances_equal_dx=True,force_bearing_contact_records=0,max_overlap_m=0.,raw_first_contact_means='First geometric r<=dx candidate, NOT first force-bearing r<dx contact',actual_nonzero_contact_observed=False,raw_stream_preserved=True))
hist=np.genfromtxt(O/'NATIVE_DYNAMIC_HISTORY.csv',delimiter=',',names=True)
assert len(hist)==201 and all(np.isfinite(hist[k]).all() for k in hist.dtype.names)
assert np.all(hist['atom_count']==3080) and np.all(hist['total_natural_failed']==0)
events=[json.loads(q) for q in (O/'EVENTS_NORMALIZED.jsonl').read_text().splitlines()];assert sum(q['kind']=='ACCEPTED' for q in events)==201
cp=json.loads((O/'full/LATEST_NATIVE_CHECKPOINT.json').read_text());assert cp['time_s']==0 and sha(cp['restart'])==cp['restart_sha256'] and sha(cp['state'])==cp['state_sha256']
s=np.load(O/'full/last_accepted.npz');assert int(s['step'])==200 and abs(float(s['time'])-50e-6)<1e-15
tree=ast.parse((A/'journal.py').read_text());fun=next(q for q in ast.walk(tree) if isinstance(q,ast.FunctionDef) and q.name=='failure');attrs={q.attr for q in ast.walk(fun) if isinstance(q,ast.Attribute)};assert not attrs.intersection({'gather','state','command','get_thermo'})
identity=json.loads((O/'SOURCE_AND_LIBRARY_IDENTITY.json').read_text());assert sha(LIB)==identity['library_sha256']==freeze['library_sha256']
build=json.loads((O/'BUILD_AUDIT.json').read_text())
for p,h in build['object_sha256'].items():assert sha(p)==h,('reused build object changed',p)
res=json.loads((O/'RUN_RESOURCES.json').read_text());res['isolated_builds']=[dict(stage='development_01_not_used_for_physics',**{k:v for k,v in json.loads((O/'development_build_01/BUILD_AUDIT.json').read_text()).items() if k in ['wall_s','CPU_s','peak_RSS_bytes','compiled_translation_units']}),dict(stage='corrected_const_reference_used_for_physics',**{k:v for k,v in build.items() if k in ['wall_s','CPU_s','peak_RSS_bytes','compiled_translation_units']})]
# Boundary-bin sensitivity is an analysis caveat; keep frozen floor assignment untouched.
ref=np.load(O/'provenance/reference.npz');q=(ref['cross_intersections'][:,:2]+SIZE[:2]/2)/DX;boundary=np.any(abs(q-np.round(q))<1e-9,axis=1)
dump(O/'PLANE_BINNING_AUDIT.json',dict(frozen_assignment='floor of FP64 original intersection coordinate, clipped to valid grid; unchanged after run',intersections_on_or_numerically_near_cell_boundaries=int(boundary.sum()),total_cross_intersections=len(q),caveat='Roundoff at exact cell boundaries can distribute tied crossings asymmetrically between adjacent cells. This is a bin-count visualization convention, not asymmetric damage.',effect_on_current_conclusion='None: zero natural damage, all initially live bridges remain alive, zero qualifying newly separated cells. Any future nonzero area claim needs an explicit robust boundary-tie convention before its formal run.',posthoc_area_or_physics_change=False))
new=f'''# 当前clot原生软化开发入口

最新任务：{TASK} / {RUN}（2026-10-14 Europe/Zurich）。

- **RESOURCE_ADMISSION_FAIL；自然裂纹/面积—耗能验证未完成。** 唯一正式路径从原始预切t=0到50µs，200接受步、0拒绝；计划2500µs仅完成2%。没有追加正式计算或修改冻结几何/载荷/材料。
- 3080实际点、144448原始参考键、3524初始预切，700×500×1100µm胞元体积尺寸；z=0预裂纹至x=−50µm，4576初始存活跨面键全保留。m保留完整参考邻域，预切D=0。
- 上下各140点仅z规定运动、x/y自由；2800内部全自由。120µm每侧半余弦2ms张开+0.5ms保持，实际到50µs。隔离fix_move::reset_dt允许成对NULL/equal；材料、接触和VV C++未变。
- 自然软化/新完全失效均0，最大自然κ=0.0031091083<s0=0.15；全部ID保留、140924存活键、一个3080点分量。无法判断裂尖/夹持区先损伤；不判定GRIP_DOMINATED_FAILURE或模型不会断裂。
- 前缀最大归一化能量余额0.0020129321%，能量/独立力学/VV一致性仅限0–50µs通过。50µs W/U/K约2.639713/0.730706/1.909028pJ，D=Uc=Jc=0，K/U=2.6126，非准静态认证。
- 主段墙钟45.700s、CPU43.929s、峰值11.2397GiB>冻结4GiB。read_data在special_bonds设置之前构建无用高阶特殊邻接；P1已约11.26GiB但正式前资源门槛漏检，资源守卫直到接受第200步后才检查。见RESOURCE_STOP_DIAGNOSIS.json。
- 最后合法FP64态50µs，最近原生restart仅t=0。异常后不做原生gather/写restart；通用失败记录candidate_not_admitted字段对此接受后停止不适用，见FAILURE_RECORD_INTERPRETATION.json。末态native端点检查NOT_PERFORMED，不是检测到端点失配。
- 先读verification/{TASK}/{RUN}/FINAL_STATUS.json、34项REVIEW与SMALL_CRACK_DESIGN_FREEZE.json。新驱动analysis/small_crack_v1/；已有run不得覆盖重跑。
- 六组真实PNG/SVG：results/{TASK}/{RUN}/OPEN_RESULTS.html。无首次断裂事件的图明确未观察到；不生成连续医学外壳。CRACK_AREA_UNRESOLVED，G_SEPARATION_PROTOCOL/GC_EFFECTIVE均UNAVAILABLE，宏观强度/V3生产均未认证。
- 新隔离库{LIB.relative_to(R)}，SHA {sha(LIB)}。旧研究库SHA 08ccf24110e2cf2d07dbc8c0f31c0f368368b19b5427e6e5de9cc108a5ae115d、生产库和env.sh均不变。
- **最新已成功的完整动力学基线仍是** lps_native_softening_domain_retention_and_failure_safe_output_v1 / 20261014T150825Z：341点、700µs、2804接受/4拒绝、最大Δd0.0978676344、最大相对余额0.372426823%、4937失效/7110存活、分量192/147/1/1。其结果不因本轮工程资源失败被撤销。
- WSL实际/home/lzy/projects/clot_rupture.md已插入本轮负向结果，before/after/diff/SHA及读回完整保存；①–⑨历史、旧G0标称值保留。
- 下一轮唯一建议：同一物理配置下把special_bonds设置移到read_data之前，验证初始拓扑排除/力/能量一致并先通过实测峰值内存准入，再经审查决定完整工况。此修复/第二正式工况本轮NOT_RUN。V3、CFD、微泡、AI、GPU/MPI扫描均NOT_RUN。
- 前一活跃入口完整保存在verification/{TASK}/{RUN}/ACTIVE_DEVELOPMENT_ENTRYPOINTS.before.md。封存后停止，等待科学审查。
'''
tmp=entry.with_suffix('.md.tmp');tmp.write_text(new);os.replace(tmp,entry);(O/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.after.md').write_text(new);(O/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.diff').write_text(''.join(difflib.unified_diff(before.read_text().splitlines(True),new.splitlines(True),fromfile='ACTIVE.before',tofile='ACTIVE.after')))
verification=dict(status='DELIVERY_VERIFIED_SCIENTIFIC_TARGET_NOT_COMPLETED',utc=datetime.now(timezone.utc).isoformat(),all_required_files_present=True,six_png_six_svg_present=True,frozen_runtime_and_reference_unchanged=True,prior_manifest_files_verified_before_entry_update=len(prior['files']),only_previous_manifest_path_advanced='ACTIVE_DEVELOPMENT_ENTRYPOINTS.md; exact old bytes saved in new run',production_library_environment_old_research_library_unchanged=True,reused_build_object_hashes_unchanged=True,WSL_actual_readback_record_verified=True,raw_failure_preserved_with_explicit_interpretation=True,formal_case_count=1,additional_physical_cases_after_stop=0,diagnostic_steps=3,no_post_exception_native_reads=True,native_final_endpoint_gather='UNAVAILABLE_NOT_A_MISMATCH',last_accepted_FP64_time_s=float(s['time']),latest_native_checkpoint_time_s=0.,all_python_files_parse=True,prior_341_full700us_status='PRESERVED_PASS')
dump(O/'DELIVERY_VERIFICATION.json',verification)
status.update(delivery_status='COMPLETE_ARCHIVED_NEGATIVE_RESULT',WSL_document_updated_and_read_back=True,source_and_library_identities_preserved=True,all_required_outputs_present=True,six_png_six_svg=True,scientific_target_completed=False);dump(O/'FINAL_STATUS.json',status)
res.update(output_inventory_scope='New drivers, full run evidence, two isolated build versions, six figure groups, WSL before/after/diff, active entry',output_bytes_before_manifest=sum(p.stat().st_size for folder in [A,O,F] for p in folder.rglob('*') if p.is_file() and '__pycache__' not in p.parts));dump(O/'RUN_RESOURCES.json',res)
files={}
for folder in [A,O,F]:
 for p in sorted(folder.rglob('*')):
  if p.is_file() and '__pycache__' not in p.parts:files[str(p.relative_to(R))]=dict(bytes=p.stat().st_size,sha256=sha(p))
files[str(entry.relative_to(R))]=dict(bytes=entry.stat().st_size,sha256=sha(entry))
dump(O/'PAYLOAD_MANIFEST.json',dict(task=TASK,run_id=RUN,scientific_result=status['status'],files=files,file_count=len(files),total_bytes=sum(z['bytes'] for z in files.values()),self_and_sidecar_excluded=True,reused_dependencies={str(OLDLIB):sha(OLDLIB),str(PREV/'PAYLOAD_MANIFEST.json'):sha(PREV/'PAYLOAD_MANIFEST.json')},portable_new_library=str(LIB.relative_to(R)),formal_cases=1))
(O/'PAYLOAD_MANIFEST.sha256').write_text(sha(O/'PAYLOAD_MANIFEST.json')+'  PAYLOAD_MANIFEST.json\n')
print('SEALED_NEGATIVE_RESULT',len(files),sum(z['bytes'] for z in files.values()),sha(O/'PAYLOAD_MANIFEST.json'))
