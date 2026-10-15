"""Verify and seal the completed 2.4--2.5 ms native continuation, once."""
from common import *
import ast,difflib,re
from datetime import datetime,timezone
assert not (O/'PAYLOAD_MANIFEST.json').exists(),'Completed archives must not be overwritten'
def read(name):return json.loads((O/name).read_text())
st=read('FINAL_STATUS.json');freeze=read('SMALL_CRACK_DESIGN_FREEZE.json')
parentfreeze=json.loads((PARENT/'SMALL_CRACK_DESIGN_FREEZE.json').read_text())
assert sha(O/'SMALL_CRACK_DESIGN_FREEZE.json')==(O/'SMALL_CRACK_DESIGN_FREEZE.sha256').read_text().split()[0]=='b1e24b8a9996f86d095ea42b717c18e4da9fe1cc83744ae0cdee5b40ff0ff007'
assert freeze['budget']['RSS_stop_GiB'] is None and freeze['budget']['main_wall_s'] is None
assert freeze['budget']['output_limit_GiB']==parentfreeze['budget']['output_limit_GiB']==2
for k in ['geometry','material','gates','dt_base_s','minimum_h_s','checkpoint_predictor_trigger','actual_damage_gate','initial_cut_energy','loading','domain','scientific_criteria','output_policy']:
 assert freeze[k]==parentfreeze[k],k
for name,h in freeze['driver_sha256'].items():assert sha(A/name)==h,name
for name,h in freeze['reference_sha256'].items():assert sha(O/'provenance'/name)==h,name
for p,h in read('PROTECTED_BASELINES.json').items():assert sha(p)==h,p
assert sha(LIB)==freeze['library_sha256']=='730e11651aee6533b4fb7d38c8341e06ad853bf3397c4908fc6ef905d4e3ac90'
prior=json.loads((PARENT/'PAYLOAD_MANIFEST.json').read_text())
assert sha(PARENT/'PAYLOAD_MANIFEST.json')==(PARENT/'PAYLOAD_MANIFEST.sha256').read_text().split()[0]
for name,z in prior['files'].items():
 p=R/name;assert p.stat().st_size==z['bytes'] and sha(p)==z['sha256'],name
for p in A.glob('*.py'):ast.parse(p.read_text())
chain=read('CONTINUATION_CHAIN_AUDIT.json');restore=read('NATIVE_RESTART_RESTORE_AUDIT.json');ledger=read('LEDGER_RESTORE_AUDIT.json')
assert chain['status']==restore['status']==ledger['status']=='PASS'
assert chain['original_steps']==9600 and chain['additional_steps']==400 and chain['total_accepted_steps']==10000
assert chain['no_fresh_t0_run'] and chain['no_duplicate_accepted_step'] and chain['strictly_monotone_time']
assert chain['energy_balance_jump_J']==ledger['R_jump_J']==0
assert restore['physical_steps_advanced']==0 and restore['original_bond_ids_exact'] and restore['reference_m_exact'] and restore['native_endpoint_graph_exact']
cp=read('full/LATEST_NATIVE_CHECKPOINT.json')
for field in ['restart','state']:assert sha(cp[field])==cp[field+'_sha256']
assert abs(cp['time_s']-END)<1e-15 and cp['step']==10000
n=read('NUMERICAL_VALIDITY.json');e=read('DYNAMIC_ENERGY_LEDGER.json');r=read('RUN_RESOURCES.json');top=read('NATIVE_COMPONENT_TOPOLOGY.json')
assert st['status']==n['status']=='FULL_PROTOCOL_NUMERICAL_PASS' and st['full_protocol_completed'] and n['complete']
assert st['accepted_steps']==10000 and st['rejected_steps']==0 and abs(st['actual_end_s']-END)<1e-15
assert n['all_original_ids_preserved'] and n['native_final_endpoint_match'] and top['native_final_endpoint_match']
assert n['max_dd']<=.10 and e['max_relative_R_percent']<=2
assert not st['crack_propagation_candidate'] and not st['unique_surface_certified']
assert st['CRACK_AREA']=='UNRESOLVED' and st['G_SEPARATION_PROTOCOL']==st['GC_EFFECTIVE']=='UNAVAILABLE'
assert st['new_formal_cases_this_authorization']==0 and st['native_continuation_segments']==1
before=O/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.before.md';entry=R/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.md'
assert sha(before)==sha(entry)==prior['files']['ACTIVE_DEVELOPMENT_ENTRYPOINTS.md']['sha256']
doc=read('WSL_DOCUMENT_UPDATE.json');assert doc['status']=='PASS' and doc['readback_exact']
assert sha(O/'wsl_document/clot_rupture.after.md')==doc['after_sha256']
assert (O/'wsl_document/clot_rupture.after.md').read_bytes().replace((O/'WSL_DOCUMENT_INSERT.md').read_bytes(),b'',1)==(O/'wsl_document/clot_rupture.before.md').read_bytes()
required=['FINAL_STATUS.json','SOURCE_AND_LIBRARY_IDENTITY.json','SMALL_CRACK_DESIGN_FREEZE.json','CONTINUATION_AUTHORIZATION.json','GEOMETRY_AND_INITIAL_CUT_AUDIT.json','NATIVE_DYNAMIC_HISTORY.csv','FIRST_SOFTENING_AND_FAILURE.json','CRACK_LOCALIZATION_AUDIT.json','CROSS_PLANE_BRIDGE_AUDIT.json','CRACK_AREA_IDENTIFIABILITY.json','NATIVE_COMPONENT_TOPOLOGY.json','DYNAMIC_ENERGY_LEDGER.json','NUMERICAL_VALIDITY.json','RUN_RESOURCES.json','REPRODUCE.txt','OPEN_RESULTS.html','REVIEW_'+TASK.upper()+'.md','PREFIX_50US_COMPARISON.json','CONTACT_STREAM_AUDIT.json','CONTINUATION_CHAIN_AUDIT.json','NATIVE_RESTART_RESTORE_AUDIT.json','LEDGER_RESTORE_AUDIT.json']
for name in required:assert (O/name).is_file(),name
assert [int(i) for i in re.findall(r'^### (\d+)\.',(O/('REVIEW_'+TASK.upper()+'.md')).read_text(),re.M)]==list(range(1,35))
assert len(list(F.glob('*.png')))==len(list(F.glob('*.svg')))==6
# Scope inherited build provenance explicitly; neither C++ nor its library changed here.
identity=read('SOURCE_AND_LIBRARY_IDENTITY.json')
identity['historical_CPP_change_from_library_build']=identity.pop('changed_CPP')
identity.update(CPP_changes_this_run=[],continuation_run=RUN,previous_run=str(PARENT),library_build_source_run=str(PRIOR_SMALL),build_inputs=str(PRIOR_SMALL/'BUILD_AUDIT.json'),parent_identity_sha256=sha(PARENT/'SOURCE_AND_LIBRARY_IDENTITY.json'))
assert Path(identity['build_inputs']).is_file()
dump(O/'SOURCE_AND_LIBRARY_IDENTITY.json',identity)
new=f'''# 当前clot原生软化开发入口

当前：{TASK} / {RUN}（取消1800秒墙钟预算后，从原生2.4ms检查点续算）。

- **FULL_PROTOCOL_NUMERICAL_PASS**：原定2.5ms时程已完整完成，累计10000接受步、0损伤拒绝；本次仅追加400步，无第二个物理工况、无t=0重跑。
- 3080点、144448原始参考键、3524预切、700×500×1100µm试样及完整m/r0/V/taper不变。上下各140点仅规定z；每侧120µm、2ms半余弦+0.5ms保持不变。
- 墙钟和RSS上限均为null，没有替代上限；2GiB输出门槛和全部损伤、能量、运动学准则保留。复用SHA730e11651aee6533b4fb7d38c8341e06ad853bf3397c4908fc6ef905d4e3ac90隔离库，无C++新编译或P1新物理运行。
- 恢复完整原生x/v/κ/d、原始ID/m及端点图，累计W/D/Jc、接触功、尺度、四项数值缺陷和区域D连续，能量R跳变量=0。旧0–2.4ms标量、节点流、损伤和接触流保持精确前缀。
- 本次末态自然软化{st['natural_softened']}、自然完全失效{st['natural_failed']}；首软化747.5µs、首完全失效799.75µs，均位于裂尖侧表面，原始bond ID68923/77580。裂尖软化领先加载端91.25µs。
- 最终存活跨面3312（长键3266），分量点数{top['component_sizes']}；主连通体仍有3076点，4个小分量是单节点，不能称为已形成两个宏观分离体。
- 最大接受Δd={n['max_dd']:.12g}、最大归一化余额{e['max_relative_R_percent']:.12g}%；原生末态端点图与独立原始图核对通过，全部3080原始ID保留。
- NOTCH_LOCALIZED_DAMAGE说明真实裂尖初始局部损伤；末态开口候选带仅50µm，未达到100µm持续前缘推进准则。后期加载端和其他自由表面大量耗能，不能由分类名认证连续裂纹传播。
- CRACK_AREA_UNRESOLVED，G_SEPARATION_PROTOCOL/GC_EFFECTIVE均UNAVAILABLE。网格代理0.025mm²、严格horizon无桥接子集0；不以总D除以代理面积。末态K/U={e['last_K_over_U']:.9g}，明显动态；最大接触重叠0.154408631dx，未认证接触物理精度。
- 最新完整原生检查点{cp['time_s']*1e3:.9g}ms：{cp['restart']}及其.state.npz。原定窗口已结束，不自动延长物理时程。
- 本次续算墙钟{r['segment_wall_s']:.3f}s、CPU{r['segment_CPU_s']:.3f}s、峰值RSS{r['segment_peak_RSS_bytes']/1024**2:.3f}MiB；原路径累计墙钟{r['protocol_total_wall_s']:.3f}s、峰值RSS{r['protocol_peak_RSS_bytes']/1024**3:.6f}GiB。1rank×1线程，无GPU。
- 先读verification/{TASK}/{RUN}/FINAL_STATUS.json、34项REVIEW及CONTINUATION_CHAIN_AUDIT.json；代码analysis/small_crack_resume_no_wall_v1/，新运行先检查冻结身份，禁止覆盖本归档。
- 六组PNG/SVG：results/{TASK}/{RUN}/OPEN_RESULTS.html；使用实际点阵/事件/曲线，无医学外壳或位移放大。
- WSL实际/home/lzy/projects/clot_rupture.md已更新并读回，before/after/diff/SHA保存；原①–⑨、旧G0标称值及全部历史文字保留。
- 先前50µs内存停止、2.4ms墙钟停止历史保留；341点700µs成功基线20261014T150825Z保留。生产库、env.sh及历史物理数据不变；仅当前入口在保存旧字节后推进。
- 当前不认证宏观强度、全局dt收敛、材料Gc或V3生产；不自动改材料/几何/载荷、运行V3/CFD/微泡/AI。已完成本轮交付，等待下一轮科学审查。
'''
assert sha(entry)==sha(before),'Concurrent active-entry edit'
tmp=entry.with_suffix('.md.tmp');tmp.write_text(new);os.replace(tmp,entry)
(O/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.after.md').write_text(new)
(O/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.diff').write_text(''.join(difflib.unified_diff(before.read_text().splitlines(True),new.splitlines(True),fromfile='ACTIVE.before',tofile='ACTIVE.after')))
audit=dict(status='DELIVERY_VERIFIED',utc=datetime.now(timezone.utc).isoformat(),new_formal_cases=0,native_continuation_segments=1,additional_accepted_steps=400,full_protocol_completed=True,previous_manifest_files_verified_before_entry_update=len(prior['files']),only_previous_manifest_path_advanced='ACTIVE_DEVELOPMENT_ENTRYPOINTS.md; previous bytes preserved in new run',frozen_geometry_material_loading_gates_unchanged=True,wall_and_RSS_caps_disabled_by_user=True,no_replacement_wall_or_RSS_cap=True,production_and_libraries_unchanged=True,no_CPP_build=True,no_new_P1_steps=True,raw_0_to_2p4ms_exact_prefix=True,WSL_actual_readback_verified=True,six_PNG_six_SVG=True,all_python_files_parse=True,scientific_status=st['status'],latest_native_checkpoint_time_s=cp['time_s'])
dump(O/'DELIVERY_VERIFICATION.json',audit)
st.update(delivery_status='COMPLETE_ARCHIVED',wall_clock_cap_s=None,WSL_document_updated_and_read_back=True,six_png_six_svg=True)
dump(O/'FINAL_STATUS.json',st)
r.update(output_bytes_before_manifest=sum(p.stat().st_size for f in [A,O,F] for p in f.rglob('*') if p.is_file() and '__pycache__' not in p.parts),output_scope='Continuation drivers, full accepted trajectory, native checkpoints, events, six figures and WSL before/after; native library reused')
assert r['output_bytes_before_manifest']<2*1024**3
dump(O/'RUN_RESOURCES.json',r)
files={}
for folder in [A,O,F]:
 for p in sorted(folder.rglob('*')):
  if p.is_file() and '__pycache__' not in p.parts:files[str(p.relative_to(R))]=dict(bytes=p.stat().st_size,sha256=sha(p))
files[str(entry.relative_to(R))]=dict(bytes=entry.stat().st_size,sha256=sha(entry))
deps={str(p):sha(p) for p in [LIB,PRIOR_SMALL/'PAYLOAD_MANIFEST.json',PARENT/'PAYLOAD_MANIFEST.json',Path(restore['source_checkpoint']['restart']),Path(restore['source_checkpoint']['state'])]}
dump(O/'PAYLOAD_MANIFEST.json',dict(task=TASK,run_id=RUN,files=files,file_count=len(files),total_bytes=sum(x['bytes'] for x in files.values()),self_and_sidecar_excluded=True,reused_dependencies=deps,scientific_status=st['status']))
(O/'PAYLOAD_MANIFEST.sha256').write_text(sha(O/'PAYLOAD_MANIFEST.json')+'  PAYLOAD_MANIFEST.json\n')
print('SEALED',len(files),sha(O/'PAYLOAD_MANIFEST.json'))
