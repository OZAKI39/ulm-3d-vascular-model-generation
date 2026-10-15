from common import *
import ast,difflib
from datetime import datetime,timezone
assert not (O/'PAYLOAD_MANIFEST.json').exists()
st=json.loads((O/'FINAL_STATUS.json').read_text());freeze=json.loads((O/'SMALL_CRACK_DESIGN_FREEZE.json').read_text());oldfreeze=json.loads((PRIOR_SMALL/'SMALL_CRACK_DESIGN_FREEZE.json').read_text())
assert freeze['budget']['RSS_stop_GiB'] is None
for k in ['geometry','material','gates','dt_base_s','minimum_h_s','checkpoint_predictor_trigger','actual_damage_gate','initial_cut_energy','loading','domain','scientific_criteria','output_policy']:assert freeze[k]==oldfreeze[k],k
for name,h in freeze['driver_sha256'].items():assert sha(A/name)==h,name
for name,h in freeze['reference_sha256'].items():assert sha(O/'provenance'/name)==h,name
for p,h in json.loads((O/'PROTECTED_BASELINES.json').read_text()).items():assert sha(p)==h,p
assert sha(LIB)==freeze['library_sha256']=='730e11651aee6533b4fb7d38c8341e06ad853bf3397c4908fc6ef905d4e3ac90'
prior=json.loads((PRIOR_SMALL/'PAYLOAD_MANIFEST.json').read_text())
for n,z in prior['files'].items():assert sha(R/n)==z['sha256'],n
for p in A.glob('*.py'):ast.parse(p.read_text())
before=O/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.before.md';entry=R/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.md';assert sha(before)==sha(entry)
doc=json.loads((O/'WSL_DOCUMENT_UPDATE.json').read_text());assert doc['status']=='PASS' and doc['readback_exact'];assert sha(O/'wsl_document/clot_rupture.after.md')==doc['after_sha256'];assert (O/'wsl_document/clot_rupture.after.md').read_bytes().replace((O/'WSL_DOCUMENT_INSERT.md').read_bytes(),b'',1)==(O/'wsl_document/clot_rupture.before.md').read_bytes()
required=['FINAL_STATUS.json','SOURCE_AND_LIBRARY_IDENTITY.json','SMALL_CRACK_DESIGN_FREEZE.json','CONTINUATION_AUTHORIZATION.json','GEOMETRY_AND_INITIAL_CUT_AUDIT.json','NATIVE_DYNAMIC_HISTORY.csv','FIRST_SOFTENING_AND_FAILURE.json','CRACK_LOCALIZATION_AUDIT.json','CROSS_PLANE_BRIDGE_AUDIT.json','CRACK_AREA_IDENTIFIABILITY.json','NATIVE_COMPONENT_TOPOLOGY.json','DYNAMIC_ENERGY_LEDGER.json','NUMERICAL_VALIDITY.json','RUN_RESOURCES.json','REPRODUCE.txt','OPEN_RESULTS.html','REVIEW_'+TASK.upper()+'.md','PREFIX_50US_COMPARISON.json','CONTACT_STREAM_AUDIT.json']
for n in required:assert (O/n).is_file(),n
assert len(list(F.glob('*.png')))==len(list(F.glob('*.svg')))==6
cp=json.loads((O/'full/LATEST_NATIVE_CHECKPOINT.json').read_text());assert sha(cp['restart'])==cp['restart_sha256'] and sha(cp['state'])==cp['state_sha256']
if st['full_protocol_completed']:assert abs(cp['time_s']-END)<1e-15
n=json.loads((O/'NUMERICAL_VALIDITY.json').read_text());e=json.loads((O/'DYNAMIC_ENERGY_LEDGER.json').read_text());r=json.loads((O/'RUN_RESOURCES.json').read_text())
new=f'''# 当前clot原生软化开发入口

当前：{TASK} / {RUN}（用户明确取消4 GiB RSS停止条件后的同工况新运行）。

- 状态 **{st['status']}**；实际到{st['actual_end_s']*1e3:.9g}ms，{st['accepted_steps']}接受/{st['rejected_steps']}损伤拒绝；计划2.5ms完成={st['full_protocol_completed']}。停止原因{n['stop_reasons']}。
- 3080点、144448原始参考键、3524预切、700×500×1100µm试样及所有原始m/r0/V/taper不变。上下各140点仅规定z、2800内部全自由；每侧120µm、2ms半余弦+0.5ms保持不变。
- 仅取消RSS停止门槛，无替代内存上限；原1800s墙钟/2GiB输出及所有数值门槛保留。复用SHA730e11651aee6533b4fb7d38c8341e06ad853bf3397c4908fc6ef905d4e3ac90隔离库，无C++新编译、无P1新物理运行。
- 旧50µs仅FP64而完整restart在0，故按用户授权从原始预切t=0重新推进。前201标量行/11节点帧与旧0–50µs逐值一致，无前缀拼接。旧资源失败及P1漏检历史保留。
- 本次自然软化{st['natural_softened']}，新增完全失效{st['natural_failed']}；首软化{st['first_softening_time_s']}s、首完全失效{st['first_failure_time_s']}s；区域分类{st['scientific_classification']}。存活跨面{st['live_cross_bonds']}，分量节点数{st['component_sizes']}。
- 最大接受Δd={n['max_dd']:.12g}，最大归一化余额{e['max_relative_R_percent']:.12g}%；全部3080原始ID保留。数值准入、损伤分类与新生裂面认证分开。
- CRACK_AREA_UNRESOLVED，G_SEPARATION_PROTOCOL/GC_EFFECTIVE均UNAVAILABLE；宏观强度、全局dt收敛及V3生产仍未认证。
- 末态仅50µm无桥接开口带，未达到100µm持续前缘推进准则；严格horizon无桥接子集为空。裂尖先软化，但末态大量失效位于端部/其他自由表面，不能把分类名理解成连续裂纹传播已认证。
- 最新完整原生检查点时间{cp['time_s']*1e3:.9g}ms：{cp['restart']}及其.state.npz，已有SHA。若下一轮补足剩余保持段，应恢复原生历史与Python能量账本，不能直接把现有fresh-t0入口称为续算；本轮未追加恢复运行。
- 主段墙钟{r['main_wall_s']:.3f}s，CPU{r['main_CPU_s']:.3f}s，峰值RSS{r['peak_RSS_bytes']/1024**3:.5f}GiB；1 rank×1线程，GPU未用。资源详情见RUN_RESOURCES.json。
- 先读verification/{TASK}/{RUN}/FINAL_STATUS.json、34项REVIEW和NUMERICAL_VALIDITY.json；代码analysis/small_crack_no_rss_v1/。已有输出禁止覆盖或自动追加工况。
- 六组PNG/SVG：results/{TASK}/{RUN}/OPEN_RESULTS.html；全部真实点阵/事件/曲线，无医学外壳或位移放大。
- WSL实际/home/lzy/projects/clot_rupture.md已更新并读回，before/after/diff/SHA完整保存；①–⑨历史与旧G0标称值保持。
- 341点700µs历史成功基线保留：lps_native_softening_domain_retention_and_failure_safe_output_v1 / 20261014T150825Z。生产库、env.sh及其他历史数据不变。
- 下一步依据本次具体事件与现有数值停止原因审查，不自动放宽门槛、不改材料/几何/载荷、不运行V3/CFD/微泡/AI。前入口备份在本run的ACTIVE_DEVELOPMENT_ENTRYPOINTS.before.md。
'''
tmp=entry.with_suffix('.md.tmp');tmp.write_text(new);os.replace(tmp,entry);(O/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.after.md').write_text(new);(O/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.diff').write_text(''.join(difflib.unified_diff(before.read_text().splitlines(True),new.splitlines(True),fromfile='ACTIVE.before',tofile='ACTIVE.after')))
audit=dict(status='DELIVERY_VERIFIED',utc=datetime.now(timezone.utc).isoformat(),new_formal_runs_authorized=1,previous_manifest_files_verified_before_entry_update=len(prior['files']),only_previous_manifest_path_advanced='ACTIVE_DEVELOPMENT_ENTRYPOINTS.md; previous bytes preserved in new run',frozen_geometry_material_loading_gates_unchanged=True,RSS_stop_disabled_by_user=True,no_replacement_RSS_cap=True,production_and_libraries_unchanged=True,no_CPP_build=True,no_new_P1_steps=True,prefix_50us_exact=True,WSL_actual_readback_verified=True,six_PNG_six_SVG=True,all_python_files_parse=True,scientific_status=st['status'],latest_native_checkpoint_time_s=cp['time_s'])
dump(O/'DELIVERY_VERIFICATION.json',audit);st.update(delivery_status='COMPLETE_ARCHIVED',WSL_document_updated_and_read_back=True,six_png_six_svg=True);dump(O/'FINAL_STATUS.json',st)
r.update(output_bytes_before_manifest=sum(p.stat().st_size for f in [A,O,F] for p in f.rglob('*') if p.is_file() and '__pycache__' not in p.parts),output_scope='New continuation drivers, accepted trajectory, checkpoints, events, analyses, six figures and WSL before/after; prior library reused');dump(O/'RUN_RESOURCES.json',r)
files={}
for folder in [A,O,F]:
 for p in sorted(folder.rglob('*')):
  if p.is_file() and '__pycache__' not in p.parts:files[str(p.relative_to(R))]=dict(bytes=p.stat().st_size,sha256=sha(p))
files[str(entry.relative_to(R))]=dict(bytes=entry.stat().st_size,sha256=sha(entry))
dump(O/'PAYLOAD_MANIFEST.json',dict(task=TASK,run_id=RUN,files=files,file_count=len(files),total_bytes=sum(x['bytes'] for x in files.values()),self_and_sidecar_excluded=True,reused_dependencies={str(LIB):sha(LIB),str(PRIOR_SMALL/'PAYLOAD_MANIFEST.json'):sha(PRIOR_SMALL/'PAYLOAD_MANIFEST.json')},scientific_status=st['status']))
(O/'PAYLOAD_MANIFEST.sha256').write_text(sha(O/'PAYLOAD_MANIFEST.json')+'  PAYLOAD_MANIFEST.json\n');print('SEALED',len(files),sha(O/'PAYLOAD_MANIFEST.json'))
