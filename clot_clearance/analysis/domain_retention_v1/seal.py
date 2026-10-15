from common import *
import difflib,shutil,ast
from datetime import datetime,timezone
assert not (O/'PAYLOAD_MANIFEST.json').exists(),'Already sealed'
status=json.loads((O/'FINAL_STATUS.json').read_text());assert status['status']=='PASS'
freeze=json.loads((O/'DOMAIN_PROTOCOL_FREEZE.json').read_text())
for n,h in freeze['driver_sha256'].items():assert sha(A/n)==h,(n,'frozen code changed')
for p,h in json.loads((O/'PROTECTED_BASELINES.json').read_text()).items():assert sha(p)==h,p
old=json.loads((PREV/'PAYLOAD_MANIFEST.json').read_text());before=O/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.before.md'
for n,m in old['files'].items():assert sha(R/n)==m['sha256'],('old artifact modified before entry update',n)
assert sha(R/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.md')==sha(before)
doc=json.loads((O/'WSL_DOCUMENT_UPDATE.json').read_text());assert doc['status']=='PASS' and doc['readback_exact']
assert sha(O/'wsl_document/clot_rupture.after.md')==doc['actual_readback_sha256']
oldtext=(O/'wsl_document/clot_rupture.before.md').read_bytes();aftertext=(O/'wsl_document/clot_rupture.after.md').read_bytes();insert=(O/'WSL_DOCUMENT_INSERT.md').read_bytes();assert aftertext.replace(insert,b'',1)==oldtext
required=['DOMAIN_CHANGE_AUDIT.json','INITIAL_STATE_PARITY.json','FAILURE_SAFE_OUTPUT_AUDIT.json','FULL_700US_DYNAMIC_AUDIT.json','ADAPTIVE_SUBSTEP_HISTORY.csv','ORIGINAL_ATOM_IDENTITY_AUDIT.json','FINAL_NATIVE_TOPOLOGY.json','INDEPENDENT_COMPONENT_RECONSTRUCTION.json','DYNAMIC_ENERGY_LEDGER.json','CONTACT_AND_OVERLAP_AUDIT.json','HISTORICAL_TRAJECTORY_COMPARISON.json','RUN_RESOURCES.json','REPRODUCE.txt','OPEN_RESULTS.html','REVIEW_'+TASK.upper()+'.md']
for n in required:assert (O/n).is_file(),n
assert len(list(F.glob('*.png')))==len(list(F.glob('*.svg')))==5
assert not (O/'full/FAILURE_RECORD.json').exists()
assert json.loads((O/'full/PROCESS_RESOURCES.json').read_text())['exit_code']==0
assert (O/'fault_test/FAILURE_RECORD.json').exists()
cp=json.loads((O/'full/LATEST_NATIVE_CHECKPOINT.json').read_text());assert sha(cp['restart'])==cp['restart_sha256'] and sha(cp['state'])==cp['state_sha256'] and abs(cp['time_s']-.0007)<1e-15
# Mechanical scope check: failure handler has no native reads or restart operations.
tree=ast.parse((A/'journal.py').read_text());fun=next(q for q in ast.walk(tree) if isinstance(q,ast.FunctionDef) and q.name=='failure');attrs={q.attr for q in ast.walk(fun) if isinstance(q,ast.Attribute)};assert not attrs.intersection({'gather','state','command','get_thermo'})
new=f'''# 当前clot原生软化开发入口

本轮：{TASK} / {RUN}（2026-10-14 Europe/Zurich）。

- **DOMAIN_RETENTION_PASS、FULL_700US_PASS、ENERGY_PASS、TOPOLOGY_PASS全部通过。** 唯一正式路径从原始无损伤t=0到700µs；未拼接旧前缀。
- 固定fff数值域由各轴±250µm扩大到±400µm，不remap、不加墙、不改周期性。依据旧真实范围、最大速度和余下时间规划。盒子不是物理血管壁。
- 全程341个原始ID保留；最小域余量145.134130µm。节点73在688µs实际y=250.054230µm，700µs为254.865870326µm，确认旧域丢点原因。
- 2804接受/4损伤拒绝，8个0.125µs接受步、2796个0.25µs步；冻结最小0.03125µs、实际损伤门槛0.10及原回滚算法保持。最大Δd=0.09786763438279739；最大相对能量余额0.372426823%。
- 正常原生末态gather核验341点、7110存活键与独立d<1原图端点集合完全一致；4937完全失效，4个实际分量192/147/1/1。节点73和105为孤立点。尚非已认证真实血栓碎块。
- W/U/K/D/Uc=4.575844829448/0.525683244103/0.674093509273/3.272865781378/0.126299310063nJ，Jc=0，R=-23.097015368980pJ。沿用完整变分能量账本和累计归一化。
- 旧0–687.75µs共2756标量时刻/18完整态比较相容；每步dt、失效数、接触数一致。旧逐步x/v缺失，只比较实际可用数据；未重跑旧域。
- 初始旧/新域物理状态逐值一致；683µs旧检查点上的一次异常输出测试通过，候选/接受态分开，失败处理不调用native gather或write_restart。
- 每步刷新FP64节点/事件/接触对流，周期两槽原生restart，沿用损伤检查点。最终有效原生restart与FP64配对见full/LATEST_NATIVE_CHECKPOINT.json。
- 先读verification/{TASK}/{RUN}/FINAL_STATUS.json、22项REVIEW及DOMAIN_PROTOCOL_FREEZE.json。新活跃驱动analysis/domain_retention_v1/，封存目录禁止重跑覆盖。
- 五组PNG/SVG：results/{TASK}/{RUN}/OPEN_RESULTS.html。真实点阵显示，不生成连续外壳。
- 复用旧隔离库{LIB.relative_to(R)}，SHA {sha(LIB)}；本轮无C++改动、无重新编译。生产库和env.sh不变。
- 主计算墙钟54.473s、CPU35.475s、峰值659.688MiB；1 MPI×1线程，无GPU。完整资源与输出规模见RUN_RESOURCES.json/PAYLOAD_MANIFEST.json。
- 1.0dx接触保持，最大重叠0.246649196dx≈12.33246µm，物理精度尚未认证。全局时间收敛、宏观强度、材料Gc和完整V3生产就绪仍未认证。
- WSL实际/home/lzy/projects/clot_rupture.md已更新并读回；before/after/diff/SHA见WSL_DOCUMENT_UPDATE.json。全部历史内容及G0标称值保留。
- 下一轮唯一建议LPS_NATIVE_SOFTENING_SMALL_NONAFFINE_CRACK_VALIDATION_V1；小型非均匀裂纹可研究，本轮NOT_RUN。未运行V3、CFD、微泡或AI，交付后停止等审查。
- 上一入口完整备份：verification/{TASK}/{RUN}/ACTIVE_DEVELOPMENT_ENTRYPOINTS.before.md；所有历史科研结果保留。
'''
entry=R/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.md';tmp=entry.with_suffix('.md.tmp');tmp.write_text(new);os.replace(tmp,entry);(O/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.after.md').write_text(new);(O/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.diff').write_text(''.join(difflib.unified_diff(before.read_text().splitlines(True),new.splitlines(True),fromfile='ACTIVE.before',tofile='ACTIVE.after')))
audit=dict(status='PASS',utc=datetime.now(timezone.utc).isoformat(),frozen_driver_hashes_unchanged=True,prior_manifest_286_files_verified_before_entry_update=True,only_old_manifest_item_changed='ACTIVE_DEVELOPMENT_ENTRYPOINTS.md (intentionally advanced; exact historical version preserved)',production_library_environment_model_unchanged=True,final_checkpoint_hashes_match=True,WSL_actual_readback_record_verified=True,five_png_and_five_svg_present=True,formal_process_exit_code=0,fault_test_count=1,formal_case_count=1,no_post_failure_native_reads_in_handler=True)
dump(O/'DELIVERY_VERIFICATION.json',audit)
status.update(WSL_document_updated_and_read_back=True,delivery_verification='PASS',source_and_library_identities_preserved=True,all_required_outputs_present=True);dump(O/'FINAL_STATUS.json',status)
resources=json.loads((O/'RUN_RESOURCES.json').read_text());resources.update(output_inventory_scope='New analysis directory, new verification run, five figure groups, active entry; reused old source/library not duplicated',output_bytes_before_manifest=sum(p.stat().st_size for folder in [A,O,F] for p in folder.rglob('*') if p.is_file() and '__pycache__' not in p.parts));dump(O/'RUN_RESOURCES.json',resources)
files={}
for folder in [A,O,F]:
 for p in sorted(folder.rglob('*')):
  if p.is_file() and '__pycache__' not in p.parts:files[str(p.relative_to(R))]=dict(bytes=p.stat().st_size,sha256=sha(p))
files[str(entry.relative_to(R))]=dict(bytes=entry.stat().st_size,sha256=sha(entry))
manifest=dict(task=TASK,run_id=RUN,files=files,file_count=len(files),total_bytes=sum(q['bytes'] for q in files.values()),self_and_sidecar_excluded=True,reused_dependencies={str(LIB):sha(LIB),str(PREV/'PAYLOAD_MANIFEST.json'):sha(PREV/'PAYLOAD_MANIFEST.json')},portable_library_existing_copy=str(PREV/'binary/liblammps.so'),scientific_cases=1)
dump(O/'PAYLOAD_MANIFEST.json',manifest);(O/'PAYLOAD_MANIFEST.sha256').write_text(sha(O/'PAYLOAD_MANIFEST.json')+'  PAYLOAD_MANIFEST.json\n')
print('SEALED',len(files),manifest['total_bytes'],sha(O/'PAYLOAD_MANIFEST.json'))
