"""P0 call-chain audit and P1 full-server baseline verification; no simulation."""
from common import *
import ast,shutil,csv,time
assert not O.exists(),'Use a new task directory, never overwrite an accepted archive'
for p in [O,F,O/'provenance',O/'derived',O/'benchmark']:p.mkdir(parents=True,exist_ok=True)
started=time.perf_counter();old=R/'analysis/v3_softening_pilot_v1'
chain={}
for p in sorted(old.glob('*.py')):
 tree=ast.parse(p.read_text());imports=[]
 for node in ast.walk(tree):
  if isinstance(node,ast.ImportFrom):imports.append(node.module)
  elif isinstance(node,ast.Import):imports.extend(x.name for x in node.names)
 chain[p.name]=dict(sha256=sha(p),imports=sorted(set(x for x in imports if x)),functions=[x.name for x in tree.body if isinstance(x,(ast.FunctionDef,ast.ClassDef))])
dump(O/'P0_CALL_CHAIN.json',chain)
(O/'P0_CLEANUP_AUDIT.md').write_text('''# P0 小范围清理与保留审计

实际运行链：launch.py → run_protocol.py → Native/native.py → 固定原始 reference、bpm/peri/variational 与 LPS/VARIATIONAL/DYNAMIC → external.py 的真实图/载荷/拖曳。Integrator 执行 trial/check/reject/subdivide；mechanics.py 独立复核完整力、能量、接触与 VV；journal.py 保存接受标量、稀疏损伤/接触和节点帧。analyze/figures/report 是归档后处理，prepare/freeze/repair_registry_build 是构建与冻结依据，archive_rejections/local_readback 是无损数据与交付核验入口。

逐文件 Python AST/导入/函数及 SHA 见 P0_CALL_CHAIN.json。四个隔离 C++ 文件与旧 endpoint-star Python 参考模型均由正式库/冻结清单引用；本轮不删除或重写。多个旧轮次中重复的 Contact、Native、图统计服务各自冻结复现；在没有证明其无历史价值时，不把“没有当前调用”误当作废弃。

本轮组织调整：新增独立目录；用 common.py 统一只读路径、SI 常量、SHA、原始状态读取，并直接导入已验证的完整 endpoint-star 能量实现。后续拓扑、候选、接触和小试样共用单个 graph_ops.py，避免重复实现材料/锚点可达性。未发现可证明安全删除且无历史价值的旧科研源码，因此旧源码删除数为 0；不为清理数量破坏正式验证或原生构建依赖。未作无关重构。
''')
manifest_path=B/'PAYLOAD_MANIFEST.json';manifest=json.loads(manifest_path.read_text());errors=[]
for name,entry in manifest['files'].items():
 p=R/name
 if not p.is_file() or p.stat().st_size!=entry['bytes'] or sha(p)!=entry['sha256']:errors.append(name)
assert not errors,errors
ref=reference();assert len(ref['X'])==3488 and int(ref['clot_count'])==1927 and len(ref['bonds'])==157637
assert np.array_equal(ref['original_bond_ids'],np.arange(1,157638))
assert np.array_equal(np.bincount(ref['bond_types'])[1:],[82270,22040,53327])
hist=np.genfromtxt(B/'V3_NATIVE_PILOT_HISTORY.csv',delimiter=',',names=True)
layout=json.loads((B/'full/STREAM_LAYOUT.json').read_text());frames=np.genfromtxt(B/'full/native_frames.csv',delimiter=',',names=True)
assert (B/'full/native_frames.bin').stat().st_size==len(frames)*np.prod(layout['nodes_per_record_shape'])*8
keys=list(keyframes());last=keys[-1][1]
assert len(hist)==10503 and int(last['step'])==10502 and abs(float(last['time'])-.00233)<1e-15
for name in ['FINAL_STATUS.json','V3_PILOT_FREEZE.json','V3_MATERIAL_CARD.json','V3_ENERGY_LEDGER.json','STREAM_LAYOUT.json']:
 src=B/name if (B/name).is_file() else B/'full'/name;shutil.copy2(src,O/'provenance'/name)
shutil.copy2(R/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.md',O/'provenance/ACTIVE_DEVELOPMENT_ENTRYPOINTS.before.md')
dump(O/'BASELINE_AUDIT.json',dict(status='PASS',source=str(B),manifest_sha256=sha(manifest_path),verified_files=len(manifest['files']),verified_bytes=manifest['total_bytes'],original_bond_identity=True,counts=dict(material=1927,anchors=1561,bulk=82270,interface=22040,scaffold=53327),saved_accepted_node_frames=len(frames),full_history_rows=len(hist),complete_typed_history_keyframes=len(keys),all_step_damage_changes=True,all_step_contact_pairs=True,complete_every_step_node_force_audit=False,missing_data=['Every-accepted-step x/v/force/material_force: NOT SAVED; raw frame every 20 accepted steps and selected keys only','Subthreshold kappa between complete keyframes: NOT SAVED; d and failure times are exact from sparse updates'],event_1943_exact_node_frame_present=bool(np.any(frames['step']==1943) or any(int(s['step'])==1943 for p,s in keys)),library_sha256=sha(LIB),baseline_recomputed=False,wall_s=time.perf_counter()-started))
dump(O/'PROTECTED_INPUTS.json',{str(p):sha(p) for p in [manifest_path,B/'V3_PILOT_FREEZE.json',LIB,R/'analysis/softening_feasibility_v1/model.py',R/'tools/lammps-native-y-wall-shared-v1-build/liblammps.so']})
print('P0/P1 PASS',len(manifest['files']),'files;',len(frames),'node frames;',len(keys),'complete keys',flush=True)
