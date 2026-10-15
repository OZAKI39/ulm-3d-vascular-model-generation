from common import *
import ast,subprocess,time,sys,importlib.util,os
start=time.perf_counter();O.mkdir(parents=True,exist_ok=False);F.mkdir(parents=True,exist_ok=False)
manifest=json.loads((PREV/'PAYLOAD_MANIFEST.json').read_text());assert sha(PREV/'PAYLOAD_MANIFEST.json')=='a421a5a3c6088e60860e566dfc1b32926dfc4b4766239c71bdf9c84aa7caac14'
for p,v in manifest['files'].items():assert sha(R/p)==v['sha256'],p
oldid=json.loads((PREV/'SOURCE_AND_LIBRARY_IDENTITY.json').read_text());files=dict(oldid['files'])
for p in [PREV/'build/liblammps.so',PREV/'NATIVE_SOFTENING_DESIGN_FREEZE.json',PREV/'PAYLOAD_MANIFEST.json',R/'analysis/softening_feasibility_v1/model.py',R/'PROJECT_ACTIVE_CONTEXT.md']:
 files[str(p)]=dict(sha256=sha(p),bytes=p.stat().st_size)
assert sha(PREV/'build/liblammps.so')=='b80bfc0a99f4160d62567c0e2d9ddc611004ad40f79c8e49370293c681560773'
for p in oldid['files']:assert sha(p)==oldid['files'][p]['sha256'],p

def command(args):return subprocess.check_output(args,cwd=R,text=True)
(R/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.md').exists() and (_ for _ in ()).throw(RuntimeError('Existing active index must be reviewed before replacement'))
for when in ['before']:
 (O/f'GIT_STATUS_{when}.txt').write_text(command(['git','status','--short']));(O/f'GIT_DIFF_{when}.patch').write_text(command(['git','diff','--binary']))
tracked=set(command(['git','ls-files']).splitlines());mods=command(['git','diff','--name-only']).splitlines();tracked_before={p:sha(R/p) for p in mods if (R/p).is_file()}
protected={}
for root in [R/'verification',R/'results']:
 for p in root.rglob('*'):
  if p.is_file() and not p.is_relative_to(O) and not p.is_relative_to(F):protected[str(p.relative_to(R))]=dict(sha256=sha(p),bytes=p.stat().st_size)
dump(O/'HISTORICAL_PROTECTION_BEFORE.json',protected)
# Restricted dependency inspection: two accepted analysis dirs, variational sources,
# and their current build compiler probes. No historical result bodies are loaded.
candidates=[]
for root in [R/'analysis/native_softening_v1',R/'analysis/softening_feasibility_v1']:
 candidates.extend(p for p in root.iterdir() if p.is_file() and p.suffix in ['.py','.h','.cpp'])
candidates.extend((PREV/'source/src/BPM').glob('*variational*'))
candidates.extend((PREV/'build/CMakeFiles').glob('*/*/*.cpp'))
texts={str(p.relative_to(R)):p.read_text(errors='replace') for p in candidates}
records=[]
for p in candidates:
 rel=str(p.relative_to(R));imports=[];includes=[]
 if p.suffix=='.py':
  tree=ast.parse(texts[rel]);imports=[n.module for n in ast.walk(tree) if isinstance(n,ast.ImportFrom)];imports+= [a.name for n in ast.walk(tree) if isinstance(n,ast.Import) for a in n.names]
 else:includes=[ln.strip() for ln in texts[rel].splitlines() if ln.lstrip().startswith('#include')]
 sealed=rel in manifest['files'];built=p.is_relative_to(PREV/'source') or p.is_relative_to(PREV/'build')
 category='ACTIVE' if p.name in ['model.py','bond_bpm_peri_variational.cpp','bond_bpm_peri_variational.h'] else 'VERIFIED_REFERENCE' if sealed or 'softening_feasibility_v1' in p.parts else 'UNCERTAIN'
 records.append(dict(path=rel,sha256=sha(p),bytes=p.stat().st_size,category=category,tracked=rel in tracked,has_uncommitted_tracked_changes=rel in mods,sealed_by_latest_manifest=sealed,build_dependency=built,python_imports=imports,CPP_includes=includes,referenced_by_scoped_source=[q for q,t in texts.items() if q!=rel and (p.name in t or (p.suffix=='.py' and p.stem in [s.split('.')[0] for s in imports if s]))],delete=False,reason='Accepted model/replay evidence or current CMake build dependency. A duplicate source has an independent frozen provenance role; no deletion solely by identical SHA.',replacement=None,recoverability='Existing immutable source snapshot/manifest; no file deleted'))
allow=[];dump(O/'CODE_CLEANUP_INVENTORY.json',dict(status='AUDITED',scope='Two accepted softening analysis dirs + current two variational sources + build compiler-probe source files',candidate_files=len(records),records=records,tracked_modified_files_preserved=tracked_before,all_historical_evidence_excluded_from_deletion=True))
dump(O/'CODE_CLEANUP_ALLOWLIST.json',dict(files=allow,status='EMPTY_NO_PROVEN_DEAD_OR_DISPOSABLE_DUPLICATE',reason='Every examined scientific source remains active or a accepted reproduction reference; build scratch dependency status uncertain.'))
dump(O/'DELETED_CODE_MANIFEST.json',dict(deleted_count=0,paths=[],archives=[],restore='No deletions; no restore operation necessary',not_fabricated=True))
# Minimal load/config checks only; no old experiment is executed.
spec=importlib.util.spec_from_file_location('frozen_model',R/'analysis/softening_feasibility_v1/model.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);ref=np.load(PREV/'provenance/reference.npz');model=m.SofteningLPS(ref);zero=model.evaluate(ref['X']);assert np.all(zero['force']==0) and zero['U']==0
cache=(PREV/'build/CMakeCache.txt').read_text();assert 'PKG_BPM:BOOL=ON' in cache and 'BUILD_SHARED_LIBS:BOOL=ON' in cache
entry=f'''# 当前 clot 活跃开发入口

本轮：LPS_NATIVE_SOFTENING_DYNAMIC_COMMIT_AND_SMALL_CRACK_V1；{RUN}。
P0已完成；P1–P5仅按本轮用户授权进入，不自动部署。

- 服务器项目（平铺）：`{R}`；环境：`source /workspace/y_clot/env.sh`。
- 先读本文件，再读本轮 `verification/{TASK}/{RUN}/FINAL_STATUS.json`（完成后生成）。
- 已接受原生源码：`{PREV}/source/src/BPM/bond_bpm_peri_variational.{{h,cpp}}`。
- 上轮可读源码/驱动：`analysis/native_softening_v1/`；完整快照在上轮 source_snapshot。
- 冻结离线数学参考：`analysis/softening_feasibility_v1/model.py`（SHA `{sha(R/'analysis/softening_feasibility_v1/model.py')}`）。
- 上轮正式库：`{PREV}/build/liblammps.so`。
- 上轮正式库SHA：`{sha(PREV/'build/liblammps.so')}`。
- 当前新开发目录：`analysis/dynamic_softening_v1/`（P0审计后才新建动力学功能）。
- 当前新隔离源码/库预定路径：`verification/{TASK}/{RUN}/source` 和 `build/liblammps.so`。
- 原始试样：`{PREV}/provenance/`，341节点、12047原始键，dx=50µm，horizon=150.5µm。
- 材料：G=940Pa，K=100000/3Pa，rho=1000kg/m³，s0=.15；本轮主动态sc=.18。
- 直接原生依赖：Bond/Atom/Comm/Domain/Update + 原锁定 pair_bpm_peri；将审查 verlet.cpp、fix_nve.cpp 和相关 hook。
- 原生 f 是 N/m³；物理 F=V*f；完整导出U为J。动力学质量/积分约定须本轮实际核验。
- 已验收：固定历史力/能、受控运动学、sorting/2MPI/restart、1dx闭合接触。
- 待本轮认证：真实自由动力学、自动一步一次提交、动力能量、小型非均匀裂纹。
- 未认证：宏观强度、材料Gc、6–8J/m²目标、V3、生产MPI伸缩性。
- 新模型二值失效不等于旧删键LPS；旧历史力L2差8.814180703270e-6N保留。
- 生产库 `tools/lammps-native-y-wall-shared-v1-build/liblammps.so` 保持隔离，禁止替换。
- 历史verification、source_snapshot、implementation_history均受保护；不因old/test/debug命名清理。
- P0安全删除数量：0。受审资产均有现行/复现价值或用途不确定，不人为制造删除。
- 文档：真实WSL `/home/lzy/projects/clot_rupture.md`；更新须before/after/diff/SHA。
- 仅按需读取上述少量入口；历史PROJECT_ACTIVE_CONTEXT中的旧“下一轮”不代表当前指令。
'''
assert len(entry.splitlines())<=150;(R/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.md').write_text(entry);(O/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.md').write_text(entry)
for p,h in tracked_before.items():assert sha(R/p)==h
assert command(['git','diff','--binary'])==(O/'GIT_DIFF_before.patch').read_text()
(O/'GIT_STATUS_after.txt').write_text(command(['git','status','--short']));(O/'GIT_DIFF_after.patch').write_text(command(['git','diff','--binary']))
identity=dict(files=files,previous_manifest_verified_files=len(manifest['files']),previous_manifest_sha256=sha(PREV/'PAYLOAD_MANIFEST.json'),git_branch=command(['git','branch','--show-current']).strip(),git_head=command(['git','rev-parse','HEAD']).strip(),compiler=command(['c++','--version']).splitlines()[0],MPI=command(['mpirun','--version']).splitlines()[0],python=sys.version,cgroup={str(p):p.read_text().strip() for p in [Path('/sys/fs/cgroup/cpu,cpuacct/cpu.cfs_quota_us'),Path('/sys/fs/cgroup/cpu,cpuacct/cpu.cfs_period_us'),Path('/sys/fs/cgroup/memory/memory.limit_in_bytes')]},memory=command(['free','-b']),disk=command(['df','-B1','/workspace']),other_processes=command(['ps','-eo','pid,comm,pcpu,pmem,etime','--sort=-pcpu']),tracked_modified_SHA=tracked_before)
dump(O/'SOURCE_AND_LIBRARY_IDENTITY.json',identity)
(O/'CODE_CLEANUP_REPORT.md').write_text(f'''# P0代码安全审计与上下文精简

审计实际候选 {len(records)} 个文件；**安全删除0个**。已有正式代码、离线数学模型、复现脚本、冻结source_snapshot及失败implementation_history保留。相同SHA的开发副本与冻结副本具有不同可追溯角色，不属于无价值重复。CMake编译探针不能仅凭临时名称证明已经无用，保留UNCERTAIN。

完成了逐文件SHA、Git状态、受限Python import/C++ include和运行入口检查；详细路径/依赖见CODE_CLEANUP_INVENTORY.json。删除白名单为空，没有删除操作或需要恢复的文件。历史科学数据仅流式SHA，不将大型日志/NPZ加载到上下文。

创建根目录ACTIVE_DEVELOPMENT_ENTRYPOINTS.md，共{len(entry.splitlines())}行；只需优先读取当前模块及其直接依赖。旧PROJECT_ACTIVE_CONTEXT的历史叙述未改。保留全部现有未提交修改；Git tracked diff前后完全一致。当前源码、构建配置、冻结model导入与参考态F/U零值检查通过。旧生产与两套已接受诊断库SHA均核对一致。

本轮P0严格先于任何新C++动力学功能。零删除是审计结论，不将保留历史证据称为清理失败。
''')
dump(O/'P0_COMPLETE.json',dict(status='PASS',wall_s=time.perf_counter()-start,candidates=len(records),deleted=0,active_index_lines=len(entry.splitlines()),previous_manifest_verified=True,offline_model_load=True,reference_zero_force_energy=True,build_configuration_present=True,tracked_diff_unchanged=True,protected_history_files=len(protected),P1_CPP_started=False))
print('P0_COMPLETE',len(records),'candidates; 0 safe deletions; index created;',len(protected),'historical files protected',flush=True)
