#!/usr/bin/env python3
"""Read-only archive verification; materialize byte-identical scratch copies."""
import csv,hashlib,json,os,posixpath,shutil,sys
from pathlib import Path
A=Path(sys.argv[1]); W=Path(sys.argv[2]); S=A/'payload'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(2**20),b''): h.update(b)
 return h.hexdigest()
assert (A/'PACKAGE_METADATA_SHA256SUMS').is_file()
records={r['relative_path']:r for r in csv.DictReader((A/'manifests/MINIMAL_BUNDLE_FILE_MANIFEST.tsv').open(),delimiter='\t')}
def resolve(rel):
 parts=rel.split('/')
 for i in range(1,len(parts)+1):
  n='/'.join(parts[:i]);r=records[n]
  if r['type']=='symlink':
   target=r['link_target']
   if target.startswith('/workspace/'): target=target[len('/workspace/'):]
   elif not target.startswith('/'): target=posixpath.normpath(posixpath.join(posixpath.dirname(n),target))
   else: raise RuntimeError('External archive link '+n)
   return resolve(posixpath.join(target,*parts[i:]))
 return rel
rows=[]; working={}
def copy(rel,dest,role):
 actual=resolve(rel); r=records[actual]; src=S/actual
 if r['type']=='directory':
  if dest: dest.mkdir(parents=True,exist_ok=True)
  for p in sorted(src.iterdir()):copy(actual+'/'+p.name,dest/p.name if dest else None,role)
 else:
  assert r['type']=='regular' and sha(src)==r['sha256'],actual
  rows.append(dict(artifact=rel,archive_path=str(src),work_path=str(dest) if dest else 'PROVENANCE_ONLY',sha256=r['sha256'],role=role,required='YES'))
  if dest:
   dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dest)
   assert sha(dest)==r['sha256']; working[str(dest.relative_to(W))]=r['sha256']
P='hemocell_gpu_poc/'; O=P+'20260913_183657/'; F=P+'20260914_stage4_batch_dispatch/'
for d in ['source','palabos','build','official_gpu_smoke','stage4_200','stage4_1000','stage4_5000','toolchain_env','slurm','logs','verification','reports','provenance','scripts']:(W/d).mkdir(parents=True,exist_ok=True)
for rel,dest,role in [
 (O+'provenance/nvhpc_download.json','provenance/nvhpc_download.json','FROZEN_OFFICIAL_NVHPC_DOWNLOAD'),
 (P+'20260914_stage3_sustained/STAGE3_EXECUTION_AND_ANALYSIS_CONTRACT.json','provenance/STAGE3_EXECUTION_AND_ANALYSIS_CONTRACT.json','FROZEN_MEMORY_AND_SUSTAINED_RULES'),
 (F+'source','source','FROZEN_STAGE4_SOURCE'),(F+'upstream/palabos','palabos','STAGE4_PATCHED_PALABOS_EXPORT'),
 (O+'upstream/palabos','official_palabos','PINNED_OFFICIAL_GPU_SOURCE_WITH_ARCHIVED_200_STEP_CAP'),
 (O+'frozen_bundle','frozen_bundle','FROZEN_GEOMETRY_PHYSICS_AND_MONITOR_INPUTS'),
 (O+'nvhpc/tbb/prefix','tbb','ARCHIVED_TBB_DEPENDENCY_NOT_NVHPC_SDK'),
 (F+'runs/batch_correctness','references/rtx200','RTX5090_STAGE4_STRONG_REFERENCE'),
 (F+'runs/batch_1000_r1','references/rtx1000','RTX5090_STAGE4_1000_REFERENCE'),
 (F+'runs/batch_5000_r1','references/rtx5000','RTX5090_STAGE4_5000_REFERENCE'),
 (P+'20260914_stage3_sustained/runs/cpu_5000_r1','references/cpu5000','CPU_MPI12_5000_REFERENCE'),
 (F+'runs/cpu_sanity','references/cpu200','CPU_MPI12_200_REFERENCE'),
 (P+'20260914_stage3_sustained/runs/cpu_1000_r1','references/cpu1000','CPU_MPI12_1000_REFERENCE'),
 (F+'verification/GPU_NUMERICAL_COMPARISON_CONTRACT.json','verification/GPU_NUMERICAL_COMPARISON_CONTRACT.json','FROZEN_TOLERANCES_AND_COMMIT_EVIDENCE'),
 (F+'scripts/evaluate_case.py','scripts/evaluate_case.py','UNMODIFIED_FROZEN_NUMERICAL_EVALUATOR'),
 (F+'scripts/runtime_support.py','provenance/original_runtime_support.py','ORIGINAL_INPUT_MATERIALIZATION_RECIPE'),
 (F+'scripts/run_case.py','provenance/original_run_case.py','ORIGINAL_LAUNCH_RECIPE'),
 (F+'STAGE4_PATCH.diff','provenance/STAGE4_PATCH.diff','FROZEN_PATCH_ALREADY_APPLIED_DO_NOT_REAPPLY'),
 (F+'STAGE4_EXECUTION_CONTRACT.json','provenance/STAGE4_EXECUTION_CONTRACT.json','STAGE4_FREEZE'),
 (F+'provenance/execution_hashes.json','provenance/original_execution_hashes.json','STAGE4_SOURCE_FREEZE'),
 (F+'provenance/STATIC_MATH_REUSE_PROOF.json','provenance/STATIC_MATH_REUSE_PROOF.json','PATCH_PROVENANCE'),
 (F+'provenance/gpu_configure.log','provenance/original_gpu_configure.log','BUILD_RECIPE'),
 (F+'provenance/gpu_build_attempt2.log','provenance/original_gpu_build.log','BUILD_RECIPE'),
 (O+'provenance/official_build_receipt.json','provenance/original_official_build_receipt.json','OFFICIAL_BUILD_RECIPE'),
 (F+'GPU_STAGE4_REPORT.md','provenance/GPU_STAGE4_REPORT.md','FINAL_CANDIDATE_AND_PERFORMANCE'),
 (F+'FINAL_GPU_SELECTION.json','provenance/FINAL_GPU_SELECTION.json','FINAL_SELECTION'),
 (P+'20260914_stage3_sustained/GPU_STAGE3_REPORT.md','provenance/GPU_STAGE3_REPORT.md','STAGE3_REFERENCE_PROVENANCE'),
 (P+'formal_step3c_gpu_recovery_20260914_183542/FORMAL_STEP3C_GPU_FREEZE.json','provenance/FORMAL_STEP3C_GPU_FREEZE.json','PRODUCTION_MULTIPLIER_AND_COMMIT_PROVENANCE')]:
  print('COPY_VERIFY',rel,flush=True);copy(rel,W/dest,role)
# Archive binary is inspected and hashed only, never copied into an execution directory.
copy(F+'build_gpu/vascular/vascularPoC',None,'OLD_RTX5090_BINARY_PROVENANCE_ONLY')
# Finalization/audit scripts are evidence, never executed against the immutable archive.
for rel in sorted(records):
 if rel.startswith(P+'formal_step3c_gpu_recovery_20260914_183542/scripts/') and records[rel]['type']=='regular':
  copy(rel,W/'provenance/formal_finalizer_scripts'/Path(rel).name,'FROZEN_FINALIZATION_AND_AUDIT_SCRIPT')
with (W/'RESTORE_INPUT_MANIFEST.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=['artifact','archive_path','work_path','sha256','role','required'],delimiter='\t');w.writeheader();w.writerows(rows)
(W/'provenance/WORKING_COPY_SHA256.json').write_text(json.dumps(working,indent=2)+'\n')
(W/'verification/RESTORE_INPUT_VERIFICATION.json').write_text(json.dumps(dict(status='PASS',critical_restore_inputs='PASS',working_copy_source_integrity='PASS',files_checked=len(rows),copied_files=len(working),archive_read_only=True,symlinks_materialized_without_content_change=True,source_diff='NONE',archive_manifest_sha256=sha(A/'PACKAGE_METADATA_SHA256SUMS')),indent=2)+'\n')
print('RESTORE_INPUTS_PASS',len(rows),flush=True)

# Verify original runtime parameter files across the frozen reference cases and reuse bytes.
params={n:(W/f'references/rtx{n}/contracts/solver_parameters.txt').read_bytes() for n in (200,1000,5000)}
assert len(set(params.values()))==1
for n in (200,1000,5000):
 d=W/f'stage4_{n}'
 for sub in ['contracts','diagnostics/field_samples','logs','provenance']:(d/sub).mkdir(parents=True,exist_ok=True)
 for p in (W/f'references/rtx{n}/contracts').iterdir():
  if p.is_file():shutil.copy2(p,d/'contracts'/p.name)
 assert (d/'contracts/solver_parameters.txt').read_bytes()==params[n]
 (d/'provenance/input_materialization.json').write_text(json.dumps(dict(status='PASS',mode='Byte-identical copy from original validated same-horizon RTX5090 run',configuration_numbers_changed=False,external_step_cap=n,hashes={p.name:sha(p) for p in (d/'contracts').iterdir() if p.is_file()}),indent=2)+'\n')
lines=params[200].decode().splitlines();nums=lines[2].split()
assert nums[6:11]==['1.9989918081065344e-07','2.0366810646671923e-09','1.0','1056.0','2.7369132390905703e-15']
assert lines[:2]==['frozen_inputs/step1_geometry_contract/geometry/cfd_surface_axis_aligned_inlet_m.stl','frozen_inputs/step2']
assert json.loads((W/'verification/GPU_NUMERICAL_COMPARISON_CONTRACT.json').read_text())['physical_inputs']['multiplier']==1.1197286861799598
(W/'verification/SCIENTIFIC_INPUT_AUDIT.json').write_text(json.dumps(dict(status='PASS',runtime_configuration_byte_identical=True,physical_values_changed=False,source_diff='NONE',source_parameter_sha256=sha(W/'references/rtx200/contracts/solver_parameters.txt')),indent=2)+'\n')
(W/'PATH_REMAP_REPORT.md').write_text('''# 路径迁移记录

归档保持只读；仅在新工作副本中解析原/workspace命名空间并复制真实目标文件。归档symlink不改写，工作副本中的source、Palabos、输入均逐文件保持SHA256一致。完整OLD_PATH→NEW_PATH见RESTORE_INPUT_MANIFEST.tsv及WORKING_COPY_SHA256.json；前者记录archive源路径，后者记录work路径及相同SHA。

| OLD_PATH | NEW_PATH（相对本work目录） | REASON |
|---|---|---|
| /workspace/hemocell_gpu_poc/20260914_stage4_batch_dispatch/source | source/ | 原冻结C++与CMake原样保存 |
| 原Stage4 upstream/palabos | palabos/ | 固定导出源码，Stage4 patch已应用，不重复打补丁 |
| 原PoC upstream/palabos | official_palabos/ | 官方cavity及200步硬限原样复制 |
| 原PoC frozen_bundle | frozen_bundle/ | 几何/物理/监控原样复制 |
| 原NVHPC toolchains/nvhpc | /workspace/hemocell_restore/toolchains/nvhpc_26_5 | 本实例独立工具链安装 |
| 原run目录 | stage4_200、stage4_1000、stage4_5000 | 每次从零，使用原对应run的合同字节 |
| 旧build_gpu | build/ | 本平台重新编译，原binary只作provenance |

runtime参数文件直接复制原已通过的RTX5090 run，不做文本或数值替换。参数首两行本来就是相对frozen_bundle的路径。历史参数文件max_steps=700000原样保留，实际Stage4 C++外部硬限仅允许200/1000/5000，本任务不会执行长程。

新独立CMake配方只把路径与GPU架构cc120→实测cc89适配；所有C++及原CMake原样封存。OMP1、MPI1、STAGE4_LBM=batch、STAGE2_PROFILE_MODE=off沿用原已验证模式。
''')
print('RUNTIME_CONTRACT_BYTES_PASS',flush=True)
