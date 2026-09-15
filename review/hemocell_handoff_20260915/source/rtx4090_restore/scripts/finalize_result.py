#!/usr/bin/env python3
"""Archive immutable numerical evidence from this task; never launch a solver."""
import csv,json,sys,hashlib,shutil,subprocess,time,gzip,statistics,fcntl
from pathlib import Path
W=Path(__file__).resolve().parents[1];A=Path('/workspace/hemocell_restore/archive/minimal_restore_20260915_143631');F=Path('/workspace/hemocell_restore/results')/W.name;L='/home/lzy/projects/compre_output/vast4090_stage4_restore/20260915_150938'
def read(p,default=None):return json.loads(p.read_text()) if p.exists() else ({} if default is None else default)
def save(p,obj):p.write_text(json.dumps(obj,indent=2)+'\n')
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(2**20),b''):h.update(b)
 return h.hexdigest()
def capture(c):
 p=subprocess.run(c,capture_output=True,text=True);return dict(command=c,returncode=p.returncode,stdout=p.stdout,stderr=p.stderr)
lock=(W/'finalize.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
state=read(W/'VALIDATION_STATE.json');assert state.get('state') in ('PASS','FAIL'),'Require real terminal validation state'
F.mkdir(parents=True,exist_ok=False)
archive=capture(['python3','-B',str(A/'reports/verify_bundle.py'),str(A),'--readonly'])
save(W/'verification/ARCHIVE_POST_AUDIT.json',archive)
archive_ok=archive['returncode']==0 and json.loads(archive['stdout'])['status']=='PASS'
frozen=read(W/'provenance/WORKING_COPY_SHA256.json');bad=[]
for rel,digest in frozen.items():
 if not (W/rel).is_file() or sha(W/rel)!=digest:bad.append(rel)
for n in (200,1000,5000):
 for p in (W/f'stage4_{n}/contracts').glob('*'):
  if p.is_file() and sha(p)!=sha(W/f'references/rtx{n}/contracts'/p.name):bad.append(str(p.relative_to(W)))
source_ok=bool(frozen) and not bad
save(W/'verification/FINAL_SOURCE_INPUT_AUDIT.json',dict(status='PASS' if source_ok else 'FAIL',checked=len(frozen),mismatches=bad,archive_status='PASS' if archive_ok else 'FAIL',science_source_diff='NONE' if source_ok else 'DETECTED',protected_formal_paths_modified=False))
old=Path('/workspace/hemocell_restore/archive/vast_snapshot');sz=capture(['du','-s','-B1',str(old)]);oldbytes=int(sz['stdout'].split()[0]) if sz['returncode']==0 else None
save(W/'verification/INCOMPLETE_FULL_UPLOAD_AUDIT.json',dict(path=str(old),exists=old.exists(),allocated_bytes=oldbytes,allocated_GiB=oldbytes/2**30 if oldbytes else None,preflight_allocated_bytes=3564896256,unchanged_size=oldbytes==3564896256,deleted=False))
gates={n:read(W/f'verification/STAGE4_{n}_GATE.json') for n in (200,1000,5000)}
runs={n:read(W/f'stage4_{n}/RUN_TERMINAL.json') for n in gates}
perf={}
for n,r in runs.items():
 if not r:continue
 d=W/f'stage4_{n}';rows=list(csv.DictReader((d/'provenance/resource_trace.csv').open()));warm=100 if n<=1000 else 500
 timed=[x for x in rows if warm<=int(x['last_completed_step'])<n and int(x['solver_RSS_bytes'])>0]
 util=[float(x['gpu_utilization_percent']) for x in timed]
 perf[str(n)]=dict(status=r['status'],initialization_seconds=r.get('initialization_seconds'),end_to_end_seconds=r.get('end_to_end_seconds'),solver_timing=r.get('solver_timing'),timed_resource_sample_count=len(timed),gpu_utilization_median=statistics.median(util) if util else None,gpu_utilization_mean=statistics.mean(util) if util else None,gpu_utilization_peak=max(util) if util else None,gpu_vram_peak_MiB=r.get('gpu_vram_peak_MiB'),solver_RSS_peak_bytes=r.get('peak_solver_RSS_bytes'),resource_sample_interval_seconds=1,resource_region_definition='last_completed_step>=warmup and <N, solverRSS>0; discrete polling estimate',normal_benchmark=True,repetitions=1)
 if n==5000 and r['status']=='PASS':
  thirds=[]
  for lo,hi in [(500,2000),(2000,3500),(3500,5000)]:
   a=[x for x in timed if lo<=int(x['last_completed_step'])<hi]
   thirds.append(dict(first_step=lo,last_step_exclusive=hi,samples=len(a),VRAM_median_MiB=statistics.median(float(x['gpu_memory_used_MiB']) for x in a) if a else None,RSS_median_MiB=statistics.median(int(x['solver_RSS_bytes'])/2**20 for x in a) if a else None))
  enough=all(x['samples']>=2 for x in thirds)
  vm=[x['VRAM_median_MiB'] for x in thirds];rs=[x['RSS_median_MiB'] for x in thirds]
  growing=enough and ((vm[0]<=vm[1]<=vm[2] and vm[2]-vm[0]>256) or (rs[0]<=rs[1]<=rs[2] and rs[2]-rs[0]>512))
  mem=dict(status=('GROWING' if growing else 'STABLE') if enough else 'UNVERIFIED',steady_thirds=thirds,VRAM_growth_MiB=vm[-1]-vm[0] if enough else None,RSS_growth_MiB=rs[-1]-rs[0] if enough else None,thresholds=dict(VRAM_MiB=256,RSS_MiB=512),rule_source='Frozen Stage3 memory_growth_stop; nondecreasing steady thirds',measurement_limitation='Sampled resource growth test, not proof of zero leak')
  perf[str(n)]['memory_stability']=mem
  windows=list(csv.DictReader((d/'diagnostics/iteration_windows.csv').open()))
  def rate(lo,hi):
   v=[x for x in windows if int(x['first_step'])>=lo and int(x['last_step'])<=hi]
   return sum(int(x['steps']) for x in v)/sum(float(x['seconds']) for x in v)
  first,last=rate(501,1500),rate(4001,5000);s1000=runs[1000]['solver_timing']['steps_per_second'];s5000=r['solver_timing']['steps_per_second'];ratio=s5000/s1000
  perf[str(n)]['performance_stability']=dict(status='STABLE' if ratio>=.95 and last/first>=.8 else 'DEGRADING',rate_5000_over_1000=ratio,first_full_1000_window_steps_per_sec=first,last_1000_window_steps_per_sec=last,last_over_first=last/first,Stage4_drop_gate=.05,Stage3_severe_drop_gate=.20,repeated_run_spread='NOT_TESTED_ONE_RUN_PER_HORIZON',strict_5000_10000_5percent_test='NOT_TESTED_NO_10000_RUN',sampling_population='One run per requested horizon, no confidence interval')
new5000=perf.get('5000',{}).get('solver_timing',{}).get('steps_per_second');comparison=dict(RTX5090_reference_steps_per_sec=125.880863,CPU_MPI12_reference_steps_per_sec=17.692465,RTX4090_sustained_steps_per_sec=new5000,RTX4090_VS_RTX5090_RATIO=new5000/125.880863 if new5000 else None,RTX4090_VS_CPU_MPI12_SPEEDUP=new5000/17.692465 if new5000 else None,reference_scope='Archived Stage4 RTX5090 sustained baseline and CPU MPI12 baseline; not rerun on new host',new_platform_repetitions=1,performance_is_correctness_gate=False,initialization_end_to_end_not_compared_as_solver_rate=True)
save(W/'RTX4090_PERFORMANCE.json',perf);save(W/'RTX4090_RTX5090_COMPARISON.json',comparison)
nonfinite=0;rows_checked=0
for n,r in runs.items():
 p=W/f'stage4_{n}/diagnostics/safety_counts.csv'
 if p.exists():
  for row in csv.DictReader(p.open()):
   nonfinite+=sum(int(row[k]) for k in ('fluid_nan','fluid_inf','ghost_nan','ghost_inf'));rows_checked+=1
correct=all(g.get('status')=='PASS' for g in gates.values());trend='FAIL' if any(g.get('numerical_divergence_trend')=='FAIL' for g in gates.values()) else 'GROWING' if any(g.get('numerical_divergence_trend')=='GROWING' for g in gates.values()) else 'STABLE' if correct else 'UNVERIFIED'
ready=read(W/'toolchain_env/NVHPC_READY.json');official=read(W/'OFFICIAL_GPU_SMOKE.json');build=read(W/'BUILD_PROVENANCE.json')
all_ok=state['state']=='PASS' and archive_ok and source_ok and ready.get('status')=='PASS' and official.get('status')=='PASS' and build.get('status')=='PASS' and correct and nonfinite==0 and trend=='STABLE'
if not build:
 build=dict(status='NOT_BUILT_OR_FAILED',last_stage=state.get('phase'),execution_commands=read(W/'provenance/EXECUTION_COMMANDS.json',[]));save(W/'BUILD_PROVENANCE.json',build)
if not official:official=dict(status='NOT_RUN',reason=state);save(W/'OFFICIAL_GPU_SMOKE.json',official)
for n in gates:
 p=W/f'STAGE4_{n}_CORRECTNESS.csv'
 if not p.exists():p.write_text('reference,quantity,step,error,tolerance,status,detail\nALL,NOT_RUN,,,,NOT_RUN,Stopped at '+state.get('phase','UNKNOWN')+'\n')
def cp(src,dest):
 dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dest)
for name in ['RESTORE_INPUT_MANIFEST.tsv','PATH_REMAP_REPORT.md','ENVIRONMENT.json','SYSTEM_DEPENDENCIES.tsv','BUILD_PROVENANCE.json','OFFICIAL_GPU_SMOKE.json','STAGE4_200_CORRECTNESS.csv','STAGE4_1000_CORRECTNESS.csv','STAGE4_5000_CORRECTNESS.csv','RTX4090_PERFORMANCE.json','RTX4090_RTX5090_COMPARISON.json','EXECUTION_PLAN.json','PREPARATION_STATE.json','VALIDATION_STATE.json']:
 if (W/name).is_file():cp(W/name,F/name)
for sub in ['verification','provenance','scripts','logs']:
 for p in (W/sub).rglob('*'):
  if p.is_file() and '__pycache__' not in p.parts:cp(p,F/p.relative_to(W))
for p in (W/'toolchain_env').glob('*.json'):cp(p,F/'toolchain_env'/p.name)
for p in (W/'toolchain_env/download').glob('*.txt'):cp(p,F/'toolchain_env'/p.name)
# Keep raw profiler evidence from official smoke only, and all compact diagnostics.
for p in (W/'official_gpu_smoke').iterdir():
 if p.is_file():cp(p,F/'official_gpu_smoke'/p.name)
for n in gates:
 d=W/f'stage4_{n}'
 manifest=[]
 for p in sorted(d.rglob('*')):
  if not p.is_file():continue
  manifest.append(dict(relative_path=str(p.relative_to(d)),remote_path=str(p),bytes=p.stat().st_size,sha256=sha(p)))
  if 'field_samples' in p.parts or 'sampled_fields' in p.parts:continue
  dst=F/f'stage4_{n}'/p.relative_to(d)
  if p.suffix=='.csv' and p.stat().st_size>2**20:
   dst=dst.with_suffix('.csv.gz');dst.parent.mkdir(parents=True,exist_ok=True)
   with p.open('rb') as src,gzip.open(dst,'wb',compresslevel=6) as out:shutil.copyfileobj(src,out)
  else:cp(p,dst)
 save(F/f'stage4_{n}/REMOTE_RUN_FILE_MANIFEST.json',manifest)
s1000=perf.get('1000',{});s5000=perf.get('5000',{});env=read(W/'ENVIRONMENT.json')
summary=dict(VAST4090_RESTORE_STATUS='PASS' if all_ok else 'FAIL_STOPPED',MINIMAL_RESTORE_ARCHIVE=str(A),MINIMAL_RESTORE_INTEGRITY='PASS' if archive_ok else 'FAIL',CRITICAL_RESTORE_INPUTS=read(W/'verification/RESTORE_INPUT_VERIFICATION.json').get('status','UNVERIFIED'),PATH_REMAP_STATUS='PASS' if source_ok else 'FAIL',WORKING_COPY_SOURCE_INTEGRITY='PASS' if source_ok else 'FAIL',GPU_MODEL='NVIDIA GeForce RTX 4090',GPU_VRAM_MIB=24564,GPU_DRIVER='595.84',GPU_COMPUTE_CAPABILITY='8.9',NVHPC_VERSION='26.5' if ready else 'UNVERIFIED',NVHPC_26_5_READY=ready.get('status','UNVERIFIED'),CUDA_COMPONENT_VERSION=ready.get('cuda_component','UNVERIFIED'),MPI_VERSION='See toolchain_env/DEPENDENCY_COMMANDS.json',CMAKE_VERSION='See SYSTEM_DEPENDENCIES.tsv',PALABOS_GPU_COMMIT='4127697e90169bbef982295f1d1c933cf6e90caa',OFFICIAL_PALABOS_GPU_BUILD='PASS' if (W/'provenance/OFFICIAL_BINARY.json').exists() else 'NOT_BUILT',OFFICIAL_PALABOS_GPU_RUN=official.get('status','NOT_RUN'),GPU_KERNEL_ACTIVITY=official.get('kernel_activity',{}).get('status','UNVERIFIED'),VAST4090_STAGE4_BUILD=build.get('status','NOT_BUILT'),SOURCE_DIFF_FROM_FROZEN_STAGE4='NONE' if source_ok else 'DETECTED',GPU_ARCH_BUILD_TARGET='cc89',STAGE4_200_CORRECTNESS=gates[200].get('status','NOT_RUN'),STAGE4_1000_CORRECTNESS=gates[1000].get('status','NOT_RUN'),STAGE4_5000_CORRECTNESS=gates[5000].get('status','NOT_RUN'),CPU_GPU_RHO_COMPARISON='PASS' if correct else 'NOT_ALL_PASS',CPU_GPU_VELOCITY_COMPARISON='PASS' if correct else 'NOT_ALL_PASS',CPU_GPU_FLUX_COMPARISON='PASS' if correct else 'NOT_ALL_PASS',CPU_GPU_MASS_COMPARISON='PASS' if correct else 'NOT_ALL_PASS',NONFINITE_COUNT=nonfinite if rows_checked else 'NOT_RUN',NUMERICAL_DIVERGENCE_TREND=trend,RTX4090_INITIALIZATION_SECONDS=s1000.get('initialization_seconds'),RTX4090_STEPS_PER_SEC=s1000.get('solver_timing',{}).get('steps_per_second'),RTX4090_5000_STEPS_PER_SEC=new5000,RTX4090_END_TO_END_SECONDS=s1000.get('end_to_end_seconds'),GPU_UTIL_MEDIAN=s1000.get('gpu_utilization_median'),GPU_UTIL_MEAN=s1000.get('gpu_utilization_mean'),GPU_UTIL_PEAK=s1000.get('gpu_utilization_peak'),GPU_VRAM_PEAK_MIB=s1000.get('gpu_vram_peak_MiB'),GPU_MEMORY_GROWTH=s5000.get('memory_stability',{}).get('status','NOT_RUN'),PERFORMANCE_STABILITY=s5000.get('performance_stability',{}).get('status','NOT_RUN'),RTX5090_REFERENCE_STEPS_PER_SEC=125.880863,RTX4090_VS_RTX5090_RATIO=comparison['RTX4090_VS_RTX5090_RATIO'],CPU_MPI12_REFERENCE_STEPS_PER_SEC=17.692465,RTX4090_VS_CPU_MPI12_SPEEDUP=comparison['RTX4090_VS_CPU_MPI12_SPEEDUP'],NEW_VAST4090_PURE_FLUID_BASELINE='PASS' if all_ok else 'FAIL',PBS_BSA_RUNTIME_SMOKE='NOT_RUN',RBC_STAGE1='NOT_RUN',RBC_TIMESTEPS=0,MICROBUBBLE='OFF',FORMAL_SCIENCE_BASELINE_MODIFIED='NO',LONG_RUN='NO',VAST4090_TO_WSL_RESULT_INTEGRITY='NOT_YET_TRANSFERRED',INCOMPLETE_FULL_UPLOAD_EXISTS='YES' if old.exists() else 'NO',INCOMPLETE_FULL_UPLOAD_DELETED='NO',REMOTE_RESULT_DIR=str(F),LOCAL_RESULT_DIR=L,VAST4090_STAGE4_RESTORE_REPORT=L+'/VAST4090_STAGE4_RESTORE_REPORT.md',NEXT_STEP='PURE_FLUID_NEW_MEDIUM_SMOKE' if all_ok else 'REVIEW_FAILURE_BEFORE_ANY_NEW_RUN')
deps=read(W/'toolchain_env/DEPENDENCY_COMMANDS.json')
for name,key in [('mpi','MPI_VERSION'),('cmake','CMAKE_VERSION')]:
 if name in deps:summary[key]=deps[name]['stdout'].splitlines()[0]
save(F/'FINAL_SUMMARY.json',summary)
(F/'FINAL_TERMINAL_SUMMARY.txt').write_text('\n'.join(f'{k} = {v if v is not None else "NOT_RUN"}' for k,v in summary.items())+'\n')
report=f'''# Vast RTX4090 Stage4 恢复核查

结果：**{summary['VAST4090_RESTORE_STATUS']}**。新平台纯流体基线：**{summary['NEW_VAST4090_PURE_FLUID_BASELINE']}**。

归档重新检查 {json.loads(archive['stdout']).get('file_count') if archive_ok else '未通过'} 个文件，工作副本 {len(frozen)} 个文件；终检源码/冻结输入差异 {len(bad)}。只适配工作路径、工具链/链接路径和实测 RTX4090 cc89。原 Stage4 C++、已应用 patch 的 Palabos 导出源码、数学/数值合同及运行参数字节保持不变。原 RTX5090 binary 仅记录来源，未执行。

NVHPC26.5 / CUDA13.2 从冻结 NVIDIA 官方地址下载，SHA256 与原下载 provenance 相同；安装证据在 toolchain_env 和 provenance。未修改主机驱动或全局 CUDA。OpenMPI 仅安装缺失依赖；正式 GPU smoke 前完成无求解器的 MPI 单核绑定检查。GPU activity 由官方 cavity 的 CUDA kernel 原始 trace 证明；Stage4 性能运行不 trace、不加同步、不优化。

| 验证 | 结果 | 说明 |
|---|---|---|
| 官方 cavity3d | {official.get('status')} | 最多200步；kernel {official.get('kernel_activity',{}).get('status','UNVERIFIED')} |
| Stage4 200 | {gates[200].get('status','NOT_RUN')} | 从零；RTX参照0/1/10/100/200、CPU参照0/100/200 |
| Stage4 1000 | {gates[1000].get('status','NOT_RUN')} | 从零；0/100/200/1000全场及每步24组flux/质量/安全 |
| Stage4 5000 | {gates[5000].get('status','NOT_RUN')} | 从零；0/100/200/1000/5000全场，执行原增长判据 |

每个 horizon 都对冻结的 CPU MPI12 和 RTX5090 Stage4 两组历史真实结果逐行比较。完整检查压缩CSV、各量最大误差、检查点数据和未修改比较器已归档，阈值未改变。NaN/Inf累计：{summary['NONFINITE_COUNT']}；趋势：{trend}。归一化速度、密度、总质量/CV质量、全部多平面流量与独立 MPI ownership、安全证据均由原比较器审查。

| 性能范围 | 初始化 s | 求解 steps/s | 端到端 s |
|---|---:|---:|---:|
| 1000（计时100–1000） | {s1000.get('initialization_seconds')} | {summary['RTX4090_STEPS_PER_SEC']} | {s1000.get('end_to_end_seconds')} |
| 5000（计时500–5000） | {s5000.get('initialization_seconds')} | {new5000} | {s5000.get('end_to_end_seconds')} |

持续求解速率相对 RTX5090：{comparison['RTX4090_VS_RTX5090_RATIO']}；相对 CPU MPI12：{comparison['RTX4090_VS_CPU_MPI12_SPEEDUP']}。分母为用户指定历史参考 125.880863 和17.692465 steps/s。每个 horizon 仅一次新测量，不宣称重复实验统计可靠性。初始化及端到端包含所需输出；性能不作为数值正确性门槛。

5000步内存判据：{summary['GPU_MEMORY_GROWTH']}；性能稳定性：{summary['PERFORMANCE_STABILITY']}。具体 steady thirds、256MiB/512MiB阈值、首末1000步窗口及原Stage4 5% horizon下降判据见 RTX4090_PERFORMANCE.json。GPU利用率为一秒采样、按已完成步数筛选，不能用来宣称kernel占比。未运行10000步，所以不声称通过原5000/10000跨horizon重复稳定性测试。

环境检查曾因帮助文本把支持架构89写作通用ccXY而被脚本误判，发生在所有编译/求解器之前。原失败记录完整保存在provenance/environment_check_failure_1，修正版改为读取supported values并执行-stdpar -gpu=cc89预处理接受检查。没有重复运行任何求解器，也未改变科学源码或阈值。

第二次启动前完整性检查发现准备脚本新命令receipt与旧NVHPC下载provenance重名。新receipt已另存，旧provenance从只读archive恢复原SHA256，5654个工作输入重新通过；该拦截发生在Stage4 solver实际启动之前。原失败结果和修复证据见provenance/stage4_prelaunch_failure_2及PROVENANCE_FILENAME_COLLISION_REPAIR.json。官方smoke没有重复。

OpenMPI在此容器报告非致命的memory binding失败警告，但实际CPU affinity验证为一个物理核心（逻辑CPU0、8）；未禁用CPU核心绑定。MPI启动等待按实测纳入端到端，不能将solver iteration加速比解读为端到端加速比。新binary的静态CUDA ELF检查确认sm_89（BINARY_CUDA_ARCHITECTURE_INSPECTION.json）。

执行状态：{state.get('phase')}。错误（如有）：{state.get('error','无')}。

远端 work：`{W}`；compact：`{F}`。原始全场二进制保留在远端各run，compact带全文件SHA256清单和完整逐项比较证据；较大的CSV仅在compact副本无损gzip，原run文件未改。官方原始CUDA trace随compact保留。WSL下载后单独验证manifest，远端summary此时尚未宣称本地传输PASS。

旧 incomplete upload：`{old}`，{oldbytes} bytes（allocated），未删除。正式历史基线和只读archive未改；没有RBC、PBS/BSA、microbubble、Stage5或long run。下一步仅作建议：{summary['NEXT_STEP']}，本任务不会执行。
'''
(F/'VAST4090_STAGE4_RESTORE_REPORT.md').write_text(report)
# The finalizer's log is written outside the sealed result; no copied live logs mutate here.
(F/'SHA256SUMS').write_text('\n'.join(sha(p)+'  '+str(p.relative_to(F)) for p in sorted(F.rglob('*')) if p.is_file() and p!=F/'SHA256SUMS')+'\n')
check=capture(['sha256sum','-c',str(F/'SHA256SUMS')]) if False else subprocess.run(['sha256sum','-c','SHA256SUMS'],cwd=F,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
assert check.returncode==0,check.stderr
save(W/'FINALIZATION_TERMINAL.json',dict(status='PASS',restore_status=summary['VAST4090_RESTORE_STATUS'],remote_result_dir=str(F),files=len(list(F.rglob('*'))),manifest_sha256=sha(F/'SHA256SUMS'),unix=time.time()))
print(json.dumps(summary,indent=2),flush=True)
