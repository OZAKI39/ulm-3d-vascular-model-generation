"""Final Stage O report from the observed windows and sole full steady run."""
import json,re,sys,xml.etree.ElementTree as ET
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from sv_validation.provenance import sha256,now
R=ROOT/'reports/sv1_3o'
def read(n):return json.loads((R/(n+'.json')).read_text())
def write(n,d):(R/(n+'.json')).write_text(json.dumps(d,indent=2,ensure_ascii=False)+'\n')
c=read('optimized_steady_candidate');w=read('winner');b=read('baseline_runtime');ref=read('reference_freeze');prof=read('profile_summary')
for n in ('optimized_steady_candidate','preservation_audit','remote_preservation','native_artifact_mirror','visuals'):assert read(n)['status']=='PASS'
cases=ET.parse(R/'pytest_full.xml').findall('.//testcase')
def counts(rows):return dict(passed=sum(not any(t.find(k) is not None for k in ('failure','error','skipped')) for t in rows),failed=sum(t.find('failure') is not None for t in rows),errors=sum(t.find('error') is not None for t in rows),skipped=sum(t.find('skipped') is not None for t in rows))
own=[t for t in cases if 'test_sv13o_' in t.get('classname','')];old=[t for t in cases if t not in own]
tests=dict(stage=counts(own),full=counts(cases),historical=counts(old),failures=[dict(test=t.get('classname')+'::'+t.get('name'),historical=t in old,message=t.find('failure').get('message')) for t in cases if t.find('failure') is not None])
prior=json.loads((ROOT/'reports/sv1_3n/test_summary.json').read_text())
tests['historical_result_matches_frozen_N']=tests['historical']==prior['full'] and {f['test'] for f in tests['failures'] if f['historical']}=={f['test'] for f in prior['failures']}
write('test_summary',tests)
assert tests['stage']['failed']==tests['stage']['errors']==0
assert tests['historical_result_matches_frozen_N'], 'Investigate changed historical test results; do not rewrite historical evidence'
status='GPU_PERFORMANCE_OPTIMIZED' if w['fixed_window_reduction']>=.10 else 'SIMPLE_GPU_TUNING_EXHAUSTED'
state=dict(completed=True,status=status,outcome='PASS',classification=c['classification'],scientific_equivalence='DEFERRED',CPU_production='CPU_EARLY_STOP_PRODUCTION',production_changed=False,supported_MPI_ranks=1,full_steady_runs=1,advanced_PC_next_stage='Suggested only; NOT STARTED' if status=='SIMPLE_GPU_TUNING_EXHAUSTED' else 'DEFERRED',timestamp=now())
runs={k:read('PERF_'+k+'_acceptance') for k in 'ABCD'}
labels={'A':'ASM overlap2 + ILU(2), restart100','B':'direct ILU(2), restart100','C':'direct ILU(1), restart100','D':'winner PC, restart200'}
tab=[]
for k,d in runs.items():
 complete=d['status']=='PASS';s=d.get('statistics',{})
 tab.append(f"| {k} | {labels[k]} | {d['status']} | "+(f"{d['wall_time_s']:.6f}" if complete else '不参与完整窗口排名')+' | '+str(s.get('total_iterations','—'))+' | '+(f"{s['mean_iterations']:.3f}" if 'mean_iterations' in s else '—')+' |')
table='\n'.join(tab)
cost_table='\n'.join(f"| {k} | {d['statistics']['KSP_solves']} | {d['statistics']['max_iterations']} | {d['statistics']['GMRES_cycles']} / {d['statistics']['GMRES_restarts']} | {d['sampled_device_memory_max_MiB']:g} |" for k,d in runs.items())
healthy=[(k,d) for k,d in runs.items() if d['status']=='PASS'];rawfast=min(healthy,key=lambda x:x[1]['wall_time_s'])
failures={}
for k,d in runs.items():
 if d['status']=='FAIL':
  e=read('remote/PERF_'+k+'_execution');failures[k]=dict(status='FAIL',stop_reason=e['monitor_stop'],timesteps_seen=d['steps'],last_step_seen=d['stop_step'],wall_before_abort_s=d['wall_time_s'],KSP_reason_values=d['KSP_reason_values'],source='PERF_'+k+'_acceptance.json',excluded_from_ranking=True)
write('candidate_failures',failures)
profile_names=list(dict.fromkeys((prof['runs']['baseline'],prof['runs']['winner'])))
pt=[]
for name in profile_names:
 events=prof['profiles'][name]['events'];pt.append('| '+name+' | '+' | '.join(f"{events.get(k,{}).get('time_s',0):.6g}" for k in ('PCSetUpOnBlocks','MatLUFactorNum','KSPSolve','PCApply','MatMult','PCSetUp'))+' |')
profile_table='\n'.join(pt)
pe=prof['profiles'][prof['runs']['baseline']]['events'];ps=prof['profiles'][prof['runs']['baseline']]['stages']
selection='\n'.join(f"- {d['candidate']}: {d['action']}；"+(f"相对当前配置减时 {d['reduction']*100:.3f}%" if d['reduction'] is not None else '未满足求解健康门槛')+'。' for d in w['decisions'])
io=read('OFFICIAL_OUTPUT_STOP_acceptance');stats=c['statistics'];m=c['measurement'];final_run=read('REAL_VASCULAR_GPU_PERF_acceptance');build=read('svmp_build')
recommend='简单调参收益不足 10%，本轮已停止；后续如继续性能研究，建议单独开展 GPU-native preconditioner study。本次没有进入 Stage SV1.3P。' if status=='SIMPLE_GPU_TUNING_EXHAUSTED' else '所选配置相对固定窗口 A 明显减时，唯一一次完整稳态验证通过；本阶段完成后停止。'
report=f'''# Stage SV1.3O — RTX 4090 单 GPU 血流求解性能优化

**PASS: {status}。保留配置 {w['candidate']}；结果为 {c['classification']}，NOT YET SCIENTIFICALLY VALIDATED。**

正式 scientific production 仍是 `CPU_EARLY_STOP_PRODUCTION`。全部 CFD 均为 1 MPI rank、1 RTX4090、OMP=1；没有 CPU 对照或多 rank CUDA 测试。

## 优化前是什么状态？

Stage N 实际完成 {b['steps']} 步，wall={b['wall_time_s']:.6f} s，{b['KSP_solves']} 次 KSP solve、{b['total_iterations']} 次 iteration，平均 {b['mean_iterations']:.6f}。首次完整稳态判据为 step{b['first_full_steady_step']}，安全停止为 step{b['stop_step']}。

这些记录为 **DEVELOPMENT_BASELINE_OBSERVED**，不是重复 benchmark median。PETSc `{ref['PETSc']['tag']}` / commit `{ref['PETSc']['commit']}`、CUDA13.2、OpenMPI4.1.6、Stage N lifecycle adapter、mesh/geometry、rho/mu/Q/BC/dt、时间积分、非线性容差及 production steady policy 均冻结。[冻结清单](reference_freeze.json)、[历史审计](preservation_audit.json)。

## 首先减少了什么无效开销？

Stage N 每步写 VTU。新的 **PERFORMANCE_OUTPUT_POLICY** 正常每 10 步写正式 VTU，继续每 10 步运行原 steady monitor；收到 STOP_SIM 后完成当前步，同时写 final VTU 和完整 native checkpoint。

源码唯一新增差异是 `main.cpp` 的输出条件：`save_vtu || (com_mod.saveVTK && reached_stop_time_step)`。原 restart 保存逻辑已经包含停止条件，未改变其二进制格式；方程、组装、积分和 lifecycle adapter 均不变。[最小补丁](../../patches/sv1_3o/final_output_on_stop.patch)、[逐文件来源](source_patch.json)、[GPU clean build](svmp_build.json)。构建沿用已冻结的 PETSc/MPI 前缀，没有重新构建 PETSc。远端原源码包缺失的 866 个受跟踪文件由 WSL 原始副本补齐，已有文件哈希全部一致。

短官方测试实际完成 {io['steps']} 步，仅保存 {io['VTU_steps']}，step11 在正常 cadence 之外仍产生完整 VTU/checkpoint；fresh reload PASS。这一机制验证没有运行额外完整 vascular steady。[官方输出测试](OFFICIAL_OUTPUT_STOP_acceptance.json)。最终完整运行保存 {c['VTU_count']} 份 VTU，步号 {c['VTU_steps']}，Stage N 为 71 份。

固定窗口从 Stage N 原生 **step60 checkpoint** 继续至 step70，实际执行61—70共10步；初始文件 SHA256 `{w['checkpoint_sha256']}`。Y_n/A_n/step/time/history 完整，PERF_A 使用新二进制的原生读取器成功加载，所有候选使用同一 SHA、XML、mesh、GPU backend 和 rank 数。继承日志的累计 elapsed 不作窗口计时，统一使用新进程单调时钟。[checkpoint 验证](checkpoint_window.json)。

## ASM 在单 rank 下有没有必要？

实际 PERF_A 的每次 PCView 均报告 **1 个 ASM block、overlap2**；子域与全局矩阵均为 **281452 × 281452**，单 rank local/global rows 相同。[实际 topology](asm_topology.json)。

因此执行了 direct ILU(2)。B 的 observed wall={runs['B']['wall_time_s']:.6f} s，相比 A 减时 {(1-runs['B']['wall_time_s']/runs['A']['wall_time_s'])*100:.3f}%。这个观测不能证明所有单 rank 场景都不需要 ASM；本次固定窗口是否保留 B 由预先冻结的收益门槛决定。

## 哪个线性配置最快？

| Candidate | 配置 | 实际状态 | 完整10步 wall / s | KSP iterations | Mean / solve |
|---|---|---|---:|---:|---:|
{table}

失败候选的 iteration 数只覆盖中止前的部分窗口，不用于性能比较。C 在 step61、D 在 step62 出现 `DIVERGED_BREAKDOWN`（reason=-5），因此均拒绝。C/D 原始失败、退出码、负 KSP reason 和日志完整保留于 [失败明细](candidate_failures.json)，未通过放松 tolerance、改变 dt 或增加参数扫描挽救。

| Candidate | KSP solves | Max iterations / solve | GMRES cycles / restarts | 整卡显存采样最大值 / MiB |
|---|---:|---:|---:|---:|
{cost_table}

GMRES cycles/restarts 由实际 iteration 与固定 restart 长度推算，未增加 trace 运行。显存每 15 秒低成本采样，不能视为精确进程峰值；C 约 10 秒即失败，1 MiB 启动采样不代表其真实显存需求。A–D 的 PCView 实际配置均已独立核对，包括 C 的 `1 level of fill`；该配置核对不会改变 C/D 的求解失败状态。[配置与成本明细](candidate_runtime_summary.json)。

完整且健康的单次原始计时最小值是 **{rawfast[0]}：{rawfast[1]['wall_time_s']:.6f} s**。按用户规则保留的 winner 是 **{w['candidate']}**，选优计时 {w['selection_wall_time_s']:.6f} s，相对 A 减时 {w['fixed_window_reduction']*100:.3f}%。低于 5% 视为无明显收益；只有 5%—10% 区间允许一次确认，不默认 median。本次没有因为失败候选的短暂运行时间而把它判为更快。

{selection}

![候选窗口](gpu_tuning_candidates.png)

## 时间改善来自哪里？

必须区分总 iteration、每次 iteration 的成本和其它开销。下图的 wall/iteration 包含启动、组装、PC、输出，不能解释为纯 GPU kernel 的单次时间。winner 的选择依据完整10步 wall 与健康门槛。

B 的总迭代从 A 的 {runs['A']['statistics']['total_iterations']} 增至 {runs['B']['statistics']['total_iterations']}，wall/iteration 从 {1000*runs['A']['wall_time_s']/runs['A']['statistics']['total_iterations']:.3f} 降至 {1000*runs['B']['wall_time_s']/runs['B']['statistics']['total_iterations']:.3f} ms；最终 wall 差异只有 {(1-runs['B']['wall_time_s']/runs['A']['wall_time_s'])*100:.3f}%，不能据单次小差异认定去掉 ASM 有稳定收益。没有为 B/C/D 追加 PC/KSP 计时 profile，相关分项标为 NOT_MEASURED；轻量诊断仅覆盖允许的 baseline/winner。

![迭代与成本](ksp_cost_comparison.png)

仅在 baseline 与最终 winner 的固定窗口进行了轻量 PETSc `log_view` / GPU event timing。两者同配置时共用一次 profile；实际 profile 次数 **{prof['runs']['runs']}**。这些运行全部排除在选优计时之外。

| Profile | PCSetUpOnBlocks / s | MatLUFactorNum / s | KSPSolve / s | PCApply / s | MatMult / s | PCSetUp / s |
|---|---:|---:|---:|---:|---:|---:|
{profile_table}

主要记录成本是 **PCSetUpOnBlocks={pe['PCSetUpOnBlocks']['time_s']:.3f} s**，其中 MatLUFactorNum={pe['MatLUFactorNum']['time_s']:.3f} s，MatILUFactorSym={pe['MatILUFactorSym']['time_s']:.4f} s；顶层 PCSetUp 不能代表全部预条件器准备成本；ASM 的块内预条件器 setup 另由 PCSetUpOnBlocks 完成。[官方说明](https://petsc.org/release/manualpages/PC/PCSetUpOnBlocks/)。KSPSolve={pe['KSPSolve']['time_s']:.3f} s 内 PCApply={pe['PCApply']['time_s']:.3f} s 是主要事件，MatMult={pe['MatMult']['time_s']:.3f} s 明显较小。应用标记的 PETSc Solve stage 共 {ps['PETSc Solve']['time_s']:.3f} s（PETSc 总计时的 {ps['PETSc Solve']['percent_total']:.1f}%）。因此现有证据支持以后优先研究 PC 的准备与应用成本，而不是继续扫描 restart。

PETSc 事件有嵌套，表中数值不能相加。同名事件按本次单 rank 日志的不同 stage 汇总，并保留逐 stage 原行。MatAssemblyBegin/End 共 {pe['MatAssemblyBegin']['time_s']+pe['MatAssemblyEnd']['time_s']:.6g} s，仅为 PETSc 矩阵最终装配，不代表完整 FEM 方程组装。FEM assembly 与 VTU/restart I/O 的独立时间均为 **NOT_SEPARATELY_MEASURED**；I/O 的已证实变化是输出文件数量减少，不能把全部 wall 差异归因于 I/O。[profile 数据](profile_summary.json)、[官方 profiling 说明](https://petsc.org/release/manual/profiling/)。

实际 Mat=`{final_run['Mat']}`、Vec=`{final_run['Vec']}`，factorization package={final_run['factor_packages']}；显式 fallback warning 记录为 {final_run['CPU_fallback_warning_lines']}。这不证明整条 PC 路径都在 GPU 上；官方 PCILU 文档对 cuSPARSE 矩阵的分解位置有 CPU 说明，端到端 residency 本阶段标为 **NOT_MEASURED**，没有追加 residency CFD。[官方 PCILU](https://petsc.org/release/manualpages/PC/PCILU/)。

![轻量时间记录](gpu_time_breakdown.png)

## 最终优化配置是什么？

最终采用 {w['candidate']}。完整 PETSc options：

```text
{w['PETSC_OPTIONS']}
```

1 MPI rank、1 RTX4090、OMP=1；VTU cadence10，并在停止时强制保存最终 VTU/restart。rtol1e-10、atol1e-24、max_it2000、right preconditioning、非线性及稳态容差保持不变。[winner](winner.json)、[完整运行 XML](../../outputs/sv1_3o/REAL_VASCULAR_GPU_PERF/solver.xml)。

## 优化后真实 vascular 到第几步稳态？

唯一一次 `REAL_VASCULAR_GPU_PERF` 从 t=0 运行。首次完整满足 production 连续 **5** 个区间是在 **step{c['first_full_steady_step']}**，安全停止 **step{c['stop_step']}**；物理时间 {c['physical_time_s']:.12g} s，wall **{c['wall_time_s']:.6f} s**。

{stats['KSP_solves']} 次 KSP solve，{stats['total_iterations']} 次 iteration，mean={stats['mean_iterations']:.6f}，max={stats['max_iterations']}；linear/nonlinear failures=0。最后正式区间 E_u={c['last_interval']['E_u']:.10g}、E_Q={c['last_interval']['E_Q']:.10g}，final mass={m['epsilon_mass']:.10g}。低成本整卡显存采样最大观测值 {c['sampled_device_memory_max_MiB']:g} MiB；不是精确进程峰值。

原 `SteadyStopMonitor` / `stop_gate` 直接复用，阈值、连续区间数、cadence 未变。第一次通过即请求 STOP_SIM，仅完成当前步并保存，没有额外余量时间步。[steady history](steady_history.json)、[停止请求](steady_stop_request.json)。

![稳态历史](optimized_steady_convergence.png)

最终 VTU SHA256 `{c['VTU_sha256']}`；checkpoint SHA256 `{c['checkpoint']['sha256']}`；solver SHA256 `{c['solver_sha256']}`。VTU 内 TimeValue 与 checkpoint 时间均为 {c['reload']['time_s']:.17g} s；fresh-process reload、有限速度/压力/流量、wall no-slip 和 mass 检查全部通过。[输出时间与重载](output_consistency.json)、[冻结 GPU 候选](optimized_steady_candidate.json)、[原生构建镜像](native_artifact_mirror.json)。

## 实际开发加速多少？

Stage N observed wall={b['wall_time_s']:.6f} s；Stage O observed wall={c['wall_time_s']:.6f} s；S_dev={c['development_speedup']:.6f}×。

**OBSERVATIONAL DEVELOPMENT SPEEDUP — SINGLE-RUN DEVELOPMENT COMPARISON。** 两者均达到相同 production steady criteria，但各只运行一次；这不是正式 benchmark，也不是 scientific production comparison。完整运行改善与固定窗口调参改善是不同记录，阶段分类按相对 low-I/O PERF_A 的调参收益判断。

![单次开发比较](gpu_runtime_before_after.png)

## 哪些工作仍然推迟？

- CPU/GPU science equivalence — **DEFERRED**。
- OLD/NEW PETSc science comparison 与 CPU benchmark — **DEFERRED**。
- Formal benchmark — **DEFERRED**。
- Multi-GPU / multi-rank CUDA 与 mpicuda reverse — **DEFERRED**。
- Advanced GPU-native PC — **DEFERRED**。

{recommend}

阶段结论 **{status}**；科学等价仍为 **DEFERRED**，正式 production 保持 `CPU_EARLY_STOP_PRODUCTION`。Stage O tests={tests['stage']}；完整 pytest={tests['full']}；历史结果={tests['historical']}。新 tests 读取已有 artifact，不启动 CFD；对 C/D 的拒绝测试不会把真实求解失败改写为 PASS。[测试记录](test_summary.json)、[本地保留审计](preservation_audit.json)、[原生保留审计](remote_preservation.json)。
'''
(R/'REPORT.md').write_text(report);write('report_manifest',dict(status='PASS',sha256=sha256(R/'REPORT.md'),timestamp=now(),stage_status=status));write('stage_result',state);print('Final report written: '+status)
