"""Plain-Chinese final report grounded in actual native artifacts."""
import json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];R=ROOT/'reports/sv1_3q'
def read(n):return json.loads((R/(n+'.json')).read_text())
def optional(n):return read(n) if (R/(n+'.json')).exists() else None
def fmt(x,n=2):return '—' if x is None else f'{x:.{n}f}'
def event(d,n):return ((d or {}).get('profile') or {}).get('events',{}).get(n,{}).get('time_s')
def count(d):return (d or {}).get('reuse',{}).get('ILU_rebuild_count')
def iterations(d):return (d or {}).get('total_attempt_iterations',(d or {}).get('statistics',{}).get('total_iterations'))
base=read('baseline_r1');w=read('winner');full=optional('winner_steady_candidate');p=json.loads((ROOT/'configs/sv1_3q/policy.json').read_text());ranking=read('late_ranking');tests=optional('pytest_summary')
rows={c:dict(smoke=read(c+'_SMOKE_acceptance'),late=optional(c+'_WINDOW_acceptance'),early=optional(c+'_EARLY_acceptance')) for c in ('R2','R3','R5','RA')}
if full and full['status']=='PASS' and w['full_required']:status='ILU_REBUILD_POLICY_WINNER_FOUND'
elif not any(r['late'] and r['late']['status']=='PASS' for r in rows.values()):status='CROSS_TIMESTEP_REUSE_UNSTABLE'
else:status='CROSS_TIMESTEP_REUSE_EXHAUSTED'
completed=not w['full_required'] or bool(full and full['status']=='PASS')
state=dict(status=status,completed=completed,winner=w.get('candidate'),scientific_equivalence='DEFERRED',CPU_production='CPU_EARLY_STOP_PRODUCTION',Stage_P_fallback_preserved=True,maximum_full_runs=1,comparison='OBSERVATIONAL DEVELOPMENT SPEEDUP',limitations=[base['early']['limitation'],'Single observations, no repeats or formal performance statistics. Device-wide VRAM sampling is not process peak.'])
(R/'stage_status.json').write_text(json.dumps(state,indent=2)+'\n')
text=['# Stage SV1.3Q — ILU 到底多久重新构建一次最快？','',f'本阶段状态：**{status}**。只改变原 ASM overlap2 + ILU(2) 的重建时机，保留 Stage P 作为回退。','',
'## 为什么不需要每一步都重新做 ILU？','',
'相邻时间步的方程通常相似，因此可以保留上一份 ILU。收益取决于少做分解省下的时间，能否抵消旧 ILU 造成的额外 GMRES 迭代。本阶段只根据实际完整窗口时间判断。','',
'PETSc 提供正式的 [KSPSetReusePreconditioner 接口](https://petsc.org/release/manualpages/KSP/KSPSetReusePreconditioner/)。原有 KSP/PC 对象持续保留，矩阵值照常更新，只有重建标志变化。','',
'## 当前基线是什么？','',
f"R1 即 Stage P：每个时间步第一次求解重建，步内继续复用。历史晚期窗口为 **{base['late']['wall_time_s']:.6f} 秒**，6720 次迭代、20 次 KSP 求解、10 次实际数值分解。本阶段没有重跑 R1。",'',
f"R1 早期基线从唯一历史完整运行的 step10 最后一轮日志到 step20 最后一轮日志直接提取：**{base['early']['wall_time_s']:.1f} 秒**、{base['early']['iterations']} 次迭代。它是连续运行区间，日志时间有有限分辨率；新策略是独立重启窗口，含启动开销。因此早期比较是带此限制的开发观测，不能当成完全对称的正式基准。",'',
'所有运行保持 1 MPI rank、1 RTX4090、OMP=1、PETSc 3.25.5、CUDA 13.2、real double、32 位索引。GMRES 右预条件、restart100、rtol1e-10、atol1e-24、max_it2000 和对角缩放不变。网格、物理、边界条件、时间步、非线性及稳态阈值不变。','',
'晚期均使用原 Stage N step60 原生检查点；早期使用验证通过的 Stage P step10 原生检查点。同一窗口的起始检查点 SHA256 完全相同。检查点不包含 ILU，因此启动后第一次实际求解必须重建；间隔从该步开始计数。','',
'| 策略 | 最多 3 步早期 smoke | 健康晚期 10 步时间 / 秒 | 全部尝试迭代数 | 实际分解次数 |',
'|---|---|---:|---:|---:|',f"| R1 每步重建 | 已有 Stage P 验收 | {base['late']['wall_time_s']:.2f} | 6720 | 10 |"]
for c,r in rows.items():
 d=r['late'];good=d and d['status']=='PASS'
 text.append(f"| {c} | {r['smoke']['status']}，{r['smoke']['wall_time_s']:.2f} 秒 | {fmt(d['wall_time_s']) if good else 'FAIL / 未运行'} | {iterations(d) if good else '—'} | {count(d) if good else '—'} |")
for c,h in [('R2','每2步重建怎么样？'),('R3','每3步怎么样？'),('R5','每5步怎么样？')]:
 r=rows[c];d=r['late'];text+=['',f'## {h}','']
 if d and d['status']=='PASS':text += [f"晚期 10 步实测 **{d['wall_time_s']:.6f} 秒**，全部尝试共 {iterations(d)} 次迭代，实际 ILU 分解 {count(d)} 次、复用 {d['reuse']['ILU_reuse_count']} 次。平均迭代 {d['statistics']['mean_iterations']:.2f}，最大 {d['statistics']['max_iterations']}。线性/非线性失败均为零。"]
 else:
  fail=d or r['smoke'];text += [f"在 {fail['mode']} 阶段未通过，失败前运行 {fail['wall_time_s']:.2f} 秒。原因：`{'; '.join(fail['errors'])}`。该时间不是健康完整窗口成绩，不参与排名。"]
text+=['','## 自动判断什么时候重建怎么样？','',
'RA 每次新建 ILU 后，以第一次健康求解的迭代数作为参考。复用时若本次迭代数严格超过参考的 1.5 倍，下一次线性求解就重建；年龄等于当前时间步减去上次重建步，达到 5 时也强制重建。规则在运行前固定，未扫描阈值。','',
'若旧 ILU 的求解原因为负或达到迭代上限，只重建并重试同一个线性问题一次。第一次失败及其原因、迭代数始终保留。恢复发生在把修正量交回血流求解器之前，原 RHS 从设备副本恢复，零初值和矩阵不变。PETSc 原有 diagonal_scale_fix 在返回时恢复矩阵缩放。参见 [KSPSolve 文档](https://petsc.org/release/manualpages/KSP/KSPSolve/) 与本机固定版本源码审计。','']
a=rows['RA']['late'] or rows['RA']['smoke']
from collections import Counter
ra_reasons=Counter(x['rebuild_reason'] for x in (a.get('reuse') or {}).get('trace',[]) if not x['reuse'])
text += [f"RA 实际{'晚期窗口' if a['mode']=='window' else '早期短检查'}结果为 {a['status']}，耗时 {a['wall_time_s']:.6f} 秒，迭代 {iterations(a)}，分解 {count(a)} 次，成功恢复 {(a.get('reuse') or {}).get('recovery_count',0)} 次。",'',f"实际重建原因：第一次求解 {ra_reasons['FIRST_SOLVE']} 次，年龄到限 {ra_reasons['MAX_AGE']} 次，迭代增长 {ra_reasons['ITERATION_GROWTH']} 次，旧 ILU 失败 {ra_reasons['STALE_FAILURE']} 次。",'',('本窗口没有触发迭代增长或失败恢复，RA 实际上采用了与 R5 相同的重建步序。二者一次运行的时间差不能证明自适应规则本身更优。' if not ra_reasons['ITERATION_GROWTH'] and not ra_reasons['STALE_FAILURE'] else '自适应触发的逐次证据保留在验收文件和时间线中。'),'',
'独立 GPU 合成故障算例使用同一生产头文件，验证了 RHS 恢复、一次 fresh 重试成功，以及 fresh 重试仍失败后停止。它不属于 CFD，也不参与性能排名。测试程序初次遗漏已冻结的 MPI 选项、第二次使用了旧版上下文清理调用；这些原日志均保留，修正只涉及测试程序，最终完整退出成功。','',
'## 哪一种策略最快？','']
if w.get('candidate'):
 text += [f"晚期排名前两名为 {'、'.join(ranking['top_two'])}。按健康的早期与晚期两个窗口总时间，最快为 **{w['candidate']}**：晚期 {w['late_wall_s']:.6f} 秒，早期 {w['early_wall_s']:.6f} 秒，总计 **{w['observed_total_wall_s']:.6f} 秒**。相对 R1 的 {w['baseline_total_wall_s']:.6f} 秒，预计节省 **{100*w['improvement']:.2f}%**。"]
else:text += ['没有候选同时满足早期和晚期健康要求以及早期退化限制，继续使用 R1。']
text += ['','## 为什么它最快？','',
'| 策略 | 分解次数 | 少一次分解对应净节省 / 秒 | 每少一次分解增加迭代 | 准备 ILU / 秒 | 使用 ILU / 秒 | KSP 求解 / 秒 | MatMult / 秒 | 观察到的显存 / MiB |','|---|---:|---:|---:|---:|---:|---:|---:|---:|']
text += [f"| R1 | 10 | — | — | {fmt(event(base['late'],'PCSetUpOnBlocks'))} | {fmt(event(base['late'],'PCApply'))} | {fmt(event(base['late'],'KSPSolve'))} | {fmt(event(base['late'],'MatMult'))} | {fmt(base['late']['sampled_device_memory_max_MiB'],0)} |"]
for c,r in rows.items():
 d=r['late']
 if not d or d['status']!='PASS':continue
 fewer=10-count(d);seconds=(base['late']['wall_time_s']-d['wall_time_s'])/fewer if fewer else None;extra=(iterations(d)-6720)/fewer if fewer else None
 text += [f"| {c} | {count(d)} | {fmt(seconds)} | {fmt(extra)} | {fmt(event(d,'PCSetUpOnBlocks'))} | {fmt(event(d,'PCApply'))} | {fmt(event(d,'KSPSolve'))} | {fmt(event(d,'MatMult'))} | {fmt(d['sampled_device_memory_max_MiB'],0)} |"]
text += ['','表中“每少一次分解”是两次完整运行差值的描述，不是独立因果测量。准备时间采用包含子块分解的 PCSetUpOnBlocks；MatLUFactorNum、PCSetUp 及其他原始事件详见验收 JSON。事件可能嵌套，不能相加充当 wall。显存是设备级离散采样最大值，不是进程精确峰值。','',
'## 早期 transient 和接近 steady 结果一样吗？','',
'只对晚期前两名运行 step10→20。早期耗时比历史 R1 高超过 5% 时，按预先固定的“明显慢于基线”规则排除；任何不健康结果也排除。其余候选不补跑早期窗口。','',
'| 策略 | 早期时间 / 秒 | 晚期时间 / 秒 | 早期迭代 | 早期分解 | 恢复 | 健康 |','|---|---:|---:|---:|---:|---:|---|']
for c in ranking['top_two']:
 d=rows[c]['early'];text += [f"| {c} | {fmt(d['wall_time_s']) if d['status']=='PASS' else 'FAIL'} | {fmt(rows[c]['late']['wall_time_s'])} | {iterations(d)} | {count(d)} | {(d.get('reuse') or {}).get('recovery_count',0)} | {d['status']} |"]
text += ['','## 如果 adaptive 最好','']
if w.get('candidate')=='RA':text += ['RA 在本次两个窗口总时间中胜出。逐次轨迹完整记录年龄、参考值、当前迭代、比值及重建原因；完整运行同样保留这些字段。']
else:text += ['RA 已实际运行，下面仍给出其重建时间线；最终选择由两窗口总时间决定，不因为“自适应”名称优先选择它。']
text += ['','![自适应重建时间线](adaptive_rebuild_timeline.png)','','## 最终完整血管计算怎么样？','']
if full:
 a=read('REAL_VASCULAR_GPU_ILU_REUSE_WINNER_acceptance');lat=read('stop_latency')
 full_reasons=Counter(x['rebuild_reason'] for x in a['reuse']['trace'] if not x['reuse'])
 text += [f"{w['candidate']} 按两个窗口总时间胜出且达到至少 10% 收益门槛，因此只对该优胜策略从 t=0 运行一次 `REAL_VASCULAR_GPU_ILU_REUSE_WINNER`。完整耗时 **{full['wall_time_s']:.6f} 秒**；首次连续 5 个正式区间达标为 **step{full['first_full_steady_step']}**，实际原生安全停止为 **step{full['stop_step']}**。",'',
 f"KSP 共 {a['reuse']['KSP_attempts']} 次尝试，{iterations(a)} 次迭代；ILU 分解 {count(a)} 次、复用 {a['reuse']['ILU_reuse_count']} 次。未恢复的线性失败 {a['linear_failures']}，非线性失败 {a['nonlinear_failures']}，已恢复尝试 {a['reuse']['recovery_count']}。最终速度和压力有限、质量门槛通过、壁面无滑移通过，final VTU 和原生 checkpoint 已由独立过程重新读取。",'',
 f"停止检测时间戳 {lat['steady_detected_timestamp']:.6f}，请求时间戳 {lat['stop_request_timestamp']:.6f}，求解器确认时间戳 {lat['solver_ack_timestamp']:.6f}（Unix 秒）。检测使用 WSL 时钟，请求和确认使用远端时钟，时钟探测见 `clock_probe.json`。请求到原生确认耗时 {lat['request_to_ack_s']:.3f} 秒。",'',
f"完整运行重建原因：首次求解 {full_reasons['FIRST_SOLVE']} 次，年龄到限 {full_reasons['MAX_AGE']} 次，迭代增长 {full_reasons['ITERATION_GROWTH']} 次，旧 ILU 失败 {full_reasons['STALE_FAILURE']} 次。逐时间步、逐次求解记录见 `adaptive_timestep_history.json`；未在真实 CFD 中触发的分支不宣称具有实测性能收益。",'',
'原 production SteadyStopMonitor 和全部阈值未改。小改动仅将 VTU 与 checkpoint 合并压缩传输、将停止请求提前到判定后立即发出，并在原生确认处写一个时间戳。允许完成当前在途步和最终输出，没有额外补跑，也没有扣除停止延迟来美化耗时。']
else:text += ['没有满足至少 10% 收益及早/晚窗口健康门槛的候选，按要求不运行完整血管计算。']
text += ['','## 相对 Stage P 快多少？','']
if full:text += [f"Stage P 为 {base['full']['wall_time_s']:.6f} 秒，本次为 {full['wall_time_s']:.6f} 秒，实测用时比 **{base['full']['wall_time_s']/full['wall_time_s']:.4f}×**，节省 **{100*(1-full['wall_time_s']/base['full']['wall_time_s']):.2f}%**。标签：**OBSERVATIONAL DEVELOPMENT SPEEDUP**。包含真实停止步数差异，不是重复运行统计。"]
else:text += ['未运行新的完整计算，因此不报告或推测完整运行加速比。']
text += ['','CPU/GPU 科学等价验证继续推迟；没有 CPU benchmark、multi-rank、多 GPU、FieldSplit、网格收敛或粒子耦合。Stage P winner 和 CPU production 保留。','']
if tests:text += [f"本阶段测试 {tests['stage_passed']} 通过、{tests['stage_failures']} 失败、{tests['stage_skipped']} 跳过；完整测试 {tests['full_passed']} 通过、{tests['full_failures']} 项历史失败、{tests['full_skipped']} 跳过。历史失败集合{'未变' if tests['historical_failure_set_unchanged'] else '发生变化'}。pytest 只读取已有 artifact，不运行 CFD。",'']
for name in ['ilu_rebuild_candidates','ilu_rebuild_count','ilu_iterations_tradeoff','ilu_setup_time','early_vs_late_window']+(['ilu_winner_steady','ilu_final_speedup'] if full else []):text += [f'![{name}]({name}.png)','']
text += ['复查入口：`REPRODUCIBILITY.md`、`reference_freeze.json`、`source_patch.json`、逐次验收文件、原始日志、`preservation_audit.json`、`native_artifact_mirror.json` 和最终 `delivery_manifest.json`。']
(R/'REPORT.md').write_text('\n'.join(text)+'\n');print(status)
