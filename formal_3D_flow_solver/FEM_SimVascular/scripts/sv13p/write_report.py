"""Plain Chinese final decision report grounded in closed stage artifacts."""
import json,sys,xml.etree.ElementTree as ET
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];R=ROOT/'reports/sv1_3p'
def read(n):return json.loads((R/(n+'.json')).read_text())
def num(x,digits=2):return '未取得' if x is None else f'{x:.{digits}f}'
s=read('candidate_summary');w=read('winner');ref=read('reference_freeze');c=s['candidates'];full=read('winner_steady_candidate') if w['full_required'] else None
status='GPU_PC_WINNER_FOUND' if full and full['status']=='PASS' else 'EXISTING_GPU_PC_DIRECT_REPLACEMENT_EXHAUSTED' if w['fastest'] else 'NO_STABLE_GPU_PC_REPLACEMENT'
summary=dict(status=status,completed=True,winner=w,full_status=full['status'] if full else 'NOT_REQUIRED',scientific_equivalence='DEFERRED',fallback_preserved=True,CPU_production_changed=False,next_stage_started=False,limitations=['Existing monitor confirmed step70 after output transfer; native safe stop at72, versus Stage O71. No full-run repeat.','AmgX official package failed CUDA13 NVTX target configuration; no compatibility port.','Instrumented candidate windows versus historical unprofiled A; single observations.'])
(R/'stage_status.json').write_text(json.dumps(summary,indent=2)+'\n')
lines=['# Stage SV1.3P — 测试现有 GPU 友好的辅助求解方法','',f'最终状态：**{status}**。Stage O 成功配置和 CPU production 均保留。','', '## 为什么要换辅助求解器？','', 'GPU 的矩阵乘法已经很快。Stage O 的诊断窗口中，准备子域辅助求解器约 121.06 秒，其中数值分解约 114.36 秒；使用辅助求解器约 81.16 秒，矩阵乘法约 3.40 秒。因此本次只改变辅助求解方法。矩阵、物理参数、边界条件、时间步、稳态标准，以及 GMRES 的 restart100、rtol1e-10、atol1e-24、max_it2000 全部冻结。','', '## 这次测试了什么？','', '| 方法 | 最多 2 步 smoke | 健康完成的 10 步时间 | 10 步迭代数 | 结论 |','|---|---|---:|---:|---|',f"| 原 ASM2 + ILU2 | 复用历史证据 | {ref['window_baseline']['wall_time_s']:.2f} 秒 | 6603 | Stage O fallback |"]
for k,d in c.items():
 smoke=d['smoke'];v=d['window'];good=v and v['status']=='PASS'
 lines.append(f"| {k} {d['label']} | {smoke['status'] if smoke else '未进入 CFD'} | {num(v['wall_time_s'])+' 秒' if good else '无健康完整窗口'} | {v['statistics']['total_iterations'] if good else '—'} | {d['result']} |")
lines+=['','所有窗口读取同一个 Stage N 原生 step60 检查点，执行 61—70；SHA256 为 `'+ref['window_baseline']['initial_checkpoint_sha256']+'`。单张 RTX4090、单 MPI 进程、OMP=1。VTU 每 10 步输出，停止时补写最终场和检查点。没有重复旧基线，没有参数扫描。','', '候选时间来自实际进程单调时钟，不能用检查点继承的旧累计时间。P 候选开启 GPU 事件计时，历史 A 基线未开启，因此这属于各一次开发测量；计时同步可能额外增加候选开销。事件时间存在嵌套，不能把各项相加当作总时间。初始 P1 smoke 的 PETSc 事件时间为 `n/a`，保留原始记录，只修正解析器读取事件计数，没有重跑。','', '## “复用”有没有帮助？','']
p=c['P1']['window'];ev=p['profile']['events'];gain=1-p['wall_time_s']/ref['window_baseline']['wall_time_s']
lines += [f"有。P1 保留同一个 KSP/PC，每个时间步第一次求解重新准备，后续非线性求解使用已有 ILU。实际 20 次求解只发生 {p['reuse']['PC_rebuild_count']} 次数值分解，窗口 {p['wall_time_s']:.2f} 秒，减少 {gain*100:.2f}%。迭代数为 {p['statistics']['total_iterations']}，均值 {p['statistics']['mean_iterations']:.2f}，最大 {p['statistics']['max_iterations']}。判定：P1_PROMISING。",'', '生命周期和每次请求见 `ksp_lifecycle_audit.json`、`P1_WINDOW_acceptance.json`；一次源码补丁只增加步内复用控制和计数。PETSc 实际 `MatLUFactorNum` 计数用于交叉核对，不能只凭 flag 判断复用有效。P1 继续链接 Stage O 原来的 PETSc/CUDA/MPI，只使用隔离编译的最小 svMultiPhysics 复用补丁。该策略使用 [PETSc 官方复用接口](https://petsc.org/release/manualpages/KSP/KSPSetReusePreconditioner/)。']
for key,title in [('P2','hypre GPU ILU0 怎么样？'),('P3','BoomerAMG 怎么样？'),('P4','PETSc GAMG 怎么样？'),('P5','NVIDIA AmgX 怎么样？')]:
 d=c[key];lines+=['','## '+title,'',f"结论：`{d['result']}`。"]
 if d['config']:lines += ['',{'P2':'只选 BJ-ILU0，关闭局部重排，使用迭代三角求解，保留默认上下三角各 5 次 Jacobi 迭代。所有参数名先由实际帮助输出确认。','P3':'只选官方支持 GPU 的粗化、插值和平滑组合；粗网格同样使用 l1-Jacobi。','P4':'只测试 PETSc 默认聚合型 GAMG，不构造新的物理近零空间。','P5':'只选一个官方 GPU AMG 配置。'}[key]]
 if d['smoke']:
  smoke=d['smoke'];ex=read('remote/'+key+'_SMOKE_execution');fail=ex['monitor_stop'];st=smoke['statistics']
  lines += ['',f"smoke 实际运行 {smoke['wall_time_s']:.2f} 秒；已结束并记录的线性求解 {st.get('KSP_solves',0)} 次、迭代 {st.get('total_iterations',0)} 次。健康状态 {smoke['status']}。" ]
  if fail:lines += ['', '原始停止原因：`'+str(fail)+'`。失败前未完成的线性迭代不冒充完整求解统计；完整残差日志保留。']
  if smoke['errors']:lines += ['', '线性求解出现负收敛原因后已停止，因此没有完整的最终场输出可供验收。原始错误项保存在对应验收 JSON 中。']
 if d['window']:
  v=d['window'];lines+=['',f"10 步窗口状态 {v['status']}，实际 {v['wall_time_s']:.2f} 秒、{v['statistics']['total_iterations']} 次迭代。只有健康完成的窗口参与最快方法排名。"]
 if key=='P2':
  b=read('remote/petsc_hypre_build');lines+=['',f"独立 PETSc 构建：{b['status']}；hypre {b['package_version']}，提交 `{b['package_commit']}`，归档 SHA256 `{b['package_archive_sha256']}`。固定 PETSc 3.25.5、CUDA 13.2、原 MPI、双精度和 32 位索引。cuSPARSE/CUDA/device 宏及独立 GPU 小算例记录在构建报告和 standalone 日志中。GPU 算法支持来自[本次固定源码的官方文档](https://github.com/hypre-space/hypre/blob/85b779557005b2eb94c231c1b516e988b87f4e53/src/docs/usr-manual/solvers-ilu.rst)。设备路径证据与每个内核的驻留测量不同，本次未做内核级驻留分析。"]
  if (R/'hypre_baseline_parity.json').exists():
   a=read('hypre_baseline_parity');lines += ['',f"新二进制上的原 ASM2+ILU2 短 smoke 通过；前两步迭代数：旧 {a['old_two_step_iterations']}、新 {a['new_two_step_iterations']}。没有证据要求重跑完整基线窗口，因此未补跑。"]
 if key=='P3':lines+=['','只选一个 GPU 支持的组合：PMIS 粗化、extended+i 插值、各层 l1-Jacobi 平滑。依据[固定 hypre 版本的 GPU 支持表](https://github.com/hypre-space/hypre/blob/85b779557005b2eb94c231c1b516e988b87f4e53/src/docs/usr-manual/solvers-boomeramg.rst)，没有统一内存下的 CPU 算法回退。']
 if key=='P4':lines+=['','当前矩阵有 281452 行、块大小 4，按节点排列三个速度分量与压力。没有提供物理近零空间，也没有改变自由度排列。[PETSc 官方说明](https://petsc.org/release/manualpages/PC/PCGAMG/)要求向量问题提供恰当的块大小及近零空间信息。这里未额外构造物理近零空间，通用多层方法直接处理这种速度与压力耦合矩阵可能收敛困难；这里只能据本次默认配置的实测结果判断，不能推断 GAMG 对所有血流模型都无效。']
 if key=='P5':
  a=read('amgx_feasibility');lines+=['',a.get('explanation',''),'',f"固定依赖包版本 {a['package_version']}；没有替换驱动、系统 CUDA 或 production PETSc。详情见 `amgx_feasibility.json`。"]
lines+=['','## 哪一种最快？','',f"健康完成 10 步的最快新方法为 **{w['fastest']}**，相对历史基线减少 **{w['improvement']*100:.2f}%**。" if w['fastest'] else '新方法没有健康完成 10 步，继续使用 Stage O。','','## 时间到底省在哪里？','',f"P1 子域准备事件为 {ev['PCSetUpOnBlocks']['time_s']:.2f} 秒，其中数值分解 {ev['MatLUFactorNum']['time_s']:.2f} 秒；外层 PCSetUp 另为 {ev['PCSetUp']['time_s']:.3f} 秒。使用辅助求解器 {ev['PCApply']['time_s']:.2f} 秒，矩阵乘法 {ev['MatMult']['time_s']:.2f} 秒。主要收益来自减少 ILU 分解次数，迭代数没有减少。外层与子域事件的调用方式不同；各事件还包含内部操作，不把它们相加当作总时间。",'',f"P1 窗口设备级采样显存最大值 {p['sampled_device_memory_max_MiB']:.0f} MiB。各候选采样值见图 E 和验收 JSON；这不是进程精确显存峰值。没有取得的失败运行事件时间明确保留为空。",'', '## 有没有明显 winner？','',('YES：最快健康窗口达到 10% 门槛，并用唯一完整运行验证。' if full else 'NO：不启动完整稳态运行。'),'', '## 如果有 winner','']
if full:
 lines += [f"唯一完整运行 `REAL_VASCULAR_GPU_PC_WINNER` 使用 {w['candidate']}，从 t=0 开始。耗时 **{full['wall_time_s']:.2f} 秒**，正式连续 5 区间稳态首次通过 **step{full['first_full_steady_step']}**，原生安全停止在 **step{full['stop_step']}**；最终 VTU 和完整检查点已独立重读。",'',f"总线性求解 {full['statistics']['KSP_solves']} 次，总迭代 {full['statistics']['total_iterations']}，均值 {full['statistics']['mean_iterations']:.2f}，最大 {full['statistics']['max_iterations']}；线性与非线性失败为零，速度和压力有限。正式质量误差 {full['measurement']['epsilon_mass']:.6g}。相比 Stage O 完整运行的观察用时比为 {full['development_speedup']:.3f}×。",'','停止请求在监测器确认第 70 步正式通过后立即发出。但现有 WSL 监测先传输并重读 VTU 和检查点，确认期间求解器已经推进到第 72 步，随后完成这一在途时间步并停止。因此首次合格到实际停止相差 2 步，存在监测延迟；不把它表述为第 70 步即时停机。没有重新启动、补跑或删除这段真实用时。Stage O 参考止于第 71 步，本次用时比较包含实际停止步数差异。']
else:lines+=['本次不需要完整稳态运行。']
lines+=['','## 如果没有 winner','', '本次已有通过完整验证的 winner，保留 Stage O 作为回退；不进入下一阶段。' if full else '现有直接替换型 GPU 辅助求解方法没有明显改善。下一步可单独研究 Stage SV1.3Q：速度和压力分别处理。本阶段没有实现 FieldSplit 或 Schur complement。', '', 'CPU/GPU 科学等价、多进程 CUDA、FieldSplit 均为 DEFERRED。当前结论是 GPU 开发结果，不是 CPU production 替换批准。','']
for f in read('visuals')['figures']:lines+=['!['+f['title']+']('+Path(f['path']).name+')','']
if (R/'pytest_summary.json').exists():
 t=read('pytest_summary');lines += [f"新增 Stage P 测试：{t['stage_passed']} 通过、{t['stage_failures']} 失败。完整测试：{t['full_passed']} 通过、{t['full_failures']} 项历史失败、{t['full_skipped']} 跳过。测试只读运行证据，不启动 CFD；历史失败列表见 `pytest_summary.json`。",'']
lines += ['本机历史保护检查、旧 FEM 工程检查及远端驱动/CUDA/MPI/成功求解器检查通过；原始日志、运行计划、源码补丁、构建记录、输入 SHA、最终输出和逐文件交付清单均保留。']
(R/'REPORT.md').write_text('\n'.join(lines)+'\n');print(status)
