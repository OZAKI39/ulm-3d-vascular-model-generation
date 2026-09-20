#!/usr/bin/env python3
"""Write the stage report only after a terminal scientific gate decision."""
import json,sys,xml.etree.ElementTree as ET
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from sv_validation.sv11 import REPORT,load
from sv_validation.provenance import write_json

def optional(name):return load(name) if (REPORT/(name+'.json')).exists() else None
def number(v):return 'N/A' if v is None else (f'{v:.16g}' if isinstance(v,float) else str(v))
def tests(path):
    root=ET.parse(path);result={'passed':0,'failed':0,'skipped':0,'failures':[]}
    for case in root.iter('testcase'):
        if case.find('failure') is not None or case.find('error') is not None:
            result['failed']+=1;result['failures'].append(case.get('classname')+'::'+case.get('name'))
        elif case.find('skipped') is not None:result['skipped']+=1
        else:result['passed']+=1
    return result

decision=load('stage_result');assert decision['status'] in {'FAIL','CONDITIONAL PASS'}
baseline=load('baseline_solver_diagnosis');comparison=load('solver_comparison');build=load('petsc_build_manifest')
smoke=load('petsc_smoke');short=load('petsc_short_gate');run=load('petsc_short_execution');qc=load('petsc_short_qc')
accepted=optional('accepted_solution');full=optional('petsc_full_qc');reload=optional('solution_reload');preservation=load('preservation_audit')
figures=load('visuals');phase_tests=tests(REPORT/'pytest_sv11.xml');all_tests=tests(REPORT/'pytest_all.xml')
test_result={'sv11':phase_tests,'all':all_tests,'historical_failures':[f for f in all_tests['failures'] if '.test_sv1_' in f],
             'scope_decision':'User selected SV1.1 tests for stage acceptance, with historical SV1 failures separately reported. No xfail, removal, or suppression.'}
write_json(REPORT/'test_results.json',test_result)
last=qc['states'][-1] if qc['states'] else None;b=comparison['baseline'];p=comparison['petsc']
field_comparison=load('field_comparison')
linear='YES（仅已执行的前 10 步；未执行区间不作结论）' if p['linear_pass'] else 'NO'
gate_names={'SHORT_RUN_MASS_NOT_IMPROVED':'前 10 步独立积分质量误差未严格优于 SV1，短程硬门槛失败。',
            'PETSC_LINEAR_SOLVE_FAIL':'PETSc 线性求解失败，按硬门槛停止。','NONLINEAR_CONVERGENCE_FAIL':'内层线性成功，但外层非线性没有收敛。',
            'STEADY_NOT_REACHED':'未通过冻结的稳态判据。','MASS_BALANCE_FAIL':'实际解的质量守恒未达最终要求。'}
lines=['# Stage SV1.1 — 调整 SimVascular 求解器并重新验证真实血流','',
       f'**STAGE SV1.1 STATUS: {decision["status"]}**  ',f'**REASON: {decision.get("reason") or "ALL_HARD_GATES_PASS_PENDING_HUMAN_REVIEW"}**','',
       gate_names.get(decision.get('reason'),'全部数值门槛通过，等待人工审核图片。'),
       f'已执行的前 10 步：LINEAR_SOLVER_PASS，非线性收敛={p["nonlinear_pass"]}；正式 full run 和稳态验收未执行。','']
def section(title,text):lines.extend(['## '+title,'',text,''])
section('SV1 为什么失败？',f'''重新读取冻结 XML、完整日志和真实 VTU 后确认：网格有效；SV1 的 {len(baseline['linear_solves'])} 次顶层线性求解中，{b['linear_failures']} 次报告未收敛，另有 {b['ill_conditioned_warnings']} 条病态 LHS 警告。问题在第 1 步已出现，并非后期才恶化。计划 {baseline['planned_steps']} 步，实际仅到 {baseline['completed_steps']} 步；第 10 步 ε_mass={number(b['epsilon_mass_step10'])}，没有达到稳态。

逐步、逐非线性迭代、逐顶层线性求解的证据见 [baseline_linear_history.csv](baseline_linear_history.csv)、[baseline_boundary_history.csv](baseline_boundary_history.csv) 和 [完整诊断](baseline_solver_diagnosis.json)。FSILS NS 内部 GM/CG 的逐次轨迹没有被旧日志记录；没有伪造这些记录，也没有重跑 FSILS。

冻结配置为 NS/FSILS，外层最大 30，rtol=1e-10，atol=1e-24，Krylov=100；内部 GM 最大 20、CG 最大 1000，二者容差 1e-3。源码显示病态分支可能把残差诊断置零，因此旧日志中的零残差不证明高精度收敛。''')
section('这次哪些东西完全没动？',f'''geometry、selected volume/surface mesh、face IDs、WALL/INLET/三个 OUTLET、rho、mu、Q_reference、入口速度规定、wall no-slip、三个 zero-Neumann outlet、初始条件、dt、保存频率和 time-policy 均保持冻结。

生产 XML 仅 LS 部分不同：[xml_diff.txt](xml_diff.txt)。前 10 步仍使用原 XML 中的 400 步计划，通过官方 STOP_SIM=10 停止；网格路径经新目录中的符号链接解析到原冻结文件。完整性审核检查 {preservation['historical_entries_checked']} 项 SV1 历史条目和 {preservation['old_fem_entries_checked']} 项旧 FEM 条目，结果 **{preservation['status']}**；旧 FEM 原有 Git 工作区差异也未变化。见 [preservation_audit.json](preservation_audit.json)。''')
section('换了什么？',r'''只更换线性代数实现：FSILS → PETSc 3.19.6，GMRES，restart=100，最大 2000 次；右侧 ASM 预条件，overlap=2，每个子域使用 PETSc 自带 ILU(2)。这一固定配置保留速度—压力耦合，适用于当前非对称系统，未引入其它线性代数包，也未进行多个 PC 或 MPI 数量扫描。

XML 中的 `petsc-jacobi` 是当前接口接受的初始 PC；现有 `KSPSetFromOptions` 在设置阶段把它覆盖为 ASM，运行日志证实实际 PC。源码中的 XML rtol/atol 写入 FSILS 参数而非 PETSc 参数，restart 调用也被注释，所以新 XML 删去这些无效字段，改由明确冻结的 `PETSC_OPTIONS` 设置并核验；没有宣称被忽略的 XML 选项已生效。

容差未放宽。两者在未约束的自由度上采用相同的对角缩放公式 `W_ii = 1 / sqrt(abs(A_ii))`，零对角按 1 处理。PETSc 使用右预条件和缩放系统的未预条件残差范数，零初始猜测；判据为 `||W(b-Ax)||₂ < max(1e-10 ||Wb||₂, 1e-24)`。固定 Dirichlet 自由度的修正及右端项为零。并行归约和递推残差并非逐位相同；日志另含 true residual。部分末次求解由绝对容差通过，不能只看相对比值。子域 `preonly` 是一次因子应用，KSP view 中的默认子域 rtol 不参与迭代停止。

详见 [实际设置与容差对应](petsc_settings.json)、[固定源码能力](linear_solver_capabilities.json)；PETSc 官方定义见 [KSPConvergedDefault](https://petsc.org/release/manualpages/KSP/KSPConvergedDefault/) 和 [KSPSetFromOptions](https://petsc.org/release/manualpages/KSP/KSPSetFromOptions/)。具体行为已核对本地 3.19.6 源码。''')
section('PETSc 真正启用了吗？',f'''**YES。Build PASS；official fluid/newtonian smoke {smoke['status']}。** 原 binary 未链接 PETSc。本机没有兼容开发库；Ubuntu 开发包的模拟安装将新增 81 个包，因此按 [官方源码构建方式](https://petsc.org/release/install/install_tutorial/) 在工程内构建最小 PETSc，复用系统 Open MPI 与 BLAS/LAPACK。未构建 Trilinos、MUMPS、HYPRE 或 CUDA。最初的 Python/库名配置问题及修正保留在 build manifest。

svMultiPhysics 仍为 `{build['commit']}`，源码未修改。新 binary SHA256：`{build['executable_sha256']}`。`ldd` 显示实际链接 PETSc；smoke 的 {len(smoke['run']['history']['linear_solves'])} 次求解均有正向收敛原因，新 VTU 含有限速度、压力。4 MPI ranks，OMP=1。

[构建清单](petsc_build_manifest.json)含版本、路径、编译器、MPI、CMake 命令、库链接及 SHA256；[smoke 结果](petsc_smoke.json)、[smoke 原始日志](../../logs/sv1_1/petsc_smoke.log)及真实案例 [原始日志](../../logs/sv1_1/petsc_short.log)提供运行证据。''')
section('前 10 步与旧 FSILS 相比如何？',f'''| 指标 | SV1 FSILS 历史 | SV1.1 PETSc |
|---|---:|---:|
| 完成时间步 | {baseline['completed_steps']} | {p['completed_steps']} |
| 线性失败次数 | {b['linear_failures']} | {p['linear_failures']} |
| 病态警告 | {b['ill_conditioned_warnings']} | {p['ill_conditioned_warnings']} |
| 平均 / 最大迭代数 | {number(b['linear_iterations']['mean'])} / {b['linear_iterations']['max']} | {number(p['linear_iterations']['mean'])} / {p['linear_iterations']['max']} |
| step10 ε_mass，独立 VTU 积分 | {number(b['epsilon_mass_step10'])} | {number(p['epsilon_mass_step10'])} |
| 单调时钟用时 / s | {number(b['elapsed_s'])} | {number(p['elapsed_s'])} |

FSILS 列是 NS 外层次数，PETSc 列是 GMRES 次数，不能按次数直接比较计算工作量。质量误差减少量（旧−新）={number(comparison.get('epsilon_mass_decrease'))}；严格改善门槛 **{short['mass_improves']}**。近舍入量级的差异不代表实质物理改善。两次第 10 步速度场的相对 L2 差={number(field_comparison['velocity_L2_relative_difference'])}，最大压力差={number(field_comparison['pressure_max_abs_difference_Pa'])} Pa；见 [field_comparison.json](field_comparison.json)。

[solver_comparison.json](solver_comparison.json)、[新逐次求解表](petsc_short_linear.csv)、[新逐步边界流量](petsc_short_boundary.csv)。

![每步线性求解是否稳定](linear_solver_before_after.png)

![线性迭代次数](linear_iterations_over_time.png)''')
section('线性问题解决了吗？',f'''**{linear}。** 已执行记录中线性失败={p['linear_failures']}，病态警告={p['ill_conditioned_warnings']}。外层非线性验收={p['nonlinear_pass']}；验收程序同时检查残差，防止源码把达到迭代上限视作结束而被误报成收敛。

![外层流体迭代](nonlinear_convergence.png)''')
section('如果 YES，完整计算跑完了吗？',f'''{('已执行正式计算，完成 '+str(full['completed_steps'])+' 步，详见 petsc_full_execution.json。') if full else '没有启动 400 步正式计算。短程硬门槛未通过，按要求停止；未从 FSILS checkpoint 启动，也未缩短 dt、改变算法或延长时间。'}''')
section('达到稳态了吗？',f'''**{'PASS' if accepted else 'NOT_REACHED / 尚无 accepted steady solution'}。** 冻结要求是最后连续 5 个保存区间 E_u≤1e-5、E_Q≤1e-6，E_Q 继续按冻结实现以 Qtarget 归一化。短程只有 {len(qc['states'])} 个已保存状态，不能形成这些区间。

![稳态可用性与判据](steady_convergence.png)''')
section('质量守恒达到要求了吗？',f'''{'正式 accepted solution 已通过质量验收。' if accepted else '尚未获得通过最终质量验收的正式解。'} 最新短程诊断的 ε_Q={number(last['epsilon_Q'] if last else None)}，ε_mass={number(last['epsilon_mass'] if last else None)}。最终阈值仍为两者≤1e-6；短程门槛只要求 ε_mass 严格优于 SV1，没有将最终阈值提前强加到第 10 步。

独立计算来自实际 VTU 的 P1 三角面通量积分，未用目标 Q 代替实测流量；原生积分与独立积分的最大差值/Qtarget={number(last.get('native_flux_difference_over_Q') if last else None)}。

质量误差轨迹接近逐步交替减半。这与固定 generalized-alpha 参数及启动入口条件一致，是源码与测量支持的**推断**，并未用解析递推值替代实测值或放宽门槛，详见 [startup_mass_analysis.json](startup_mass_analysis.json)。

![完整瞬态质量误差](mass_balance_over_time.png)''')
section('正式速度场和压力场',f'''{'正式解见 accepted_solution.json。' if accepted else '没有 accepted steady solution，所以未生成正式 velocity_global / velocity_slices / pressure_global / pressure_sections 图。没有复用 SV1 图片。'}

短程真实字段的 finite velocity={last['velocity_finite'] if last else 'N/A'}，finite pressure={last['pressure_finite'] if last else 'N/A'}，wall no-slip={last['wall_noslip_pass'] if last else 'N/A'}。这些是瞬态诊断，不能替代正式流场验收。''')
section('三个出口最终各流多少？',('正式出口份额见 accepted_solution.json。' if accepted else '最终份额 **N/A**；稳态未通过，未发布正式 outlet fractions，也未生成 outlet_flow_split.png。')+('\n\n![瞬态诊断流量，非最终份额](flux_balance.png)' if last else ''))
section('资源消耗',f'''短程用时={number(run['elapsed_s'])} s；最大单进程 RSS={number(run['peak_rss_kib'])} KiB，不能当作四个 MPI 进程的内存总和；MPI=4，OMP=1。正式 full runtime={'见完整执行记录' if full else 'N/A'}。完整资源原文见 [petsc_short_resources.txt](../../logs/sv1_1/petsc_short_resources.txt)。

SV1 的 Python 单调时钟与 GNU time wall 记录不一致，二者均保留，仅作描述性比较，不推断加速比。

![资源消耗](solver_resource_usage.png)''')
section('结果能否重新读取？',f'''{('**'+reload['status']+'（'+reload['artifact_scope']+'）**。新 Python 进程重新读取同一 VTU，速度 L2、压力范围、Qin、三个 Qout、质量闭合和 wall no-slip 均逐项比较；见 [solution_reload.json](solution_reload.json)。') if reload else 'N/A：没有可复读的真实血管已保存 VTU。'}

正式 accepted solution reload={'PASS' if accepted and reload and reload['status']=='PASS' else 'N/A，尚无正式 accepted solution'}。''')
failures='\n'.join('- `'+f+'`' for f in all_tests['failures']) or '- 无'
section('下一步还剩什么？',f'''本阶段停止于 **{decision.get('reason') or decision['status']}**。{decision.get('next_step','后续工作需另立阶段；本阶段不继续调整网格、边界、物性、dt 或算法。')}

SV1.1 测试：{phase_tests['passed']} passed，{phase_tests['failed']} failed，{phase_tests['skipped']} skipped。完整 pytest：{all_tests['passed']} passed，{all_tests['failed']} failed，{all_tests['skipped']} skipped。按用户确认，SV1 的历史失败另列，旧测试及证据未改；没有用 xfail 或删除测试隐藏失败。未执行的下游验收明确 skipped，不视作通过。

失败项目：

{failures}

[SV1.1 JUnit](pytest_sv11.xml)、[完整 JUnit](pytest_all.xml)、[测试口径与分类](test_results.json)。''')
lines+=['本阶段新图（诊断图不等于正式结果）：','']+['- ['+Path(f['path']).name+']('+Path(f['path']).name+') — '+f['scope'] for f in figures['files']]
lines+=['','未生成的正式图：'+', '.join(figures['withheld_until_accepted'])+'；原因：缺少 accepted steady solution。' if not accepted else '','']
(REPORT/'REPORT.md').write_text('\n'.join(lines))
print('Wrote reports/sv1_1/REPORT.md')
