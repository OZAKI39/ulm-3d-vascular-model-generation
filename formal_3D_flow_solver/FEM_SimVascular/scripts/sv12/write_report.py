#!/usr/bin/env python3
"""Write the requested stage report and terminal summary from final evidence."""
import json
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from sv_validation.sv12 import REPORT, CONFIG, load
from sv_validation.provenance import write_json, now

flow, qc, history = load('flow_execution'), load('saved_state_qc'), load('solver_history')
resources, audit, reload = load('solver_resource_usage'), load('preservation_audit'), load('solution_reload')
tests, manifest, visuals = load('test_results'), load('frozen_input_manifest'), load('visuals')
policy = json.loads((CONFIG/'policy.json').read_text())
model = ET.parse(CONFIG/'sv_flow.xml')
rho = float(model.findtext('.//Density'))
mu = float(model.findtext('.//Viscosity/Value'))
steady_assessment = load('steady_assessment')
accepted = flow['accepted_solution_available'] and reload['status'] == 'PASS'
q = load('accepted_solution') if accepted else qc['states'][-1]
latest = qc['intervals'][-1] if qc['intervals'] else {'E_u': None, 'E_Q': None}
linear, nonlinear = history['new_linear_iteration_distribution'], history['new_nonlinear_iteration_distribution']
reason = flow['reason']
if audit['status'] != 'PASS': reason = 'FROZEN_INPUT_CHANGED'
elif flow['accepted_solution_available'] and reload['status'] != 'PASS': reason = 'RESULT_RELOAD_MISMATCH'
current_ok = tests['sv12']['failed'] == 0 and tests['sv12']['errors'] == 0 and tests['sv12']['skipped'] == 0
regression_ok = not tests['unexpected_regression_failures'] and tests['all']['errors'] == 0
complete = accepted and not reason and current_ok and regression_ok
status = 'CONDITIONAL PASS' if complete else 'FAIL'
if not reason and not current_ok: reason = 'CURRENT_ACCEPTANCE_TEST_FAILURE'
if not reason and not regression_ok: reason = 'UNEXPECTED_REGRESSION_FAILURE'
result = {'stage': 'SV1.2', 'status': status, 'reason': reason, 'timestamp': now(),
          'numerical_accepted_solution_available': accepted, 'current_acceptance_tests_pass': current_ok,
          'historical_regression_failures_classified': tests['historical_expected_failures'],
          'history_and_old_FEM_unchanged': audit['status'] == 'PASS',
          'pending': 'User visual review of formal accepted field figures' if complete else None,
          'next_stage_started': False}
write_json(REPORT/'stage_result.json', result)

def num(value): return 'N/A' if value is None else f'{value:.12g}'
def pic(name):
    return f'![{name}]({name})\n' if (REPORT/name).exists() else '未生成：未获得 accepted steady solution。\n'
def table(rows):
    return '| 项目 | 实测值 |\n|---|---:|\n'+''.join(f'| {key} | {value} |\n' for key, value in rows)+'\n'
def distribution(d):
    return f"总数 {d['count']}，总迭代 {d['total']}，mean={num(d['mean'])}，median={num(d['median'])}，P95={num(d['P95'])}，max={num(d['max'])}"

text = '# Stage SV1.2 — 使用稳定 PETSc 求解器完成真实三维血管稳态流场\n\n'
text += f"**状态：{status}**"+(f"；原因：`{reason}`。\n\n" if reason else '；数值验收完成，等待用户人工审核正式场图。\n\n')
text += '## 为什么现在可以继续完整计算？\n\n'
text += '只重查 SV1.1 已有证据，未重复短程 benchmark。前 10 步的 30 次 PETSc 求解全部收敛，0 linear failures、0 ill-conditioned warnings，外层非线性收敛、速度和压力有限、实际壁面零速。SV1 与 SV1.1 的历史 FAIL 状态保留。\n\n'
text += '[前置证据](prerequisites.json)、[原生续算审计](restart_decision.json)。checkpoint 含完整 Y_n、A_n、时间步、物理时间和方程初始范数；同一网格、二进制与四个 MPI ranks。VTU 仅用于独立 QC，未用作 restart。\n\n'
text += '## 这次改了什么？\n\n'
text += '没有改变科学模型。SV1.2 policy 取消 step10 质量误差必须严格优于 FSILS 的短程排名门槛；没有修改 SV1.1 历史报告。本阶段继续原定完整时间块，移除短程 STOP_SIM=10，XML 只启用合法 native restart；若批准规则触发唯一一次 extension，则仅把总步数改为 800。\n\n'
text += '## 哪些内容完全冻结？\n\n'
text += table([(key, str(value)) for key, value in manifest['unchanged'].items()])
text += f"dt={policy['dt_s']:.17g} s，每 {policy['save_interval_steps']} 步保存；rho={rho:g} kg/m³，mu={mu:.12g} Pa·s，Qtarget={q['Q_target_m3_s']:.17g} m³/s。上述值读取自冻结 artifact。原始入口流量、wall no-slip、三个 zero-Neumann outlets、积分参数及初值历史均保留。PETSc 3.19.6，right GMRES(100)，max 2000，ASM overlap 2 + ILU(2)，4 MPI ranks，OMP=1。逐项 SHA 与实际运行 OPTIONS 见 [冻结清单](frozen_input_manifest.json)。\n\n"
text += f"最终保护审计 **{audit['status']}**：{audit['historical_entries_checked']} 个历史条目、旧 FEM {audit['old_fem_entries_checked']} 个条目及 Git 状态。没有历史文件、旧测试或旧失败状态被改写。\n\n"
text += '## 计算到底跑了多少步？\n\n'
text += f"从合法原生 step10 checkpoint 开始；绝对时间步达到 **{flow['last_step']}**，SV1.2 新完成 {max(0, flow['last_step']-10)} 步；extension 使用 {flow['extension_used']} 次。累计物理时间 {flow['physical_time_s']:.12e} s，即 {flow['physical_time_s']/policy['viscous_time_s']:.6g} 个粘性时间尺度。最终保存步为 {q['step']}。前 10 步属于已验证继承证据，未重复计算。\n\n"
if flow['reason']:
    text += f"本次真实终止原因为 `{flow['reason']}`，不把未完成的时间块称为完成 400/800 步。详见 [执行记录](flow_execution.json) 与原始日志。\n\n"
text += '## PETSc 全程稳定吗？\n\n'
text += f"新增 linear failures={history['linear_failures']}，ill-conditioned warnings={history['ill_conditioned_warnings']}。逐次 KSP reason、报告残差和 true residual monitor 均保留，不以退出码代替收敛证据。\n\n新增求解迭代分布：{distribution(linear)}。\n\n包含继承 step1–10 的完整历史：{distribution(history['complete_linear_iteration_distribution'])}。\n\n"
text += pic('linear_iterations_over_time.png')+'\n'
text += '## 外层非线性求解稳定吗？\n\n'
text += f"新增失败步数 {history['nonlinear_failures']}；{distribution(nonlinear)}。每步记录初始和最终缩放残差以及 Ri/R0、Ri/R1；最少 2 次、最多 12 次，残差门限沿用 1e-10。达到迭代上限不能单独视为收敛。\n\n"
text += pic('nonlinear_convergence.png')+'\n'
text += '## 质量误差如何随时间变化？\n\n'
text += f"所有 {len(qc['states'])} 个保存状态均对原生 P1 速度做真实端面积分。step10 是复制的 SV1.1 QC 起点；其 εmass={qc['states'][0]['epsilon_mass']:.12g}。最新保存状态 εmass={q['epsilon_mass']:.12g}。早期 startup transient 不能替代最终质量验收。完整序列见 [mass_history.csv](../../outputs/sv1_2/qc/mass_history.csv)。\n\n"
text += pic('mass_balance_over_time.png')+'\n'
text += '## 是否达到稳态？\n\n'
text += f"首个通过质量门限的保存场：step {steady_assessment['first_saved_state_passing_mass_step']}；首次完整满足连续 5 个联合稳态区间：step {steady_assessment['first_five_joint_intervals_end_step']}。这些记录不触发提前停止，最终验收仍使用完整时间块的最后保存场。\n\n"
text += f"联合 steady 判据：**{qc['steady_reached']}**。最后 E_u={num(latest['E_u'])}，E_Q={num(latest['E_Q'])}；连续共同通过 {qc['consecutive_passing_intervals']} 个保存区间。E_u 由四面体 P1 速度差的精确体积 L2 积分计算，E_Q 为四端口通量差最大值/Qtarget；要求最后连续 5 区间 E_u≤1e-5 且 E_Q≤1e-6，至少需要 6 个连续保存状态。\n\n"
text += '| 保存区间 | E_u | E_Q | 联合通过 |\n|---|---:|---:|---|\n'
for row in qc['intervals'][-5:]:
    text += f"| {row['previous_step']} → {row['step']} | {row['E_u']:.8e} | {row['E_Q']:.8e} | {row['E_u']<=1e-5 and row['E_Q']<=1e-6} |\n"
text += '\n'+pic('steady_convergence.png')+'\n'
text += '## 最终质量守恒如何？\n\n'
text += ('下表来自 accepted 最后保存状态。\n\n' if accepted else '未获得 accepted 最终状态。下表仅保留最新瞬态 QC，不能作为正式稳态结论。\n\n')
text += table([('Qtarget / m³·s⁻¹', num(q['Q_target_m3_s'])), ('Qin / m³·s⁻¹', num(q['Q_in_m3_s'])),
    *[(role+' / m³·s⁻¹', num(q['outlet_flows_m3_s'][role])) for role in ('OUTLET_01', 'OUTLET_02', 'OUTLET_03')],
    ('Qout_total / m³·s⁻¹', num(q['Q_out_total_m3_s'])), ('epsilon_Q', num(q['epsilon_Q'])), ('epsilon_mass', num(q['epsilon_mass']))])
text += '两项相对误差均以 Qtarget 归一化，验收门限均为 1e-6。\n\n'+pic('flux_balance.png')+'\n'
text += '## 三个出口最终各流多少？\n\n'
if accepted:
    text += table([(role, f'{fraction:.10%}') for role, fraction in q['outlet_fractions'].items()])
    text += '上述比例为实际 Qout_i / Qout_total，属于 accepted steady solution 的 final outlet flow split。\n\n'
else:
    text += '没有 accepted steady solution，不报告正式出口分流比例。\n\n'
text += pic('outlet_flow_split.png')+'\n'
text += '## 最终速度场\n\n'
if accepted:
    text += f"原生 [VTU](../../{q['path']})，SHA256 `{q['sha256']}`；velocity finite={q['velocity_finite']}，max |u|={q['velocity_max_m_s']:.12g} m/s，volume L2={q['velocity_L2']:.12g}。没有平滑、裁剪、过滤或覆写保存场；全局图显示全部原生节点，截面只用于显示原生 P1 场。\n\n"
text += pic('velocity_global.png')+'\n'+pic('velocity_slices.png')+'\n'
text += '## 最终压力场\n\n'
if accepted:
    text += f"原生全局范围 [{q['pressure_range_pa'][0]:.12g}, {q['pressure_range_pa'][1]:.12g}] Pa；pressure finite={q['pressure_finite']}。保持 zero-Neumann outlet reference，未减出口压力、未 pin pressure。\n\n"
    text += table([(role+' 面积平均 / Pa', num(value)) for role, value in q['area_average_pressure_pa'].items()])
    text += '12 个几何主轴截面使用三角面面积加权 P1 压力。平面可能跨多个分支，不能解释为单根血管中心线压降。\n\n'
text += pic('pressure_global.png')+'\n'+pic('pressure_sections.png')+'\n'
text += '## 壁面条件\n\n'
text += f"{'最终 accepted 场' if accepted else '最新瞬态场'} WALL max |u|={q['wall_velocity_max_m_s']:.12g} m/s，P95={q['wall_velocity_P95_m_s']:.12g} m/s；冻结门限 {q['wall_noslip_tolerance_m_s']:.12g} m/s，实际场检查={q['wall_noslip_pass']}。\n\n"
text += '## 结果重读\n\n'
text += f"**{reload['status']}**。"+('求解器退出后，独立新 Python 进程重读原生 VTU，重算速度 L2、压力范围、所有端口流量、质量误差、壁面速度和分流比例。比较采用 rtol=1e-12 与量纲对应绝对容差，未使用会淹没微小 SI 流量的默认绝对容差。' if accepted else '没有 accepted steady solution，因此不以瞬态重读替代正式验收。')+'\n\n'
text += '[独立重读证据](solution_reload.json)。\n\n'
text += '## 资源\n\n'
text += table([('SV1.2 solver wall / s', num(resources['total_solver_wall_s'])), ('Python monotonic / s', num(resources['total_python_monotonic_s'])),
    ('peak single-process RSS / KiB', resources['peak_single_process_RSS_KiB']), ('MPI ranks / OMP', '4 / 1'),
    ('结果文件数', resources['result_file_count']), ('结果逻辑字节数', resources['result_logical_size_bytes'])])
text += 'GNU time 的 RSS 是最大单进程常驻量，不是所有 MPI rank 内存之和。原生日志 elapsed 继承 checkpoint 计时，GNU/Python 新时钟分别保留；没有推算加速比或比较 FSILS 性能。PETSc 和非线性迭代分布见前述完整记录。\n\n'+pic('solver_resource_usage.png')+'\n'
text += f"SV1.2 当前测试：{tests['sv12']['passed']} passed，{tests['sv12']['failed']} failed，{tests['sv12']['skipped']} skipped，{tests['sv12']['errors']} errors。完整回归：{tests['all']['passed']} passed，{tests['all']['failed']} failed，{tests['all']['skipped']} skipped，{tests['all']['errors']} errors。\n\n"
text += '历史阶段预期失败按原样单独列出，未删除、修改或 xfail：\n\n'
text += ''.join('- `'+name+'`\n' for name in tests['historical_expected_failures'])+'\n'
if tests['unexpected_regression_failures']:
    text += '额外失败：'+', '.join(tests['unexpected_regression_failures'])+'。\n\n'
text += '## 尚未完成什么？\n\n- experimental pump Q\n- mesh convergence\n- WSS validation\n- RBC/microbubble\n'
(REPORT/'REPORT.md').write_text(text)

lines = ['Stage SV1.2 completed.', '', 'frozen inputs:']
lines += [f'    {key} unchanged = {value and audit["frozen_inputs_unchanged"]}' for key, value in manifest['unchanged'].items()]
lines += ['', 'run:', '    started from = native restart, step10', f"    steps completed = {flow['last_step']} ({flow['last_step']-10} newly computed)",
          f"    extension used = {flow['extension_used']}", f"    physical time = {flow['physical_time_s']:.12e} s", '', 'linear solver:',
          f"    failures = {history['linear_failures']}"]
lines += [f'    {key} iterations = {num(linear[key])}' for key in ('mean', 'median', 'P95', 'max')]
lines += ['', 'nonlinear:', f"    failures = {history['nonlinear_failures']}", f"    mean iterations = {num(nonlinear['mean'])}", f"    max iterations = {num(nonlinear['max'])}",
          '', 'steady:', f"    E_u final = {num(latest['E_u'])}", f"    E_Q final = {num(latest['E_Q'])}",
          f"    consecutive passing intervals = {qc['consecutive_passing_intervals']}", f"    reached = {qc['steady_reached']}", '', 'flow:']
lines += [f'    {name} = {num(value)}' for name, value in [('Qtarget', q['Q_target_m3_s']), ('Qin', q['Q_in_m3_s']),
          *q['outlet_flows_m3_s'].items(), ('Qout_total', q['Q_out_total_m3_s']), ('epsilon_Q', q['epsilon_Q']), ('epsilon_mass', q['epsilon_mass'])]]
lines += ['', 'outlet fractions:']
lines += [f'    {role} = {num(q["outlet_fractions"][role]) if accepted else "N/A: no accepted state"}' for role in ('OUTLET_01', 'OUTLET_02', 'OUTLET_03')]
lines += ['', 'fields:', f"    velocity finite = {q['velocity_finite']}", f"    pressure finite = {q['pressure_finite']}", f"    wall max velocity = {num(q['wall_velocity_max_m_s'])} m/s", '', 'resources:',
          f"    runtime = {resources['total_solver_wall_s']} s GNU; {resources['total_python_monotonic_s']} s Python", f"    peak RSS = {resources['peak_single_process_RSS_KiB']} KiB (single process)",
          '', 'reload:', '    '+reload['status'], '', 'SV1.2 tests:']
lines += [f'    {key} = {tests["sv12"][key]}' for key in ('passed', 'failed', 'skipped', 'errors')]
lines += ['', 'full regression:']+[f'    {key} = {tests["all"][key]}' for key in ('passed', 'failed', 'skipped', 'errors')]
lines += [f"    historical expected failures = {len(tests['historical_expected_failures'])}", '', 'report:', '    reports/sv1_2/REPORT.md', '', 'human review:']
for name in ('mass_balance_over_time', 'linear_iterations_over_time', 'nonlinear_convergence', 'steady_convergence',
             'velocity_global', 'velocity_slices', 'pressure_global', 'pressure_sections', 'flux_balance', 'outlet_flow_split', 'solver_resource_usage'):
    lines.append(f'    {name}.png = '+('available' if (REPORT/(name+'.png')).exists() else 'not generated: no accepted solution'))
lines += ['', 'STAGE SV1.2 STATUS:', '    '+status+(f' / {reason}' if reason else ''), '']
summary = '\n'.join(lines)
(REPORT/'terminal_summary.txt').write_text(summary)
print(summary, flush=True)
