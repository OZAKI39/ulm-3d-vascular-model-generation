#!/usr/bin/env python3
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from sv_validation.sv11 import load,REPORT
def optional(name):return load(name) if (REPORT/(name+'.json')).exists() else None
stage=load('stage_result');compare=load('solver_comparison');test=load('test_results');build=load('petsc_build_manifest')
qc=load('petsc_short_qc');last=qc['states'][-1] if qc['states'] else {};accepted=optional('accepted_solution') or {}
full=optional('petsc_full_qc');reload=optional('solution_reload');na='N/A — no accepted steady solution'
out=['Stage SV1.1 completed.','','baseline FSILS:',
     '    linear failures = '+str(compare['baseline']['linear_failures']),
     '    ill-conditioned warnings = '+str(compare['baseline']['ill_conditioned_warnings']),
     '    step10 epsilon_mass = '+str(compare['baseline']['epsilon_mass_step10']),
     '', 'PETSc:', '    enabled = true','    version = '+build['version'],
     '    KSP = GMRES, restart 100, max 2000', '    PC = right ASM overlap 2 + local ILU(2)',
     '    tolerance = rtol 1e-10, atol 1e-24, diagonally scaled unpreconditioned norm',
     '', '10-step comparison:',
     '    PETSc linear failures = '+str(compare['petsc']['linear_failures']),
     '    mean iterations = '+str(compare['petsc']['linear_iterations']['mean']),
     '    max iterations = '+str(compare['petsc']['linear_iterations']['max']),
     '    epsilon_mass = '+str(compare['petsc']['epsilon_mass_step10']),
     '    improvement vs FSILS = '+str(compare['strict_mass_gate']),
     '    epsilon_mass decrease (old - new) = '+str(compare.get('epsilon_mass_decrease')),
     '', 'full run:',
     '    steps completed = '+(str(full['completed_steps']) if full else '0 — short gate did not pass'),
     '    linear failures = '+('see full execution' if full else 'N/A — not run'),
     '    nonlinear convergence = '+('see full execution' if full else 'N/A — not run; short run = '+str(compare['petsc']['nonlinear_pass'])),
     '    steady reached = '+str(bool(accepted)),
     '', 'final flow:', '    Qtarget = '+str(load('flow_qc','sv1')['Q_target_m3_s'])]
for key,label in [('Q_in_m3_s','Qin'),('Q_out_total_m3_s','Qout_total'),('epsilon_Q','epsilon_Q'),('epsilon_mass','epsilon_mass')]:out.append('    '+label+' = '+str(accepted.get(key,na)))
for i in range(1,4):out.append(f'    Qout0{i} = '+str(accepted.get('outlet_flows_m3_s',{}).get(f'OUTLET_0{i}',na)))
out+=['', 'fields:']
for key,label in [('velocity_finite','velocity finite'),('pressure_finite','pressure finite'),('wall_noslip_pass','wall no-slip')]:out.append('    '+label+' = '+str(accepted.get(key,na))+'; short diagnostic = '+str(last.get(key,'N/A')))
out+=['','resources:','    runtime = '+str(compare['petsc']['elapsed_s'])+' s (short run); full = '+('see full execution' if full else 'N/A'),
      '    peak RSS = '+str(compare['petsc']['peak_rss_kib'])+' KiB, maximum single process, not MPI sum',
      '    MPI ranks = 4; OMP_NUM_THREADS = 1','','reload:',
      '    '+(reload['status']+' — '+reload['artifact_scope'] if reload else 'N/A'),'', 'tests (SV1.1):']
for key in ('passed','failed','skipped'):out.append('    '+key+' = '+str(test['sv11'][key]))
out+=['','full pytest:']
for key in ('passed','failed','skipped'):out.append('    '+key+' = '+str(test['all'][key]))
out+=['    historical SV1 failures = '+str(len(test['historical_failures'])),'','report:','    reports/sv1_1/REPORT.md','','human review:']
for name in ['linear_solver_before_after.png','linear_iterations_over_time.png','mass_balance_over_time.png','nonlinear_convergence.png','steady_convergence.png','velocity_global.png','velocity_slices.png','pressure_global.png','pressure_sections.png','flux_balance.png','outlet_flow_split.png','solver_resource_usage.png']:
    out.append('    '+name+(' — generated' if (REPORT/name).exists() else ' — NOT GENERATED: no accepted solution'))
out+=['','STAGE SV1.1 STATUS:','    '+stage['status'],'REASON:','    '+str(stage.get('reason')),'']
text='\n'.join(out);(REPORT/'terminal_summary.txt').write_text(text);print(text)
