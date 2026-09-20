#!/usr/bin/env python3
"""Reparse complete run logs and form a comparison without trusting exit status."""
import json,statistics,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from sv_validation.sv11 import REPORT,load,parse_solver_log,linear_gate,nonlinear_gate,write_csv
from sv_validation.provenance import write_json,sha256

def stats(rows):
    values=[r['linear_iterations'] for r in rows]
    return {'count':len(values),'mean':statistics.mean(values) if values else None,'min':min(values,default=None),
            'max':max(values,default=None),'median':statistics.median(values) if values else None,
            'note':'FSILS values count NS outer iterations; PETSc values count GMRES iterations, not equivalent work units'}

base=load('baseline_solver_diagnosis');baseline=load('flow_qc','sv1')
policy=json.loads((ROOT/'configs/time_policy.json').read_text())
run=load('petsc_short_execution');history=parse_solver_log((ROOT/run['log']).read_text(),policy['dt_s'])
assert not history['unparsed_rows']
run['history']=history;write_json(REPORT/'petsc_short_execution.json',run)
qc=load('petsc_short_qc');last=qc['states'][-1] if qc['states'] else None
try:linear_gate(run);linear=True;linear_error=None
except ValueError as exc:linear=False;linear_error=str(exc)
try:nonlinear_gate(history);nonlinear=True;nonlinear_error=None
except ValueError as exc:nonlinear=False;nonlinear_error=str(exc)
assert linear==load('petsc_short_gate')['linear_pass']
comparison={'baseline':{'linear_failures':base['failed_linear_solves'],'ill_conditioned_warnings':base['ill_conditioned_warnings'],
                       'linear_iterations':stats(base['linear_solves']),'epsilon_mass_step10':baseline['epsilon_mass'],
                       'elapsed_s':base['execution']['runs'][0]['elapsed_s'] if 'runs' in base['execution'] else load('flow_execution','sv1')['runs'][0]['elapsed_s']},
            'petsc':{'linear_pass':linear,'linear_error':linear_error,'nonlinear_pass':nonlinear,'nonlinear_error':nonlinear_error,
                      'linear_failures':history['failed_linear_solves'],'ill_conditioned_warnings':history['ill_conditioned_warnings'],
                      'linear_iterations':stats(history['linear_solves']),'elapsed_s':run['elapsed_s'],
                      'peak_rss_kib':run['peak_rss_kib'],'completed_steps':qc['completed_steps'],
                      'epsilon_mass_step10':last['epsilon_mass'] if last and last['step']==10 else None},
            'step_comparison':[],'strict_mass_gate':load('petsc_short_gate')['mass_improves']}
for step in sorted({r['step'] for r in base['linear_solves']} | {r['step'] for r in history['linear_solves']}):
    item={'step':step}
    for label,rows in [('FSILS',base['linear_solves']),('PETSc',history['linear_solves'])]:
        selected=[r for r in rows if r['step']==step]
        item[label]={'nonlinear_iterations':max((r['nonlinear_iteration'] for r in selected),default=0),
                     'linear_failures':sum(not r['linear_converged'] for r in selected),
                     'linear_iterations':[r['linear_iterations'] for r in selected],
                     'final_nonlinear_Ri_over_R0':selected[-1]['nonlinear_Ri_over_R0'] if selected else None}
    comparison['step_comparison'].append(item)
if last and last['step']==10:
    comparison['epsilon_mass_decrease']=baseline['epsilon_mass']-last['epsilon_mass']
    comparison['epsilon_mass_relative_decrease']=(baseline['epsilon_mass']-last['epsilon_mass'])/baseline['epsilon_mass']
    comparison['mass_change_interpretation']='Strict numerical comparison required by stage; a change near roundoff is not evidence of material conservation improvement.'
for label,rows in [('FSILS',base['boundary_history']),('PETSc',qc['boundary_history'])]:
    comparison[label+'_startup_closure_diagnostic']={'maximum_deviation_from_alternating_halving':max((abs(r['signed_closure_over_Q']+(-.5)**r['step']) for r in rows),default=None),
      'interpretation':'Observed recurrence is compatible with alpha_f=1/(1+rho_infinity)=2/3 and generalized-alpha interpolation. This is an inference, not an independent proof of its cause. No time parameters were changed.'}
write_json(REPORT/'solver_comparison.json',comparison)
write_csv(REPORT/'petsc_short_linear.csv',history['linear_solves'])
print(json.dumps(comparison['petsc'],indent=2))
