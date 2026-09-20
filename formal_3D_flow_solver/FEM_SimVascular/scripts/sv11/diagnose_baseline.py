#!/usr/bin/env python3
import json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from sv_validation.sv11 import load,parse_solver_log,parse_boundary,linear_settings,write_csv,REPORT
from sv_validation.provenance import write_json,sha256
policy=json.loads((ROOT/'configs/time_policy.json').read_text())
qc=load('flow_qc','sv1');execution=load('flow_execution','sv1');mesh=load('mesh_validity','sv1');quality=load('mesh_quality','sv1')
log=ROOT/execution['runs'][0]['log'];parsed=parse_solver_log(log.read_text(),policy['dt_s'])
assert not parsed['unparsed_rows'] and not parsed['unassigned_warnings']
assert parsed['failed_linear_solves']==qc['linear_nonconvergence_warnings']
assert parsed['ill_conditioned_warnings']==qc['ill_conditioned_warnings']
flows=parse_boundary(ROOT/'outputs/sv1/vascular_flow/4-procs/B_NS_Velocity_flux.txt',qc['Q_target_m3_s'])
by_step={r['step']:r for r in flows}
for row in parsed['linear_solves']:row['boundary_at_end_of_step']=by_step.get(row['step'])
result={**parsed,'status':'DIAGNOSED','configuration':linear_settings(ROOT/'configs/sv_flow.xml'),
        'planned_steps':policy['first_block_steps'],'completed_steps':max(r['step'] for r in parsed['linear_solves']),
        'step10_independent_vtu_qc':qc,'boundary_history':flows,'mesh_validity':mesh,'mesh_quality':quality,
        'execution':execution,'steady':load('steady_state','sv1'),
        'residual_semantics':{'Ri/R1':'Nonlinear current initial residual / initial residual of this time step, as printed by output.cpp',
           'Ri/R0':'FSILS.RI.iNorm / eq.iNorm, as printed by output.cpp; not an absolute physical norm',
           'R/Ri':'FSILS.RI.fNorm / FSILS.RI.iNorm, linear residual ratio',
           'absolute_residual':'Not logged by this historical run; unavailable, never reconstructed from rounded ratios',
           'ill_conditioned_zero':'ns_solver.cpp resets fNorm and dB to zero on negative Resc/Resm; zero is not independent evidence of convergence',
           'inner_subsolves':'Only one top-level solve per nonlinear iteration is logged; inner GM/CG subsolve history is unavailable'},
        'diagnosis':'Failures begin at step 1, nonlinear iteration 2; ill-conditioned warning already accompanies step 1 iteration 1. Failure is present from startup, not only at late time.',
        'artifacts':{str(p.relative_to(ROOT)):sha256(p) for p in [ROOT/'configs/sv_flow.xml',ROOT/'configs/time_policy.json',log,ROOT/'reports/sv1/flow_execution.json',ROOT/'reports/sv1/flow_qc.json']}}
write_json(REPORT/'baseline_solver_diagnosis.json',result)
write_csv(REPORT/'baseline_linear_history.csv',parsed['linear_solves'])
write_csv(REPORT/'baseline_boundary_history.csv',flows)
print('Parsed baseline:',len(parsed['linear_solves']),'solves;',parsed['failed_linear_solves'],'failures;',parsed['ill_conditioned_warnings'],'ill-conditioned warnings; first failed step',parsed['first_failed_step'])
