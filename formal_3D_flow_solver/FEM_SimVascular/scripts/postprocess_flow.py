#!/usr/bin/env python3
"""Measure every saved native state, then evaluate the frozen final-state gates."""
import json
import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from sv_validation.postprocess import SolutionMeasurements
from sv_validation.provenance import sha256, write_json
from sv_validation.validation import steady_state, require

R = ROOT / 'reports/sv1'
run = json.loads((R / 'flow_execution.json').read_text())
require(run['runs'][-1]['exit_code'] == 0, 'Native process did not exit successfully')
policy = json.loads((ROOT / 'configs/time_policy.json').read_text())
Q = json.loads((R / 'inlet_normalization.json').read_text())['Q_target_m3_s']
measure = SolutionMeasurements(ROOT/'outputs/sv1/SV_MESH/mesh_arrays.npz', Q, policy['Umean_m_s'])
directory = ROOT / 'outputs/sv1/vascular_flow/4-procs'
paths = sorted(directory.glob('result_*.vtu'), key=lambda path: int(path.stem.rsplit('_', 1)[1]))
require(paths, 'No actual native saved field exists; do not substitute prescribed fields')
states, velocity_errors, flow_errors = [], [], []
previous_velocity = previous_summary = None
for path in paths:
    u, p = measure.read(path)
    summary = measure.measure(u, p)
    step = int(path.stem.rsplit('_', 1)[1])
    summary.update(step=step, time_s=step*policy['dt_s'], path=str(path.relative_to(ROOT)), sha256=sha256(path),
                   size_bytes=path.stat().st_size)
    if previous_velocity is not None:
        error_u = measure.velocity_l2(u-previous_velocity)/max(summary['velocity_L2'], np.finfo(float).tiny)
        error_q = max(abs(summary['signed_outward_boundary_flows_m3_s'][name]-previous_summary['signed_outward_boundary_flows_m3_s'][name])/Q
                      for name in summary['signed_outward_boundary_flows_m3_s'])
        velocity_errors.append(error_u)
        flow_errors.append(error_q)
        summary.update(velocity_relative_change=error_u, boundary_flow_relative_change=error_q)
    states.append(summary)
    previous_velocity, previous_summary = u, summary
target_step = policy['first_block_steps'] + run['extensions_used']*policy['extension_steps']
scheduled_block_complete = states[-1]['step'] == target_step
log_text = (ROOT / run['runs'][-1]['log']).read_text()
linear_failures = log_text.count('The linear system solution has not converged')
stopped_for_failure = (R / 'flow_stop_decision.json').exists()
is_steady = len(velocity_errors) >= 5 and steady_state(velocity_errors, flow_errors)
steady = {'status': 'PASS' if is_steady else 'NOT_REACHED', 'saved_states': len(states),
          'required_consecutive_intervals': 5, 'velocity_change_limit': 1e-5, 'flow_change_limit': 1e-6,
          'velocity_errors': velocity_errors, 'flow_errors': flow_errors,
          'last_five_velocity_errors': velocity_errors[-5:], 'last_five_flow_errors': flow_errors[-5:],
          'extensions_used': run['extensions_used'], 'scheduled_block_complete': scheduled_block_complete,
          'scheduled_final_step': target_step, 'actual_final_step': states[-1]['step'],
          'reason': 'Insufficient saved intervals after failure stop' if len(velocity_errors) < 5 else 'Frozen final interval criteria'}
write_json(R / 'steady_state.json', steady)
write_json(R / 'saved_states.json', {'states': states})
final = states[-1]
final['mass_balance_pass'] = final['epsilon_Q'] <= 1e-6 and final['epsilon_mass'] <= 1e-6
final['steady_pass'] = is_steady
final['solver_success'] = scheduled_block_complete and not stopped_for_failure and linear_failures == 0
final['linear_nonconvergence_warnings'] = linear_failures
final['ill_conditioned_warnings'] = log_text.count('ill-conditioned LHS')
final['artifact_kind'] = 'accepted_solution_candidate' if final['solver_success'] else 'unaccepted_transient_diagnostic_state'
final['status'] = 'PASS' if all(final[k] for k in ('mass_balance_pass', 'steady_pass', 'wall_noslip_pass', 'solver_success', 'velocity_finite', 'pressure_finite')) else 'FAIL'
final['failure_reasons'] = ([] if final['solver_success'] and is_steady else ['FLOW_SOLVE_FAIL']) + ([] if final['mass_balance_pass'] else ['MASS_BALANCE_FAIL'])
final['integration_method'] = 'Exact P1 triangle flux and pressure integrals; exact P1 tetrahedral velocity mass matrix'
# Independently compare the native integral table with the VTU triangle integrals.
lines = (directory / 'B_NS_Velocity_flux.txt').read_text().splitlines()
roles = next(line.split()[2:] for line in lines if line.startswith('step'))
row = next(line.split() for line in lines if line.split() and line.split()[0] == str(final['step']))
native_flux = dict(zip(roles, map(float, row[2:])))
final['native_boundary_fluxes_m3_s'] = native_flux
final['native_vs_vtu_max_difference_over_Q'] = max(abs(native_flux[n]-v)/Q for n, v in final['signed_outward_boundary_flows_m3_s'].items())
write_json(R / 'flow_qc.json', final)
print(json.dumps({'final': final, 'steady': steady}, indent=2))
