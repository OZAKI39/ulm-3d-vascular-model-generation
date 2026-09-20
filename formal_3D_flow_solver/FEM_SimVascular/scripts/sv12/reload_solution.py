#!/usr/bin/env python3
"""Run only after the native solver exits; independently reread its final VTU."""
import json
import os
import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from sv_validation.sv12 import REPORT, CONFIG, load
from sv_validation.provenance import sha256, write_json, now
from sv_validation.postprocess import SolutionMeasurements

run = load('flow_execution')
if not run['accepted_solution_available']:
    write_json(REPORT/'solution_reload.json', {
        'status': 'NOT_APPLICABLE', 'fresh_process': True, 'pid': os.getpid(),
        'checks': {}, 'reason': 'No accepted steady solution; transient fields are not promoted.'})
    raise SystemExit(0)

accepted = load('accepted_solution')
policy = json.loads((CONFIG/'policy.json').read_text())
measure = SolutionMeasurements(ROOT/'outputs/sv1/SV_MESH/mesh_arrays.npz',
                               accepted['Q_target_m3_s'], policy['Umean_m_s'])
u, p = measure.read(ROOT/accepted['path'])
actual = measure.measure(u, p)
# Unit-aware absolute tolerances prevent small SI flows from passing a default
# 1e-8 absolute comparison. Relative tolerance is applied to nonzero observables.
checks = {'file_sha256': sha256(ROOT/accepted['path']) == accepted['sha256']}
differences = {}
scales = {
    'velocity_L2': max(accepted['velocity_L2'], np.finfo(float).tiny),
    'pressure_range_pa': max(abs(v) for v in accepted['pressure_range_pa']),
    'Q_in_m3_s': measure.Q, 'Q_out_total_m3_s': measure.Q,
    'epsilon_Q': 1., 'epsilon_mass': 1.,
    'wall_velocity_max_m_s': policy['Umean_m_s'],
    'wall_velocity_P95_m_s': policy['Umean_m_s'],
}
def compare(name, before, after, scale):
    a, b = np.asarray(before), np.asarray(after)
    checks[name] = bool(np.allclose(a, b, rtol=1e-12, atol=1e-14*scale))
    differences[name] = float(np.max(np.abs(a-b)))
for name, scale in scales.items():
    compare(name, accepted[name], actual[name], scale)
for name, scale in [('outlet_flows_m3_s', measure.Q), ('outlet_fractions', 1.),
                    ('area_average_pressure_pa', scales['pressure_range_pa'])]:
    for role, value in accepted[name].items():
        compare(name+'.'+role, value, actual[name][role], scale)
for name in ('velocity_finite', 'pressure_finite', 'wall_noslip_pass'):
    checks[name] = bool(actual[name] and accepted[name])
status = 'PASS' if all(checks.values()) else 'FAIL'
write_json(REPORT/'solution_reload.json', {
    'status': status, 'fresh_process': True, 'pid': os.getpid(), 'timestamp': now(),
    'native_solver_already_exited': True, 'file': accepted['path'],
    'sha256': sha256(ROOT/accepted['path']), 'checks': checks,
    'absolute_differences': differences, 'relative_tolerance': 1e-12,
    'absolute_tolerance': '1e-14 times the recorded unit-aware scale',
    'absolute_tolerance_scales': scales, 'measurements': actual})
print('SV1.2 independent fresh-process reload:', status, flush=True)
raise SystemExit(0 if status == 'PASS' else 1)
