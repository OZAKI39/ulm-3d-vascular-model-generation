#!/usr/bin/env python3
"""Read the final result in a fresh process and recompute its physical integrals."""
import json
import os
import sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from sv_validation.postprocess import SolutionMeasurements
from sv_validation.provenance import write_json, sha256

R = ROOT/'reports/sv1'
expected = json.loads((R/'flow_qc.json').read_text())
policy = json.loads((ROOT/'configs/time_policy.json').read_text())
path = ROOT/expected['path']
assert sha256(path) == expected['sha256']
measure = SolutionMeasurements(ROOT/'outputs/sv1/SV_MESH/mesh_arrays.npz', expected['Q_target_m3_s'], policy['Umean_m_s'])
u, p = measure.read(path)
actual = measure.measure(u, p)
Q = expected['Q_target_m3_s']
flow_difference = max(abs(actual['signed_outward_boundary_flows_m3_s'][name]-expected['signed_outward_boundary_flows_m3_s'][name])/Q
                      for name in actual['signed_outward_boundary_flows_m3_s'])
norm_difference = abs(actual['velocity_L2']-expected['velocity_L2'])/expected['velocity_L2']
pressure_equal = np.allclose(actual['pressure_range_pa'], expected['pressure_range_pa'], rtol=1e-12, atol=1e-14)
passed = flow_difference <= 1e-12 and norm_difference <= 1e-12 and pressure_equal
write_json(R/'solution_reload.json', {'status': 'PASS' if passed else 'FAIL', 'fresh_process': True,
           'process_id': os.getpid(), 'artifact_kind': expected['artifact_kind'],
           'physical_acceptance': expected['status'],
           'result_sha256': sha256(path), 'measurements': actual,
           'flow_difference_over_Q': flow_difference, 'velocity_norm_relative_difference': norm_difference,
           'pressure_range_matches': bool(pressure_equal)})
print('Solution reload:', 'PASS' if passed else 'FAIL')
raise SystemExit(0 if passed else 1)
