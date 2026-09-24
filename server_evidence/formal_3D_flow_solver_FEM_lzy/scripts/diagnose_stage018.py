#!/usr/bin/env python3
"""Complete physical diagnostics before any repair; WSL-only assembly."""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from fem3d.audit import sha256, timestamp, write_json
from fem3d.residual_geometry import diagnose, tetra_shape

inp = ROOT / 'inputs/stage01_8'
out = ROOT / 'outputs/stage01_8/diagnostics'
if (out / 'residual_cells.json').exists():
    raise RuntimeError('Diagnosis already exists')
policy = json.loads((inp / 'acceptance_policy.json').read_text())
lock = json.loads((inp / 'freeze_lock.json').read_text())
assert sha256(inp / 'acceptance_policy.json') == lock['policy_sha256']
baseline = json.loads((ROOT / 'outputs/stage01_8/baseline/baseline_recomputed.json').read_text())
assert baseline['status'] == 'PASS'
data = dict(np.load(inp / 'stage017_selected/mesh/volume_mesh.npz'))
assert sha256(inp / 'stage017_selected/mesh/volume_mesh.npz') == lock['sole_input_mesh_sha256']
data['min_sicn'] = np.load(ROOT / 'outputs/stage01_8/baseline/min_sicn_recomputed.npy')
result = diagnose(data, policy)
assert result['residual_count'] == baseline['quality']['total_below_0_1']
result.update(timestamp=timestamp(), policy_sha256=lock['policy_sha256'], mesh_sha256=lock['sole_input_mesh_sha256'], repairs_started=False)
write_json(out / 'residual_cells.json', result)
# Reproducible healthy comparison sample, spread over the actual tetra ordering.
ordinary = np.flatnonzero(data['min_sicn'] >= baseline['quality']['min_sicn']['median'])
sample = ordinary[np.unique(np.linspace(0, len(ordinary) - 1, min(256, len(ordinary))).astype(int))]
records = [{'cell_index': int(i), 'min_sicn': float(data['min_sicn'][i]), **tetra_shape(data['points_m'][data['tetra'][i]])} for i in sample]
write_json(out / 'ordinary_tetra_sample.json', {'method': '256 deterministic equally spaced indices among cells at or above baseline median quality; descriptive comparison, not random sampling', 'cells': records})
write_json(ROOT / 'reports/stage01_8/baseline_recomputed.json', baseline)
print(json.dumps({'residual_count': result['residual_count'], 'cells': [{
    'cell_index': c['cell_index'], 'q': c['min_sicn'], 'reason': c['reason'],
    'frozen_vertices': c['boundary_constraints']['boundary_vertex_count'],
    'dihedral': [c['minimum_dihedral_degrees'], c['maximum_dihedral_degrees']],
    'wall_q': c['wall_diagnostics']['nearest_wall_triangle']['q_tri'],
    'normal_variation': c['wall_diagnostics']['normal_variation']['maximum_degrees'],
    'rings': c['neighborhoods']} for c in result['cells']]}, indent=2))
