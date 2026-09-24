#!/usr/bin/env python3
"""Choose among actual safe repairs and the unmodified Stage 1.7 original."""
import json
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from fem3d.audit import sha256, timestamp, write_json
from fem3d.residual_repair import acceptance, selection

I, O = ROOT / 'inputs/stage01_8', ROOT / 'outputs/stage01_8'
read = lambda p: json.loads(p.read_text())
policy = read(I / 'acceptance_policy.json'); lock = read(I / 'freeze_lock.json')
assert sha256(I / 'acceptance_policy.json') == lock['policy_sha256']
assert not (O / 'selection.json').exists(), 'Selection evidence already exists'
baseline = read(O / 'baseline/baseline_recomputed.json')
records = [read(O / name / 'qc/result.json') for name in policy['repair_methods'] if (O / name / 'qc/result.json').exists()]
assert records and records[0]['name'] == 'repair_00'
last = records[-1]
if last['name'] == 'repair_00':
    assert last['acceptance'].get('low_cells_eliminated') or last['acceptance'].get('significant_q_min_improvement'), 'repair_01 assessment still required'
if last['name'] == 'repair_01':
    assert last['acceptance'].get('PASS_A') or last.get('boundary_lock', {}).get('status') != 'PASS', 'repair_02 assessment still required'
chosen = selection(records, baseline, policy)
if chosen == 'KEEP_STAGE017':
    source = I / 'stage017_selected'
    measured = baseline
    residuals = read(I / 'residual_cells.json')
    reason = 'The original is at least equivalent under the frozen low-count, minimum-quality and cost ordering; preserve accepted evidence.'
else:
    source = O / chosen
    record = next(r for r in records if r['name'] == chosen)
    measured, residuals = record['volume'], record['residual_diagnostics']
    reason = 'Best measured hard-safety-passing mesh by minimum residual count, maximum q_min, minimum P2 proxy, minimum tetra count.'
dest = O / 'chosen'
if dest.exists():
    raise RuntimeError('Chosen directory already exists')
(dest / 'mesh').mkdir(parents=True)
for name in ('volume_mesh.npz', 'fluid.msh'):
    shutil.copyfile(source / 'mesh' / name, dest / 'mesh' / name)
(dest / 'surface').mkdir()
shutil.copyfile(I / 'stage017_selected/surface/tagged_surface_si.npz', dest / 'surface/tagged_surface_si.npz')
shutil.copyfile(I / 'stage017_selected/planar_port_contract_v2.json', dest / 'planar_port_contract_v2.json')
write_json(dest / 'qc/volume_quality.json', measured)
write_json(dest / 'qc/residual_cells.json', residuals)
acc = acceptance(measured, baseline, policy)
assert acc['status'] == 'PASS'
causes = {c['reason']['category'] for c in residuals['cells']}
pass_b = bool(acc['significant_q_min_improvement'] and causes <= {'WALL_TRIANGLE_CONSTRAINED', 'RAPID_GEOMETRIC_VARIATION'} and last['name'] == 'repair_02')
assessment = 'PASS-A' if acc['PASS_A'] else 'PASS-B' if pass_b else 'SAFE_KEEP' if chosen == 'KEEP_STAGE017' else 'SAFE_MEASURED_IMPROVEMENT'
result = {'timestamp': timestamp(), 'chosen': chosen, 'recommendation': chosen if chosen == 'KEEP_STAGE017' else 'USE_' + chosen.upper(),
          'reason': reason, 'repair_success_class': assessment, 'acceptance': acc,
          'chosen_mesh_sha256': sha256(dest / 'mesh/volume_mesh.npz'), 'source_mesh_sha256': sha256(source / 'mesh/volume_mesh.npz'),
          'source_mesh_path': 'outputs/stage01_7/selected/mesh/volume_mesh.npz' if chosen == 'KEEP_STAGE017' else str((source / 'mesh/volume_mesh.npz').relative_to(ROOT)),
          'cap_surface_bytes_identical': sha256(dest / 'surface/tagged_surface_si.npz') == sha256(I / 'stage017_selected/surface/tagged_surface_si.npz'),
          'records': records, 'not_run': [n for n in policy['repair_methods'] if n not in {r['name'] for r in records}],
          'policy_sha256': lock['policy_sha256'], 'stage3_started': False, 'fem_solved': False}
write_json(O / 'selection.json', result); write_json(dest / 'metadata/selection.json', {k: v for k, v in result.items() if k != 'records'})
print(json.dumps({k: v for k, v in result.items() if k != 'records'}, indent=2))
