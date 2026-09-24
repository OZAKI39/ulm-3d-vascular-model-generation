#!/usr/bin/env python3
"""One predeclared repair attempt per invocation; no retries or new methods."""
import argparse
import json
import shutil
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from fem3d.audit import sha256, timestamp, write_json
from fem3d.adaptive_qc import measure_volume
from fem3d.adaptive_volume import dolfinx_loadability
from fem3d.residual_geometry import diagnose
from fem3d.residual_meshing import run_mesher
from fem3d.residual_repair import acceptance

ap = argparse.ArgumentParser(); ap.add_argument('name', choices=['repair_00', 'repair_01', 'repair_02'])
ap.add_argument('--recover-validation-error', action='store_true')
args = ap.parse_args(); name = args.name
I, O = ROOT / 'inputs/stage01_8', ROOT / 'outputs/stage01_8'
read = lambda p: json.loads(p.read_text())
policy = read(I / 'acceptance_policy.json'); lock = read(I / 'freeze_lock.json')
assert sha256(I / 'acceptance_policy.json') == lock['policy_sha256']
diagnosis = read(I / 'diagnostic_summary.json'); residuals = read(I / 'residual_cells.json')
assert diagnosis['status'] == 'DIAGNOSIS_COMPLETE' and not diagnosis['repairs_started']
assert sha256(I / 'residual_cells.json') == diagnosis['residual_cells_sha256']
assert diagnosis['policy_sha256'] == lock['policy_sha256']
baseline = read(O / 'baseline/baseline_recomputed.json')
base = I / 'stage017_selected'; input_msh = base / 'mesh/fluid.msh'
if name == 'repair_01':
    previous = read(O / 'repair_00/qc/result.json')
    acc = previous['acceptance']
    assert not (acc.get('low_cells_eliminated') or acc.get('significant_q_min_improvement')), 'repair_00 already improves/eliminates residuals'
if name == 'repair_02':
    previous = read(O / 'repair_01/qc/result.json')
    assert not previous['acceptance'].get('PASS_A'), 'repair_01 already PASS-A'
    assert previous.get('boundary_lock', {}).get('status') == 'PASS', 'Cannot optimize a boundary-unsafe predecessor'
    input_msh = O / 'repair_01/mesh/fluid.msh'
out = O / name
if args.recover_validation_error:
    previous = read(out / 'qc/result.json')
    known_errors = {'repair_00': "ValueError('Frozen boundary vertices changed: maximum displacement=0.0')",
                    'repair_01': "ValueError('Boundary vertex IDs/coordinates changed')"}
    assert name in known_errors and previous.get('error') == known_errors[name]
    archive = O / ('failed_attempts/' + name + '_boundary_validator')
    assert not archive.exists()
    archive.parent.mkdir(parents=True, exist_ok=True)
    out.rename(archive)
    recovery = {
        'timestamp': timestamp(), 'method': name, 'failed_attempt_path': str(archive.relative_to(ROOT)),
        'reason': ('Boundary checker incorrectly supplied all volume vertices to a surface-only bijection; corrected to compact boundary vertices.' if name == 'repair_00' else
                   'Boundary identity checker conflated Gmsh transient numbering with physical vertex/connectivity identity; corrected to compare canonical physical bytes and independently validate exact tagged boundary.'),
        'scientific_policy_changed': False, 'policy_sha256': lock['policy_sha256'],
        'administrative_policy_deviation': 'One extra execution of the same method for implementation-error recovery, beyond the self-imposed per_method_runs=1; user limit of three distinct repair methods unchanged.',
        'failed_result_counted_as_pass': False, 'all_failed_logs_retained': True}
    if (O / 'technical_recovery.json').exists():
        first = read(O / 'technical_recovery.json'); first.setdefault('additional_recoveries', []).append(recovery)
        write_json(O / 'technical_recovery.json', first)
    else:
        write_json(O / 'technical_recovery.json', recovery)
if out.exists():
    raise RuntimeError('Repair attempt already exists; no overwrite or hidden retry')
out.mkdir()
attempt = {'timestamp': timestamp(), 'name': name, 'policy_sha256': lock['policy_sha256'],
           'diagnostic_summary_sha256': sha256(I / 'diagnostic_summary.json'), 'diagnosis_completed_utc': diagnosis['timestamp'],
           'input_mesh': str(input_msh.relative_to(ROOT)), 'surface_sha256': sha256(base / 'surface/tagged_surface_si.npz')}
write_json(out / 'metadata/attempt.json', attempt)
source = dict(np.load(base / 'surface/tagged_surface_si.npz')); original_mesh = dict(np.load(base / 'mesh/volume_mesh.npz'))
assert sha256(base / 'mesh/volume_mesh.npz') == lock['sole_input_mesh_sha256']
record = {'name': name, 'attempt': attempt}
try:
    data, boundary = run_mesher(name, source, original_mesh, input_msh, out, policy, residuals)
    record['boundary_lock'] = boundary
    contract = read(base / 'planar_port_contract_v2.json')
    original = dict(np.load(ROOT / 'inputs/stage01/tagged_surface_si.npz'))
    measured = measure_volume(data, source, original, contract, policy)
    load = dolfinx_loadability(out / 'mesh/fluid.msh')
    acc = acceptance(measured, baseline, policy, boundary['status'] == 'PASS', load['status'] == 'PASS')
    residual = diagnose(data, policy)
    record.update(volume=measured, acceptance=acc, residual_diagnostics=residual,
                  cap_surface_sha256=sha256(base / 'surface/tagged_surface_si.npz'))
    write_json(out / 'qc/volume_quality.json', measured); write_json(out / 'qc/dolfinx_loadable.json', load)
    write_json(out / 'qc/residual_cells.json', residual)
    (out / 'surface').mkdir(); shutil.copyfile(base / 'surface/tagged_surface_si.npz', out / 'surface/tagged_surface_si.npz')
    assert sha256(out / 'surface/tagged_surface_si.npz') == attempt['surface_sha256']
    shutil.copyfile(base / 'planar_port_contract_v2.json', out / 'planar_port_contract_v2.json')
except Exception as exc:
    record.update(acceptance={'status': 'FAIL', 'PASS_A': False, 'error': repr(exc)}, error=repr(exc))
finally:
    record['finished_utc'] = timestamp(); write_json(out / 'qc/result.json', record)
print(json.dumps({k: v for k, v in record.items() if k not in ('volume', 'residual_diagnostics', 'boundary_lock')}, indent=2))
if 'volume' in record:
    print(json.dumps({'proxy': record['volume']['proxy'], 'quality': record['volume']['quality']['min_sicn'],
                      'low_count': record['volume']['quality']['total_below_0_1']}, indent=2))
