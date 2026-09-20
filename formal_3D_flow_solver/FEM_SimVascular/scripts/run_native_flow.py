#!/usr/bin/env python3
"""Run the fixed native vascular case; at most one checkpoint extension is allowed."""
import argparse
import json
import os
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from sv_validation.execution import run_logged
from sv_validation.provenance import write_json, sha256

parser = argparse.ArgumentParser()
parser.add_argument('--extend', action='store_true')
args = parser.parse_args()
R = ROOT / 'reports/sv1'
native = json.loads((R / 'native_solver.json').read_text())
time_policy = json.loads((ROOT / 'configs/time_policy.json').read_text())
frozen = json.loads((R / 'flow_config_provenance.json').read_text())
assert sha256(ROOT / 'configs/time_policy.json') == frozen['time_policy_sha256']
assert sha256(ROOT / 'configs/sv_reference.yaml') == frozen['reference_physics_sha256']
assert sha256(ROOT / 'configs/sv_flow.xml') == frozen['xml_sha256']
for name in ('official_smoke', 'mesh_validity', 'inlet_normalization'):
    assert json.loads((R / (name + '.json')).read_text())['status'] == 'PASS'
case = ROOT / 'outputs/sv1/vascular_flow'
previous = json.loads((R / 'flow_execution.json').read_text()) if (R / 'flow_execution.json').exists() else {'runs': []}
if args.extend:
    assert len(previous['runs']) == 1 and previous['runs'][0]['exit_code'] == 0
    steady = json.loads((R / 'steady_state.json').read_text())
    qc = json.loads((R / 'flow_qc.json').read_text())
    assert steady['status'] == 'NOT_REACHED' and steady['scheduled_block_complete']
    assert qc['solver_success'], 'A failed numerical solve must stop; extension is only for a completed, converged time block'
    assert not (case / 'STOP_SIM').exists(), 'An explicit failure stop cannot be resumed automatically'
    assert (case / '4-procs/stFile_last.bin').exists()
    tree = ET.parse(ROOT / 'configs/sv_flow.xml')
    tree.find('.//Continue_previous_simulation').text = 'true'
    tree.find('.//Number_of_time_steps').text = str(time_policy['first_block_steps'] + time_policy['extension_steps'])
    input_name = 'solver_extension.xml'
    tree.write(case / input_name, encoding='utf-8', xml_declaration=True)
    label = 'vascular_flow_extension'
else:
    assert not previous['runs'] and not list(case.glob('*-procs/result_*.vtu')), 'Do not reuse prior solver results'
    assert sha256(case / 'solver.xml') == frozen['xml_sha256']
    input_name = 'solver.xml'
    label = 'vascular_flow_first_block'
env = dict(os.environ, PATH='/usr/bin:/bin:/usr/local/bin',
           LD_LIBRARY_PATH=native['runtime_library_path'], OMP_NUM_THREADS='1')
result = run_logged(['/usr/bin/mpiexec', '--bind-to', 'none', '-n', '4', native['executable'], input_name],
                    label, cwd=case, env=env, timeout=3600)
result['solver_input_sha256'] = sha256(case / input_name)
result['mpi_ranks'] = 4
result['omp_num_threads'] = 1
result['extension'] = args.extend
previous['runs'].append(result)
previous['status'] = 'EXECUTED' if result['exit_code'] == 0 else 'FAIL'
previous['mode'] = 'transient-to-steady'
previous['solver_commit'] = native['commit']
previous['executable_sha256'] = native['executable_sha256']
previous['extensions_used'] = int(args.extend)
write_json(R / 'flow_execution.json', previous)
print(json.dumps(result, indent=2))
raise SystemExit(0 if result['exit_code'] == 0 else 1)
