#!/usr/bin/env python3
"""Run a project script with the already available official embedded Python."""
import argparse
import json
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'mesh_generate/src'))
from vascular_validation.execution import run_logged, meshing_environment
from vascular_validation.provenance import write_json, sha256

parser = argparse.ArgumentParser()
parser.add_argument('script', nargs='?', default='probe_meshing_api.py')
parser.add_argument('--label', default='simvascular_api')
args = parser.parse_args()
distribution = json.loads((ROOT / 'reports/mesh_and_flow/simvascular_distribution.json').read_text())
launcher = ROOT / distribution['path'] / 'simvascular'
script = Path(__file__).resolve().parent / args.script
result = run_logged([launcher, '-python', '--', script], args.label, env=meshing_environment())
if args.script == 'probe_meshing_api.py':
    text = (ROOT / result['log']).read_text()
    lines = [line.split('=', 1)[1] for line in text.splitlines() if line.startswith('MESH_API_PROBE_JSON=')]
    data = json.loads(lines[-1]) if lines else {}
    result.update(status='PASS' if result['exit_code'] == 0 and data.get('meshing_api_present') else 'FAIL',
                  result=data, launcher=str(launcher), launcher_sha256=sha256(launcher),
                  package_version='2023.05.31')
    write_json(ROOT / 'reports/mesh_and_flow/simvascular_api.json', result)
    print(json.dumps(result, indent=2))
raise SystemExit(result['exit_code'] if result['exit_code'] is not None else 1)
