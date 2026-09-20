#!/usr/bin/env python3
"""Run one unmodified official 3D Newtonian fluid case with the native executable."""
import hashlib
import json
import os
import shutil
import sys
import urllib.request
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from sv_validation.execution import run_logged
from sv_validation.provenance import sha256, write_json
from sv_validation.validation import parse_result_vtu

native = json.loads((ROOT / 'reports/sv1/native_solver.json').read_text())
assert native['status'] == 'PASS'
commit = native['commit']
relative = Path('tests/cases/fluid/newtonian')
source = ROOT / 'external/svMultiPhysics' / relative
case = ROOT / 'outputs/sv1/official_fluid_smoke'
case.mkdir(parents=True, exist_ok=True)
assert not list(case.glob('*-procs/*.vtu')), 'Do not reuse existing solution files as a fresh smoke result'
inputs = [source / 'solver.xml', source / 'lumen_inlet.flow']
inputs += sorted(p for p in (source / 'mesh').rglob('*') if p.is_file())
records = []
for original in inputs:
    name = original.relative_to(source)
    target = case / name
    target.parent.mkdir(parents=True, exist_ok=True)
    data = original.read_bytes()
    url = None
    if data.startswith(b'version https://git-lfs.github.com/spec/v1'):
        lines = data.decode().splitlines()
        expected = next(line.split('sha256:', 1)[1] for line in lines if line.startswith('oid '))
        size = int(next(line.split()[1] for line in lines if line.startswith('size ')))
        url = 'https://media.githubusercontent.com/media/SimVascular/svMultiPhysics/' + commit + '/' + str(relative / name)
        print('Fetching official case data:', name, flush=True)
        with urllib.request.urlopen(url, timeout=180) as response:
            data = response.read()
        assert len(data) == size and hashlib.sha256(data).hexdigest() == expected
    target.write_bytes(data)
    records.append({'path': str(target.relative_to(ROOT)), 'sha256': sha256(target),
                    'official_relative_path': str(relative / name), 'url': url})
write_json(ROOT / 'reports/sv1/official_smoke_inputs.json', {'commit': commit, 'files': records})
env = dict(os.environ, PATH='/usr/bin:/bin:/usr/local/bin',
           LD_LIBRARY_PATH=native['runtime_library_path'], OMP_NUM_THREADS='1')
run = run_logged(['/usr/bin/mpiexec', '--bind-to', 'none', '-n', '4', native['executable'], 'solver.xml'],
                 'official_fluid_smoke', cwd=case, env=env, timeout=1800)
files = sorted(case.glob('4-procs/result_*.vtu'))
result = dict(run, status='FAIL', case='fluid/newtonian', commit=commit, mpi_ranks=4,
              solver_xml_unmodified=True, result_files=[], velocity_finite=False, pressure_finite=False)
for path in files:
    mesh, u, p = parse_result_vtu(path)
    result['result_files'].append({'path': str(path.relative_to(ROOT)), 'sha256': sha256(path),
                                  'points': mesh.n_points, 'cells': mesh.n_cells,
                                  'velocity_max': float(np.linalg.norm(u, axis=1).max()),
                                  'pressure_range': [float(p.min()), float(p.max())]})
if run['exit_code'] == 0 and files:
    result.update(status='PASS', velocity_finite=True, pressure_finite=True)
write_json(ROOT / 'reports/sv1/official_smoke.json', result)
print(json.dumps(result, indent=2))
raise SystemExit(0 if result['status'] == 'PASS' else 1)
