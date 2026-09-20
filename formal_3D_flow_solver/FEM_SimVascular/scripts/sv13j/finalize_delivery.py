"""Verify native mirrors, summarize real test outcomes, and seal J-only evidence."""
import ast
import json
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'src'))
from sv_validation.provenance import now, sha256, write_json

R = ROOT / 'reports/sv1_3j'


def summarize(filename):
    tree = ET.parse(R / filename)
    result = dict(passed=0, failed=0, errors=0, skipped=0, failures=[])
    for case in tree.iter('testcase'):
        identity = case.attrib['classname'] + '::' + case.attrib['name']
        if case.find('failure') is not None:
            result['failed'] += 1
            result['failures'].append(identity)
        elif case.find('error') is not None:
            result['errors'] += 1
            result['failures'].append(identity)
        elif case.find('skipped') is not None:
            result['skipped'] += 1
        else:
            result['passed'] += 1
    return result


stage = summarize('pytest_stage.xml')
full = summarize('pytest_full.xml')
write_json(R / 'test_summary.json', {
    'timestamp': now(), 'stage': stage, 'full_suite': full,
    'stage_test_files': len(list((ROOT / 'tests').glob('test_sv13j_*.py'))),
    'historical_failures': [name for name in full['failures'] if '.test_sv13j_' not in name],
    'actual_stage_failure': 'Official fluid GPU smoke: MPI_ERR_TYPE; retained as a real failing acceptance test.',
    'skip_policy': 'Only unexecuted downstream runtime gates skip; B/C NOT_REQUIRED assertions pass.',
})

native = json.loads((R / 'native_artifacts.json').read_text())
entries = []
for item in native['artifacts']:
    path = ROOT / 'outputs/sv1_3j/native' / Path(item['remote_path']).name
    assert path.stat().st_size == item['bytes'], path
    assert sha256(path) == item['sha256'], path
    entries.append({**item, 'local_path': str(path.relative_to(ROOT)), 'verified': True})
write_json(R / 'native_mirror_verification.json', {
    'timestamp': now(), 'status': 'PASS', 'artifacts': entries,
    'CUDA_installer_mirrored': False,
    'CUDA_installer_recovery': 'Official URL, verified official MD5 and local SHA256 in cuda123_source_manifest.json',
    'portability': 'Native binaries retain the recorded remote library paths; this mirror is an archive, not a WSL installation.',
})

terminal = f'''Stage SV1.3J completed.

toolchain:
    default GCC = 13.3.0 (unchanged)
    compatibility GCC = Ubuntu gcc/g++-12 12.4.0; C++17; sm_89
    MPI = existing project Open MPI 4.1.6 (unchanged)

compatibility matrix:
CUDA 12.3.2:
    Thrust = 2.2.0 (THRUST_VERSION=200200)
    kernel = 3/3 PASS
    configure = PASS
    make = PASS
    error = NONE
CUDA 12.2.2:
    status = NOT_REQUIRED (stop on first success)
    Thrust = NOT RUN
    kernel = NOT RUN
    configure = NOT RUN
    make = NOT RUN
    error = NOT RUN
CUDA 12.1.1:
    status = NOT_REQUIRED (stop on first success)
    Thrust = NOT RUN
    kernel = NOT RUN
    configure = NOT RUN
    make = NOT RUN
    error = NOT RUN

compatibility winner:
    CUDA 12.3.2 / NVCC 12.3.107 / GCC 12.4.0 / unpatched PETSc 3.19.6

PETSc GPU:
    smoke = 3/3 PASS; CPU/MPI/CUDA self-tests PASS
    Mat = seqaijcusparse
    Vec = seqcuda
    KSP = GMRES, reason=2, iterations=21, true relative residual=3.1767506586176323e-13

svMultiPhysics GPU:
    build = PASS (fixed commit, source unchanged)
    smoke = FAIL; exit 3; MPI_Bcast MPI_ERR_TYPE; no VTU
    linked PETSc = sv1_3j/external/compat_cuda/petsc-cuda123/lib/libpetsc.so.3.19.6
    runtime CUDA backend = NOT REACHED

vascular proof:
    GPU steps = NOT RUN
    CPU steps = NOT RUN
    GPU linear failures = NOT RUN
    GPU nonlinear failures = NOT RUN

science equivalence:
    velocity L2 = NOT RUN
    pressure L2 = NOT RUN
    Qin difference = NOT RUN
    max Qout difference = NOT RUN

residency:
    matrix resident = NOT RUN (vascular ASM/ILU2)
    vectors resident = NOT RUN (vascular ASM/ILU2)
    transfer bound = NOT RUN
    peak VRAM = NOT RUN

benchmark:
    CPU1 = NOT RUN
    CPU4 = NOT RUN
    GPU1 = NOT RUN
    speedup = NOT RUN

classification:
    SVMULTIPHYSICS_GPU_FAIL / MPI_FORTRAN_PREDEFINED_DATATYPES_UNAVAILABLE
    performance classification = BLOCKED; no speed claim
    production = CPU_EARLY_STOP_PRODUCTION (unchanged)

tests:
    passed = {stage['passed']}
    failed = {stage['failed']}
    skipped = {stage['skipped']}
    full suite = {full['passed']} passed, {full['failed']} failed, {full['skipped']} skipped, {full['errors']} errors
    full failures = {len(full['failures']) - len(stage['failures'])} historical + {len(stage['failures'])} current-stage acceptance failure

report:
    reports/sv1_3j/REPORT.md

STAGE SV1.3J STATUS:
    FAIL: SVMULTIPHYSICS_GPU_FAIL
'''
(R / 'terminal_summary.txt').write_text(terminal)

paths = set()
for folder in ('reports/sv1_3j', 'outputs/sv1_3j', 'logs/sv1_3j',
               'configs/sv1_3j', 'benchmarks/sv1_3j', 'scripts/sv13j',
               'external/compat_cuda'):
    paths.update(p for p in (ROOT / folder).rglob('*') if p.is_file()
                 and '__pycache__' not in p.parts and 'plot_cache' not in p.parts)
paths.update((ROOT / 'tests').glob('test_sv13j_*.py'))
paths.update(ROOT / p for p in ('tests/sv13j_support.py',
             'src/sv_validation/sv13j.py', 'scripts/use_cuda123_env.sh'))
paths.discard(R / 'delivery_manifest.json')
for path in paths:
    if path.suffix == '.py':
        ast.parse(path.read_text(), filename=str(path))
write_json(R / 'delivery_manifest.json', {
    'timestamp': now(), 'status': 'EVIDENCE_COMPLETE_STAGE_FAILED',
    'stage_status': 'FAIL', 'reason': 'SVMULTIPHYSICS_GPU_FAIL',
    'entries': len(paths),
    'files': [{'path': str(p.relative_to(ROOT)), 'sha256': sha256(p),
               'size': p.stat().st_size} for p in sorted(paths)],
    'scope': 'Only SV1.3J additions; historical files covered by preservation_audit.json',
    'exclusions': ['delivery_manifest.json itself', 'plot_cache', '__pycache__'],
})
report = (R / 'REPORT.md').read_text()
assert len(re.findall(r'^## ', report, re.M)) == 11
for link in re.findall(r'\]\(([^)]+)\)', report):
    if not link.startswith('https://'):
        assert (R / link).exists(), link
for item in json.loads((R / 'delivery_manifest.json').read_text())['files']:
    assert sha256(ROOT / item['path']) == item['sha256']
print(terminal)
print(f'Verified {len(paths)} delivery hashes, 5 native mirrors, 11 report sections and all local links.')
