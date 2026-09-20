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

R = ROOT / 'reports/sv1_3l'


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
    'stage_test_files': len(list((ROOT / 'tests').glob('test_sv13l_*.py'))),
    'historical_failures': [name for name in full['failures'] if '.test_sv13l_' not in name],
    'actual_stage_failure': 'Official fluid GPU smoke: Vector is not ghosted; retained as a real failing acceptance test.',
    'skip_policy': 'Only unexecuted downstream runtime gates skip; observed official failure remains a real failure.',
})

native = json.loads((R / 'native_artifacts.json').read_text())
entries = []
for item in native['artifacts']:
    path = ROOT / 'outputs/sv1_3l/native' / Path(item['remote_path']).name
    assert path.stat().st_size == item['bytes'], path
    assert sha256(path) == item['sha256'], path
    entries.append({**item, 'local_path': str(path.relative_to(ROOT)), 'verified': True})
write_json(R / 'native_mirror_verification.json', {
    'timestamp': now(), 'status': 'PASS', 'artifacts': entries,
    'CUDA_installer_mirrored': False,
    'CUDA_installer_recovery': 'Unchanged Stage J CUDA12.3.2 official source manifest',
    'portability': 'Native binaries retain the recorded remote library paths; this mirror is an archive, not a WSL installation.',
})

terminal = f'''Stage SV1.3L completed.

MPI:
    version = Open MPI 4.1.6 (new gpu_mpi_fortran prefix)
    C compiler = gcc-12 12.4.0
    C++ compiler = g++-12 12.4.0
    Fortran compiler = gfortran-12 12.4.0
    Fortran bindings = all: mpif.h, use mpi, use mpi_f08
    rank1 = 5/5 PASS
    rank2 = PASS
    Fortran rank1/rank2 = PASS
    C++ datatype complete rank1/rank2 repetitions = 3/3 PASS

Fortran datatype:
    MPI_INTEGER size = 4 bytes
    bcast = PASS; content verified
    MPI_DOUBLE_PRECISION size = 8 bytes
    bcast = PASS; content verified
    MPI_CHARACTER size = 1 byte
    bcast = PASS; content verified
    MPI_LOGICAL size = 4 bytes
    bcast = PASS; content verified

PETSc GPU:
    rebuild = PASS; CPU/MPI/CUDA self-tests PASS
    Mat = seqaijcusparse
    Vec = seqcuda
    smoke = 3/3 PASS (GMRES, reason=2, iterations=21)

svMultiPhysics:
    linked MPI = sv1_3l/external/gpu_mpi_fortran/lib/libmpi.so.40
    linked PETSc = sv1_3l/external/compat_cuda/petsc-cuda123-mpif/lib/libpetsc.so.3.19.6
    official GPU smoke = FAIL: Vector is not ghosted
    original MPI_ERR_TYPE = DISAPPEARED
    runtime Mat = seqaijcusparse
    runtime Vec = seqcuda
    runtime KSP/PC = GMRES / ASM overlap2 / ILU(2)
    observed KSP convergence = one record, 18 iterations; does not override prior PETSc error
    first error = VecGhostUpdateBegin, commonmpvec.c:216, petsc_impl.cpp:669
    process = SEGV / MPI_ABORT 59; launcher return -11; error monitor triggered
    VTU count = 0; field/reload validation NOT RUN

vascular GPU proof:
    steps = NOT RUN
    linear failures = NOT RUN
    nonlinear failures = NOT RUN
    velocity finite = NOT RUN
    pressure finite = NOT RUN
    mass = NOT RUN
    reload = NOT RUN

science equivalence:
    velocity L2 = NOT RUN
    pressure L2 = NOT RUN
    Qin difference = NOT RUN
    max Qout difference = NOT RUN

residency:
    matrix resident = NOT RUN (vascular profiling)
    vectors resident = NOT RUN (vascular profiling)
    PC location = NOT RUN (complete execution/transfer audit)
    transfer pattern = NOT RUN
    peak VRAM = NOT RUN

benchmark:
    CPU1 median = NOT RUN
    CPU4 median = NOT RUN
    GPU1 median = NOT RUN
    speedup = NOT RUN

classification:
    FAIL: SVMULTIPHYSICS_GPU_SMOKE_FAIL
    detail: CUDA_VECTOR_GHOST_LAYOUT_INCOMPATIBILITY
    MPI repair PASS; CUDA backend active; performance NOT MEASURED
    production = CPU_EARLY_STOP_PRODUCTION (unchanged)

tests:
    passed = {stage['passed']}
    failed = {stage['failed']}
    skipped = {stage['skipped']}
    full suite = {full['passed']} passed, {full['failed']} failed, {full['skipped']} skipped, {full['errors']} errors
    historical failures = {len(full['failures']) - len(stage['failures'])}

report:
    reports/sv1_3l/REPORT.md

human review:
    mpi_datatype_before_after.png
    gpu_stack_pipeline.png
    petsc_gpu_revalidation.png
    svmultiphysics_gpu_smoke.png
    gpu_cpu_solution_difference.png = NOT RUN; not generated (official smoke gate failed)
    gpu_residency.png = NOT RUN; not generated
    gpu_transfer_cost.png = NOT RUN; not generated
    gpu_memory.png = NOT RUN; not generated
    cpu_gpu_runtime.png = NOT RUN; not generated

STAGE SV1.3L STATUS:
    FAIL
'''
(R / 'terminal_summary.txt').write_text(terminal)

paths = set()
for folder in ('reports/sv1_3l', 'outputs/sv1_3l', 'logs/sv1_3l',
               'configs/sv1_3l', 'benchmarks/sv1_3l', 'scripts/sv13l',
               'external/gpu_mpi_fortran', 'external/compat_cuda/petsc-cuda123-mpif'):
    paths.update(p for p in (ROOT / folder).rglob('*') if p.is_file()
                 and '__pycache__' not in p.parts and 'plot_cache' not in p.parts)
paths.update((ROOT / 'tests').glob('test_sv13l_*.py'))
paths.update(ROOT / p for p in ('tests/sv13l_support.py',
             'src/sv_validation/sv13l.py', 'scripts/run_gpu_mpi_fortran.sh'))
paths.discard(R / 'delivery_manifest.json')
for path in paths:
    if path.suffix == '.py':
        ast.parse(path.read_text(), filename=str(path))
write_json(R / 'delivery_manifest.json', {
    'timestamp': now(), 'status': 'EVIDENCE_COMPLETE_STAGE_FAILED',
    'stage_status': 'FAIL', 'reason': 'SVMULTIPHYSICS_GPU_SMOKE_FAIL',
    'entries': len(paths),
    'files': [{'path': str(p.relative_to(ROOT)), 'sha256': sha256(p),
               'size': p.stat().st_size} for p in sorted(paths)],
    'scope': 'Only SV1.3L additions; historical files covered by preservation_audit.json',
    'exclusions': ['delivery_manifest.json itself', 'plot_cache', '__pycache__'],
})
report = (R / 'REPORT.md').read_text()
assert len(re.findall(r'^## ', report, re.M)) == 12
for link in re.findall(r'\]\(([^)]+)\)', report):
    if not link.startswith('https://'):
        assert (R / link).exists(), link
for item in json.loads((R / 'delivery_manifest.json').read_text())['files']:
    assert sha256(ROOT / item['path']) == item['sha256']
print(terminal)
print(f'Verified {len(paths)} delivery hashes, 7 native mirrors, 12 report sections and all local links.')
