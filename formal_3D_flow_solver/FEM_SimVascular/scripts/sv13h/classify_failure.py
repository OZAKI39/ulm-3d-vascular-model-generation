"""Publish measured failure and explicit unexecuted dependants, preserving raw evidence."""
import json, re, shutil, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'src'))
from sv_validation.provenance import sha256, write_json, now
R = ROOT / 'reports/sv1_3h'
remote = R / 'remote'
for p in remote.glob('*.json'):
    shutil.copyfile(p, R / p.name)
def load(name): return json.loads((R / (name + '.json')).read_text())
build = load('petsc_cuda12_build')
audit = load('petsc_failure_audit')
assert build['status'] == 'FAIL' and build['make_exit'] != 0 and build['attempt'] == 2
assert audit['status'] == 'PASS' and audit['source_unmodified']
make = build['steps'][-1]
log = ROOT / 'logs/sv1_3h/remote' / Path(make['log']).name
txt = log.read_text()
assert sha256(log) == make['sha256']
errors = [{'line': n, 'text': x} for n, x in enumerate(txt.splitlines(), 1) if 'error:' in x]
assert any('tuple' in x['text'] and 'has no member "get"' in x['text'] for x in errors)
old_patterns = ('clockRate', 'memoryClockRate', 'unary_function')
repeated = [x for x in errors if any(p in x['text'] for p in old_patterns)]
variables = (remote / 'build_failure_evidence/petscvariables').read_text()
compiler_flags = {key: next(x.split(' = ', 1)[1] for x in variables.splitlines() if x.startswith(key + ' = '))
                  for key in ('CC', 'CXX', 'CUDAC', 'CXX_FLAGS', 'CUDAC_FLAGS', 'CUDA_INCLUDE', 'CUDA_LIB')}
assert all('-std=c++17' in compiler_flags[k] for k in ('CXX_FLAGS', 'CUDAC_FLAGS'))
commands = [x for x in txt.splitlines() if compiler_flags['CUDAC'] in x and ' -c ' in x and 'cupm.cu' in x]
assert commands, 'Missing full compiler invocation'
write_json(R / 'compiler_commands.json', {'commands': commands, 'flags': compiler_flags, 'log_sha256': sha256(log)})
reason = 'UNEXPECTED_CUDA12_API_FAILURE' if repeated else 'CUDA12_THRUST_TUPLE_API_INCOMPATIBILITY'
build.update(reason='PETSC_CUDA12_BUILD_FAIL', failure_class=reason, source_unmodified=True,
             source_files_verified=audit['source_files_verified'], source_patch_applied=False, PETSc_upgraded=False,
             compiler=compiler_flags, library_sha256=None, self_test='NOT_RUN', install='NOT_RUN',
             build_time_s=sum(s['wall_time_s'] for s in build['steps']),
             attempted_build_time_s=sum(s['wall_time_s'] for s in load('petsc_cuda12_build_attempt1')['steps']) + sum(s['wall_time_s'] for s in build['steps']))
write_json(R / 'petsc_cuda12_build.json', build)
write_json(R / 'failure_classification.json', {
    'status': 'FAIL', 'reason': 'PETSC_CUDA12_BUILD_FAIL', 'failure_class': reason,
    'previous_API_errors_recurred': bool(repeated), 'compiler_errors': errors,
    'error_log': str(log.relative_to(ROOT)), 'error_log_sha256': sha256(log),
    'upstream_explanation': 'PETSc has recorded the same tuple member-access incompatibility since CUDA 12.4. No patch was applied.',
    'official_sources': [
        {'url': 'https://gitlab.com/petsc/petsc/-/merge_requests/7354', 'supports': 'PETSc tuple.get incompatibility with CUDA 12.4'},
        {'url': 'https://github.com/NVIDIA/cccl/releases/tag/v2.3.0', 'supports': 'Thrust tuple implementation transition'},
        {'url': 'https://developer.nvidia.com/cuda-12-6-3-download-archive', 'supports': 'Official CUDA 12.6.3 release'},
        {'url': 'https://petsc.org/release/changes/319/', 'supports': 'CUDA 12 support claim does not guarantee compatibility with later minor releases'}]})
names = ['petsc_self_test', 'petsc_cuda_types', 'petsc_gpu_smoke', 'svmp_gpu_build', 'svmp_cuda_link',
         'svmp_gpu_smoke', 'gpu_proof', 'cpu_proof', 'science_equivalence', 'gpu_residency',
         'gpu_transfer', 'gpu_memory', 'benchmark', 'speedup']
for name in names:
    d = {'status': 'NOT_RUN', 'executed': False, 'reason': 'PETSC_CUDA12_BUILD_FAIL: prerequisite build gate failed', 'measurements': None}
    if name in ('petsc_cuda_types', 'svmp_gpu_smoke'): d.update(mat_type=None, vec_type=None, KSP=None, PC=None)
    if name in ('gpu_proof', 'cpu_proof'): d.update(steps=None, initial_state='t=0', linear_failures=None, nonlinear_failures=None, velocity_finite=None, pressure_finite=None, mass_error=None, reload_pass=None)
    if name == 'science_equivalence': d.update(velocity_relative_L2=None, pressure_relative_L2=None, Qin_relative=None, Qout_relative=None, mass_error_difference=None, max_velocity_relative=None)
    if name == 'gpu_residency': d.update(measured=False, matrix_resident=None, vectors_resident=None, pc_apply_device=None)
    if name == 'gpu_transfer': d.update(measured=False, H2D_count=None, D2H_count=None, H2D_bytes=None, D2H_bytes=None, large_transfer_per_iteration=None)
    if name == 'gpu_memory': d.update(measured=False, peak_bytes=None, device_capacity_MiB=24564, device_capacity_is_solver_peak=False)
    if name in ('benchmark', 'speedup'): d.update(CPU_1R=None, CPU_4R=None, GPU_1R=None, speedup=None)
    write_json(R / (name + '.json'), d)
write_json(ROOT / 'benchmarks/sv1_3h/status.json', load('benchmark'))
write_json(ROOT / 'outputs/sv1_3h/status.json', {'native_CUDA_kernel': 'PASS', 'PETSc_and_vascular_fields': 'NOT_RUN', 'reason': 'PETSC_CUDA12_BUILD_FAIL'})
write_json(R / 'native_artifact_mirror.json', {'files': [{'path': p, 'sha256': sha256(ROOT / p), 'size': (ROOT / p).stat().st_size} for p in ('outputs/sv1_3h/native/cuda_kernel_smoke', 'scripts/use_cuda12_gpu_env.sh', 'benchmarks/sv1_3h/cuda_kernel_smoke.cu')], 'scope': 'Native remote kernel and wrapper; no CUDA PETSc library produced. Official 4.45 GB installer not mirrored; verified URL and hashes retained.'})
write_json(R / 'stage_result.json', {'timestamp': now(), 'status': 'FAIL', 'reason': 'PETSC_CUDA12_BUILD_FAIL',
    'classification': 'BLOCKED', 'CUDA12_runtime': 'PASS', 'MPI': 'PASS', 'PETSc_configure': 'PASS', 'PETSc_build': 'FAIL',
    'GPU_runtime': 'NOT_RUN', 'production': 'CPU_EARLY_STOP_PRODUCTION', 'production_changed': False,
    'mesh_convergence_started': False, 'GPU_full_production_started': False, 'recommend_SV1_3I': False})
print('FAIL: PETSC_CUDA12_BUILD_FAIL; CUDA runtime PASS; MPI PASS; downstream NOT_RUN')
