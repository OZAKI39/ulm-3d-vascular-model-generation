"""Correct CUDA dialect using supported configure option in a fresh PETSC_ARCH."""
import shutil
from pathlib import Path
from runner_remote import *

old = load('petsc_cuda12_build')
assert old['status'] == 'FAIL'
assert load('petsc_failure_audit')['source_unmodified']
assert not (R / 'petsc_cuda12_build_attempt1.json').exists()
shutil.copyfile(R / 'petsc_cuda12_build.json', R / 'petsc_cuda12_build_attempt1.json')
shutil.copyfile(R / 'petsc_failure_audit.json', R / 'petsc_failure_audit_attempt1.json')
shutil.copytree(R / 'build_failure_evidence', R / 'build_failure_evidence_attempt1')
source = Path(old['PETSC_DIR'])
arch = 'arch-sv13h-cuda12-cxx17'
assert not (source / arch).exists()
command = [('PETSC_ARCH=' + arch) if x.startswith('PETSC_ARCH=') else x for x in old['configure_command']]
command.append('--with-cuda-dialect=C++17')
d = dict(old, status='BUILDING', attempt=2, PETSC_ARCH=arch, configure_command=command,
         steps=[], configure_exit=None, make_exit=None, cuda_enabled=None, source_unmodified=None,
         reason=None, failed_step=None,
         correction='Host C++17 was insufficient: set CUDA dialect explicitly to C++17; fresh architecture; original sources unchanged')
write('petsc_cuda12_build', d)
def step(args, name):
    r = run(args, name, cwd=source, timeout=3600, extra_env={'PETSC_OPTIONS': '-use_gpu_aware_mpi 0'})
    d['steps'].append(r)
    if 'configure' in name: d['configure_exit'] = r['exit_code']
    if 'make' in name: d['make_exit'] = r['exit_code']
    if not okay(r):
        d.update(status='FAIL', reason='PETSC_CUDA12_BUILD_FAIL', failed_step=name)
        write('petsc_cuda12_build', d)
        print(text(r)[-3500:], flush=True)
        raise SystemExit(1)
    write('petsc_cuda12_build', d)
step(command, 'petsc_cuda12_cxx17_configure')
variables = (source / arch / 'lib/petsc/conf/petscvariables').read_text()
for key in ('CUDAC_FLAGS', 'CXX_FLAGS'):
    line = next(x for x in variables.splitlines() if x.startswith(key + ' ='))
    assert '-std=c++17' in line and '-std=c++20' not in line
assert '#define PETSC_HAVE_CUDA 1' in (source / arch / 'include/petscconf.h').read_text()
d['cuda_enabled'] = True
step(['make', '-j8', 'V=1', 'PETSC_DIR=' + str(source), 'PETSC_ARCH=' + arch, 'all'], 'petsc_cuda12_cxx17_make')
step(['make', 'PETSC_DIR=' + str(source), 'PETSC_ARCH=' + arch, 'install'], 'petsc_cuda12_cxx17_install')
step(['make', 'PETSC_DIR=' + d['prefix'], 'PETSC_ARCH=', 'check'], 'petsc_cuda12_cxx17_check')
d.update(status='PASS', self_test='PASS', library_sha256=digest(Path(d['prefix']) / 'lib/libpetsc.so'))
write('petsc_cuda12_build', d)
