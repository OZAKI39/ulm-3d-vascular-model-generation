"""Preserve the failed build and verify original sources; do not retry or patch."""
import json, shutil
from pathlib import Path
from runner_remote import BASE, R, L, load, write, digest
from environment_remote import snapshot

build = load('petsc_cuda12_build')
assert build['status'] == 'FAIL' and build['make_exit'] != 0
source = Path(build['PETSC_DIR'])
original = load('petsc_original_source_files')['files']
changed = [name for name, sha in original.items()
           if not (source / name).is_file() or digest(source / name) != sha]
evidence = R / 'build_failure_evidence'
evidence.mkdir(exist_ok=True)
paths = [source / build['PETSC_ARCH'] / 'lib/petsc/conf/make.log',
         source / build['PETSC_ARCH'] / 'lib/petsc/conf/configure.log',
         source / build['PETSC_ARCH'] / 'lib/petsc/conf/petscvariables',
         source / build['PETSC_ARCH'] / 'include/petscconf.h',
         source / 'src/vec/vec/impls/seq/cupm/vecseqcupm.hpp',
         Path(build['CUDA_prefix']) / 'include/thrust/tuple.h',
         Path(build['CUDA_prefix']) / 'include/thrust/version.h']
mirrored = []
for p in paths:
    assert p.is_file(), str(p)
    dest = evidence / p.name
    shutil.copyfile(p, dest)
    mirrored.append({'remote_path': str(p), 'evidence': str(dest.relative_to(BASE)),
                     'sha256': digest(p), 'size': p.stat().st_size})
prior = json.loads((BASE / 'configs/mpi_native_artifacts.json').read_text())
mpi_changes = [a['path'] for a in prior if digest(a['realpath']) != a['sha256']]
resolution = json.loads((BASE / 'configs/mpi_resolution.json').read_text())
mpi_wrapper_unchanged = digest(resolution['working_launcher']) == resolution['wrapper_sha256']
write('petsc_failure_audit', {
    'status': 'PASS' if not changed and not mpi_changes and mpi_wrapper_unchanged else 'FAIL',
    'source_unmodified': not changed, 'source_files_verified': len(original),
    'changed_original_files': changed, 'mirrored': mirrored,
    'MPI_changes': mpi_changes, 'MPI_wrapper_unchanged': mpi_wrapper_unchanged,
    'MPI_artifacts_verified': len(prior), 'source_patch_applied': False,
    'PETSc_upgraded': False, 'usable_library_exists': (Path(build['prefix']) / 'lib/libpetsc.so').exists()})
after = snapshot()
write('final_environment', after)
before = load('pre_install_environment')
keys = ('driver_files', 'cuda_default_link', 'cuda13_prefix', 'cuda13_files', 'ld_configuration')
checks = {key: before[key] == after[key] for key in keys}
write('final_cuda13_preservation', {'status': 'PASS' if all(checks.values()) else 'FAIL', 'checks': checks})
print(json.dumps({'source_files_verified': len(original), 'changed_sources': changed,
                  'MPI_changes': mpi_changes, 'final_environment_checks': checks}), flush=True)
