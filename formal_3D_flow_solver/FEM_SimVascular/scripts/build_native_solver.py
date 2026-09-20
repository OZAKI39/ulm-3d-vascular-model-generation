#!/usr/bin/env python3
"""Build only the pinned official fluid solver with the existing native MPI."""
import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from sv_validation.execution import run_logged
from sv_validation.provenance import sha256, write_json

source = ROOT / 'external/svMultiPhysics'
commit = (ROOT / 'external/svMultiPhysics_commit.txt').read_text().strip()
assert subprocess.check_output(['git', '-C', str(source), 'rev-parse', 'HEAD'], text=True).strip() == commit
sv = ROOT / 'external/SimVascularDistribution/usr/local/sv/simvascular/2023-05-31'
vtk = ROOT / 'external/install/vtk'
build = ROOT / 'external/build/native'
env = dict(os.environ, PATH='/usr/bin:/bin:/usr/local/bin',
           LD_LIBRARY_PATH=str(vtk / 'lib'), PKG_CONFIG_PATH='',
           CMAKE_BUILD_PARALLEL_LEVEL='4', OMP_NUM_THREADS='1')
extra = ['-DSV_USE_TETGEN:BOOL=OFF',
         '-DTETGEN_LIBRARY_NAME:FILEPATH=' + str(sv / 'lib/lib_simvascular_thirdparty_tetgen.a'),
         '-DCMAKE_CXX_FLAGS:STRING=-I' + str(sv / 'include/thirdparty/tetgen'),
         '-DBLAS_LIBRARIES:FILEPATH=/lib/x86_64-linux-gnu/libblas.so.3',
         '-DLAPACK_LIBRARIES:FILEPATH=/lib/x86_64-linux-gnu/liblapack.so.3']
command = ['/usr/bin/cmake', '-S', source, '-B', build,
           '-DCMAKE_C_COMPILER=/usr/bin/gcc', '-DCMAKE_CXX_COMPILER=/usr/bin/g++',
           '-DCMAKE_BUILD_TYPE=Release', '-DSV_LOCAL_VTK_PATH=' + str(vtk),
           '-DSV_ADDITIONAL_CMAKE_ARGS=' + ';'.join(extra)]
configured = run_logged(command, 'native_configure', env=env)
built = run_logged(['/usr/bin/cmake', '--build', build, '--parallel', '4'], 'native_build', env=env) if configured['exit_code'] == 0 else None
with (ROOT / 'reports/sv1/native_build_log.txt').open('a') as stream:
    for attempt in (configured, built):
        if attempt:
            stream.write((ROOT / attempt['log']).read_text())
executable = build / 'svMultiPhysics-build/bin/svmultiphysics'
passed = bool(built and built['exit_code'] == 0 and executable.is_file())
result = {'status': 'PASS' if passed else 'FAIL', 'commit': commit,
          'executable': str(executable) if passed else None,
          'configure': configured, 'build': built, 'mpi': '/usr/bin/mpiexec',
          'extra_linear_algebra_packages': False,
          'official_existing_library_reused': str(sv / 'lib/lib_simvascular_thirdparty_tetgen.a'),
          'mesh_library_recompiled': False, 'runtime_library_path': str(vtk / 'lib')}
if passed:
    result['executable_sha256'] = sha256(executable)
    linked = subprocess.run(['/usr/bin/ldd', str(executable)], env=env, text=True, capture_output=True)
    result['linked_libraries'] = linked.stdout
    result['linked_libraries_complete'] = linked.returncode == 0 and 'not found' not in linked.stdout
    assert result['linked_libraries_complete']
write_json(ROOT / 'reports/sv1/native_solver.json', result)
print(json.dumps(result, indent=2))
raise SystemExit(0 if passed else 1)
