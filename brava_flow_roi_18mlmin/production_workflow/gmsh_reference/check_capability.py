"""Record the isolated installed runtime; never install or mesh on import."""
from pathlib import Path
import datetime
import hashlib
import json
import sys

HERE = Path(__file__).resolve().parent
RUNTIME = HERE.parent / 'runtime/gmsh'


def main():
    if (HERE/'CAPABILITY.json').exists():
        raise SystemExit('Capability receipt is frozen; do not overwrite this one-candidate experiment.')
    sys.path.insert(0, str(RUNTIME))
    import gmsh
    gmsh.initialize([], readConfigFiles=False)
    try:
        files = [Path(gmsh.__file__), Path(gmsh.libpath), *RUNTIME.glob('*.whl'), *RUNTIME.glob('*.deb')]
        start = datetime.datetime.fromisoformat('2026-10-03T19:59:05.495020+00:00')
        now = datetime.datetime.now(datetime.timezone.utc)
        result = dict(version=gmsh.__version__, build_info=gmsh.option.getString('General.BuildInfo'),
                      build_options=gmsh.option.getString('General.BuildOptions'),
                      capability_start_UTC=start.isoformat(), checked_UTC=now.isoformat(),
                      capability_elapsed_seconds=(now-start).total_seconds(),
                      existing_before_install=dict(executable=None, system_module=None, main_venv_module=None,
                                                   population_venv_module=None),
                      install_method='pip install --no-deps --target runtime/gmsh gmsh==4.15.2; extract original libgmsh from the same pinned wheel to repair pip --target layout; apt-get download and dpkg-deb -x libglu1-mesa/libxinerama1/libopengl0 into runtime/gmsh/dependencies only',
                      global_environment_modified=False, cpu_quota=7.68, threads=7,
                      algorithm_documentation='https://gmsh.info/doc/texinfo/#Mesh-options',
                      algorithm_3d=dict(name='HXT', id=10, version='4.15.2'),
                      binary_path=gmsh.libpath,
                      files_sha256={str(p.relative_to(RUNTIME)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files})
        assert all(s in result['build_options'] for s in ['Hxt', 'OpenCASCADE', 'OpenMP'])
        result['status'] = 'PASS' if result['capability_elapsed_seconds'] <= 600 else 'FAIL_CAPABILITY_TIME_LIMIT'
        (HERE/'CAPABILITY.json').write_text(json.dumps(result, indent=2)+'\n')
        print(json.dumps(result, indent=2))
    finally:
        gmsh.finalize()


if __name__ == '__main__':
    main()
