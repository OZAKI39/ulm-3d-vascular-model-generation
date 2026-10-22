"""Native TetGen interface preserving the official embedded-Python execution route."""
from pathlib import Path
import os,subprocess
from .mesh_tools.tetgen_native import generate as generate_tetra_mesh, validate_root


def launch_embedded(launcher,prepared_root,*,factor=1.0):
    """Explicit future meshing with installed SimVascular, no old installation borrowed.

    The child script remains Python 3.5 compatible. Its native sv extension is
    not loaded into the modern project interpreter's incompatible Python ABI.
    """
    root=validate_root(prepared_root);binary=Path(launcher).expanduser().resolve()
    if not binary.is_file():raise FileNotFoundError(binary)
    if factor not in (1.,.8):raise ValueError('Only original primary/fallback factors allowed')
    script=Path(__file__).parent/'mesh_tools/tetgen_native.py'
    env=dict(os.environ,MOUSE_MESH_ROOT=str(root),MESH_EDGE_FACTOR=str(factor),
        PYTHONDONTWRITEBYTECODE='1',PYTHONNOUSERSITE='1',QT_QPA_PLATFORM='offscreen',OMP_NUM_THREADS='1')
    env.pop('PYTHONPATH',None)
    with (root/'meshing_execution.log').open('x') as stream:
        return subprocess.run([str(binary),'-python','--',str(script)],cwd=root,env=env,
                              stdout=stream,stderr=subprocess.STDOUT,check=True)
