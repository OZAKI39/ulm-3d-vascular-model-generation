#!/usr/bin/env python3
"""Read CPU, memory limits and the installed FEM API; no factorization."""
import inspect
import json
import os
import platform
import shutil
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from fem3d.audit import command,sha256,timestamp,write_json
import dolfinx
from dolfinx.fem import petsc
from petsc4py import PETSc
from mpi4py import MPI

O=ROOT/'outputs/stage03/resources';O.mkdir(parents=True,exist_ok=True)
assert not (O/'initial.json').exists()
mem={line.split(':')[0]:int(line.split()[1])*1024 for line in Path('/proc/meminfo').read_text().splitlines() if len(line.split())>=3}
cgroup={}
for name in ('memory.max','memory.current','memory.peak','memory.events','cpu.max','cpuset.cpus.effective'):
    path=Path('/sys/fs/cgroup')/name
    if path.exists():cgroup[name]=path.read_text().strip()
guide=Path('/etc/vast-agents-guide.md')
if guide.exists():(O/'instance_guide.txt').write_text(guide.read_text())
packages={json.loads(p.read_text())['name']:json.loads(p.read_text())['version'] for p in (ROOT/'remote/.env/conda-meta').glob('*.json')}
metadata={'timestamp':timestamp(),'hostname':platform.node(),'cpu':command(['lscpu']),
          'affinity_cpu_count':len(os.sched_getaffinity(0)),'memory_bytes':mem,'cgroup':cgroup,
          'disk':dict(zip(('total','used','free'),shutil.disk_usage(ROOT))),
          'dolfinx_version':dolfinx.__version__,'petsc_version':list(PETSc.Sys.getVersion()),
          'mumps_available':PETSc.Sys.hasExternalPackage('mumps'),'mpi_library':MPI.Get_library_version(),
          'packages':{n:v for n,v in packages.items() if any(k in n for k in ('mumps','petsc','mpi','dolfinx','adios'))},
          'gpu_used':False,'time_binary':shutil.which('time'),'python':sys.executable,
          'fem_python_sha256':sha256(Path(sys.executable))}
write_json(O/'initial.json',metadata)
(O/'installed_linear_problem_solve.py.txt').write_text(inspect.getsource(petsc.LinearProblem.solve))
(O/'assembly_signatures.txt').write_text('\n'.join(f'{n}: {inspect.signature(getattr(petsc,n))}' for n in ('assemble_matrix','assemble_vector','apply_lifting','set_bc')))
print(json.dumps(metadata,indent=2))
print(inspect.getsource(petsc.LinearProblem.solve))
