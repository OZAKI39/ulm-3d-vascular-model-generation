#!/usr/bin/env python3
"""Use the MPI and Python from exactly one project environment; 1/2/4 ranks."""
import json
import os
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from fem3d.audit import sha256, write_json

prefix = ROOT/"remote/.env"
python = prefix/"bin/python"
mpiexec = prefix/"bin/mpiexec"
results = []
for ranks in (1,2,4):
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1",
               PYTHONDONTWRITEBYTECODE="1",FEM3D_EXPECTED_RANKS=str(ranks),
               XDG_CACHE_HOME=str(ROOT/"outputs/stage00/cache"),
               PATH=str(prefix/"bin")+os.pathsep+env.get("PATH",""))
    output = f"outputs/stage00/smoke_r{ranks}.json"
    command = [str(python),"-B","scripts/record_run.py","--label",f"dolfinx_smoke_r{ranks}",
               "--ranks",str(ranks),"--remote-dir",str(ROOT),"--",str(python),"-B",
               "scripts/launch_mpi.py","--ranks",str(ranks),"--output",output]
    start = time.perf_counter()
    process = subprocess.run(command,cwd=ROOT,env=env)
    results.append({"mpi_ranks":ranks,"returncode":process.returncode,
                    "wrapper_wall_time_s":time.perf_counter()-start,"result":output,
                    "result_sha256":sha256(ROOT/output) if (ROOT/output).exists() else None})
    write_json(ROOT/"outputs/stage00/smoke_matrix.json",{"runs":results,"status":"PASS" if len(results)==3 and all(r['returncode']==0 for r in results) else "INCOMPLETE_OR_FAIL"})
    # All three small cases are attempted; each failed command retains its exit status.
sys.exit(any(r["returncode"] != 0 for r in results))
