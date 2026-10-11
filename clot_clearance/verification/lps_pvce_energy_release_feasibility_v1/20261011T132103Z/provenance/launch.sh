#!/bin/bash
set -euo pipefail
source /workspace/y_clot/env.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
cd /workspace/y_clot/project
python - <<'RUN'
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import subprocess,json,time
b=Path('verification/lps_pvce_energy_release_feasibility_v1/20261011T132103Z');t=time.perf_counter()
def one(c):
 with (b/c/'stdout.log').open('w') as f:r=subprocess.run(['python','-m','analysis.lps_pvce_v1.run',c],stdout=f,stderr=subprocess.STDOUT)
 return dict(case=c,exit_code=r.returncode)
with ThreadPoolExecutor(max_workers=3) as pool:rows=list(pool.map(one,['A0','A1','A2']))
(b/'provenance/COMPUTE_EXECUTION.json').write_text(json.dumps(dict(cases=rows,concurrent_wall_seconds=time.perf_counter()-t,formal_cases=3),indent=2))
assert all(r['exit_code']==0 for r in rows),rows
RUN
