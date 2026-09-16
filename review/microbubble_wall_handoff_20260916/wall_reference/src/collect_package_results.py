"""First formal output collector. Never retries or repairs a package result."""
from pathlib import Path
import json,sys,warnings,time,os,hashlib
import numpy as np
R=Path(__file__).resolve().parents[1]
plan=json.loads((R/'WALL_REFERENCE_CONVENTION.json').read_text());label=sys.argv[1]
if label=='PYSTOKES':from pystokes_reference_adapter import get_wall_mobility
else:
    from rmbw_reference_adapter import get_wall_mobility as _get
    def get_wall_mobility(a,h,mu):return _get(a,h,mu,'lubrication' if label=='RMBW_LUBRICATION' else 'sphere_semianalytical')
file=R/'raw'/(label+'.jsonl');assert not file.exists(),'Formal outputs already exist; no retry/overwrite'
with file.open('x') as out:
    for size,a in plan['radii_m'].items():
        for epsilon in plan['gaps_h_over_a']:
            for mu in plan['viscosities_Pa_s']:
                row={'implementation':label,'size':size,'radius_m':a,'gap_m':a*epsilon,'epsilon':epsilon,'viscosity_Pa_s':mu,'first_formal_output':True}
                with warnings.catch_warnings(record=True) as seen:
                    warnings.simplefilter('always')
                    try:
                        t=time.perf_counter();M,meta=get_wall_mobility(a,a*epsilon,mu);row.update(matrix_SI=M.tolist(),metadata=meta,finite=bool(np.all(np.isfinite(M))),status='RETURNED',wall_seconds=time.perf_counter()-t)
                    except Exception as e:row.update(status='EXCEPTION',error=repr(e),matrix_SI=None,finite=False)
                    row['warnings']=[str(w.message) for w in seen]
                out.write(json.dumps(row,allow_nan=True)+'\n');out.flush()
print(label,'FORMAL_ROWS',len(plan['radii_m'])*len(plan['gaps_h_over_a'])*len(plan['viscosities_Pa_s']),flush=True)
