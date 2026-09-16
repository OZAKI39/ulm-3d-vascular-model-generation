from pathlib import Path
import shutil,json,subprocess,os,concurrent.futures,time
R=Path(__file__).resolve().parents[1];W=Path(json.loads((R/'provenance/TASK_PATHS.json').read_text())['local_work']);D=W/'scalar_work_beta025'
assert not D.exists(),'NO_GENERATOR_OVERWRITE'
shutil.copytree(W/'upstream/stokesian-dynamics/find_resistance_scalars',D)
shutil.copy2(W/'scalar_work/helen_fortran/lamb.exe',D/'helen_fortran/lamb.exe');(D/'values_of_lambda.txt').write_text('0.25\n')
(R/'contracts/SCALAR_GENERATION_STAGE2.json').write_text(json.dumps(dict(status='FROZEN_BEFORE_GENERATION',beta=.25,lambdas=[.25,4],grid='unmodified upstream 38 nodes',Fortran_globalerror=1e-8,source_math_changes=False,time_bound_seconds_each=1200),indent=2))
def generate(which):
 t=time.monotonic()
 with (R/f'logs/GENERATE_{which}_0P25.log').open('w') as f:p=subprocess.run([str(W/'env/bin/python'),'-B',str(D/f'generate_scalars_{which}field.py')],cwd=D,env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1'),stdout=f,stderr=subprocess.STDOUT,timeout=1200)
 return dict(stage=which,exit_code=p.returncode,seconds=time.monotonic()-t)
with concurrent.futures.ThreadPoolExecutor(2) as pool:results=list(pool.map(generate,['mid','near']))
(R/'validation/SCALAR_GENERATOR_BETA025_RECEIPT.json').write_text(json.dumps(results,indent=2));print(results,flush=True);assert all(x['exit_code']==0 for x in results)
