from pathlib import Path
import json,subprocess,time,concurrent.futures,sys
b=Path(__file__).parent
cases=['M2P','M5B']
rows=[dict(case=c,command=['python','-u','analysis/shear_boundary_notch_v2/run.py',c,str(b/c)],state='REGISTERED') for c in cases]
(b/'RUN_LEDGER.json').write_text(json.dumps(rows,indent=2))
def execute(row):
 start=time.perf_counter()
 with (b/(row['case']+'.process.log')).open('w') as f:p=subprocess.run(row['command'],stdout=f,stderr=subprocess.STDOUT)
 row.update(returncode=p.returncode,wall_seconds=time.perf_counter()-start,state='EXITED')
 (b/(row['case']+'.process.json')).write_text(json.dumps(row,indent=2));return row
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:results=list(pool.map(execute,rows))
(b/'CAMPAIGN_EXECUTION.json').write_text(json.dumps(results,indent=2));print(json.dumps(results,indent=2));sys.exit(0 if all(r['returncode']==0 for r in results) else 1)
