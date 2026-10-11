from concurrent.futures import ThreadPoolExecutor
import subprocess,json,time
from pathlib import Path
out=Path(__file__).resolve().parent
cases=[('M1',[]),('M2',[]),('M3',[]),('M5',[]),('M4',[]),('M6',[]),('T1',['--dt','2.5e-7']),('D1',['--dx','4e-5'])]
(out/'RUN_LEDGER.json').write_text(json.dumps([{'case':c,'command':['python','-m','pd_clot.lammps_fragmentation.macro_lab_v1','--output',str(out/c),'--case',c]+a,'uncertainty':c+' registered macro response / numerical sensitivity','decision':'Support or limit mapping; no tuning','status':'REGISTERED'} for c,a in cases],indent=2))
def work(item):
 c,args=item;command=['python','-m','pd_clot.lammps_fragmentation.macro_lab_v1','--output',str(out/c),'--case',c]+args
 with (out/(c+'.process.log')).open('w') as f:
  t=time.monotonic();p=subprocess.run(command,stdout=f,stderr=subprocess.STDOUT)
 rec={'case':c,'returncode':p.returncode,'wall_seconds':time.monotonic()-t,'command':command}
 (out/(c+'.process.json')).write_text(json.dumps(rec,indent=2));print(json.dumps(rec),flush=True)
 return rec
with ThreadPoolExecutor(max_workers=4) as pool: results=list(pool.map(work,cases))
(out/'CAMPAIGN_EXECUTION.json').write_text(json.dumps(results,indent=2))
