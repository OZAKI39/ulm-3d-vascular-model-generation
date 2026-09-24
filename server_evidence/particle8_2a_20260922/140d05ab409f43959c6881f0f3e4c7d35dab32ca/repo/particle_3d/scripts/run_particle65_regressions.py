#!/usr/bin/env python3
"""Re-run preserved Frozen/P0-P6 commands plus new V1 suite; save logs and JUnit."""
from pathlib import Path
import concurrent.futures,json,subprocess,sys
REPO=Path(__file__).resolve().parents[2];LOGS=REPO/'particle_3d/reports/particle6_5/logs'

def main():
 rows=json.loads((LOGS/'baseline_commands.json').read_text())
 commands=[]
 for row in rows:
  commands.append(dict(name=row['name'].replace('baseline_','final_'),cwd=row['cwd'],
   command=[s.replace('/baseline_','/final_') for s in row['command']]))
 commands.append(dict(name='final_particle6_5',cwd=str(REPO),command=[str(REPO/'.venv/bin/python'),'-B','-m','pytest','-q','-p','no:cacheprovider',
  '--junitxml='+str(LOGS/'final_particle6_5.xml'),'particle_3d/tests/particle6_5']))
 def run(row):
  with (LOGS/(row['name']+'.log')).open('w') as out:result=subprocess.run(row['command'],cwd=row['cwd'],stdout=out,stderr=subprocess.STDOUT)
  print(row['name'],result.returncode,flush=True);return dict(row,returncode=result.returncode)
 with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:results=list(pool.map(run,commands))
 (LOGS/'final_commands.json').write_text(json.dumps(results,indent=2)+'\n')
 if any(r['returncode'] for r in results):raise SystemExit('REGRESSION_FAILURE; SEE LOGS')
if __name__=='__main__':main()
