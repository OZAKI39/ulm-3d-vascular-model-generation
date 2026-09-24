#!/usr/bin/env python3
"""Run the frozen legacy suites and/or P7, retaining exact commands and JUnit."""
from pathlib import Path
import argparse,json,subprocess,concurrent.futures
repo=Path(__file__).resolve().parents[2]; logs=repo/'particle_3d/reports/particle7/logs'
p=argparse.ArgumentParser(); g=p.add_mutually_exclusive_group(); g.add_argument('--legacy-only',action='store_true'); g.add_argument('--p7-only',action='store_true'); a=p.parse_args()
rows=[]
if not a.p7_only:
 for r in json.loads((logs/'baseline_commands.json').read_text()):
  rows.append(dict(name=r['name'].replace('baseline_','final_'),command=[x.replace('baseline_','final_') for x in r['command']],cwd=r['cwd']))
if not a.legacy_only:
 rows.append(dict(name='final_particle7',command=[str(repo/'.venv/bin/python'),'-B','-m','pytest','-q','-p','no:cacheprovider','--junitxml='+str(logs/'final_particle7.xml'),'particle_3d/tests/particle7'],cwd=str(repo)))
def run(row):
 with (logs/(row['name']+'.log')).open('w') as stream: result=subprocess.run(row['command'],cwd=row['cwd'],stdout=stream,stderr=subprocess.STDOUT)
 record=dict(row,returncode=result.returncode); print(record['name'],record['returncode'],flush=True); return record
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool: results=list(pool.map(run,rows))
path=logs/'final_commands.json'; old=json.loads(path.read_text()) if path.exists() else []
all_records={r['name']:r for r in old}; all_records.update({r['name']:r for r in results}); path.write_text(json.dumps(list(all_records.values()),indent=2)+'\n')
if any(r['returncode'] for r in results): raise SystemExit(1)
