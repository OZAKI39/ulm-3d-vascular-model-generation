#!/usr/bin/env python3
"""Run upstream regression without rewriting historical P7 logs, plus P8 tests."""
from pathlib import Path
import argparse,json,subprocess,concurrent.futures
repo=Path(__file__).resolve().parents[2];logs=repo/'particle_3d/reports/particle8/logs';logs.mkdir(exist_ok=True)
p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group();g.add_argument('--upstream-only',action='store_true');g.add_argument('--p8-only',action='store_true');a=p.parse_args()
rows=[]
if not a.p8_only:
 old=json.loads((repo/'particle_3d/reports/particle7/logs/final_commands.json').read_text())
 for r in old:
  rows.append(dict(name=r['name'],command=[x.replace('reports/particle7/logs/','reports/particle8/logs/') for x in r['command']],cwd=r['cwd']))
if not a.upstream_only:
 rows.append(dict(name='final_particle8',command=[str(repo/'.venv/bin/python'),'-B','-m','pytest','-q','-p','no:cacheprovider','--junitxml='+str(logs/'final_particle8.xml'),'particle_3d/tests/particle8'],cwd=str(repo)))
def run(row):
 with (logs/(row['name']+'.log')).open('w') as f:r=subprocess.run(row['command'],cwd=row['cwd'],stdout=f,stderr=subprocess.STDOUT)
 print(row['name'],r.returncode,flush=True);return dict(row,returncode=r.returncode)
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:results=list(pool.map(run,rows))
path=logs/'test_commands.json';previous=json.loads(path.read_text()) if path.exists() else []
all_rows={r['name']:r for r in previous};all_rows.update({r['name']:r for r in results});path.write_text(json.dumps(list(all_rows.values()),indent=2)+'\n')
if any(r['returncode'] for r in results):raise SystemExit(1)
