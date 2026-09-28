"""Supervised sequential execution gate, with persistent case reports."""
from pathlib import Path
import argparse,json,time,subprocess,sys
p=argparse.ArgumentParser();p.add_argument('--wait-for',type=Path,required=True);p.add_argument('--case',type=Path,required=True);a=p.parse_args()
def cancelled():
 if (a.case/'reports/user_cancellation.json').exists() or (a.case.parent/'user_cancellation.json').exists():
  print('SKIPPED_BY_USER: dependent case cancelled');raise SystemExit(3)
cancelled()
while not a.wait_for.exists():cancelled();time.sleep(15)
cancelled()
r=json.loads(a.wait_for.read_text());assert r['status']=='PASS','Upstream actual solver did not finish successfully; do not launch dependent run'
raise SystemExit(subprocess.run([sys.executable,'-B',str(Path(__file__).with_name('run_solver.py')),'--case',str(a.case)]).returncode)
