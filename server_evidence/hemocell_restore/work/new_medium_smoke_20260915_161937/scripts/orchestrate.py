from pathlib import Path
import subprocess,json,time,traceback
R=Path(__file__).resolve().parents[1]
def state(**v):(R/'SMOKE_STATE.json').write_text(json.dumps(dict(unix=time.time(),**v),indent=2)+'\n')
with (R/'SMOKE_STARTED.json').open('x') as f:json.dump(dict(unix=time.time(),runs=[500,5000],from_zero=True,retries=0),f)
try:
 for n in (500,5000):
  state(state='RUNNING',phase='SOLVER',horizon=n)
  rc=subprocess.call(['/usr/bin/python3','-u','-B',str(R/'scripts/run_smoke.py'),str(n)])
  state(state='RUNNING',phase='AUDIT',horizon=n,solver_returncode=rc)
  audit=subprocess.call(['/usr/bin/python3','-u','-B',str(R/'scripts/evaluate_smoke.py'),str(n)])
  if rc or audit:raise RuntimeError(f'{n} gate failed; downstream launch prohibited')
 state(state='PASS',phase='NUMERICAL_RUNS_AND_INDEPENDENT_AUDITS_COMPLETE',completed_runs=[500,5000],total_steps=5500,RBC_steps=0)
except BaseException as exc:
 state(state='FAIL',phase='STOPPED',error=str(exc),traceback=traceback.format_exc());raise
