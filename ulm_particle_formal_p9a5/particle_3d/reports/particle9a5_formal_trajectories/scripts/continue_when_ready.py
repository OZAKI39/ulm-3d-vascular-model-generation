"""Finish the fixed benchmark before launching authorized, test-gated production."""
from pathlib import Path
import json,subprocess,sys,time,datetime
ROOT=Path(__file__).resolve().parents[4];R=ROOT/'particle_3d/reports/particle9a5_formal_trajectories';c=json.loads((R/'data/run_context.json').read_text())
probe='''from pathlib import Path
import json
root=Path(REMOTE);done=root/'particle_3d/reports/particle9a5_formal_trajectories/data/worker_scaling_results.json';pid=json.loads((root/'logs/benchmark_pid.json').read_text())['pid'];proc=Path('/proc')/str(pid)
alive=proc.exists() and str(root).encode() in (proc/'cmdline').read_bytes()
print(json.dumps(dict(complete=done.exists(),alive=alive)))
'''.replace('REMOTE',repr(c['remote']))
failures=0
while True:
 try:
  result=subprocess.run(['ssh','-o','BatchMode=yes','-o','ConnectTimeout=15','vast4090','python3','-'],input=probe,text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=40,check=True)
  state=json.loads(result.stdout);failures=0
 except Exception as e:
  failures+=1;print('CONNECTION_RETRY',failures,str(e),flush=True)
  if failures>=5:raise
  time.sleep(45);continue
 if state['complete']:break
 if not state['alive']:raise RuntimeError('Benchmark stopped without completed gate; preserve evidence and inspect logs')
 time.sleep(45)
print('BENCHMARK_COMPLETE_STARTING_TEST_GATED_PRODUCTION',datetime.datetime.now(datetime.timezone.utc).isoformat(),flush=True)
subprocess.run([sys.executable,str(R/'scripts/start_production.py')],cwd=ROOT,check=True)
print('PRODUCTION_AND_POSTHOC_POINTS_LAUNCHED',flush=True)
