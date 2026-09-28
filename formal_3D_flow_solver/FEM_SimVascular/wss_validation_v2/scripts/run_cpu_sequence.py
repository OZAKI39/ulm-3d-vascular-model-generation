from pathlib import Path
import sys,subprocess,json,time
root=Path(__file__).resolve().parents[1]
state=root/'logs/cpu_sequence_status.json'
for rel in sys.argv[1:]:
 case=root/rel;out=case/'reports/execution.json'
 if (case/'reports/user_cancellation.json').exists() or (case.parent/'user_cancellation.json').exists():
  print(rel,'SKIPPED_BY_USER',flush=True);continue
 if out.exists() and json.loads(out.read_text()).get('status')=='PASS':continue
 state.write_text(json.dumps(dict(current=rel,status='RUNNING',unix=time.time()))+'\n')
 r=subprocess.run([sys.executable,'-B',str(root/'scripts/run_solver.py'),'--case',str(case)])
 if r.returncode:
  state.write_text(json.dumps(dict(current=rel,status='FAILED',returncode=r.returncode,unix=time.time()))+'\n');raise SystemExit(r.returncode)
state.write_text(json.dumps(dict(status='COMPLETE',unix=time.time(),cases=sys.argv[1:]))+'\n')
