"""Wait for actual CFD AND independent medium audit, initialize fine, then actually solve."""
from pathlib import Path
import json,time,sys,subprocess,hashlib
V=Path(__file__).resolve().parents[1];source=V/'stage3/vessel_medium';case=V/'stage3/vessel_fine'
assert V==Path('/workspace/wss_validation_v2_20260927T1230Z')
def cancelled():
 if (case/'reports/user_cancellation.json').exists():
  print('SKIPPED_BY_USER: fine vascular CFD and initialization disabled',flush=True);raise SystemExit(3)
cancelled()
while not (source/'reports/execution.json').exists():cancelled();time.sleep(15)
assert json.loads((source/'reports/execution.json').read_text())['status']=='PASS'
while not (source/'reports/independent_acceptance.json').exists():cancelled();time.sleep(15)
cancelled()
gate=json.loads((source/'reports/independent_acceptance.json').read_text());assert gate['accepted']
assert hashlib.sha256((source/'frozen_flow/flow_arrays_si.npz').read_bytes()).hexdigest()==gate['source_flow_sha256']
subprocess.run([sys.executable,'-B',str(V/'scripts/prepare_initial_guess.py'),'--source',str(source),'--case',str(case)],check=True)
(case/'reports/fine_initialization_gate.json').write_text(json.dumps(gate,indent=2)+'\n')
raise SystemExit(subprocess.run([sys.executable,'-B',str(V/'scripts/run_solver.py'),'--case',str(case)]).returncode)
