"""Read-only status, safe after a disconnected client. Never starts a run."""
from pathlib import Path
import json, subprocess, datetime
R=Path(__file__).resolve().parents[3]
O=R/'verification/y_clot_irregular_current_method_v4/20261020T041100Z'
def read(p):
 try:return json.loads(p.read_text())
 except FileNotFoundError:return None
services=subprocess.run(['supervisorctl','status'],capture_output=True,text=True).stdout.splitlines()
rot=read(O/'rotation/VISUALIZATION_VALIDATION.json')
result=read(O/'NATIVE_V4_RESULT.json')
print(json.dumps(dict(checked_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),
 supervisor=[line for line in services if '20261020T041100Z' in line],
 progress=read(O/'cases/H_Y_IRREGULAR_LPS_V4/PROGRESS.json'),
 result={k:result[k] for k in ['status','accepted_time_s','accepted_steps','max_N20_200','max_relative_R','exception']} if result else None,
 last_legal_checkpoint=read(O/'cases/H_Y_IRREGULAR_LPS_V4/LAST_LEGAL_CHECKPOINT.json'),
 rotation_validated=rot['status'] if rot else None,
 postprocess=read(O/'POSTPROCESS_STATUS.json')),indent=2))
