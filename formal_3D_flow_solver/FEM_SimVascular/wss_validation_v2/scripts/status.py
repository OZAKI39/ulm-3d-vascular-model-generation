"""Report accepted, running, pending and intentionally stopped cases separately."""
from pathlib import Path
import json,time
v=Path(__file__).resolve().parents[1];rows=[]
for c in sorted(v.glob('stage[234]/*')):
 if not (c/'policy.json').exists():continue
 ex=c/'reports/execution.json';pr=c/'reports/progress.json';cancel=c/'reports/intentional_cancellation.json';d=json.loads(ex.read_text()) if ex.exists() else json.loads(pr.read_text()) if pr.exists() else {};s=d.get('states',[]);last=s[-1] if s else {};status=d.get('status','RUNNING' if pr.exists() else 'NOT_STARTED')
 if cancel.exists():status='INTENTIONALLY_STOPPED_NOT_VALIDATION_RESULT'
 if (c/'reports/baseline_reuse.json').exists():status='REUSED_CFD_CHECKED'
 if (c/'reports/user_cancellation.json').exists():status='SKIPPED_BY_USER'
 rows.append(dict(case=str(c.relative_to(v)),status=status,last_step=last.get('step'),elapsed_s=round(d.get('elapsed_s',0)),mass=last.get('epsilon_mass'),du=last.get('velocity_relative_change'),dp=last.get('pressure_relative_change'),dw=last.get('wss_area_L2_relative_change'),read_errors=len(d.get('read_errors',[]))))
out=dict(observed_unix=time.time(),cases=rows);(v/'logs/current_case_status.json').write_text(json.dumps(out,indent=2)+'\n')
for row in rows:print(json.dumps(row))
