"""Offline fresh-process verification of actual VTU/checkpoint timestamps; no CFD."""
import json,math,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from sv_validation.sv13n import checkpoint_one_rank
R=ROOT/'reports/sv1_3o';cases={}
for path in sorted(R.glob('*_acceptance.json')):
 a=json.loads(path.read_text())
 if a['reload'] is None:continue
 v=json.loads(subprocess.check_output([sys.executable,'-B',ROOT/'scripts/sv13o/reload_field.py',a['reload']['path']],text=True))
 e=json.loads((R/'remote'/(a['name']+'_execution.json')).read_text())
 cp=checkpoint_one_rank(a['checkpoint']['path'],a['stop_step'],e['dt_s'])
 assert v['sha256']==a['reload']['sha256'] and cp['sha256']==a['checkpoint']['sha256']
 assert math.isclose(v['time_s'],cp['time_s'],rel_tol=1e-12)
 assert v['points']==cp['nodes'] and v['velocity_finite'] and v['pressure_finite']
 cases[a['name']]=dict(output_consistency='PASS',solver_status=a['status'],VTU=v,checkpoint=cp,
                      interpretation='File integrity/time agreement only; does not override solver failure')
assert 'OFFICIAL_OUTPUT_STOP' in cases and 'REAL_VASCULAR_GPU_PERF' in cases
(R/'output_consistency.json').write_text(json.dumps(dict(status='PASS',cases=cases),indent=2)+'\n')
print('Actual VTU TimeValue and native checkpoint step/time agree; fresh-process reloads PASS.')
