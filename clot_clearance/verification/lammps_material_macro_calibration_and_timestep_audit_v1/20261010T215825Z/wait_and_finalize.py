import time,json,subprocess
from pathlib import Path
base=Path(__file__).parent
for _ in range(720):
 if (base/'CAMPAIGN_EXECUTION.json').exists() and (base/'M3P/EXECUTION.json').exists():break
 time.sleep(5)
else:raise RuntimeError('Bounded wait for current campaign exceeded')
results=json.loads((base/'CAMPAIGN_EXECUTION.json').read_text())
assert len(results)==8 and all(r['returncode']==0 for r in results),results
assert json.loads((base/'M3P/EXECUTION.json').read_text())['status']=='COMPLETED'
subprocess.run(['python','analysis/material_macro_audit_v1/finalize.py',str(base),'/workspace/y_clot/project/results/lammps_material_macro_calibration_and_timestep_audit_v1/20261010T215825Z'],check=True)
