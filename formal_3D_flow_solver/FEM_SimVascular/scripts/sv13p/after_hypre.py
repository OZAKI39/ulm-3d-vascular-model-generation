"""Run P4 only after P3 closed; then return for AmgX feasibility work."""
import json,subprocess,sys
from remote import ROOT,upload
S=ROOT/'scripts/sv13p';R=ROOT/'reports/sv1_3p'
assert (R/'P3_SMOKE_acceptance.json').exists()
p3=json.loads((R/'P3_SMOKE_acceptance.json').read_text())
assert p3['status']=='FAIL' or (R/'P3_WINDOW_acceptance.json').exists()
ref=json.loads((R/'reference_freeze.json').read_text());p=R/'baseline_solver_build.json';p.write_text(json.dumps(dict(ref['svmp'],PETSc_build=ref['PETSc_build']),indent=2)+'\n');upload(p,'reports/svmp_baseline_build.json')
assert subprocess.run([sys.executable,'-B',S/'make_candidates.py','P4']).returncode==0
code=subprocess.run([sys.executable,'-B',S/'run_case.py','P4_SMOKE','smoke','P4']).returncode
if code==0:subprocess.run([sys.executable,'-B',S/'run_case.py','P4_WINDOW','window','P4'])
print('P4 bounded evaluation finished',flush=True)
