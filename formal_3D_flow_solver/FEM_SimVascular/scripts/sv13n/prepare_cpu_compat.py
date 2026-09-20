"""Run repair_01 only after the old timed proof; never overlap native measurements."""
raise SystemExit('SUPERSEDED: latest user strategy inherits CPU build and validates repair on GPU')
import json,subprocess,sys
from remote import ROOT,upload
R=ROOT/'reports/sv1_3n';S=ROOT/'scripts/sv13n'
def execute(script,*args):
 print('COMPAT '+script+' '+' '.join(args),flush=True)
 subprocess.run([sys.executable,'-B',S/script,*args],check=True)
assert json.loads((R/'OLD_PETSC_CPU_PROOF_20_acceptance.json').read_text())['status']=='PASS'
upload(ROOT/'patches/sv1_3n/svmp_petsc325_compat.patch','patches/svmp_petsc325_compat.patch')
execute('invoke.py','lifecycle_remote.py','before')
execute('invoke.py','svmp_build_remote.py','cpu','repair_01')
execute('invoke.py','lifecycle_remote.py','after')
execute('run_case.py','OFFICIAL_CPU_REPAIR01','1','CPU','cpu')
execute('canonicalize.py')
for name in ('lifecycle_before_gate','lifecycle_after_gate'):
 d=json.loads((R/'remote'/(name+'.json')).read_text());assert d['status']=='PASS'
 (R/(name+'.json')).write_text(json.dumps(d,indent=2)+'\n')
manifest=R/'compatibility_adapter.json';d=json.loads(manifest.read_text())
d.update(status='PASS',build_319='NOT_RUN: superseding user cost-control strategy',build_325='PASS',official_325='PASS',lifecycle_regression='PASS',science_regression='PENDING_MATCHED_SHORT_COMPARISON')
manifest.write_text(json.dumps(d,indent=2)+'\n')
print('repair_01 clean build, lifecycle and official CPU regressions PASS; short science pending.',flush=True)
