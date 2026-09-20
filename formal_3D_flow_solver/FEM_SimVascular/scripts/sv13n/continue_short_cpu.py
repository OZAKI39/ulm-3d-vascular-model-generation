"""Updated strategy: safe old stop -> repair regression -> one matched CPU short run."""
raise SystemExit('DEFERRED by USER_GPU_PRIORITY_UPDATE.txt: direct GPU development')
import json,subprocess,sys,time
from remote import ROOT
R=ROOT/'reports/sv1_3n';S=ROOT/'scripts/sv13n'
def execute(script,*args):
 print('SHORT '+script+' '+' '.join(args),flush=True)
 subprocess.run([sys.executable,'-B',S/script,*args],check=True)
while True:
 p=R/'OLD_PETSC_CPU_PROOF_20_acceptance.json'
 if p.exists():
  try:d=json.loads(p.read_text())
  except json.JSONDecodeError:time.sleep(1);continue
  assert d['status']=='PASS',d.get('acceptance_errors');break
 time.sleep(3)
execute('freeze_short_reference.py')
execute('prepare_cpu_compat.py')
execute('run_case.py','NEW_CPU_SHORT','1','CPU','cpu')
execute('compare_science.py','old_new_cpu_science','OLD_CPU_SHORT_REFERENCE','NEW_CPU_SHORT')
print('Matched short CPU version science complete; no benchmark queued.',flush=True)
