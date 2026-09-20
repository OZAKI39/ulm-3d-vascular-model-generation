"""After the one real GPU run: archive, audit, figures and fast artifact tests only."""
import json,subprocess,sys,time
from remote import ROOT
R=ROOT/'reports/sv1_3n';S=ROOT/'scripts/sv13n'
while True:
 p=R/'gpu_steady_candidate.json'
 if p.exists():
  try:d=json.loads(p.read_text())
  except json.JSONDecodeError:time.sleep(1);continue
  assert d['status']=='PASS';break
 p=R/'gpu_monitor_failure.json'
 if p.exists():raise SystemExit('Monitor failure needs diagnosis; no closure or new CFD')
 p=R/'REAL_VASCULAR_GPU_acceptance.json'
 if p.exists() and json.loads(p.read_text())['status']=='FAIL':raise SystemExit('Vascular failure needs diagnosis; no closure')
 time.sleep(5)
def execute(script,*args):
 print('CLOSE '+script,flush=True)
 subprocess.run([sys.executable,'-B',S/script,*args],check=True)
execute('invoke.py','final_audit_remote.py')
execute('sync_remote.py');execute('canonicalize.py','gpu13')
execute('mirror_native.py');execute('audit_preservation.py');execute('figures.py')
log=ROOT/'logs/sv1_3n/pytest_full.log';assert not log.exists()
with log.open('x') as f:
 r=subprocess.run([sys.executable,'-B','-m','pytest','-q','-p','no:cacheprovider','--junitxml='+str(R/'pytest_full.xml')],cwd=ROOT,stdout=f,stderr=subprocess.STDOUT)
(R/'pytest_execution.json').write_text(json.dumps(dict(exit_code=r.returncode,log=str(log.relative_to(ROOT)),historical_failures_retained=True),indent=2)+'\n')
execute('audit_preservation.py')
print('Candidate archived; figures/full tests ready for final human-readable review. No new CFD started.',flush=True)
