"""Sequential gated continuation while long CPU proofs finish; no overlapping native jobs."""
raise SystemExit('SUPERSEDED by convergence-driven user strategy; use revised short/steady pipeline')
import json,subprocess,sys,time
from pathlib import Path
from remote import ROOT,upload
R=ROOT/'reports/sv1_3n';S=ROOT/'scripts/sv13n'
while True:
 p=R/'old_new_cpu_science.json'
 if p.exists():
  try:d=json.loads(p.read_text())
  except json.JSONDecodeError:time.sleep(1);continue
  assert d['status']=='PASS','CPU_SCIENCE_GATE_FAILED'
  break
 for name in ('OLD_PETSC_CPU_PROOF_20','NEW_PETSC_CPU_PROOF_20'):
  p=R/(name+'_acceptance.json')
  if p.exists():
   try:d=json.loads(p.read_text())
   except json.JSONDecodeError:continue
   assert d['status']=='PASS','CPU_PROOF_GATE_FAILED: '+name
 time.sleep(5)
def execute(script,*args):
 print('START '+script+' '+' '.join(args),flush=True)
 subprocess.run([sys.executable,'-B',S/script,*args],check=True)
execute('benchmark_cpu.py')
upload(R/'old_new_cpu_science.json','reports/old_new_cpu_science.json')
upload(R/'petsc325_source.json','configs/petsc325_source.json')
execute('invoke.py','petsc_build_remote.py','gpu13')
execute('invoke.py','petsc_smoke_remote.py','gpu13')
execute('invoke.py','ghost_probe_remote.py','old')
execute('invoke.py','ghost_probe_remote.py','gpu13')
execute('invoke.py','coherence_remote.py','gpu13')
execute('sync_remote.py');execute('canonicalize.py','gpu13')
print('G1 standalone and ghost observations ready for review; no unvalidated solver adoption.',flush=True)
