"""One explicit smoke/window; fixed input, no implicit retry or CFD in tests."""
import json,subprocess,sys
from remote import ROOT,upload
R=ROOT/'reports/sv1_3p';C=ROOT/'configs/sv1_3p';S=ROOT/'scripts/sv13p'
name,mode,candidate=sys.argv[1:4];assert mode in ('smoke','window')
config=json.loads((C/'candidates'/f'{candidate}.json').read_text())
if mode=='window':assert json.loads((R/(candidate+'_SMOKE_acceptance.json')).read_text())['status']=='PASS'
plan=dict(name=name,mode=mode,candidate=candidate,MPI_ranks=1,GPUs=1,OMP_NUM_THREADS=1,start_step=60,end_step=62 if mode=='smoke' else 70,PETSC_OPTIONS=config['PETSC_OPTIONS'],build_report=config['build_report'])
p=C/'runplans'/(name+'.json');p.parent.mkdir(exist_ok=True);assert not p.exists();p.write_text(json.dumps(plan,indent=2)+'\n');upload(p,'configs/runplans/'+p.name)
native=subprocess.run([sys.executable,'-B',S/'invoke.py','flow_remote.py',name])
subprocess.run([sys.executable,'-B',S/'fetch_case.py',name],check=True)
accept=subprocess.run([sys.executable,'-B',S/'accept_case.py',name])
raise SystemExit(0 if native.returncode==accept.returncode==0 else 1)
