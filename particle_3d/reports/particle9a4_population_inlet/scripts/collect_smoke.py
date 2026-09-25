from pathlib import Path
import hashlib,json,subprocess,shlex
ROOT=Path(__file__).resolve().parents[4];R=ROOT/'particle_3d/reports/particle9a4_population_inlet'
remote=json.loads((R/'data/server_deployment.json').read_text())['server_root'];remote_out=remote+'/particle_3d/reports/particle9a4_population_inlet/outputs'
assert not (R/'outputs').exists()
code='from pathlib import Path;import json,hashlib;p=Path('+repr(remote_out)+');assert (p/"smoke30/point_completed.json").is_file();print(json.dumps({str(f.relative_to(p)):hashlib.sha256(f.read_bytes()).hexdigest() for f in p.rglob("*") if f.is_file()}))'
s=subprocess.run(['ssh','-o','BatchMode=yes','vast4090','python3 -c '+shlex.quote(code)],capture_output=True,text=True,check=True);manifest=json.loads(s.stdout)
(R/'data/server_output_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
subprocess.run(['scp','-q','-r','vast4090:'+remote_out,str(R)+'/'],check=True)
for name,expected in manifest.items():
 p=R/'outputs'/name
 assert hashlib.sha256(p.read_bytes()).hexdigest()==expected,name
with (R/'logs/server_smoke30.txt').open('x') as f:subprocess.run(['ssh','-o','BatchMode=yes','vast4090','cat '+shlex.quote(remote+'/smoke30.log')],stdout=f,check=True)
print('COLLECTED_AND_SHA_VERIFIED',len(manifest),flush=True)
