"""Mirror only one new stage's case; every solver result is checked against its runtime SHA."""
import json,re,shlex,subprocess,sys
from pathlib import Path
from remote import ROOT,REMOTE,ssh_prefix
name=sys.argv[1];assert re.fullmatch(r'[A-Za-z0-9_]+',name)
subprocess.run([sys.executable,'-B',str(ROOT/'scripts/sv13m/sync_remote.py')],check=True)
ssh=ssh_prefix();target=ROOT/'outputs/sv1_3m'/name;target.mkdir(parents=True,exist_ok=True)
subprocess.run(['rsync','-a','-e',shlex.join(ssh[:-1]),ssh[-1]+':'+REMOTE+'/outputs/'+name+'/',str(target)+'/'],check=True)
sys.path.insert(0,str(ROOT/'src'))
from sv_validation.provenance import sha256
record=json.loads((ROOT/'reports/sv1_3m/remote'/(name+'_execution.json')).read_text())
for f in record['results']:
    p=ROOT/'outputs/sv1_3m'/Path(f['path']).relative_to('outputs')
    assert sha256(p)==f['sha256'],p
print(name+': native result hashes verified')
