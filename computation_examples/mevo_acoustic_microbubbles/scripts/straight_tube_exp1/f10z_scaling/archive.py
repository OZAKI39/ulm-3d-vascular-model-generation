"""Pull completed campaigns and verify every transferred byte."""
import subprocess,json
from .remote import SSH,SERVER
from .common import *
def main():
 remote=SERVER+'/results/straight_tube_exp1/f10z_scaling'
 code='root='+repr(remote)+'''\nfrom pathlib import Path
import json,hashlib
r=Path(root);rows=[]
for p in r.rglob('*'):
 if p.is_file() and p.name!='REMOTE_MANIFEST.json':rows.append({'path':str(p.relative_to(r)),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
(r/'REMOTE_MANIFEST.json').write_text(json.dumps(rows,indent=2))
print(len(rows))
'''
 print(subprocess.check_output(SSH+['python3 -'],input=code,text=True))
 subprocess.run(['rsync','-a','-e','ssh -i /home/lzy/.ssh/vast_step3b_ed25519 -o BatchMode=yes -p 15283','root@95.3.33.46:'+remote+'/',str(BASE)+'/'],check=True)
 rows=read(BASE/'REMOTE_MANIFEST.json')
 for row in rows:assert sha(BASE/row['path'])==row['sha256'],row['path']
 print('Verified transferred files:',len(rows))
if __name__=='__main__':main()
