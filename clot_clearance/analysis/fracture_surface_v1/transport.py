import subprocess,sys,shlex,json
from pathlib import Path
C=json.loads(Path(__file__).with_name('config.json').read_text());L=Path(C['local']);REL=f"verification/{C['task']}/{C['run']}";OUT=L/REL
SSH=['ssh','-p','15283','-i','/home/lzy/.ssh/vast_step3b_ed25519','-o','IdentitiesOnly=yes','-o','BatchMode=yes','root@95.3.33.46']
def remote(s):return subprocess.run(SSH+['bash -lc '+shlex.quote('source /workspace/y_clot/env.sh\ncd /workspace/y_clot/project\n'+s)],check=True)
def sync():
 remote('mkdir -p '+shlex.quote(REL)+' analysis/fracture_surface_v1 results/'+C['task']+'/'+C['run'])
 for rel in ['analysis/fracture_surface_v1',REL]:subprocess.run(['rsync','-a',*(['--ignore-existing'] if rel==REL else []),'-e',shlex.join(SSH[:-1]),str(L/rel)+'/',SSH[-1]+':'+C['server']+'/'+rel+'/'],check=True)
def get():
 for rel in [REL,'results/'+C['task']+'/'+C['run']]:subprocess.run(['rsync','-a','-e',shlex.join(SSH[:-1]),SSH[-1]+':'+C['server']+'/'+rel+'/',str(L/rel)+'/'],check=True)
if __name__=='__main__':
 if sys.argv[1]=='sync':sync()
 elif sys.argv[1]=='get':get()
 else:remote(sys.argv[2])
