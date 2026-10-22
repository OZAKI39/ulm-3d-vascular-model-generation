"""Immutable commit deployments; isolated environment and bounded Supervisor jobs."""
import argparse,subprocess,tarfile,io,re
from .common import *
SSH=['ssh','-i','/home/lzy/.ssh/vast_step3b_ed25519','-o','BatchMode=yes','-p','15283','root@95.3.33.46']
SERVER='/workspace/mevo_acoustic_microbubbles/f10aa_simulation_20261018'
def call(code):return subprocess.check_output(SSH+['python3 -'],input=code,text=True)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--module',required=True);ap.add_argument('--label',required=True);ap.add_argument('--seconds',type=int,default=600);a=ap.parse_args()
 assert re.fullmatch('[a-z0-9_]+',a.label) and re.fullmatch('[a-z0-9_]+',a.module)
 paths=['src','scripts/straight_tube_exp1','configs/straight_tube_exp1','scripts/fast_pzt/python.sh','tests/f10aa']
 assert not subprocess.check_output(['git','status','--porcelain','--',*paths],cwd=ROOT,text=True).strip()
 commit=identity()['commit'];prefix='computation_examples/mevo_acoustic_microbubbles/'
 raw=subprocess.check_output(['git','archive',commit,*[prefix+p for p in paths if (ROOT/p).exists()]],cwd=ROOT.parents[1])
 b=io.BytesIO()
 with tarfile.open(fileobj=io.BytesIO(raw)) as src,tarfile.open(fileobj=b,mode='w:gz') as dst:
  for item in src:
   item.name=item.name[len(prefix):]
   if item.name:dst.addfile(item,src.extractfile(item) if item.isfile() else None)
  for name in ['PREREGISTRATION.json','STATES.h5','STATE_RULES.json']:
   dst.add(BASE/name,arcname='results/straight_tube_exp1/f10aa/'+name)
 archive=BASE/'evidence'/(commit+'_'+a.label+'.tar.gz');archive.write_bytes(b.getvalue())
 ex=SERVER+'/execution/'+commit+'_'+a.label
 print(call('from pathlib import Path\nPath('+repr(ex)+').mkdir(parents=True,exist_ok=False)'))
 subprocess.run(['scp','-q','-i',SSH[2],'-P','15283',str(archive),'root@95.3.33.46:'+ex+'/source.tar.gz'],check=True)
 d=dict(server=SERVER,execution=ex,commit=commit,label=a.label,module=a.module,seconds=a.seconds,sha256=sha(archive),name='mevo_f10aa_'+a.label)
 code='d='+repr(d)+'''\nfrom pathlib import Path
import tarfile,hashlib,json,subprocess
ex=Path(d['execution']);base=Path(d['server'])/'results/straight_tube_exp1/f10aa'
base.mkdir(parents=True,exist_ok=True)
archive=ex/'source.tar.gz';assert hashlib.sha256(archive.read_bytes()).hexdigest()==d['sha256']
with tarfile.open(archive) as f:f.extractall(ex,filter='data')
local=ex/'results/straight_tube_exp1'
for p in (local/'f10aa').iterdir():
 q=base/p.name
 if q.exists():assert q.read_bytes()==p.read_bytes()
 else:q.write_bytes(p.read_bytes())
 p.unlink()
(local/'f10aa').rmdir();(local/'f10aa').symlink_to(base,target_is_directory=True)
old=Path('/workspace/mevo_acoustic_microbubbles/f10s_control_validation/results/straight_tube_exp1')
for name in ['f10n','f10q','f10r','f10s']:(local/name).symlink_to(old/name,target_is_directory=True)
(local/'f10z_scaling').symlink_to('/workspace/mevo_acoustic_microbubbles/f10z_scaling_20261018/results/straight_tube_exp1/f10z_scaling',target_is_directory=True)
(local/'f10z').symlink_to('/workspace/mevo_acoustic_microbubbles/f10z_simulation_20261018/results/straight_tube_exp1/f10z',target_is_directory=True)
(ex/'.venv').symlink_to(Path(d['server'])/'.venv',target_is_directory=True)
(ex/'SOURCE_F10AA_COMMIT.txt').write_text(d['commit']+'\\n')
reg=json.loads((base/'PREREGISTRATION.json').read_text())
for k in ['input_sha256','frozen_source_sha256']:
 for rel,h in reg[k].items():assert hashlib.sha256((ex/rel).read_bytes()).hexdigest()==h,rel
assert not subprocess.check_output(['nvidia-smi','--query-compute-apps=pid','--format=csv,noheader'],text=True).strip(),'GPU occupied; no unrelated process will be stopped'
logs=base/'logs';logs.mkdir(exist_ok=True)
cmd='/usr/bin/timeout --signal=TERM --kill-after=15s '+str(d['seconds'])+'s /bin/bash scripts/fast_pzt/python.sh -m scripts.straight_tube_exp1.f10aa.'+d['module']
conf=Path('/etc/supervisor/conf.d')/(d['name']+'.conf');assert not conf.exists()
conf.write_text('[program:'+d['name']+']\\ncommand='+cmd+'\\ndirectory='+str(ex)+'\\nautostart=false\\nautorestart=false\\nstartsecs=0\\nstopasgroup=true\\nkillasgroup=true\\nstdout_logfile='+str(logs/(d['label']+'.log'))+'\\nstderr_logfile='+str(logs/(d['label']+'.stderr.log'))+'\\nstdout_logfile_maxbytes=0\\nstderr_logfile_maxbytes=0\\nenvironment=PYTHONUNBUFFERED="1",WARP_CACHE_PATH="/workspace/mevo_acoustic_microbubbles/.cache/warp"\\n')
(logs/(d['label']+'_supervisor.conf')).write_text(conf.read_text());(logs/(d['label']+'_EXECUTION.json')).write_text(json.dumps(d,indent=2))
for args in [['reread'],['update',d['name']],['start',d['name']]]:subprocess.run(['supervisorctl',*args],check=True)
print(json.dumps(d))
'''
 print(call(code));write(BASE/'evidence'/(a.label+'_EXECUTION.json'),d)
if __name__=='__main__':main()
