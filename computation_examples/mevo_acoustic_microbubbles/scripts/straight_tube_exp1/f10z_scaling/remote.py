"""Immutable Git deployment, existing environment, unique Supervisor jobs."""
import argparse,subprocess,io,tarfile,shlex,json
from .common import *
SSH=['ssh','-i','/home/lzy/.ssh/vast_step3b_ed25519','-o','BatchMode=yes','-p','15283','root@95.3.33.46']
SCP=['scp','-q','-i','/home/lzy/.ssh/vast_step3b_ed25519','-P','15283']
SERVER='/workspace/mevo_acoustic_microbubbles/f10z_scaling_20261018'
def call(code):return subprocess.check_output(SSH+['python3 -'],input=code,text=True)
def main():
 p=argparse.ArgumentParser();p.add_argument('--label',required=True);a=p.parse_args()
 commit=identity()['execution_commit'];prefix='computation_examples/mevo_acoustic_microbubbles/';repo=ROOT.parents[1]
 paths=['src','scripts/straight_tube_exp1','configs/straight_tube_exp1','scripts/fast_pzt/python.sh','tests/f10z_scaling']
 assert not subprocess.check_output(['git','status','--porcelain','--',*paths],cwd=ROOT,text=True).strip()
 raw=subprocess.check_output(['git','archive',commit,*[prefix+k for k in paths]],cwd=repo);buffer=io.BytesIO()
 with tarfile.open(fileobj=io.BytesIO(raw)) as src,tarfile.open(fileobj=buffer,mode='w:gz') as dst:
  for item in src:
   item.name=item.name[len(prefix):]
   if item.name:dst.addfile(item,src.extractfile(item) if item.isfile() else None)
  for n in ['PREREGISTRATION.json','INITIAL_STATES.h5','INITIAL_STATE_AUDIT.json','RUN_SCHEDULES.json','N32_REGRESSION.h5']:
   dst.add(BASE/n,arcname='results/straight_tube_exp1/f10z_scaling/'+n)
 (BASE/'evidence').mkdir(exist_ok=True)
 ex=f'{SERVER}/execution/{commit}_{a.label}';arc=BASE/'evidence'/f'{commit}_{a.label}.tar.gz';arc.write_bytes(buffer.getvalue())
 print(call('from pathlib import Path\np=Path('+repr(ex)+')\np.mkdir(parents=True,exist_ok=False)'))
 subprocess.run(SCP+[str(arc),'root@95.3.33.46:'+ex+'/source.tar.gz'],check=True)
 data=dict(server=SERVER,ex=ex,sha=commit,archive_sha256=sha(arc),label=a.label,name='mevo_f10z_scaling_'+a.label, historical_protected=PROTECTED)
 code='d='+repr(data)+'''\nfrom pathlib import Path
import subprocess,json,hashlib,tarfile
ex=Path(d['ex']);root=Path(d['server']);archive=ex/'source.tar.gz'
assert hashlib.sha256(archive.read_bytes()).hexdigest()==d['archive_sha256']
with tarfile.open(archive) as f:f.extractall(ex,filter='data')
# Separate F10Z_SCALING output root, read-only historical inputs by path and hash.
local=ex/'results/straight_tube_exp1';out=root/'results/straight_tube_exp1/f10z_scaling';out.mkdir(parents=True,exist_ok=True)
for p in (local/'f10z_scaling').iterdir():
 q=out/p.name
 if q.exists():assert q.read_bytes()==p.read_bytes()
 else:q.write_bytes(p.read_bytes())
for p in (local/'f10z_scaling').iterdir():p.unlink()
(local/'f10z_scaling').rmdir();(local/'f10z_scaling').symlink_to(out,target_is_directory=True)
old=Path('/workspace/mevo_acoustic_microbubbles/f10s_control_validation/results/straight_tube_exp1')
for name in ['f10n','f10q','f10r','f10s']:(local/name).symlink_to(old/name,target_is_directory=True)
(ex/'.venv').symlink_to('/workspace/mevo_acoustic_microbubbles/.venv',target_is_directory=True)
(ex/'SOURCE_F10Z_SCALING_COMMIT.txt').write_text(d['sha']+'\\n')
reg=json.loads((out/'PREREGISTRATION.json').read_text())
for rel,digest in reg['input_sha256'].items():assert hashlib.sha256((ex/rel).read_bytes()).hexdigest()==digest,rel
latest=Path('/workspace/mevo_acoustic_microbubbles/f10s_control_validation/.execution/4d12e514e8ed342c2fe7eaa7083aa1487c9bff6c')
mismatches=[]
for rel,digest in reg['protected_source_sha256'].items():
 assert hashlib.sha256((ex/rel).read_bytes()).hexdigest()==digest,rel
 if rel in d['historical_protected'] and (latest/rel).exists() and hashlib.sha256((latest/rel).read_bytes()).hexdigest()!=digest:mismatches.append(rel)
assert not mismatches,mismatches
assert not subprocess.check_output(['nvidia-smi','--query-compute-apps=pid','--format=csv,noheader'],text=True).strip(),'GPU occupied; no unrelated process will be stopped'
log=out/'logs';log.mkdir(exist_ok=True)
command='/usr/bin/timeout --signal=TERM --kill-after=15s 1250s /bin/bash scripts/fast_pzt/python.sh -m scripts.straight_tube_exp1.f10z_scaling.run --label '+d['label']
conf=Path('/etc/supervisor/conf.d')/(d['name']+'.conf');assert not conf.exists()
conf.write_text('[program:'+d['name']+']\\ncommand='+command+'\\ndirectory='+str(ex)+'\\nautostart=false\\nautorestart=false\\nstartsecs=0\\nstopasgroup=true\\nkillasgroup=true\\nstdout_logfile='+str(log/(d['label']+'.log'))+'\\nstderr_logfile='+str(log/(d['label']+'.stderr.log'))+'\\nstdout_logfile_maxbytes=0\\nstderr_logfile_maxbytes=0\\nenvironment=PYTHONUNBUFFERED="1",WARP_CACHE_PATH="/workspace/mevo_acoustic_microbubbles/.cache/warp",XDG_CACHE_HOME="/workspace/mevo_acoustic_microbubbles/.cache"\\n')
(log/(d['label']+'_supervisor.conf')).write_text(conf.read_text())
(log/(d['label']+'_EXECUTION.json')).write_text(json.dumps(d,indent=2))
for args in [['reread'],['update',d['name']],['start',d['name']]]:subprocess.run(['supervisorctl',*args],check=True)
print(json.dumps(d))
'''
 print(call(code));write(BASE/'evidence'/f'{a.label}_EXECUTION.json',data)
if __name__=='__main__':main()
