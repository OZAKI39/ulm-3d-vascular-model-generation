"""Bounded immutable execution snapshots; never overwrite F10AA or stop other jobs."""
import argparse,subprocess,tarfile,io,re
from .common import *
SSH=['ssh','-i','/home/lzy/.ssh/vast_step3b_ed25519','-o','BatchMode=yes','-p','15283','root@95.3.33.46']
SERVER='/workspace/mevo_acoustic_microbubbles/f10ab_shape_20261018'
AA='/workspace/mevo_acoustic_microbubbles/f10aa_simulation_20261018'
def call(code):return subprocess.check_output(SSH+['python3 -'],input=code,text=True)
def audit():
 local=read(BASE/'BASELINE_AUDIT.json');paths=local['source_sha256'];inputs=read(BASE/'PREREGISTRATION.json')['input_sha256']
 code='paths='+repr(list(paths))+'\nhistory='+repr(local['historical_trajectories'])+'''\nfrom pathlib import Path
import hashlib,json,subprocess,os
base=Path('/workspace/mevo_acoustic_microbubbles/f10aa_simulation_20261018/delivery/6028598d0012912cc19aa766fa786c4afb0a91a5')
roots={'f10s':'f10s_control_validation','f10z':'f10z_simulation_20261018','f10z_scaling':'f10z_scaling_20261018','f10aa':'f10aa_simulation_20261018'}
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
files={p:sha(base/p) for p in paths}
checks=[]
for item in history:
 stage=item['path'].split('/')[2];raw=Path('/workspace/mevo_acoustic_microbubbles')/roots[stage]/item['path']
 checks.append(dict(path=item['path'],exists=raw.exists(),matches=raw.exists() and sha(raw)==item['sha256']))
repo=Path('/workspace/mevo_acoustic_microbubbles/repo')
print(json.dumps(dict(delivery=json.loads((base/'DELIVERY.json').read_text()),source_sha256=files,history=checks,
 repo_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip(),repo_status=subprocess.check_output(['git','status','--porcelain','-uno'],cwd=repo,text=True),
 gpu=subprocess.check_output(['nvidia-smi','--query-gpu=name,memory.used,utilization.gpu','--format=csv,noheader'],text=True),
 cpus=os.cpu_count(),load=os.getloadavg(),memory=Path('/proc/meminfo').read_text(),
 disk=subprocess.check_output(['df','-h','/workspace'],text=True))))
'''
 result=json.loads(call(code));result['all_sources_match']=result['source_sha256']==paths
 assert result['all_sources_match'];assert all(r['matches'] for r in result['history'])
 write(BASE/'evidence/REMOTE_BASELINE_AUDIT.json',result);print('Verified server delivery source and',len(result['history']),'historical trajectories')
def deploy(stage,label,seconds,post=False):
 assert re.fullmatch('[a-z0-9_]+',label) and re.fullmatch('[a-z0-9_]+',stage)
 paths=['src','scripts/straight_tube_exp1','scripts/fast_pzt/python.sh','configs/straight_tube_exp1','tests/f10ab']
 assert not subprocess.check_output(['git','status','--porcelain','--',*paths],cwd=ROOT,text=True).strip()
 commit=identity()['source_commit'];prefix='computation_examples/mevo_acoustic_microbubbles/'
 raw=subprocess.check_output(['git','archive',commit,*[prefix+p for p in paths]],cwd=ROOT.parents[1]);b=io.BytesIO()
 with tarfile.open(fileobj=io.BytesIO(raw)) as src,tarfile.open(fileobj=b,mode='w:gz') as dst:
  for item in src:
   item.name=item.name[len(prefix):]
   if item.name:dst.addfile(item,src.extractfile(item) if item.isfile() else None)
  for name in ['PREREGISTRATION.json','INITIAL_STATES.h5','INITIAL_STATE_AUDIT.json','ACTIONS.json','BASELINE_AUDIT.json','ACTION_SELECTION.json','TARGET_PREREGISTRATION.json','EXTENSION_PREREGISTRATION.json','TARGET_SELECTION.json']:
   if (BASE/name).exists():dst.add(BASE/name,arcname='results/straight_tube_exp1/f10ab/'+name)
  if (BASE/'plans'/f'{stage}.json').exists():dst.add(BASE/'plans'/f'{stage}.json',arcname='results/straight_tube_exp1/f10ab/plans/'+stage+'.json')
 archive=BASE/'evidence'/(commit+'_'+label+'.tar.gz');archive.write_bytes(b.getvalue());ex=SERVER+'/execution/'+commit+'_'+label
 call('from pathlib import Path\nPath('+repr(ex)+').mkdir(parents=True,exist_ok=False)')
 subprocess.run(['scp','-q','-i',SSH[2],'-P','15283',str(archive),'root@95.3.33.46:'+ex+'/source.tar.gz'],check=True)
 d=dict(server=SERVER,execution=ex,commit=commit,label=label,stage=stage,seconds=seconds,sha256=sha(archive),name='mevo_f10ab_'+label,post=post)
 print(call('d='+repr(d)+'''\nfrom pathlib import Path
import tarfile,hashlib,json,subprocess
ex=Path(d['execution']);server=Path(d['server']);base=server/'results/straight_tube_exp1/f10ab';base.mkdir(parents=True,exist_ok=True)
a=ex/'source.tar.gz';assert hashlib.sha256(a.read_bytes()).hexdigest()==d['sha256']
with tarfile.open(a) as f:f.extractall(ex,filter='data')
local=ex/'results/straight_tube_exp1'
for p in (local/'f10ab').rglob('*'):
 if not p.is_file():continue
 q=base/p.relative_to(local/'f10ab');q.parent.mkdir(parents=True,exist_ok=True)
 if q.exists():assert q.read_bytes()==p.read_bytes()
 else:q.write_bytes(p.read_bytes())
import shutil
shutil.rmtree(local/'f10ab');(local/'f10ab').symlink_to(base,target_is_directory=True)
old=Path('/workspace/mevo_acoustic_microbubbles/f10s_control_validation/results/straight_tube_exp1')
for n in ['f10n','f10q','f10r','f10s']:(local/n).symlink_to(old/n,target_is_directory=True)
for n,folder in [('f10z','f10z_simulation_20261018'),('f10z_scaling','f10z_scaling_20261018'),('f10aa','f10aa_simulation_20261018')]:
 (local/n).symlink_to('/workspace/mevo_acoustic_microbubbles/'+folder+'/results/straight_tube_exp1/'+n,target_is_directory=True)
if not (server/'.venv').exists():(server/'.venv').symlink_to('/workspace/mevo_acoustic_microbubbles/f10aa_simulation_20261018/.venv',target_is_directory=True)
(ex/'.venv').symlink_to(server/'.venv',target_is_directory=True)
(ex/'SOURCE_F10AB_COMMIT.txt').write_text(d['commit']+'\\n')
logs=base/'logs';logs.mkdir(exist_ok=True)
script=ex/'run.sh';script.write_text('#!/bin/bash\\nset -euo pipefail\\nexport OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg\\nexport PYTHONPATH=src:.\\nexec '+(str(server/'.postvenv/bin/python') if d.get('post') else '.venv/bin/python')+' -m scripts.straight_tube_exp1.f10ab.'+('post' if d.get('post') else 'run')+' '+d['stage']+'\\n')
conf=Path('/etc/supervisor/conf.d')/(d['name']+'.conf');assert not conf.exists()
conf.write_text('[program:'+d['name']+']\\ncommand=/usr/bin/timeout --signal=TERM --kill-after=15s '+str(d['seconds'])+' /bin/bash '+str(script)+'\\ndirectory='+str(ex)+'\\nautostart=false\\nautorestart=false\\nstartsecs=0\\nstopasgroup=true\\nkillasgroup=true\\nstdout_logfile='+str(logs/(d['label']+'.log'))+'\\nstderr_logfile='+str(logs/(d['label']+'.stderr.log'))+'\\nstdout_logfile_maxbytes=0\\nstderr_logfile_maxbytes=0\\n')
(logs/(d['label']+'_EXECUTION.json')).write_text(json.dumps(d,indent=2));(logs/(d['label']+'_supervisor.conf')).write_text(conf.read_text())
for a in [['reread'],['update',d['name']],['start',d['name']]]:subprocess.run(['supervisorctl',*a],check=True)
print(json.dumps(d))
'''))
 write(BASE/'evidence'/(label+'_EXECUTION.json'),d)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('stage');ap.add_argument('--label');ap.add_argument('--seconds',type=int,default=600);ap.add_argument('--post',action='store_true');a=ap.parse_args()
 if a.stage=='audit':audit()
 else:deploy(a.stage,a.label,a.seconds,a.post)

def reuse(stage,label,snapshot_label,seconds=600,post=False):
 """Reuse an unchanged execution snapshot; launch a new external immutable wrapper."""
 assert re.fullmatch('[a-z0-9_]+',label) and re.fullmatch('[a-z0-9_]+',stage)
 d=read(BASE/'evidence'/(snapshot_label+'_EXECUTION.json')).copy()
 d.update(label=label,stage=stage,seconds=seconds,post=post,name='mevo_f10ab_'+label,reused_snapshot=snapshot_label)
 payload={}
 for p in [BASE/'plans'/f'{stage}.json',BASE/'TARGET_PREREGISTRATION.json',BASE/'TARGET_SELECTION.json',BASE/'EXTENSION_PREREGISTRATION.json']:
  if p.exists():payload[str(p.relative_to(BASE))]=p.read_text()
 print(call('d='+repr(d)+'\npayload='+repr(payload)+'''\nfrom pathlib import Path
import json,subprocess,hashlib
ex=Path(d['execution']);server=Path(d['server']);base=server/'results/straight_tube_exp1/f10ab';logs=base/'logs'
assert (ex/'SOURCE_F10AB_COMMIT.txt').read_text().strip()==d['commit']
for name,content in payload.items():
 p=base/name;p.parent.mkdir(parents=True,exist_ok=True)
 if p.exists():assert p.read_text()==content
 else:p.write_text(content)
script=logs/(d['label']+'_run.sh');assert not script.exists()
python=str(server/('.postvenv' if d['post'] else '.venv')/'bin/python')
module='scripts.straight_tube_exp1.f10ab.'+('post' if d['post'] else 'run')
script.write_text('#!/bin/bash\\nset -euo pipefail\\nexport OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg PYTHONPATH=src:.\\nexec '+python+' -m '+module+' '+d['stage']+'\\n')
d['wrapper_sha256']=hashlib.sha256(script.read_bytes()).hexdigest()
conf=Path('/etc/supervisor/conf.d')/(d['name']+'.conf');assert not conf.exists()
conf.write_text('[program:'+d['name']+']\\ncommand=/usr/bin/timeout --signal=TERM --kill-after=15s '+str(d['seconds'])+' /bin/bash '+str(script)+'\\ndirectory='+str(ex)+'\\nautostart=false\\nautorestart=false\\nstartsecs=0\\nstopasgroup=true\\nkillasgroup=true\\nstdout_logfile='+str(logs/(d['label']+'.log'))+'\\nstderr_logfile='+str(logs/(d['label']+'.stderr.log'))+'\\nstdout_logfile_maxbytes=0\\nstderr_logfile_maxbytes=0\\n')
(logs/(d['label']+'_EXECUTION.json')).write_text(json.dumps(d,indent=2));(logs/(d['label']+'_supervisor.conf')).write_text(conf.read_text())
for a in [['reread'],['update',d['name']],['start',d['name']]]:subprocess.run(['supervisorctl',*a],check=True)
print(json.dumps(d))
'''))
 write(BASE/'evidence'/(label+'_EXECUTION.json'),d)
if __name__=='__main__':main()
