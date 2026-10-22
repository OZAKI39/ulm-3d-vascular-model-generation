"""Git snapshots and bounded Supervisor jobs, no unrelated services touched."""
import subprocess,tarfile,io,re,json
from .common import *
SSH=['ssh','-i','/home/lzy/.ssh/vast_step3b_ed25519','-o','BatchMode=yes','-p','15283','root@95.3.33.46']
SERVER='/workspace/mevo_acoustic_microbubbles/f10ae_target_20261019'
def call(code):return subprocess.check_output(SSH+['python3 -'],input=code,text=True)
def deploy(module,label,seconds=900):
 assert re.fullmatch('[a-z0-9_]+',label)
 commit=identity()['source_commit'];paths=['src','scripts/straight_tube_exp1','configs','tests/f10aa','tests/f10ab','tests/f10ac','tests/f10ae']
 paths=[x for x in paths if (ROOT/x).exists()]
 assert not subprocess.check_output(['git','status','--porcelain','--',*paths],cwd=ROOT,text=True).strip()
 prefix='computation_examples/mevo_acoustic_microbubbles/'
 raw=subprocess.check_output(['git','archive',commit,*[prefix+x for x in paths]],cwd=ROOT.parents[1]);b=io.BytesIO()
 with tarfile.open(fileobj=io.BytesIO(raw)) as src,tarfile.open(fileobj=b,mode='w:gz') as dst:
  for item in src:
   item.name=item.name[len(prefix):]
   if item.name:dst.addfile(item,src.extractfile(item) if item.isfile() else None)
 archive=BASE/'evidence'/(label+'.tar.gz');archive.parent.mkdir(parents=True,exist_ok=True);archive.write_bytes(b.getvalue())
 ex=SERVER+'/execution/'+commit+'_'+label
 call('from pathlib import Path\nPath('+repr(ex)+').mkdir(parents=True,exist_ok=False)')
 reuse=json.loads(call('from pathlib import Path\nimport json,hashlib,shutil\na=list(Path('+repr(SERVER+'/execution')+').glob('+repr(commit+'*/source.tar.gz')+'))\na=[p for p in a if str(p.parent)!='+repr(ex)+']\np=a[0] if a else None\nif p:shutil.copy2(p,'+repr(ex+'/source.tar.gz')+')\nprint(json.dumps(dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) if p else None))'))
 if not reuse:subprocess.run(['scp','-q','-i',SSH[2],'-P','15283',str(archive),'root@95.3.33.46:'+ex+'/source.tar.gz'],check=True)
 d=dict(ex=ex,server=SERVER,module=module,label=label,seconds=seconds,commit=commit,sha256=reuse['sha256'] if reuse else sha(archive),source_archive_reuse=reuse,name='mevo_f10ae_'+label)
 print(call('d='+repr(d)+r"""
from pathlib import Path
import tarfile,hashlib,json,subprocess
ex=Path(d['ex']);root=Path(d['server']);a=ex/'source.tar.gz';assert hashlib.sha256(a.read_bytes()).hexdigest()==d['sha256']
with tarfile.open(a) as f:f.extractall(ex,filter='data')
r=ex/'results/straight_tube_exp1';r.mkdir(parents=True,exist_ok=True)
for name in ['f10ae']:
 target=root/'results/straight_tube_exp1'/name;target.mkdir(parents=True,exist_ok=True);(r/name).symlink_to(target,target_is_directory=True)
for name in ['f10ac_qd','f10ad_multifrequency']:(r/name).symlink_to('/workspace/mevo_acoustic_microbubbles/f10ac_ad_20261018/results/straight_tube_exp1/'+name,target_is_directory=True)
for name in ['f10n','f10q','f10r','f10s']:(r/name).symlink_to('/workspace/mevo_acoustic_microbubbles/f10s_control_validation/results/straight_tube_exp1/'+name,target_is_directory=True)
for name,folder in [('f10aa','f10aa_simulation_20261018'),('f10ab','f10ab_shape_20261018'),('f10z','f10z_simulation_20261018'),('f10z_scaling','f10z_scaling_20261018')]:
 (r/name).symlink_to('/workspace/mevo_acoustic_microbubbles/'+folder+'/results/straight_tube_exp1/'+name,target_is_directory=True)
(ex/'.venv').symlink_to(root/'.venv',target_is_directory=True)
(ex/'third_party').symlink_to(root/'third_party',target_is_directory=True)
(ex/'SOURCE_F10AE_COMMIT.txt').write_text(d['commit']+'\n')
logs=root/'logs';script=ex/'run.sh'
script.write_text('#!/bin/bash\nset -euo pipefail\nexport OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg PYTHONPATH=src:.\nexec .venv/bin/python -m scripts.straight_tube_exp1.'+d['module']+'\n')
conf=Path('/etc/supervisor/conf.d')/(d['name']+'.conf');assert not conf.exists()
conf.write_text('[program:'+d['name']+']\ncommand=/usr/bin/timeout --signal=TERM --kill-after=15s '+str(d['seconds'])+' /bin/bash '+str(script)+'\ndirectory='+str(ex)+'\nautostart=false\nautorestart=false\nstartsecs=0\nstopasgroup=true\nkillasgroup=true\nstdout_logfile='+str(logs/(d['label']+'.log'))+'\nredirect_stderr=true\nstdout_logfile_maxbytes=0\n')
(logs/(d['label']+'_EXECUTION.json')).write_text(json.dumps(d,indent=2))
(logs/(d['label']+'_supervisor.conf')).write_text(conf.read_text())
for args in [['reread'],['update',d['name']],['start',d['name']]]:subprocess.run(['supervisorctl',*args],check=True)
print(json.dumps(d))
"""))
 write(BASE/'evidence'/(label+'_EXECUTION.json'),d)
if __name__=='__main__':
 import sys
 deploy(sys.argv[1],sys.argv[2],int(sys.argv[3]) if len(sys.argv)>3 else 900)
