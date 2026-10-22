"""Publish an immutable local Git snapshot to the user's server, no GitHub push."""
import subprocess,json,argparse,shlex
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];SERVER='/workspace/f10w_vision_20261017'
SSH=['ssh','-i','/home/lzy/.ssh/vast_step3b_ed25519','-p','15283','root@95.3.33.46'];SCP=['scp','-q','-i','/home/lzy/.ssh/vast_step3b_ed25519','-P','15283']
p=argparse.ArgumentParser();p.add_argument('--jobs',type=Path,required=True);a=p.parse_args();jobs=json.loads(a.jobs.read_text())
repo=subprocess.check_output(['git','rev-parse','--show-toplevel'],cwd=ROOT,text=True).strip();head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip();short=head[:9];base='external_connect/vision_tracking_research/'
paths=[base+x for x in ['scripts','protocol.json','evidence/input_identity.json','collective_motion_v1/src','collective_motion_v1/configs','collective_motion_v2/src','collective_motion_v2/scripts','collective_motion_v2/configs','collective_motion_v3/src','collective_motion_v3/scripts','collective_motion_v3/configs','collective_motion_v3/tests']]+['external_connect/CHRONOS_control']
folder=ROOT/'evidence/deployment'/short;folder.mkdir(parents=True,exist_ok=True);archive=folder/'source.tar'
with archive.open('wb') as f:subprocess.run(['git','archive',head,*paths],cwd=repo,stdout=f,check=True)
identity=folder/'EXECUTION.json';identity.write_text(json.dumps(dict(git_head=head,git_branch='dev/f10w-pretrained-full-video-20261017'),indent=2)+'\n')
subprocess.run(SSH+[f'mkdir -p {SERVER}/setup/{short} {SERVER}/snapshots/{short}'],check=True)
subprocess.run(SCP+[str(archive),str(identity),f'root@95.3.33.46:{SERVER}/setup/{short}/'],check=True)
module=f'{SERVER}/snapshots/{short}/'+base+'collective_motion_v3'
subprocess.run(SSH+[f'test ! -f {module}/EXECUTION.json && tar -xf {SERVER}/setup/{short}/source.tar -C {SERVER}/snapshots/{short} && cp {SERVER}/setup/{short}/EXECUTION.json {module}/EXECUTION.json'],check=True)
for job in jobs:
    name=job['name'];assert all(c.isalnum() or c=='_' for c in name)
    run=folder/(name+'.sh');conf=folder/(name+'.conf')
    run.write_text('#!/bin/bash\nset -euo pipefail\nsource '+SERVER+'/.venv/bin/activate\nexport OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=2\nexport F10W_VIDEO=/workspace/f10v_vision_20261017/setup/f10v_source.mp4.complete\nexport F10W_OUTPUT_ROOT='+shlex.quote(job.get('output',SERVER+'/runs/models_v1'))+'\nexport F10W_THIRD_PARTY='+SERVER+'/third_party\nexport F10W_MODELS='+SERVER+'/models\ncd '+module+'\n'+job['command']+'\n')
    conf.write_text(f'[program:{name}]\ncommand=/usr/bin/timeout {job.get("timeout",1800)} /bin/bash {SERVER}/setup/{short}/{name}.sh\nautostart=false\nautorestart=false\nstartsecs=0\nstopasgroup=true\nkillasgroup=true\nstdout_logfile={SERVER}/logs/{name}_stdout.log\nstderr_logfile={SERVER}/logs/{name}_stderr.log\n')
    subprocess.run(SCP+[str(run),f'root@95.3.33.46:{SERVER}/setup/{short}/'],check=True)
    subprocess.run(SCP+[str(conf),'root@95.3.33.46:/etc/supervisor/conf.d/'],check=True)
subprocess.run(SSH+['supervisorctl reread'],check=True)
for job in jobs:
    subprocess.run(SSH+['supervisorctl update '+job['name']],check=True);subprocess.run(SSH+['supervisorctl start '+job['name']],check=True)
print('Deployed',head,[j['name'] for j in jobs])
