"""Immutable Git execution snapshots and bounded, individually named Supervisor jobs."""
import argparse,json,subprocess,shlex
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
SERVER='/workspace/f10x_piv_20261017'
SSH=['ssh','-i','/home/lzy/.ssh/vast_step3b_ed25519','-o','BatchMode=yes','-p','15283','root@95.3.33.46']
SCP=['scp','-q','-i','/home/lzy/.ssh/vast_step3b_ed25519','-P','15283']
def main():
 p=argparse.ArgumentParser();p.add_argument('name');p.add_argument('command');p.add_argument('--timeout',type=int,default=1200);p.add_argument('--env',default='.venv');p.add_argument('--output',default='benchmark_v1');a=p.parse_args()
 assert a.name.startswith('f10x_') and all(c.isalnum() or c=='_' for c in a.name)
 repo=subprocess.check_output(['git','rev-parse','--show-toplevel'],cwd=ROOT,text=True).strip();head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip();short=head[:9]
 folder=ROOT/'evidence/deployment'/short;folder.mkdir(parents=True,exist_ok=True)
 archive=folder/'source.tar';base='external_connect/vision_tracking_research/'
 paths=[base+x for x in ['scripts','protocol.json','evidence/input_identity.json','collective_motion_v1/src','collective_motion_v1/configs','collective_motion_v2/src','collective_motion_v2/scripts','collective_motion_v2/configs','collective_motion_v3/src','collective_motion_v3/scripts','collective_motion_v3/configs','collective_motion_v4']]+['external_connect/CHRONOS_control']
 with archive.open('wb') as f:subprocess.run(['git','archive',head,*paths],cwd=repo,stdout=f,check=True)
 module=f'{SERVER}/snapshots/{short}/{base}collective_motion_v4'
 subprocess.run(SSH+[f'mkdir -p {SERVER}/setup/{short} {SERVER}/snapshots/{short} {SERVER}/logs'],check=True)
 subprocess.run(SCP+[str(archive),f'root@95.3.33.46:{SERVER}/setup/{short}/'],check=True)
 ident=json.dumps({'git_head':head,'branch':subprocess.check_output(['git','branch','--show-current'],cwd=repo,text=True).strip()})
 subprocess.run(SSH+[f'if test ! -f {module}/EXECUTION.json; then tar -xf {SERVER}/setup/{short}/source.tar -C {SERVER}/snapshots/{short}; printf %s {shlex.quote(ident)} > {module}/EXECUTION.json; fi'],check=True)
 script=folder/(a.name+'.sh');conf=folder/(a.name+'.conf')
 script.write_text('#!/bin/bash\nset -euo pipefail\nexport OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=2 RAYON_NUM_THREADS=4 MPLBACKEND=Agg\nexport F10X_VIDEO=/workspace/f10v_vision_20261017/setup/f10v_source.mp4.complete\nexport F10X_SERVER='+SERVER+'\nexport F10X_OUTPUT='+SERVER+'/runs/'+a.output+'\nsource '+SERVER+'/'+a.env+'/bin/activate\ncd '+module+'\n'+a.command+'\n')
 if a.command=='SETUP':script.write_text('#!/bin/bash\nset -euo pipefail\nbash '+module+'/scripts/setup_remote.sh\n')
 conf.write_text(f'[program:{a.name}]\ncommand=/usr/bin/timeout {a.timeout} /bin/bash {SERVER}/setup/{short}/{a.name}.sh\nautostart=false\nautorestart=false\nstartsecs=0\nstopasgroup=true\nkillasgroup=true\nstdout_logfile={SERVER}/logs/{a.name}_stdout.log\nstderr_logfile={SERVER}/logs/{a.name}_stderr.log\n')
 subprocess.run(SCP+[str(script),f'root@95.3.33.46:{SERVER}/setup/{short}/'],check=True)
 subprocess.run(SSH+[f'test ! -e /etc/supervisor/conf.d/{a.name}.conf'],check=True)
 subprocess.run(SCP+[str(conf),'root@95.3.33.46:/etc/supervisor/conf.d/'],check=True)
 for cmd in ['supervisorctl reread','supervisorctl update '+a.name,'supervisorctl start '+a.name]:subprocess.run(SSH+[cmd],check=True)
 print('EXECUTION',head,'JOB',a.name,'MODULE',module)
if __name__=='__main__':main()
