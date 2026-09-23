from pathlib import Path
import datetime,json,shlex,subprocess,time
SSH=['ssh','-p','4159','-i','/home/lzy/.ssh/vast_step3b_ed25519','-o','IdentitiesOnly=yes','-o','BatchMode=yes','-o','ConnectTimeout=15','root@50.115.148.16']
CODE='''from pathlib import Path
import json,subprocess
root=Path('/workspace/particle8_2a_20260922')
compute=(root/'results/COMPUTE_AND_MEDIA_COMPLETE').exists()
quality=(root/'results/FINAL_QUALITY_COMPLETE').exists()
log=root/('bdae7ed73ff51049768fbcaee371de6f76a6543c/final_quality.log' if compute else 'fbc4a143212a46dd74069c44dd1207c9b6d29978/remaining.log')
lines=log.read_text().splitlines() if log.exists() else []
session='particle82a_final_quality' if compute else 'particle82a_remaining'
active=subprocess.run(['tmux','has-session','-t',session],capture_output=True).returncode==0
print(json.dumps(dict(progress=lines[-1] if lines else 'WAITING',compute_complete=compute,quality_complete=quality,active=active)))
'''
while True:
    t=datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')
    p=subprocess.run(SSH+['python3 -c '+shlex.quote(CODE)],capture_output=True,text=True,timeout=30)
    if p.returncode:print(t,'SSH_CHECK_FAILED',p.returncode,flush=True)
    else:
        value=json.loads(p.stdout);print(t,json.dumps(value),flush=True)
        if value['quality_complete'] or not value['active']:break
    time.sleep(50)
