"""Incremental off-server copies; no task generation or solver invocation."""
from pathlib import Path
import subprocess,time,json
HERE=Path(__file__).resolve().parents[1]
REMOTE='/workspace/microbubble_roi_only_dt0p5ms_n1500_20260928/'
while True:
    result=subprocess.run(['rsync','-a','--partial','--delay-updates','--exclude=scripts/','--exclude=assets/','--exclude=config.json','--exclude=supervisord.conf','--exclude=*.so','--exclude=__pycache__','--exclude=supervisor.sock','--exclude=supervisord.pid',
        'vast4090:'+REMOTE,str(HERE)+'/'],capture_output=True,text=True)
    with (HERE/'logs/collection.log').open('a') as stream:
        stream.write(json.dumps(dict(unix_s=time.time(),returncode=result.returncode,stderr=result.stderr))+'\n')
    if (HERE/'DELIVERY_COMPLETE.json').exists():
        print('COLLECTED_DELIVERY_COMPLETE',flush=True);break
    if (HERE/'data/DELIVERY_BLOCKED.json').exists():
        print('REMOTE_DELIVERY_BLOCKED',flush=True);break
    if (HERE/'data/progress.json').exists():
        print((HERE/'data/progress.json').read_text(),flush=True)
    time.sleep(45)
