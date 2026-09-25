#!/usr/bin/env python3
import os,json,time,subprocess,fcntl,sys
from pathlib import Path
W=Path(__file__).resolve().parents[1]
f=(W/'continue.lock').open('a');fcntl.flock(f,fcntl.LOCK_EX|fcntl.LOCK_NB)
with (W/'CONTINUE_STARTED.json').open('x') as out:json.dump(dict(pid=os.getpid(),unix=time.time(),automatic_retries=0),out)
while True:
 s=json.loads((W/'PREPARATION_STATE.json').read_text())
 if s['state'] in ('PASS','FAIL'):break
 time.sleep(10)
if s['state']!='PASS':
 (W/'VALIDATION_STATE.json').write_text(json.dumps(dict(state='FAIL',phase='PREPARATION_FAILED',error=s.get('error'),preparation_state=s),indent=2))
else:
 p=subprocess.run(['/usr/bin/python3','-u','-B',str(W/'scripts/validation_controller.py')])
print('FINALIZING',flush=True)
p=subprocess.run(['/usr/bin/python3','-u','-B',str(W/'scripts/finalize_result.py')])
sys.exit(p.returncode)
