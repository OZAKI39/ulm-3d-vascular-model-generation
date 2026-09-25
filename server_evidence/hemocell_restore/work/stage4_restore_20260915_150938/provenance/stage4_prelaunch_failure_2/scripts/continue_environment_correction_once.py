#!/usr/bin/env python3
from pathlib import Path
import subprocess,sys,json,time,os,fcntl
W=Path(__file__).resolve().parents[1]
f=(W/'continue.lock').open('a');fcntl.flock(f,fcntl.LOCK_EX|fcntl.LOCK_NB)
with (W/'ENVIRONMENT_CORRECTION_CONTINUE_STARTED.json').open('x') as out:json.dump(dict(pid=os.getpid(),unix=time.time(),reason='Correct help parser: NVHPC lists ccXY supported values89 rather than literalcc89; no solver/build ran before repair',automatic_solver_retries=0),out)
p=subprocess.run(['/usr/bin/python3','-u','-B',str(W/'scripts/validation_controller_v2.py')])
p=subprocess.run(['/usr/bin/python3','-u','-B',str(W/'scripts/finalize_result.py')])
sys.exit(p.returncode)
