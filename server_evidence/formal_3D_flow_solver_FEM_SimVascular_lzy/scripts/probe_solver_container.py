#!/usr/bin/env python3
"""Attempt the requested official image pull and record unsupported runtime evidence."""
import json
import platform
import shutil
import subprocess
import sys
import time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from sv_validation.provenance import now,write_json
site=sys.argv[1]
args=['docker','pull','simvascular/solver:latest'];start=time.monotonic()
try:
    p=subprocess.run(args,capture_output=True,text=True,timeout=45)
    result={'returncode':p.returncode,'stdout':p.stdout,'stderr':p.stderr}
except (OSError,subprocess.TimeoutExpired) as exc:result={'returncode':None,'error':str(exc)}
result.update(timestamp=now(),site=site,hostname=platform.node(),command=args,elapsed_s=time.monotonic()-start,
              docker_path=shutil.which('docker'),image_pulled=False,container_started=False)
write_json(ROOT/f'outputs/sv0/environment/{site}_solver_pull.json',result)
print(json.dumps(result,indent=2))
