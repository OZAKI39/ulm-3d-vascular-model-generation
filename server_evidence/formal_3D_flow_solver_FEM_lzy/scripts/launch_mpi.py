#!/usr/bin/env python3
import argparse
import json
import os
import subprocess
import time
from pathlib import Path

p=argparse.ArgumentParser()
p.add_argument("--ranks",type=int,required=True)
p.add_argument("--output",required=True)
a=p.parse_args()
root=Path(__file__).resolve().parents[1]
prefix=root/"remote/.env"
command=[str(prefix/"bin/mpiexec"),"-n",str(a.ranks),str(prefix/"bin/python"),"-B",
         "tests/test_dolfinx_remote_smoke.py","--output",a.output]
env=os.environ.copy()
env["FEM3D_LAUNCH_EPOCH"]=str(time.time())
start=time.perf_counter()
completed=subprocess.run(command,env=env,timeout=180)
elapsed=time.perf_counter()-start
result={"command":command,"wall_time_s":elapsed,"returncode":completed.returncode,"ranks":a.ranks}
out=Path(a.output).with_name(f"smoke_r{a.ranks}_launch.json")
out.parent.mkdir(parents=True,exist_ok=True)
out.write_text(json.dumps(result,indent=2)+"\n")
raise SystemExit(completed.returncode)
