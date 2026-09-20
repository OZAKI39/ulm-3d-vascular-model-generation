"""Explicit SSH transport using the existing user-owned instance connection."""
import json,shlex,subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
REMOTE_ROOT='/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3'
def ssh_prefix():
    prior=json.loads(subprocess.check_output(['git','-C',str(ROOT),'show','HEAD:reports/sv0/evidence/remote_init.json'],text=True))
    return prior['command'][:-1]
def run_python(code,timeout=60):
    return subprocess.run(ssh_prefix()+['python3 -B -'],input=code,capture_output=True,text=True,timeout=timeout)
def run_command(args,**kwargs):
    return subprocess.run(ssh_prefix()+[shlex.join(args)],**kwargs)
