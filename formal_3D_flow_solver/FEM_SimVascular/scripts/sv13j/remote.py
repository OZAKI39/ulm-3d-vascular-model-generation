"""Explicit native SSH transport; no container engines or inherited credentials."""
import json,shlex,subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
REMOTE='/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3j'
def ssh_prefix():
    old=json.loads(subprocess.check_output(['git','-C',str(ROOT),'show','HEAD:reports/sv0/evidence/remote_init.json'],text=True))
    return old['command'][:-1]
def python(code,timeout=60):
    return subprocess.run(ssh_prefix()+['/usr/bin/python3 -B -'],input=code,capture_output=True,text=True,timeout=timeout)
def upload(local,relative):
    content=Path(local).read_bytes()
    import base64
    code=f"from pathlib import Path\nimport base64\np=Path({(REMOTE+'/'+relative)!r});p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(base64.b64decode({base64.b64encode(content).decode()!r}))"
    r=python(code);r.check_returncode()
def fetch(relative,local):
    r=subprocess.run(ssh_prefix()+[shlex.join(['cat',REMOTE+'/'+relative])],capture_output=True,timeout=60)
    r.check_returncode();Path(local).parent.mkdir(parents=True,exist_ok=True);Path(local).write_bytes(r.stdout)

