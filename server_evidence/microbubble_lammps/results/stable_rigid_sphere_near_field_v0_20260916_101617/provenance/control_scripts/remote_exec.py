from pathlib import Path
import json, subprocess, sys
STATE=json.loads((Path(__file__).parent/'TASK_PATHS.json').read_text())
SSH=['ssh','-o','BatchMode=yes','-o','ConnectTimeout=15','-o','ControlMaster=auto','-o','ControlPersist=300','-o','ControlPath=/tmp/rigid_20260916_101617',STATE['ssh_alias']]
def run(script,timeout=90):
    return subprocess.run(SSH+['python3 -'],input=script,text=True,capture_output=True,timeout=timeout)
if __name__=='__main__':
    result=run(Path(sys.argv[1]).read_text(),timeout=int(sys.argv[2]) if len(sys.argv)>2 else 90)
    print(result.stdout,end='')
    print(result.stderr,end='',file=sys.stderr)
    raise SystemExit(result.returncode)
