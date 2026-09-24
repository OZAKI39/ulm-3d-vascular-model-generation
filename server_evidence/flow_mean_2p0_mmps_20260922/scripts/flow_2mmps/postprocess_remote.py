"""Run dependent export/render/audit only after the native solve is accepted."""
from pathlib import Path
import sys,time,json,subprocess,traceback
ROOT=Path(__file__).resolve().parents[2]
REPORTS=ROOT/'flow_cases/mean-2p0-mmps/reports'

def state(status,phase,**extra):
    p=REPORTS/'postprocess_status.json';tmp=p.with_suffix('.pending')
    tmp.write_text(json.dumps(dict(status=status,phase=phase,unix_time=time.time(),**extra),indent=2)+'\n');tmp.replace(p)

def main():
    state('WAITING','native solve');deadline=time.time()+14400
    while not (REPORTS/'execution.json').exists():
        assert time.time()<deadline,'Native solve did not finish within waiting budget'
        time.sleep(5)
    execution=json.loads((REPORTS/'execution.json').read_text())
    assert execution['status']=='PASS','Native solver health gate failed'
    for script in ('validate.py','render.py','audit_outputs.py'):
        state('RUNNING',script)
        with (REPORTS/(Path(script).stem+'_remote.log')).open('w') as log:
            result=subprocess.run([sys.executable,str(ROOT/'scripts/flow_2mmps'/script)],cwd=ROOT,
                stdout=log,stderr=subprocess.STDOUT)
        assert result.returncode==0,f'{script} failed; inspect retained log'
    state('PASS','complete')

if __name__=='__main__':
    try:main()
    except Exception as exc:
        state('FAIL','failed',error=str(exc));traceback.print_exc();raise
