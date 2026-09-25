from pathlib import Path
import subprocess,sys,concurrent.futures
S=Path(__file__).resolve().parents[1]
assert str(S).startswith('/workspace/') and 'surfacecontact_' in str(S)
names=sorted(p.stem for p in (S/'configs').glob('SENS_*.cfg'))
def run(name):
 p=subprocess.run([sys.executable,str(S/'scripts/run_case.py'),name,'1'],capture_output=True,text=True)
 print(p.stdout,p.stderr,flush=True);assert p.returncode==0
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:list(ex.map(run,names))
