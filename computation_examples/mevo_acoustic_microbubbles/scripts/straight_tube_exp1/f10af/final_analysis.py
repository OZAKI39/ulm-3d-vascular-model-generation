"""Wait for bounded physics, verify history, run regression tests and analyze."""
import time,subprocess,importlib.metadata as md
from .common import *
def main():
    start=time.perf_counter()
    while not (BASE/'PHYSICS_WORKFLOW_COMPLETE.json').exists():
        if list(BASE.glob('*_WORKFLOW_FAILURE.json')):raise RuntimeError('Physics workflow reported a failure; inspect preserved evidence')
        if time.perf_counter()-start>3400:raise RuntimeError('Analysis wait budget exhausted')
        time.sleep(3)
    from .analyze import main as analyze
    analyze();verify();rows=[];missing=[]
    # Every previously sealed historical output is checked if present on this host.
    for rel,d in read(BASE/'HISTORICAL_AE_SHA256.json')['files'].items():
        p=AE/rel
        if p.exists():rows.append(dict(path=rel,matched=sha(p)==d['sha256']))
        else:missing.append(rel)
    write(BASE/'HISTORICAL_FINAL_AUDIT.json',dict(checked=len(rows),all_present_match=all(r['matched'] for r in rows),mismatches=[r for r in rows if not r['matched']],not_on_this_host=missing,identity=identity()))
    result=subprocess.run([str(ROOT/'.venv/bin/python'),'-m','pytest','tests/f10aa','tests/f10ab','tests/f10ae','tests/f10af','-q'],cwd=ROOT,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    (BASE/'evidence/REGRESSION_TESTS.txt').write_text(result.stdout);write(BASE/'TEST_RESULTS.json',dict(exit_code=result.returncode,stdout_path='evidence/REGRESSION_TESTS.txt',identity=identity()));assert result.returncode==0,result.stdout
    lock=subprocess.check_output(['/usr/local/bin/uv','pip','freeze','--python',str(ROOT/'.venv/bin/python')],text=True);(BASE/'evidence/remote_environment.lock').write_text(lock)
    write(BASE/'ANALYSIS_COMPLETE.json',dict(identity=identity()))
if __name__=='__main__':main()
