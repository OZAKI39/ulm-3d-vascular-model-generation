"""Finish archives, visual evidence and artifact-only tests after the unique full run."""
import json,os,subprocess,sys,time
from remote import ROOT
R=ROOT/'reports/sv1_3p';S=ROOT/'scripts/sv13p'
while not (R/'winner_steady_candidate.json').exists():
 if (R/'full_monitor_failure.json').exists():raise SystemExit('Monitor failure: preserve evidence and diagnose; no automatic repeat')
 time.sleep(5)
assert json.loads((R/'winner_steady_candidate.json').read_text())['status']=='PASS'
def call(script,*args):
 print('CLOSE '+script,flush=True);subprocess.run([sys.executable,'-B',S/script,*args],check=True)
call('invoke.py','build_evidence_close_remote.py');call('invoke.py','final_audit_remote.py');call('sync_remote.py');call('mirror_native.py');call('audit_preservation.py');call('figures.py')
log=ROOT/'logs/sv1_3p/pytest_full.log';assert not log.exists()
with log.open('x') as out:
 p=subprocess.run([sys.executable,'-B','-m','pytest','-q','-p','no:cacheprovider','--junitxml='+str(R/'pytest_full.xml')],cwd=ROOT,env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1'),stdout=out,stderr=subprocess.STDOUT)
(R/'pytest_execution.json').write_text(json.dumps(dict(exit_code=p.returncode,log=str(log.relative_to(ROOT)),historical_failures_preserved=True),indent=2)+'\n')
call('test_summary.py');call('audit_preservation.py');call('write_report.py');call('terminal_summary.py')
print('Closure prepared. Review plot rendering and report, then freeze delivery. No further CFD.',flush=True)
