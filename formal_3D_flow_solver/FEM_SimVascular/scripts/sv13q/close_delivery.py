"""Complete read-only audit, native artifact mirror, report and tests after CFD is finished."""
import json,subprocess,sys,time
from pathlib import Path
from remote import ROOT
S=ROOT/'scripts/sv13q';R=ROOT/'reports/sv1_3q';L=ROOT/'logs/sv1_3q'
def run(script,*args):subprocess.run([sys.executable,'-B',S/script,*args],check=True)
w=json.loads((R/'winner.json').read_text());assert not w['full_required'] or json.loads((R/'winner_steady_candidate.json').read_text())['status']=='PASS'
if not w['full_required']:(R/'full_steady_decision.json').write_text(json.dumps(dict(status='NOT_REQUIRED',reason='No >=10% eligible winner; no full CFD started'),indent=2)+'\n')
else:(R/'full_steady_decision.json').write_text(json.dumps(dict(status='PASS',full_runs=1,reason='Healthy early/late winner meets >=10% threshold'),indent=2)+'\n')
run('summarize_cases.py');run('invoke.py','final_audit_remote.py');run('sync_remote.py');run('audit_preservation.py');run('mirror_native.py');run('figures.py');run('write_report.py')
for scope,cmd in [('stage',[sys.executable,'-B','-m','pytest','-q','-p','no:cacheprovider',*map(str,sorted((ROOT/'tests').glob('test_sv13q_*.py'))),'--junitxml='+str(R/'pytest_stage.xml')]),('full',[sys.executable,'-B','-m','pytest','-q','-p','no:cacheprovider','--junitxml='+str(R/'pytest_full.xml')])]:
 start=time.monotonic()
 with (L/('pytest_'+scope+'.log')).open('w') as out:p=subprocess.run(cmd,cwd=ROOT,stdout=out,stderr=subprocess.STDOUT)
 (R/('pytest_'+scope+'_execution.json')).write_text(json.dumps(dict(command=cmd,exit_code=p.returncode,wall_time_s=time.monotonic()-start,no_CFD=True),indent=2)+'\n')
 print('pytest',scope,'exit',p.returncode,flush=True)
 if scope=='stage':assert p.returncode==0,'Fix real acceptance or implementation error; never hide it'
run('test_summary.py');run('write_report.py');run('terminal_summary.py')
print('Ready for visual inspection and final delivery freeze.',flush=True)
