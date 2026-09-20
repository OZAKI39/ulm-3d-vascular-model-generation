"""Only integrity, archives, figures and artifact tests after the unique CFD run."""
import hashlib,json,shlex,subprocess,sys,time
from pathlib import Path
from remote import ROOT,REMOTE,ssh_prefix
R=ROOT/'reports/sv1_3o';S=ROOT/'scripts/sv13o'
while not (R/'optimized_steady_candidate.json').exists():
 if (R/'full_monitor_failure.json').exists():raise SystemExit('Full monitor failure requires diagnosis; no closure')
 time.sleep(5)
assert json.loads((R/'optimized_steady_candidate.json').read_text())['status']=='PASS'
def call(script,*args):
 print('CLOSE '+script,flush=True);subprocess.run([sys.executable,'-B',S/script,*args],check=True)
call('invoke.py','final_audit_remote.py');call('sync_remote.py')
for name in ('remote_preservation','final_input_integrity','svmp_build'):
 (R/(name+'.json')).write_bytes((R/'remote'/(name+'.json')).read_bytes())
d=json.loads((R/'remote/native_artifacts.json').read_text());local=ROOT/'outputs/sv1_3o/native';local.mkdir(exist_ok=True);ssh=ssh_prefix()
for f in d['files']:
 p=local/Path(f['remote_path']).name;subprocess.run(['rsync','-a','-e',shlex.join(ssh[:-1]),ssh[-1]+':'+f['remote_path'],str(p)],check=True)
 assert hashlib.sha256(p.read_bytes()).hexdigest()==f['sha256'];f['local_path']=str(p.relative_to(ROOT))
(R/'native_artifact_mirror.json').write_text(json.dumps(d,indent=2)+'\n')
call('audit_preservation.py');call('output_consistency.py');call('profile_summary.py');call('figures.py')
log=ROOT/'logs/sv1_3o/pytest_full.log';assert not log.exists()
with log.open('x') as out:
 p=subprocess.run([sys.executable,'-B','-m','pytest','-q','-p','no:cacheprovider','--junitxml='+str(R/'pytest_full.xml')],cwd=ROOT,stdout=out,stderr=subprocess.STDOUT)
(R/'pytest_execution.json').write_text(json.dumps(dict(exit_code=p.returncode,log=str(log.relative_to(ROOT)),historical_failures_preserved=True),indent=2)+'\n')
call('audit_preservation.py');print('Full result and provenance closed. Review figures/tests, then write report. No more CFD.',flush=True)
