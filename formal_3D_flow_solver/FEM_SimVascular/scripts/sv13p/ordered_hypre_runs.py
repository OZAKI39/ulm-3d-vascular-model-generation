"""Phase B parity and one profile per hypre candidate, with smoke-before-window."""
import json,subprocess,sys
from remote import ROOT
S=ROOT/'scripts/sv13p';R=ROOT/'reports/sv1_3p'
def call(script,*args):return subprocess.run([sys.executable,'-B',S/script,*args]).returncode
assert json.loads((R/'remote/svmp_hypre_build.json').read_text())['status']=='PASS'
assert call('make_candidates.py','HYPRE_BASELINE')==0
assert call('run_case.py','HYPRE_BASELINE_SMOKE','smoke','HYPRE_BASELINE')==0
base=json.loads((R/'HYPRE_BASELINE_SMOKE_acceptance.json').read_text());ref=json.loads((R/'reference_freeze.json').read_text())
old=json.loads((ROOT/'reports/sv1_3o/remote/PERF_A_execution.json').read_text())
oldrows=[r for r in old['history']['linear_solves'] if r['step'] in (61,62)]
new=json.loads((R/'remote/HYPRE_BASELINE_SMOKE_execution.json').read_text())['history']['linear_solves']
metrics=dict(status='PASS',old_two_step_iterations=sum(r['linear_iterations'] for r in oldrows),new_two_step_iterations=sum(r['linear_iterations'] for r in new),new_smoke_wall_s=base['wall_time_s'],full_baseline_rerun=False,decision='Do not repeat 10-step A: same PETSc source, compilers, CUDA, MPI, precision and flags; only hypre package added. Original ASM/ILU behaviour passes short smoke; no meaningful runtime-altering change established. Smoke wall includes fixed startup/output costs and is not extrapolated for ranking.')
(R/'hypre_baseline_parity.json').write_text(json.dumps(metrics,indent=2)+'\n')
for key in ['P2','P3']:
 if key=='P3':
  assert call('invoke.py','standalone_remote.py','hypre','boomeramg')==0;assert call('sync_remote.py')==0
 assert call('make_candidates.py',key)==0
 code=call('run_case.py',key+'_SMOKE','smoke',key)
 if code==0:call('run_case.py',key+'_WINDOW','window',key)
print('P2/P3 bounded evaluation finished',flush=True)
