"""Continue the ordered evaluation after the single P1 window has closed."""
import json,subprocess,sys,time
from remote import ROOT,upload
R=ROOT/'reports/sv1_3p';S=ROOT/'scripts/sv13p'
while not (R/'P1_WINDOW_acceptance.json').exists():time.sleep(3)
assert (R/'remote/P1_WINDOW_execution.json').exists()
upload(S/'gpu_sparse.c','gpu_sparse.c')
def call(script,*args):return subprocess.run([sys.executable,'-B',S/script,*args]).returncode
code=call('invoke.py','package_build_remote.py','hypre');call('sync_remote.py')
if code:raise SystemExit(code)
code=call('invoke.py','standalone_remote.py','hypre','ilu');call('sync_remote.py')
if code:raise SystemExit(code)
code=call('invoke.py','package_svmp_remote.py','hypre');call('sync_remote.py')
raise SystemExit(code)
