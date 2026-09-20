"""Execute, mirror, independently reload and validate a fresh flow case."""
import subprocess,sys
from remote import ROOT
S=ROOT/'scripts/sv13n'
assert sys.argv[1]!='NEW_PETSC_CPU_PROOF_20','Fixed20 strategy superseded by USER_STRATEGY_UPDATE.txt'
assert not sys.argv[1].startswith(('BENCH_','PROFILE_')),'User cancelled repeated benchmarks and dedicated profiling'
assert sys.argv[1] not in ('NEW_CPU_SHORT','NEW_GPU_SHORT'),'Short science comparisons deferred by latest user decision'
run=subprocess.run([sys.executable,'-B',S/'invoke.py','flow_remote.py',*sys.argv[1:]])
subprocess.run([sys.executable,'-B',S/'fetch_case.py',sys.argv[1]],check=True)
raise SystemExit(subprocess.call([sys.executable,'-B',S/'accept_case.py',sys.argv[1]]))
