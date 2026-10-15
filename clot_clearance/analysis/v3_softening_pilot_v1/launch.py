from common import *
import sys,subprocess,time
name=sys.argv[1];assert name in ['legacy','preflight','ledger','pilot']
env=os.environ.copy();lib=OLDLIB if name=='legacy' else LIB
env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',LD_LIBRARY_PATH=str(lib.parent)+os.pathsep+env.get('LD_LIBRARY_PATH',''),PYTHONPATH=str(NATIVE/'source/python')+os.pathsep+str(R)+os.pathsep+env.get('PYTHONPATH',''))
script={'legacy':'legacy_check.py','preflight':'preflight.py','ledger':'ledger_check.py','pilot':'run_protocol.py'}[name]
start=time.perf_counter();proc=subprocess.run([sys.executable,'-u',str(A/script)],env=env,cwd=A)
dump(O/(name.upper()+'_PROCESS.json'),dict(returncode=proc.returncode,wall_s=time.perf_counter()-start,script=script,library=str(lib),library_sha256=sha(lib)))
sys.exit(proc.returncode)
