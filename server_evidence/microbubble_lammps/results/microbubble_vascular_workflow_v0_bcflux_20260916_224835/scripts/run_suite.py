from pathlib import Path
import subprocess,sys,json,hashlib
S=Path(__file__).resolve().parents[1];ranks=int(sys.argv[1]);variant=sys.argv[2] if len(sys.argv)>2 else 'release';exe=Path(str(S).replace('/results/','/work/'))/'build/workflow_lmp'
cases=['LOW','MEDIUM','CAPACITY_STRESS_LIMIT','FLUX_WEIGHTED']
for case in cases:
    label=f'{case}_{variant}_mpi{ranks}';out=S/'runs'/label
    if out.exists():
        receipt=json.loads((out/'EXECUTION_RECEIPT.json').read_text())
        assert receipt['returncode']==0 and receipt['binary_sha256']==hashlib.sha256(exe.read_bytes()).hexdigest()
        print('REUSE_COMPLETED',label,flush=True)
    else:subprocess.run([sys.executable,str(S/'scripts/run_case.py'),case,str(ranks),label],check=True)
print('SUITE_EXECUTED_MPI',ranks,flush=True)
