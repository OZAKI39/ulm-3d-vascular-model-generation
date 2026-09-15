"""CPU controller service: evaluator subprocesses are never forked from CUDA solver."""
import subprocess,threading,os
from pathlib import Path
from remote_common import write
def serve_evaluations(root,done):
    root=Path(root);last=0
    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1',TMPDIR=str(root/'tmp'))
    while not done.is_set():
        req=root/'diagnostics/evaluation_request.txt'
        try:step=int(req.read_text().strip()) if req.exists() else 0
        except (OSError,ValueError):step=0
        if step>last:
            try:
                if step!=last+5000:raise RuntimeError('Nonconsecutive evaluator request')
                with (root/'logs/online_evaluator.log').open('a') as log:
                    p=subprocess.run(['/usr/bin/python3','-B',str(root/'scripts/monitor_online.py'),str(root),str(root),str(step)],env=env,stdout=log,stderr=subprocess.STDOUT,timeout=280)
                if p.returncode!=0:raise RuntimeError('Evaluator returned '+str(p.returncode))
                last=step
            except Exception as e:
                write(root/'diagnostics/EXTERNAL_EVALUATOR_FAILURE.json',dict(step=step,error=str(e)))
                (root/'diagnostics/evaluation_error.txt').write_text(str(step)+'\n')
                return
        done.wait(.05)
