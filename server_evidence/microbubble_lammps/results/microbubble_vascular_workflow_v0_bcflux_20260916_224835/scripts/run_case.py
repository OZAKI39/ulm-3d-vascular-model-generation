#!/usr/bin/env python3
"""Run only in this new results stage, using the isolated work binary."""
from pathlib import Path
import sys,subprocess,json,hashlib,time,os,shutil
S=Path(__file__).resolve().parents[1]
assert str(S).startswith('/workspace/microbubble_lammps/results/microbubble_vascular_workflow_v0_bcflux_')
name=sys.argv[1];ranks=int(sys.argv[2]);label=sys.argv[3] if len(sys.argv)>3 else name+'_mpi'+str(ranks)
cfg=S/'configs'/f'{name}.cfg';out=S/'runs'/label;out.mkdir(parents=True)
work=Path(str(S).replace('/results/','/work/'));exe=work/'build/workflow_lmp'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
cmd=['mpirun','--allow-run-as-root','--bind-to','none','--oversubscribe','-n',str(ranks),str(exe),str(cfg)]
shutil.copy2(cfg,out/'case.cfg')
receipt=dict(command=cmd,cwd=str(out),config_sha256=sha(cfg),binary_sha256=sha(exe),library_sha256=sha(work/'build/librigid_math.so'),launch_environment={'HWLOC_COMPONENTS':'-gl','OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1'},started_epoch=time.time(),disclaimer='NOT EXPERIMENTAL CONCENTRATION')
(out/'RUN_STARTED.json').write_text(json.dumps(receipt,indent=2)+'\n')
with (out/'engine.log').open('w') as log:
    p=subprocess.run(cmd,cwd=out,stdout=log,stderr=subprocess.STDOUT,env={**os.environ,'HWLOC_COMPONENTS':'-gl','OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1'},timeout=600)
receipt.update(returncode=p.returncode,finished_epoch=time.time());(out/'EXECUTION_RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(label,'exit',p.returncode);print((out/'engine.log').read_text()[-1200:]);assert p.returncode==0
