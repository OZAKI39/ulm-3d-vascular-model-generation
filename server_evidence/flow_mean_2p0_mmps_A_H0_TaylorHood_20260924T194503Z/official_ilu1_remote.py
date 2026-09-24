"""Solver-only RCM ordering check on the official pinned P2/P1 smoke."""
from pathlib import Path
import json,subprocess,os,time,sys,re,argparse
parser=argparse.ArgumentParser();parser.add_argument('root');parser.add_argument('--label',default='N004_shift_rcm');parser.add_argument('--factor-levels',default=2,type=int);a=parser.parse_args()
r=Path(a.root);b=json.loads((r/'audit/remote_protection_manifest.json').read_text())['build']
out=r/'official_smoke_lfs'/a.label;out.mkdir()
(out/'solver.xml').write_bytes((r/'official_smoke_lfs/N004_shift/solver.xml').read_bytes())
extra=' -sub_pc_factor_shift_type nonzero -sub_pc_factor_shift_amount 1e-6 -sub_pc_factor_mat_ordering_type rcm'
extra+=' -sub_pc_factor_levels '+str(a.factor_levels)
opts=(r/'baseline_PETSC_OPTIONS.txt').read_text().strip()+extra;(out/'PETSC_OPTIONS.txt').write_text(opts+'\n')
env={k:os.environ[k] for k in ('HOME','USER','LOGNAME','LANG') if k in os.environ}
env.update(PATH='/usr/bin:/bin:/usr/local/cuda/bin',LC_ALL='C',OMP_NUM_THREADS='1',CUDA_VISIBLE_DEVICES='0',LD_LIBRARY_PATH=b['runtime_library_path'],PETSC_OPTIONS=opts)
w=b['PETSc_build'];cmd=[w['candidate_wrapper'],w['launcher'],'-n','1',w['candidate_wrapper'],b['executable'],'solver.xml'];start=time.time()
with (out/'solver.log').open('w') as f:p=subprocess.run(cmd,cwd=out,env=env,stdout=f,stderr=subprocess.STDOUT,timeout=180)
text=(out/'solver.log').read_text();solves=len(re.findall('SV13Q_END',text));failed=len(re.findall(r'SV13Q_END[^\n]*reason=-',text))
d=dict(exit_code=p.returncode,elapsed_seconds=time.time()-start,linear_solves=solves,linear_failures=failed,
    gpu_matrix='seqaijcusparse' in text,gpu_vector='seqcuda' in text,ordering_verified='matrix ordering: rcm' in text,solver_only_added=extra)
d['status']='PASS' if p.returncode==0 and solves==200 and not failed and d['gpu_matrix'] and d['gpu_vector'] and d['ordering_verified'] else 'FAIL'
(out/'execution.json').write_text(json.dumps(d,indent=2)+'\n');print(json.dumps(d),flush=True)
