"""Independent CPU8 -> GPU8 execution experiment; identical CFD and KSP inputs."""
import argparse,shutil,json
from case_common import *
p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--dest',type=Path,required=True);p.add_argument('--ranks',type=int);a=p.parse_args()
src=a.source.resolve();dst=a.dest.resolve();assert dst.is_relative_to(V) and not dst.exists();dst.mkdir(parents=True)
shutil.copytree(src/'SV_MESH',dst/'SV_MESH');(dst/'run').mkdir();shutil.copy2(src/'run/solver.xml',dst/'run/solver.xml')
old=(src/'run/PETSC_OPTIONS.txt').read_text();assert '-mat_type aij -vec_type standard' in old
new=old.replace('-mat_type aij -vec_type standard','-mat_type aijcusparse -vec_type cuda');(dst/'run/PETSC_OPTIONS.txt').write_text(new)
pol=json.loads((src/'policy.json').read_text());pol.update(case=dst.name,linear_algebra_backend='GPU PETSc aijcusparse/cuda; same KSP/PC and physical method',backend_reference_case=src.name)
if a.ranks is not None:pol['MPI_ranks']=a.ranks
dump(dst/'policy.json',pol);lock_case(dst)
dump(dst/'reports/backend_only_change.json',dict(reference=str(src),same_solver_xml=sha(src/'run/solver.xml')==sha(dst/'run/solver.xml'),same_mesh={str(p.relative_to(src/'SV_MESH')):sha(p)==sha(dst/'SV_MESH'/p.relative_to(src/'SV_MESH')) for p in (src/'SV_MESH').rglob('*') if p.is_file()},options_before=old,options_after=new,numerical_criteria_unchanged=True))
print(dst)
