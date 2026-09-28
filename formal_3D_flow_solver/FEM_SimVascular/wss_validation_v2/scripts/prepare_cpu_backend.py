"""Identical CFD input, PETSc CPU matrix/vector backend only; explicit comparison."""
from pathlib import Path
import argparse,json,shutil
from case_common import *
p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--dest',type=Path,required=True);a=p.parse_args();src=a.source.resolve();dst=a.dest.resolve();assert dst.is_relative_to(V) and not dst.exists();dst.mkdir(parents=True)
shutil.copytree(src/'SV_MESH',dst/'SV_MESH');(dst/'run').mkdir()
for name in ['solver.xml','inlet_profile.txt','outlet_traction.vtp']:
 if (src/'run'/name).exists():shutil.copy2(src/'run'/name,dst/'run'/name)
if (src/'analytical_nodal.npz').exists():shutil.copy2(src/'analytical_nodal.npz',dst/'analytical_nodal.npz')
old=(src/'run/PETSC_OPTIONS.txt').read_text();assert '-mat_type aijcusparse -vec_type cuda' in old;new=old.replace('-mat_type aijcusparse -vec_type cuda','-mat_type aij -vec_type standard');(dst/'run/PETSC_OPTIONS.txt').write_text(new)
pol=json.loads((src/'policy.json').read_text());pol.update(case=dst.name,linear_algebra_backend='CPU PETSc aij/standard; same KSP/PC tolerances and physical method',backend_reference_case=src.name);dump(dst/'policy.json',pol);lock_case(dst)
dump(dst/'reports/backend_only_change.json',dict(reference=str(src),same_solver_xml=sha(src/'run/solver.xml')==sha(dst/'run/solver.xml'),same_mesh={str(p.relative_to(src/'SV_MESH')):sha(p)==sha(dst/'SV_MESH'/p.relative_to(src/'SV_MESH')) for p in (src/'SV_MESH').rglob('*') if p.is_file()},options_before=old,options_after=new,numerical_criteria_unchanged=True))
print(dst)
