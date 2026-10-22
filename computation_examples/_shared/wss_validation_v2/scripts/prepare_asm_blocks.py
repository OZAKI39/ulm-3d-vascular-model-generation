"""Execution variant: fixed PDE, smaller ASM subdomains, unchanged residual gates."""
import argparse,shutil,json
from case_common import *
p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--dest',type=Path,required=True);p.add_argument('--blocks',type=int,default=16);a=p.parse_args();src=a.source.resolve();dst=a.dest.resolve();assert not dst.exists();dst.mkdir(parents=True);shutil.copytree(src/'SV_MESH',dst/'SV_MESH');(dst/'run').mkdir()
for name in ['solver.xml','inlet_profile.txt','outlet_traction.vtp','PETSC_OPTIONS.txt']:
 if (src/'run'/name).exists():shutil.copy2(src/'run'/name,dst/'run'/name)
if (src/'analytical_nodal.npz').exists():shutil.copy2(src/'analytical_nodal.npz',dst/'analytical_nodal.npz')
opts=(dst/'run/PETSC_OPTIONS.txt').read_text().strip();assert '-pc_asm_local_blocks' not in opts;opts+=' -pc_asm_local_blocks '+str(a.blocks);(dst/'run/PETSC_OPTIONS.txt').write_text(opts+'\n');pol=json.loads((src/'policy.json').read_text());pol.update(case=dst.name,asm_local_blocks=a.blocks,execution_variant_reference=src.name);dump(dst/'policy.json',pol);lock_case(dst);dump(dst/'reports/execution_variant.json',dict(reference=str(src),only_change='PETSc ASM local block count; residual tolerances unchanged',blocks=a.blocks,same_solver_xml=sha(src/'run/solver.xml')==sha(dst/'run/solver.xml'),solver_and_discretization_unchanged=True))
