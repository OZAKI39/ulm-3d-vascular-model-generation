"""Independent execution experiment, only replace ASM subdomain ILU with exact LU."""
from case_common import *
import shutil,time
source=V/'stage3/vessel_baseline_cpu_mpi8';dest=V/'stage3/execution_checks/vessel_baseline_cpu8_subLU'
assert not dest.exists();dest.mkdir(parents=True);shutil.copytree(source/'SV_MESH',dest/'SV_MESH');(dest/'run').mkdir();shutil.copy2(source/'run/solver.xml',dest/'run/solver.xml')
s=(source/'run/PETSC_OPTIONS.txt').read_text();assert '-sub_pc_type ilu -sub_pc_factor_levels 2' in s
s=s.replace('-sub_pc_type ilu -sub_pc_factor_levels 2','-sub_pc_type lu');(dest/'run/PETSC_OPTIONS.txt').write_text(s)
policy=json.loads((source/'policy.json').read_text());policy.update(case=dest.name,kind='linear_algebra_performance_and_matched_field_check',maximum_wall_time_s=1200,linear_algebra_backend='CPU8 GMRES ASM overlap2; exact sequential subdomain LU instead of ILU2',scientific_comparison_result=False);dump(dest/'policy.json',policy);lock_case(dest)
dump(dest/'reports/preregistration.json',dict(registered_unix=time.time(),reference=str(source),same_mesh_and_XML=all(sha(p)==sha(dest/'SV_MESH'/p.relative_to(source/'SV_MESH')) for p in (source/'SV_MESH').rglob('*') if p.is_file()) and sha(source/'run/solver.xml')==sha(dest/'run/solver.xml'),change='ASM subdomain ILU2 to LU only; not a PDE/discretization change',matched_saved_step=5,maximum_relative_L2_difference_u_p_WSS=1e-9,criteria='All accepted linear corrections converged at original tolerance; compare actual velocity,pressure,WSS with original CPU8 at same step. May stop experiment after matched snapshot; not a steady CFD validation case.',concurrent_run='CPU8 medium continues; raw wall-clock timing is not isolated benchmarking'))
print(dest)
