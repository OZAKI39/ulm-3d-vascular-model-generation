"""Bounded second ASM experiment: higher incomplete factorization fill, same equations."""
from case_common import *
import shutil,time
source=V/'stage3/vessel_baseline_cpu_mpi8';dest=V/'stage3/execution_checks/vessel_baseline_cpu8_ILU6'
assert not dest.exists();dest.mkdir(parents=True);shutil.copytree(source/'SV_MESH',dest/'SV_MESH');(dest/'run').mkdir();shutil.copy2(source/'run/solver.xml',dest/'run/solver.xml')
s=(source/'run/PETSC_OPTIONS.txt').read_text();assert '-sub_pc_factor_levels 2' in s
(dest/'run/PETSC_OPTIONS.txt').write_text(s.replace('-sub_pc_factor_levels 2','-sub_pc_factor_levels 6'))
policy=json.loads((source/'policy.json').read_text());policy.update(case=dest.name,kind='linear_algebra_performance_and_matched_field_check',maximum_wall_time_s=1200,linear_algebra_backend='CPU8 GMRES ASM overlap2; subdomain ILU6 instead of ILU2',scientific_comparison_result=False);dump(dest/'policy.json',policy);lock_case(dest)
dump(dest/'reports/preregistration.json',dict(registered_unix=time.time(),reference=str(source),change='Only ASM subdomain ILU fill level 2 to 6',matched_saved_step=5,maximum_relative_L2_difference_u_p_WSS=1e-9,criteria='Original residual target and same physical configuration. Only an execution experiment; not an accepted steady flow. Compare field and factor memory/iterations before adoption.',reason='Exact LU reduced iterations but substantial fill observed; bounded alternative checks less memory-intensive stronger incomplete factorization',concurrent_medium_case=True))
print(dest)
