from case_common import *
import shutil,time
src=V/'stage3/execution_checks/vessel_baseline_cpu8_subLU';dst=V/'stage3/execution_checks/vessel_baseline_cpu16_subLU';assert not dst.exists();dst.mkdir(parents=True);shutil.copytree(src/'SV_MESH',dst/'SV_MESH');(dst/'run').mkdir()
for name in ['solver.xml','PETSC_OPTIONS.txt']:shutil.copy2(src/'run'/name,dst/'run'/name)
pol=json.loads((src/'policy.json').read_text());pol.update(case=dst.name,MPI_ranks=16,linear_algebra_backend='CPU16 GMRES ASM overlap2; native subdomain LU',maximum_solver_RSS_bytes=50*1024**3);dump(dst/'policy.json',pol);lock_case(dst)
dump(dst/'reports/preregistration.json',dict(registered_unix=time.time(),reference=str(src),change='MPI ranks8 to16 only; same LU/options/XML/mesh',matched_saved_step=5,maximum_relative_L2_difference_u_p_WSS=1e-9,reason='Verify fine-grid execution partition before fine CFD starts; reduce per-subdomain factor storage',not_a_steady_validation_result=True));print(dst)
