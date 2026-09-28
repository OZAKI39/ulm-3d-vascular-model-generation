"""Preserve unfinished ILU2 run; switch execution PC only after matched field test."""
from pathlib import Path
import shutil,json,hashlib,time
V=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def put(p,d):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(d,indent=2)+'\n')
proof=json.loads((V/'data/vessel_subdomain_LU_step5.json').read_text());assert proof['backend_equivalence_work_gate'] and max(proof[k] for k in ['velocity_relative_L2','pressure_relative_L2','wss_area_relative_L2'])<1e-9
case=V/'stage3/vessel_medium';old=V/'stage3/execution_checks/vessel_medium_cpu8_ILU2_stopped';assert not old.exists() and (case/'reports/intentional_cancellation.json').exists()
# No detached solver may survive a rename of its case directory.
for p in Path('/proc').iterdir():
 if not p.name.isdigit():continue
 try:cwd=(p/'cwd').resolve();cmd=(p/'cmdline').read_bytes().decode(errors='replace')
 except OSError:continue
 assert not(cwd==case/'run' and any(n in cmd for n in ['svmultiphysics','mpiexec','mpirun'])),'Solver process still active'
case.rename(old);case.mkdir();shutil.copytree(old/'SV_MESH',case/'SV_MESH');shutil.copytree(old/'reports',case/'reports');(case/'run').mkdir()
for name in ['launch.json','progress.json','execution.json','intentional_cancellation.json','linear_attempt_acceptance.json','external_interruption.json']:
 p=case/'reports'/name
 if p.exists():p.unlink()
for name in ['solver.xml','PETSC_OPTIONS.txt']:shutil.copy2(old/'run'/name,case/'run'/name)
shutil.copy2(old/'policy.json',case/'policy.json')
for name in ['mesh_request.json','mesh_launch.json']:
 if (old/name).exists():shutil.copy2(old/name,case/name)
for name,ranks in [('vessel_medium',8),('vessel_fine',16)]:
 c=V/'stage3'/name;assert not (c/'run/solver.log').exists();a=c/'reports/before_subdomain_LU';a.mkdir()
 for rel in ['policy.json','input_hashes.json','run/PETSC_OPTIONS.txt']:
  source=(old/rel) if name=='vessel_medium' else c/rel;target=a/rel;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,target)
 s=(c/'run/PETSC_OPTIONS.txt').read_text();assert '-sub_pc_type ilu -sub_pc_factor_levels 2' in s;(c/'run/PETSC_OPTIONS.txt').write_text(s.replace('-sub_pc_type ilu -sub_pc_factor_levels 2','-sub_pc_type lu'))
 pol=json.loads((c/'policy.json').read_text());pol.update(MPI_ranks=ranks,linear_algebra_backend=f'CPU{ranks} PETSc aij/standard GMRES ASM overlap2; subdomain exact LU (native PETSc)',maximum_solver_RSS_bytes=50*1024**3,execution_selection_reason='Matched CPU8 ILU2/LU step5 fields below 5e-15 relative; subdomain LU reduces iterations. Fine uses16MPI to reduce local factor storage; its matched-step verification required before launch.');put(c/'policy.json',pol)
 files=[p for root in [c/'SV_MESH',c/'run'] for p in root.rglob('*') if p.is_file()]+[c/'policy.json'];put(c/'input_hashes.json',{str(p.relative_to(c)):sha(p) for p in files})
 put(c/'reports/subdomain_LU_selection.json',dict(unix=time.time(),matched_field_proof=proof,mesh_XML_dt_physics_and_residual_targets_unchanged=True,MPI_ranks=ranks,zero_initial_state=True,not_a_new_physical_model=True))
 print(name,'LU adopted',ranks)
