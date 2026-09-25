from pathlib import Path
import shutil,json
R=Path('/workspace/microbubble_lammps/results/stable_rigid_sphere_near_field_v0_20260916_101617');D=R/'failed_attempts/initial_finalizer_serialization';D.mkdir(parents=True,exist_ok=True)
for s,n in [('scripts/validate_cases.py','validate_cases_before_fix.py'),('scripts/run_phase.py','run_phase_before_fix.py'),('logs/PHASE_INITIAL.log','PHASE_INITIAL.log'),('PHASE_INITIAL_STATE.json','PHASE_INITIAL_STATE.json'),('cases/A_normal_cg04/INDEPENDENT_VALIDATION.log','A_VALIDATION.log')]:shutil.copy2(R/s,D/n)
(D/'RECOVERY.json').write_text(json.dumps({'cause':'numpy.bool_ is not JSON serializable in validation report writer','physics_changed':False,'solver_retry':False,'recovery':'convert report scalars to Python bool/int; validate completed A and proceed only after PASS'},indent=2)+'\n')
print((R/'cases/A_normal_cg04/RUN_STATE.json').read_text())
