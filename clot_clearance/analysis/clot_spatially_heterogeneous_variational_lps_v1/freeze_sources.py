from common import *
import subprocess,difflib
assert not (O/'SOURCE_FREEZE.json').exists()
for p in ['KEY_CHECKS.json','VARIATIONAL_ENERGY_AUDIT.json','PRESCRIBED_PATH_ENERGY.json']:assert json.loads((O/p).read_text())['status']=='PASS'
names=['common.py','native.py','run_cases.py','current_surface.py','shear_field_library.py','multifragment_shear_objective.py','bond_bpm_peri_variational.cpp','bond_bpm_peri_variational.h','fix_lps_variational_dynamic.cpp','fix_lps_variational_dynamic.h']
deps=[R/'analysis/v3_softening_pilot_v1/native.py',R/'analysis/softening_feasibility_v1/model.py',R/'analysis/v3_chunk_fragmentation_mechanism_v1/graph_ops.py',LIB]
inputs=['HETEROGENEITY_FIELD.npz','HETEROGENEITY_SPEC.json','SHEAR_DRIVE.npz','SHEAR_DRIVE_FREEZE.json','CURRENT_EXPOSED_SURFACE_SPEC.json','MATERIAL_AND_SPECIMEN_FREEZE.json','KEY_CHECKS.json','VARIATIONAL_ENERGY_AUDIT.json','PRESCRIBED_PATH_ENERGY.json','provenance/reference.npz','provenance/reference.native','provenance/body.data']
patch=''
for name in names:
 base=(HA if name.startswith('bond_') else R/'analysis/clot_multifragment_breakup_v1' if name.endswith(('.cpp','.h')) else MV_A)/name
 if base.exists():patch+=''.join(difflib.unified_diff(base.read_text().splitlines(True),(A/name).read_text().splitlines(True),fromfile=str(base.relative_to(R)),tofile=str((A/name).relative_to(R))))
(O/'SOURCE_CHANGES.diff').write_text(patch)
dump(O/'SOURCE_FREEZE.json',dict(git_sha=subprocess.check_output(['git','-c',f'safe.directory={R}','rev-parse','HEAD'],cwd=R,text=True).strip(),task=TASK,run=RUN,formal_cases=CASES,maximum_formal_runs=2,implementation_sha256={str(p.relative_to(R)):sha(p) for p in [*[A/name for name in names],*deps]},input_sha256={p:sha(O/p) for p in inputs},library_sha256=sha(LIB),historical_library_sha256=sha(OLDLIB),native_energy='Unmodified endpoint-star minimized series compliance U(x,d); only original-ID bulk completion threshold is new',material='Frozen reference and typed card; G,K,rho,s0,interface,original m and graph unchanged',surface='Byte-identical current_surface.py from previous round',load='Fixed old accepted checkerboard amplitude curve, no new work controller',output='Every accepted scalar and damage event; FP64 full states every 25 us and named events; no interpolation of physical trajectories'))
assert sha(A/'current_surface.py')==sha(MV_A/'current_surface.py')
print('SOURCE_FROZEN',sha(O/'SOURCE_FREEZE.json'),flush=True)
