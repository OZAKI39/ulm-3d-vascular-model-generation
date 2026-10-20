from common import *
from current_surface import CurrentSurface
from shear_field_library import ShearField
import subprocess,difflib
assert not (O/'SOURCE_FREEZE.json').exists()
for p in ['KEY_CHECKS.json','VARIATIONAL_ENERGY_AUDIT.json','PRESCRIBED_PATH_ENERGY.json']:assert json.loads((O/p).read_text())['status']=='PASS'
assert sha(A/'current_surface.py')==sha(PREVIOUS_A/'current_surface.py')
ref=load_npz(O/'provenance/reference.npz');surface=CurrentSurface(ref,'uniform');spec=json.loads((PREVIOUS/'CURRENT_EXPOSED_SURFACE_SPEC.json').read_text());spec.update(tetra_count=len(surface.tet),unique_triangle_count=len(surface.tri));dump(O/'CURRENT_EXPOSED_SURFACE_SPEC.json',spec)
field=ShearField('uniform');t=np.linspace(0,END,3601);tau=1000*field.envelope(t);np.savez_compressed(O/'SHEAR_DRIVE.npz',time_s=t,nominal_shear_Pa=tau)
assert np.allclose(field.envelope(np.array([0,RISE/2,RISE,HOLD,(HOLD+OFF)/2,OFF,END])),[0,.5,1,1,.5,0,0])
dump(O/'FROZEN_SHEAR_DRIVE.json',dict(mode='UNIFORM_SURFACE_SHEAR_X',label=LABEL,direction=[1,0,0],peak_Pa=1000.,smoothstep='S(u)=3u^2-2u^3',rise_s=[0,RISE],hold_s=[RISE,HOLD],down_s=[HOLD,OFF],off_observation_s=[OFF,END],spatial_window=False,projection='eX - n*(n dot eX), not renormalized',traction='beta*(1000/beta*g(t)*P_t(eX) - P_t(v_surface))',passive_drag_after_off=True,beta_Pa_s_m=BETA,current_surface_sha256=sha(A/'current_surface.py'),same_frozen_drive_both_cases=True,work_feedback=False,applicability=field.applicability()))
from multifragment_shear_objective import freeze
freeze()
names=['common.py','native.py','run_cases.py','current_surface.py','shear_field_library.py','multifragment_shear_objective.py','bond_bpm_peri_variational.cpp','bond_bpm_peri_variational.h','fix_lps_variational_dynamic.cpp','fix_lps_variational_dynamic.h']
deps=[R/'analysis/v3_softening_pilot_v1/native.py',R/'analysis/softening_feasibility_v1/model.py',R/'analysis/v3_chunk_fragmentation_mechanism_v1/graph_ops.py',LIB]
inputs=['BIPHASIC_FIELD.npz','BIPHASIC_MATERIAL_SPEC.json','GEOMETRY_FREEZE.json','PREFLIGHT_STRUCTURE_AUDIT.json','SHEAR_DRIVE.npz','FROZEN_SHEAR_DRIVE.json','CURRENT_EXPOSED_SURFACE_SPEC.json','MULTIFRAGMENT_OBJECTIVE.json','KEY_CHECKS.json','VARIATIONAL_ENERGY_AUDIT.json','PRESCRIBED_PATH_ENERGY.json','provenance/reference.npz','provenance/reference.native','provenance/body.data']
patch=''
for name in names:
 base=PREVIOUS_A/name
 if base.exists():patch+=''.join(difflib.unified_diff(base.read_text().splitlines(True),(A/name).read_text().splitlines(True),fromfile=str(base.relative_to(R)),tofile=str((A/name).relative_to(R))))
(O/'SOURCE_CHANGES.diff').write_text(patch)
dump(O/'SOURCE_FREEZE.json',dict(git_sha=subprocess.check_output(['git','-c',f'safe.directory={R}','rev-parse','HEAD'],cwd=R,text=True).strip(),task=TASK,run=RUN,formal_cases=CASES,maximum_formal_runs=2,execution_order='uniform then biphasic, each once',implementation_sha256={str(p.relative_to(R)):sha(p) for p in [*[A/name for name in names],*deps]},input_sha256={p:sha(O/p) for p in inputs},library_sha256=sha(LIB),historical_library_sha256=sha(OLDLIB),native_energy='Unmodified endpoint-star minimized series compliance U(x,d); atomic original-ID threshold setter and reference/restart capacity guards are the only C++ changes',material='Frozen binary phase field; exact user threshold ratios, no energy renormalization; unchanged G,K,rho,interface,original m and graph',surface='Byte-identical previous current_surface.py',load='1000 Pa uniform +X tangential drive, specified cubic smoothstep; no spatial window or work controller',output='Every accepted scalar/graph and damage event; FP64 full states every 25 us and named events; no trajectory interpolation'))
print('SOURCE_FROZEN',sha(O/'SOURCE_FREEZE.json'),flush=True)
