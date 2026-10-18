"""Once-only formal H launch, after completed P0-P4 evidence."""
from common import *
import subprocess,shutil,datetime

def main():
    assert not (O/'FORMAL_REQUEST.json').exists(),'One formal run only; no outcome retry'
    assert not (O/'cases/biphasic/history.csv').exists()
    for name in ['KEY_CHECKS.json','SURFACE_GEOMETRY_AND_MAPPING_AUDIT.json','NATIVE_SURFACE_HISTORY_RESTART_AUDIT.json']:
        assert json.loads((O/name).read_text())['status']=='PASS',name
    assert json.loads((O/'DOMAIN_SCALE_PROVENANCE.json').read_text())['formal_material']=='EXACT_FROZEN_V2_10_DOMAIN_BASELINE'
    assert len(json.loads((O/'P0_TARGET_SURFACE_DIAGNOSIS.json').read_text())['targets'])==7
    shutil.copy2(A/'SURFACE_METHOD_SPEC.md',O/'SURFACE_METHOD_SPEC.md')
    spec=dict(label='CURRENT_SURFACE_APPROXIMATE_COUPLING',L1='Frozen V2 supported-simplex path retained',L2='Conservative free-component material-dual-face quadrature with bounded affine mapping',
        reference_volume_m3=V,beta_Pa_s_m=BETA,opening_m=[1e-6,5e-6],reference_and_current_condition_limit=100,minimum_affine_weight=-.25,maximum_affine_absolute_weight_sum=2,
        source_spec_sha256=sha(A/'SURFACE_METHOD_SPEC.md'),true_unrepresented_surface_area='UNKNOWN',fluid_loading_unresolved_is_not_zero_physical_drag=True)
    dump(O/'CURRENT_EXPOSED_SURFACE_SPEC.json',spec)
    name='clot_surface_coupling_v3_'+RUN;wrapper=A/'run_formal.sh'
    wrapper.write_text(f'#!/bin/bash\nset -euo pipefail\nsource /workspace/y_clot/env.sh\nexport OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1\nexport LD_LIBRARY_PATH={V2_O}/build:${{LD_LIBRARY_PATH:-}}\ncd {R}\nexec python -u {A}/run_cases.py --name biphasic\n');wrapper.chmod(0o755)
    runtime=['common.py','native.py','current_surface.py','ray_cells.cpp','run_cases.py','fragment_recording.py','shear_field_library.py','multifragment_shear_objective.py','run_formal.sh']
    dependencies=[V2_A/'current_surface.py',R/'analysis/v3_softening_pilot_v1/native.py',R/'analysis/softening_feasibility_v1/model.py',R/'analysis/v3_chunk_fragmentation_mechanism_v1/graph_ops.py',LIB,O/'build/libsurface_ray.so']
    inputs=['BIPHASIC_FIELD.npz','BIPHASIC_MATERIAL_SPEC.json','GEOMETRY_FREEZE.json','FROZEN_SHEAR_DRIVE.json','SHEAR_DRIVE.npz','MULTIFRAGMENT_OBJECTIVE.json','CURRENT_EXPOSED_SURFACE_SPEC.json','SURFACE_METHOD_SPEC.md','KEY_CHECKS.json','SURFACE_GEOMETRY_AND_MAPPING_AUDIT.json','NATIVE_SURFACE_HISTORY_RESTART_AUDIT.json','DOMAIN_SCALE_PROVENANCE.json','provenance/reference.npz','provenance/reference.native','provenance/body.data','provenance/dual_face_crossings.npz']
    old=json.loads((O/'V2_DEPENDENCY_FREEZE.json').read_text())
    for p,want in old['implementation_sha256'].items():assert sha(R/p)==want,p
    for p,want in old['input_sha256'].items():assert sha(V2_O/p)==want,p
    freeze=dict(task=TASK,run=RUN,formal_cases=['H_BIPHASIC_SURFACE_COUPLING_V3'],maximum_formal_runs=1,formal_case_argument='biphasic',V2_baseline_read_only=True,
        implementation_sha256={str(p.relative_to(R)):sha(p) for p in [*[A/q for q in runtime],*dependencies]},input_sha256={q:sha(O/q) for q in inputs},native_library_sha256=sha(LIB),
        numerical_gates=dict(base_dt_s=BASE,min_dt_s=MIN_H,max_accepted_delta_d=.1,max_relative_energy_residual=.01),load='Frozen V2 +X uniform tangential 1000Pa; 0-3.6ms; drive off at2.8ms; passive relative drag remains',
        only_formal_changed_factor='Current exposed fragment surface quadrature and resulting relative fluid forces',physical_scale_status='MICROSTRUCTURE_SCALE_UNCALIBRATED')
    dump(O/'SOURCE_FREEZE.json',freeze)
    config=f'[program:{name}]\ncommand={wrapper}\nautostart=false\nautorestart=false\nstartsecs=0\nredirect_stderr=true\nstdout_logfile={O}/formal.supervisor.log\nstdout_logfile_maxbytes=0\nstopasgroup=true\nkillasgroup=true\n'
    cfg=Path('/etc/supervisor/conf.d')/(name+'.conf');assert not cfg.exists();cfg.write_text(config);(O/(name+'.supervisor.conf')).write_text(config)
    with (O/'FORMAL_REQUEST.json').open('x') as f:
        json.dump(dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),name=name,formal_runs=1,source_freeze_sha256=sha(O/'SOURCE_FREEZE.json'),wall_time_budget=None,memory_budget=None,automatic_restart=False),f,indent=2)
    subprocess.run(['supervisorctl','reread'],check=True);subprocess.run(['supervisorctl','update',name],check=True);subprocess.run(['supervisorctl','start',name],check=True)
    print('ONE_FORMAL_H_STARTED',name,flush=True)

if __name__=='__main__':main()
