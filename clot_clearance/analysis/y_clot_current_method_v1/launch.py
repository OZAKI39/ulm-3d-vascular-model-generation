from common import *
import subprocess,datetime

def main():
    assert json.loads((O/'KEY_CHECKS.json').read_text())['status']=='PASS'
    assert not (O/'FORMAL_REQUEST.json').exists()
    runtime=['common.py','native.py','current_surface.py','core_surface.py','supported_surface.py','ray_cells.cpp',
        'run_cases.py','fragment_recording.py','shear_field_library.py','multifragment_shear_objective.py']
    dependencies=[R/'analysis/v3_softening_pilot_v1/native.py',R/'analysis/softening_feasibility_v1/model.py',
        R/'analysis/v3_chunk_fragmentation_mechanism_v1/graph_ops.py',LIB,O/'build/libsurface_ray.so',
        R/'pd_clot/lammps_fragmentation/native_y_wall_query.py',R/'pd_clot/lammps_fragmentation/native_y_wall_contact.py']
    inputs=['BIPHASIC_FIELD.npz','BIPHASIC_MATERIAL_SPEC.json','GEOMETRY_FREEZE.json','FROZEN_SHEAR_DRIVE.json',
        'CURRENT_EXPOSED_SURFACE_SPEC.json','KEY_CHECKS.json','provenance/reference.npz','provenance/reference.native',
        'provenance/body.data','provenance/dual_face_crossings.npz']
    dump(O/'SOURCE_FREEZE.json',dict(task=TASK,run=RUN,maximum_formal_runs=1,
        implementation_sha256={str(p.relative_to(R)):sha(p) for p in [*[A/q for q in runtime],*dependencies]},
        input_sha256={p:sha(O/p) for p in inputs},native_library_sha256=sha(LIB),
        numerical_gates=dict(base_dt_s=BASE,min_dt_s=MIN_H,max_accepted_delta_d=.1,max_relative_energy_residual=.01)))
    wrapper=A/'run_formal.sh'
    wrapper.write_text(f'#!/bin/bash\nset -euo pipefail\nsource /workspace/y_clot/env.sh\nexport LD_LIBRARY_PATH={LIB.parent}:${{LD_LIBRARY_PATH:-}}\ncd {R}\nexec python -u {A}/run_cases.py --name biphasic\n');wrapper.chmod(0o755)
    name='y_clot_current_method_'+RUN
    config=f'[program:{name}]\ncommand={wrapper}\nautostart=false\nautorestart=false\nstartsecs=0\nredirect_stderr=true\nstdout_logfile={O}/formal.supervisor.log\nstdout_logfile_maxbytes=0\nstopasgroup=true\nkillasgroup=true\n'
    path=Path('/etc/supervisor/conf.d')/(name+'.conf');assert not path.exists();path.write_text(config)
    (O/'supervisor.conf').write_text(config)
    with (O/'FORMAL_REQUEST.json').open('x') as f:
        json.dump(dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),formal_runs=1,supervisor_name=name,
            wall_clock_budget=None,memory_budget=None,automatic_restart=False,source_freeze_sha256=sha(O/'SOURCE_FREEZE.json')),f,indent=2)
    for cmd in [['supervisorctl','reread'],['supervisorctl','update',name],['supervisorctl','start',name]]:subprocess.run(cmd,check=True)
    print('ONE_Y_FORMAL_STARTED',name,flush=True)

if __name__=='__main__':main()
