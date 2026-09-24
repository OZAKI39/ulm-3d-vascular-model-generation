#!/usr/bin/env python3
"""Freeze the environment stop before any vascular import, remesh or solver run."""
import json
import shutil
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from sv_validation.provenance import now,sha256,write_json
O=ROOT/'outputs/sv0';R=ROOT/'reports/sv0'
read=lambda p:json.loads(p.read_text())
for name in ('remote_probe.json','remote_solver_pull.json'):
    source=O/'remote_return/outputs/sv0/environment'/name;target=O/'environment'/name
    if target.exists():assert sha256(target)==sha256(source)
    shutil.copy2(source,target)
api=read(O/'environment/simvascular_headless_attempt_01.json');line=next(s for s in api['stdout'].splitlines() if s.startswith('SV0_API_PROBE_JSON='))
api_result=json.loads(line.split('=',1)[1]);assert api['returncode']==0 and api_result['meshing_api_present'] and not api_result['options_use_mmg']
write_json(R/'simvascular_api_probe.json',{'status':'HEADLESS_API_PASS','scope':'Embedded Python imports and actual SimVascular mesher/model constructors; no real geometry import or mesh generation yet',
    'result':api_result,'raw_attempt':'outputs/sv0/environment/simvascular_headless_attempt_01.json',
    'raw_attempt_sha256':sha256(O/'environment/simvascular_headless_attempt_01.json'),
    'warnings':'Some unused MITK DICOM modules could not auto-load libicuuc.so.66; headless Python meshing API returned successfully.',
    'bundled_library_note':'The official distribution includes its internal mesh libraries. No external mesher package/environment was installed; use_mmg is explicitly false; no mesh generator has been executed.'})
probe=read(R/'environment_probe.json');probe['post_install_wsl_api']=read(R/'simvascular_api_probe.json');probe['initial_probe_only']=False
write_json(R/'environment_probe.json',probe)
registry=read(R/'solver_registry_resolution.json')
pulls={site:read(O/f'environment/{site}_solver_pull.json') for site in ('wsl','remote')}
assert not any(row['image_pulled'] or row['container_started'] for row in pulls.values())
container={'status':'BLOCKED','reason':'CONTAINER_RUNTIME_UNAVAILABLE','requested_image':'simvascular/solver:latest',
           'registry_metadata':registry,'pull_attempts':pulls,'local_image_id':None,'local_image_digest':None,
           'container_inspect':None,'container_created':False,'case_mounts':[],
           'pin_for_future_actual_pull':registry['resolved_immutable_reference'],
           'image_layers_pulled':False,'digest_is_registry_resolution_only':True,
           'remote_restrictions':'Unprivileged Docker instance: no nested Docker engine; no installed Docker/Podman/Apptainer/Singularity; unshare -Ur fails EPERM; no Docker socket or /dev/fuse',
           'standalone_solver_substitution':False,'old_fem_environment_modified':False}
write_json(R/'solver_container_manifest.json',container)
write_json(R/'official_smoke.json',{'status':'NOT_EXECUTED_ENVIRONMENT_BLOCKED','case_requested':'fluid/pipe_RCR_3d',
    'repository':'https://github.com/SimVascular/svMultiPhysics','compatible_revision':None,'solver_exit_code':None,
    'result_vtu':None,'velocity_finite':None,'pressure_finite':None,'official_expectation_checked':False,
    'reason':'Required official container could not be pulled or started. No real vascular debugging is allowed before smoke PASS.'})
status={'stage':'SV0','status':'BLOCKED','reason':'ENVIRONMENT','detail':'CONTAINER_RUNTIME_UNAVAILABLE',
        'classification':'BLOCKED','feasibility_classification':'SOLVER_ENVIRONMENT_FAIL',
        'stopped_at_step':5,'last_completed_step':4,'simvascular_headless_api':'PASS',
        'meshing_generation_capability':'NOT_YET_TESTED','official_fluid_smoke':'NOT_EXECUTED_ENVIRONMENT_BLOCKED',
        'vascular_geometry_import':'NOT_EXECUTED','meshing_policy_frozen':False,'candidate_M0':'NOT_GENERATED','candidate_M1':'NOT_GENERATED',
        'new_pressure_support_audit':'NOT_EXECUTED','real_solver':'NOT_EXECUTED','steady':'NOT_EVALUATED',
        'flow_gates':'NOT_EVALUATED','solution_reload':'NOT_EXECUTED','profile_sensitivity':'NOT_AVAILABLE',
        'experimental':False,'warning':'NOT_EXPERIMENTAL_PUMP_FLOW','gpu_used':False,
        'next_stage_started':False,'timestamp':now()}
write_json(R/'execution_status.json',status)
names=['source_geometry_and_faces','surface_before_after','surface_geometry_error','simvascular_mesh_cutaway','mesh_quality_comparison',
       'pressure_support_comparison','real_geometry_and_bc','steady_convergence','velocity_global','velocity_slices','pressure_global',
       'pressure_sections','flux_balance','outlet_flow_split','profile_sensitivity','solver_resource_usage']
write_json(R/'visualization_manifest.json',{'status':'NOT_GENERATED_ENVIRONMENT_BLOCKED','generated_figures':0,
    'policy':'Stop at Step 5; no artificial flow or mesh plots and no state pages substituted for scientific data.',
    'figures':{name+'.png':{'generated':False,'reason':'ENVIRONMENT_BLOCKED_BEFORE_SMOKE','optional':name=='profile_sensitivity'} for name in names}})
print(json.dumps(status,indent=2))
