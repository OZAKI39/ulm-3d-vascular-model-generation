#!/usr/bin/env python3
"""Freeze historical evidence before any performance experiment."""
import gzip,json,shutil,sys,xml.etree.ElementTree as ET
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from sv_validation.sv13 import *
from sv_validation.provenance import inventory,git_state,write_json,now
from sv_validation.postprocess import SolutionMeasurements
from sv_validation.sv12 import checkpoint_audit

assert not (REPORT/'reference_freeze.json').exists(), 'Do not replace frozen evidence'
for p in (REPORT,CONFIG,OUTPUT,LOG,ROOT/'benchmarks/sv1_3',OUTPUT/'qc'):p.mkdir(parents=True,exist_ok=True)
scopes=[f'{kind}/{stage}' for kind in ('reports','outputs','logs','configs') for stage in ('sv1','sv1_1','sv1_2') if (ROOT/kind/stage).exists()]+['inputs']
history={'timestamp':now(),'files':{},'scopes':scopes}
for name in scopes:history['files'].update({name+'/'+k:v for k,v in inventory(ROOT/name)['files'].items()})
for name in ('configs','src','scripts','tests'):
    for k,v in inventory(ROOT/name)['files'].items():
        if 'sv13' in k or 'sv1_3' in k or '__pycache__' in k:continue
        history['files'][name+'/'+k]=v
(REPORT/'history_baseline.json.gz').write_bytes(gzip.compress(json.dumps(history).encode(),mtime=0))
old=inventory(ROOT.parent/'FEM');old['git']=git_state(ROOT.parent/'FEM')
(REPORT/'old_fem_baseline.json.gz').write_bytes(gzip.compress(json.dumps(old).encode(),mtime=0))

accepted=load('accepted_solution','sv1_2');p12=json.loads((ROOT/'configs/sv1_2/policy.json').read_text())
opts=json.loads((ROOT/'configs/sv1_2/petsc_options.json').read_text())
build=load('petsc_build_manifest','sv1_1')
manifest=load('frozen_input_manifest','sv1_2')['files']
paths={Path(f['path']) for f in manifest}
paths.update(ROOT/'inputs'/p for p in ('MANIFEST.json','fem_reference/exterior_surface.npz','fem_reference/volume_mesh.npz'))
paths.update(ROOT/'reports/sv1_2'/p for p in ('accepted_solution.json','steady_assessment.json','solver_resource_usage.json','solver_history.json','saved_state_qc.json'))
paths.update((ROOT/accepted['path'],ROOT/'configs/sv1_2/sv_flow.xml',ROOT/'configs/sv1_2/policy.json',ROOT/'configs/sv1_2/petsc_options.json'))
files=[{'path':str(p.relative_to(ROOT)),'sha256':sha256(p)} for p in sorted(paths)]
for f in manifest:assert sha256(f['path'])==f['sha256']
assert sha256(ROOT/accepted['path'])==accepted['sha256']
policy=dict(p12,mass_limit=1e-6,
    equivalence={'velocity_relative_volume_L2':1e-5,'pressure_relative_volume_L2':1e-5,
      'relative_Qin':1e-6,'relative_Qout':1e-5,'absolute_fraction_difference':1e-5,
      'port_pressure_scale_normalized':1e-5,'relative_max_velocity':1e-5},
    stop_policy='After first qualifying saved state, atomically request next saved step via native STOP_SIM; revalidate final state and checkpoint; no signals',
    gpu_minimum_speedup=1.25,benchmark_window_steps=20,benchmark_repetitions=2,
    third_repeat_relative_difference=.1,maximum_gpu_B_candidates=1,
    CPU_initial_state='same legal native step10 checkpoint as SV1.2; inherited first10 timing separately reported',
    pressure_scale='max(abs(reference pressure range)); common scale across all ports; no pressure shift')
write_json(CONFIG/'policy.json',policy);write_json(CONFIG/'petsc_options.json',opts)
xml=ET.parse(ROOT/'configs/sv1_2/sv_flow.xml')
ref={'status':'PASS','timestamp':now(),'designation':'VALIDATION_REFERENCE','files':files,
     'source_stage':'sv1_2','accepted_solution':accepted,'physics':{'rho':float(xml.findtext('.//Density')),'mu':float(xml.findtext('.//Viscosity/Value'))},
     'BC':[ET.tostring(b,encoding='unicode') for b in xml.findall('.//Add_BC')],
     'dt':float(xml.findtext('.//Time_step_size')),'time_integrator':{'name':'native generalized alpha','spectral_radius':float(xml.findtext('.//Spectral_radius_of_infinite_time_step'))},
     'PETSc':opts,'build':{'commit':build['commit'],'binary_sha256':build['executable_sha256'],'PETSc_version':build['version']},
     'policy':policy,'baseline_resources':load('solver_resource_usage','sv1_2'),
     'baseline_steady':load('steady_assessment','sv1_2'),
     'user_review':'SV1.3 request explicitly confirms SV1.2 numerical and manual acceptance; historical report remains unchanged'}
write_json(REPORT/'reference_freeze.json',ref);check_reference()
case=OUTPUT/'cpu_early_stop';(case/'4-procs').mkdir(parents=True)
seed=ROOT/'outputs/sv1_2/vascular_flow/4-procs/stFile_010.bin'
audit=checkpoint_audit(seed,10,policy['dt_s'])
shutil.copyfile(seed,case/'4-procs/stFile_last.bin')
shutil.copyfile(seed,case/'4-procs/stFile_010.bin')
shutil.copyfile(ROOT/'outputs/sv1_2/vascular_flow/4-procs/result_010.vtu',case/'4-procs/result_010.vtu')
shutil.copyfile(ROOT/'configs/sv1_2/sv_flow.xml',case/'solver.xml')
shutil.copyfile(case/'solver.xml',CONFIG/'cpu_early_stop.xml')
(OUTPUT/'SV_MESH').symlink_to(ROOT/'outputs/sv1/SV_MESH',target_is_directory=True)
write_json(REPORT/'cpu_initial_state.json',audit)
write_json(REPORT/'termination_source_audit.json',{'status':'PASS','source':'external/svMultiPhysics/Code/Source/solver/main.cpp',
    'sha256':sha256(ROOT/'external/svMultiPhysics/Code/Source/solver/main.cpp'),
    'lines':{'read_STOP_SIM':431,'write_restart':471,'write_VTU':483,'break_after_output':531},
    'no_scientific_source_change':True,'final_log_flush':'normal native process exit; Popen.wait before evidence acceptance',
    'scheduling':policy['stop_policy']})
print('SV1.3 reference frozen; legal step10 copy; policy frozen before replay or run',flush=True)
