#!/usr/bin/env python3
import gzip,json,shutil,subprocess,sys,xml.etree.ElementTree as ET
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from sv_validation.sv12 import REPORT,OUTPUT,CONFIG,LOG,load,checkpoint_audit,validate_xml,validate_frozen
from sv_validation.sv11 import parse_solver_log,linear_gate,nonlinear_gate
from sv_validation.provenance import inventory,git_state,sha256,write_json,now
from sv_validation.postprocess import SolutionMeasurements
for p in (REPORT,OUTPUT,CONFIG,LOG,OUTPUT/'qc'):p.mkdir(parents=True,exist_ok=True)
assert not (REPORT/'history_baseline.json.gz').exists(),'Do not replace a frozen baseline'
scopes=[f'{kind}/{stage}' for kind in ('reports','outputs','logs') for stage in ('sv1','sv1_1')]+['inputs']
history={'root':str(ROOT),'timestamp':now(),'files':{},'scopes':scopes}
for name in scopes:
    history['files'].update({name+'/'+k:v for k,v in inventory(ROOT/name)['files'].items()})
for directory,glob in [('configs','*'),('configs/sv1_1','*'),('tests','*.py'),('src/sv_validation','*.py'),('scripts','*.py'),('scripts/sv11','*.py')]:
    for p in (ROOT/directory).glob(glob):
        if p.is_file() and p.name!='sv12.py':
            s=p.stat();history['files'][str(p.relative_to(ROOT))]={'sha256':sha256(p),'size':s.st_size,'mtime_ns':s.st_mtime_ns,'mode':s.st_mode}
(REPORT/'history_baseline.json.gz').write_bytes(gzip.compress(json.dumps(history).encode(),mtime=0))
old=inventory(ROOT.parent/'FEM');old['git']=git_state(ROOT.parent/'FEM')
(REPORT/'old_fem_baseline.json.gz').write_bytes(gzip.compress(json.dumps(old).encode(),mtime=0))

build=load('petsc_build_manifest','sv1_1');prior=load('petsc_short_execution','sv1_1')
policy=json.loads((ROOT/'configs/time_policy.json').read_text());dt=policy['dt_s']
prior['history']=parse_solver_log((ROOT/prior['log']).read_text(),dt)
linear_gate(prior);nonlinear_gate(prior['history'])
q=load('petsc_short_qc','sv1_1')['states'][-1]
measurement=SolutionMeasurements(ROOT/'outputs/sv1/SV_MESH/mesh_arrays.npz',q['Q_target_m3_s'],policy['Umean_m_s'])
u,p=measurement.read(ROOT/q['path']);fresh=measurement.measure(u,p)
assert fresh['velocity_finite'] and fresh['pressure_finite'] and fresh['wall_noslip_pass']
assert build['version']=='3.19.6' and sha256(build['executable'])==build['executable_sha256']
assert subprocess.check_output(['git','-C',str(ROOT/'external/svMultiPhysics'),'rev-parse','HEAD'],text=True).strip()==build['commit']
assert not subprocess.check_output(['git','-C',str(ROOT/'external/svMultiPhysics'),'status','--porcelain'],text=True).strip()
write_json(REPORT/'prerequisites.json',{'status':'PASS','source_stage_status_preserved':load('stage_result','sv1_1')['status'],
           'linear_solves':len(prior['history']['linear_solves']),'linear_failures':0,'ill_conditioned_warnings':0,
           'nonlinear_pass':True,'finite_velocity':True,'finite_pressure':True,'wall_noslip':True,
           'repeat_short_benchmark':False,'old_step10_mass_ranking_gate_applied':False,'actual_vtu_sha256':sha256(ROOT/q['path'])})
files=[]
def record(path,role):
    path=Path(path);files.append({'path':str(path),'sha256':sha256(path),'semantic_role':role})
for path in (ROOT/'outputs/sv1/SV_MESH').rglob('*'):
    if path.is_file():record(path,'frozen volume/surface mesh, boundary identity or integration topology')
for name,role in [('configs/sv1_1/sv_flow_petsc.xml','SV1.1 production model, BC, time and LS'),('configs/sv_reference.yaml','reference rho mu Q'),('configs/face_map.json','face IDs and semantic roles'),('configs/time_policy.json','frozen time policy'),('outputs/sv1/vascular_flow/prescribed_inlet.vtp','frozen prescribed inlet diagnostic'),('reports/sv1_1/petsc_settings.json','verified runtime PETSc configuration'),('reports/sv1_1/petsc_short_execution.json','actual PETSC_OPTIONS and MPI/OMP'),('reports/sv1_1/petsc_build_manifest.json','pinned executable and libraries')]:record(ROOT/name,role)
record(build['executable'],'identical PETSc-enabled solver binary')
record(Path(build['petsc_path'])/'lib/libpetsc.so','identical PETSc 3.19.6 shared library')
record(Path(build['petsc_path'])/'include/petscversion.h','PETSc version header')
options={'PETSC_OPTIONS':prior['PETSC_OPTIONS'],'options':load('petsc_settings','sv1_1')['options'],'mpi_ranks':prior['mpi_ranks'],'OMP_NUM_THREADS':prior['omp_num_threads']}
write_json(CONFIG/'petsc_options.json',options)
write_json(CONFIG/'policy.json',dict(policy,early_transient_mass_solver_ranking=False,maximum_total_steps=800,
           steady_flow_normalization='Qtarget',extension_eligibility='No failures; all finite; last-five error maxima <= max(threshold, preceding-five error maxima)',
           native_restart_allowed=True,VTU_restart_forbidden=True))
write_json(REPORT/'frozen_input_manifest.json',{'status':'PASS','files':files,'timestamp':now(),'source_stage':'sv1_1',
           'unchanged':{'geometry':True,'mesh':True,'physics':True,'BC':True,'dt':True,'PETSc_settings':True,'MPI_OMP':True},
           'source_binary_sha256':build['executable_sha256'],'original_initial_condition_retained_in_native_history':True})
validate_frozen()

checkpoint=ROOT/'outputs/sv1_1/vascular_short/4-procs/stFile_010.bin'
restart=checkpoint_audit(checkpoint,10,dt)
assert sha256(checkpoint)==sha256(checkpoint.with_name('stFile_last.bin'))
sources=['Code/Source/solver/output.cpp','Code/Source/solver/initialize.cpp','Code/Source/solver/Integrator.cpp','Code/Source/solver/SPLIT.c']
restart.update(started_from='native restart',first_new_step=11,source_checkpoint_sha256=sha256(checkpoint),
               source_code={s:sha256(ROOT/'external/svMultiPhysics'/s) for s in sources},
               partition_evidence='Same binary, mesh, four ranks; pinned SPLIT.c sets ParMETIS seed=10. Native headers include local node counts.',
               history_evidence='write_restart saves Y_n and A_n; init_from_bin reads into old state; Integrator constructor copies old into current. First-order generalized-alpha needs these plus step/time. Rigid fluid dFlag=0; no displacement or coupled boundary state required.')
case=OUTPUT/'vascular_flow';(case/'4-procs').mkdir(parents=True,exist_ok=True)
shutil.copyfile(checkpoint,case/'4-procs/stFile_last.bin')
shutil.copyfile(checkpoint,case/'4-procs/stFile_010.bin')
shutil.copyfile(ROOT/q['path'],case/'4-procs/result_010.vtu')
assert sha256(case/'4-procs/stFile_last.bin')==restart['sha256']
restart['copied_checkpoint_path']=str(case/'4-procs/stFile_last.bin');restart['seed_state_origin']='immutable SV1.1 step10; copied only as initial QC state, not a new solve'
write_json(REPORT/'restart_decision.json',restart)
(OUTPUT/'SV_MESH').symlink_to(ROOT/'outputs/sv1/SV_MESH',target_is_directory=True)
tree=ET.parse(ROOT/'configs/sv1_1/sv_flow_petsc.xml');tree.find('.//Continue_previous_simulation').text='true'
ET.indent(tree);tree.write(CONFIG/'sv_flow.xml',encoding='utf-8',xml_declaration=True)
validate_xml(CONFIG/'sv_flow.xml');shutil.copyfile(CONFIG/'sv_flow.xml',case/'solver.xml')
assert not (case/'STOP_SIM').exists()
write_json(REPORT/'xml_invariance.json',{'status':'PASS','only_difference':'Continue_previous_simulation false -> true, explicitly authorized native checkpoint continuation',
           'baseline_xml_sha256':sha256(ROOT/'configs/sv1_1/sv_flow_petsc.xml'),'new_xml_sha256':sha256(CONFIG/'sv_flow.xml'),
           'target_step':400,'STOP_SIM_present':False})
print('SV1.2 prepared: verified native step10 checkpoint, frozen PETSc, target400, no short-run STOP_SIM')
