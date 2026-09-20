"""Classify the observed official failure without confusing it with the repaired MPI issue."""
import json,re,shutil,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from sv_validation.provenance import now,sha256,write_json
R=ROOT/'reports/sv1_3l'
for p in (R/'remote').glob('*.json'):shutil.copyfile(p,R/p.name)
def load(n):return json.loads((R/(n+'.json')).read_text())
execution=load('official_fluid_gpu_smoke_execution')
log=ROOT/'logs/sv1_3l/remote/official_fluid_gpu_smoke.log';text=log.read_text()
assert sha256(log)==execution['log_sha256']
assert load('mpi_application_gate')['status']=='PASS'
assert execution['exit_code']!=0 and 'Vector is not ghosted' in text and 'MPI_ERR_TYPE' not in text
assert 'seqaijcusparse' in execution['matrix_types'] and 'seqcuda' in execution['vector_types']
ksp=re.search(r'^KSP Object:[^\n]*\n\s*type:\s*(\S+)',text,re.M)[1]
pc=re.search(r'^PC Object:[^\n]*\n\s*type:\s*(\S+)',text,re.M)[1]
assert ksp=='gmres' and pc=='asm' and '2 levels of fill' in text and 'amount of overlap = 2' in text
reason='SVMULTIPHYSICS_GPU_SMOKE_FAIL';detail='CUDA_VECTOR_GHOST_LAYOUT_INCOMPATIBILITY'
backend={'status':'PASS','executed':True,'mat_type':'seqaijcusparse','vec_type':'seqcuda','KSP':ksp,'PC':pc,
    'ASM_overlap':2,'sub_PC':'ilu','ILU_levels':2,'factorization_package_view':'cusparse',
    'matrix_shape':[12428,12428],'matrix_nonzeros':600272,'vector_length':12428,
    'evidence':'Actual MatView, VecView, KSPView and PCView in failed official run; backend activation does not imply scientific correctness.',
    'log_sha256':execution['log_sha256'],'scientific_acceptance':False}
write_json(R/'runtime_backend.json',backend)
smoke={'status':'FAIL','reason':reason,'failure_class':detail,'executed':True,
    'startup_pass':True,'startup_evidence':'MPI_ERR_TYPE absent; actual PETSc CUDA objects and KSP reached',
    'exit_code':execution['exit_code'],'wall_time_s':execution['wall_time_s'],'monitor_stop':execution['monitor_stop'],
    'linear_solves':len(execution['history']['petsc_reasons']),
    'completed_solver_nonlinear_rows':len(execution['history']['linear_solves']),
    'linear_failures':None,'nonlinear_failures':None,'solver_status':'FAIL',
    'KSP_reasons':execution['history']['petsc_reasons'],'KSP_last_monitor':execution['history']['unassigned_petsc_monitor'][-1],
    'VTU_count':len(execution['results']),'velocity_finite':None,'pressure_finite':None,'reload_pass':None,
    'backend_status':'PASS','mat_type':'seqaijcusparse','vec_type':'seqcuda','KSP':ksp,'PC':pc,
    'PETSc_error_detected':True,'first_error':'Vector is not ghosted','first_error_call':'VecGhostUpdateBegin',
    'first_error_PETSc_source':'src/vec/vec/impls/mpi/commonmpvec.c:216',
    'first_error_solver_source':'Code/Source/solver/petsc_impl.cpp:669',
    'note':'One KSP reports convergence after an earlier ghost-vector error; no completed nonlinear row or VTU, so no solver success or zero-failure claim.'}
write_json(R/'svmp_gpu_smoke.json',smoke)
write_json(R/'official_gpu_smoke_resolution.json',{'status':'FAIL','MPI_ERR_TYPE_disappeared':True,
    'MPI_datatype_repair':'PASS','official_smoke_pass':False,'stage_J_official_smoke_blocker_formally_closed':False,
    'new_failure':detail,'source_modified':False,'PETSc_upgraded':False})
for name in ('gpu_proof','cpu_proof','science_equivalence','gpu_residency','gpu_transfer','gpu_memory','benchmark','speedup'):
    d={'status':'NOT_RUN','executed':False,'reason':reason+': '+detail,'measurements':None}
    if name in ('gpu_proof','cpu_proof'):d.update(steps=None,initial_state='t=0',linear_failures=None,nonlinear_failures=None,velocity_finite=None,pressure_finite=None,mass_error=None,reload_pass=None)
    if name=='science_equivalence':d.update(velocity_relative_L2=None,pressure_relative_L2=None,Qin_relative=None,Qout_relative=None,mass_error_difference=None,max_velocity_relative=None)
    if name=='gpu_residency':d.update(measured=False,matrix_resident=None,vectors_resident=None,MatMult_location=None,ASM_location=None,ILU2_location=None)
    if name=='gpu_transfer':d.update(measured=False,H2D_count=None,D2H_count=None,H2D_bytes=None,D2H_bytes=None,transfer_bound=None)
    if name=='gpu_memory':
        probe=next(p for p in load('pre_install_environment')['probes'] if any('--query-gpu=name,driver_version,memory.total,compute_cap'==arg for arg in p['command']))
        capacity=int(re.search(r'(\d+)\s+MiB',probe['stdout'])[1])
        d.update(measured=False,peak_bytes=None,device_capacity_MiB=capacity,capacity_source='pre_install_environment.json: nvidia-smi device capacity, not peak allocation')
    if name in ('benchmark','speedup'):d.update(CPU_1R=None,CPU_4R=None,GPU_1R=None,speedup=None)
    write_json(R/(name+'.json'),d)
write_json(ROOT/'benchmarks/sv1_3l/status.json',load('benchmark'))
write_json(R/'stage_result.json',{'timestamp':now(),'status':'FAIL','reason':reason,'failure_class':detail,
    'MPI_datatype_repair':'PASS','MPI_ERR_TYPE_disappeared':True,'PETSc_GPU_revalidation':'PASS',
    'svMultiPhysics_build':'PASS','svMultiPhysics_link':'PASS','runtime_CUDA_types':'PASS',
    'official_smoke':'FAIL','classification':'FAIL','performance_classification':'NOT_MEASURED',
    'production':'CPU_EARLY_STOP_PRODUCTION','production_changed':False,
    'GPU_full_production_started':False,'mesh_convergence_started':False,
    'classification_note':'The observed active-CUDA ghost-layout failure is not covered by the request’s named MPI-failure or backend-inactive cases; report its own signature instead of misclassifying it.'})
print('MPI repair PASS; CUDA backend active; official smoke FAIL: '+detail)
