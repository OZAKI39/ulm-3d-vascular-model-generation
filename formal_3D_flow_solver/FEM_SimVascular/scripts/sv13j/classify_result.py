"""Publish the winner and observed svMultiPhysics failure without hiding raw evidence."""
import json,shutil,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from sv_validation.provenance import sha256,write_json,now
R=ROOT/'reports/sv1_3j'
for p in (R/'remote').glob('*.json'):shutil.copyfile(p,R/p.name)
def load(n):return json.loads((R/(n+'.json')).read_text())
matrix=load('compatibility_matrix');assert matrix['winner']=='cuda123' and len(matrix['candidates'])==1
for key,version in [('cuda122','12.2.2'),('cuda121','12.1.1')]:
    row={'key':key,'CUDA_version':version,'status':'NOT_REQUIRED','executed':False,'reason':'STOP_ON_FIRST_SUCCESS: CUDA 12.3.2 make PASS',
         'nvcc_version':None,'Thrust_version':None,'host_compiler':None,'kernel':'NOT_RUN','MPI':'NOT_RUN','configure':'NOT_RUN','make':'NOT_RUN','first_error':None,'classification':'NOT_REQUIRED'}
    matrix['candidates'].append(row)
    write_json(R/(key+'_runtime.json'),row)
write_json(R/'compatibility_matrix.json',matrix)
diagnosis=load('mpi_datatype_diagnosis');rows=[row for run in diagnosis['runs'] for row in run['rows']]
fnames={'MPI_INTEGER','MPI_DOUBLE_PRECISION','MPI_CHARACTER','MPI_LOGICAL'}
assert len(rows)==24
assert all(row['error_class']==3 and row['size']==0 for row in rows if row['name'] in fnames)
assert all(row['bcast_rc']==0 and row['size']>0 for row in rows if row['name'] not in fnames)
diagnosis.update(raw_classification=diagnosis['status'],status='CONFIRMED',failed_Fortran_datatypes=sorted(fnames),
    classification_correction='Initial analysis expected MPI_Type_size to return an error. It actually returns success with size 0; MPI_Bcast returns MPI_ERR_TYPE (3). Classification corrected from the same raw runs, without rerunning or changing MPI.')
write_json(R/'mpi_datatype_diagnosis.json',diagnosis)
execution=load('official_fluid_gpu_smoke_execution');log=ROOT/'logs/sv1_3j/remote/official_fluid_gpu_smoke.log'
assert sha256(log)==execution['log_sha256'] and 'MPI_ERR_TYPE: invalid datatype' in log.read_text()
assert execution['exit_code']!=0 and not execution['results'] and not execution['history']['linear_solves']
write_json(R/'svmp_gpu_smoke.json',{'status':'FAIL','reason':'SVMULTIPHYSICS_GPU_FAIL',
    'failure_class':'MPI_FORTRAN_PREDEFINED_DATATYPES_UNAVAILABLE','executed':True,'execution':execution,
    'linear_solves_reached':False,'nonlinear_solves_reached':False,'runtime_cuda_types':'NOT_REACHED',
    'mat_type':None,'vec_type':None,'KSP':None,'PC':None,'VTU_count':0,'solver_source_modified':False,
    'note':'Failure occurs before runtime CUDA backend can be inspected; not evidence that PETSc options were ignored.'})
for name in ('gpu_proof','cpu_proof','science_equivalence','gpu_residency','gpu_transfer','gpu_memory','benchmark','speedup'):
    d={'status':'NOT_RUN','executed':False,'reason':'SVMULTIPHYSICS_GPU_FAIL: official fluid smoke failed at MPI_Bcast','measurements':None}
    if name in ('gpu_proof','cpu_proof'):d.update(steps=None,initial_state='t=0',linear_failures=None,nonlinear_failures=None,velocity_finite=None,pressure_finite=None,mass_error=None,reload_pass=None)
    if name=='science_equivalence':d.update(velocity_relative_L2=None,pressure_relative_L2=None,Qin_relative=None,Qout_relative=None,mass_error_difference=None,max_velocity_relative=None)
    if name=='gpu_residency':d.update(measured=False,matrix_resident=None,vectors_resident=None,MatMult_location=None,ASM_location=None,ILU2_location=None)
    if name=='gpu_transfer':d.update(measured=False,H2D_count=None,D2H_count=None,H2D_bytes=None,D2H_bytes=None,transfer_bound=None)
    if name=='gpu_memory':d.update(measured=False,peak_bytes=None,device_capacity_MiB=24564)
    if name in ('benchmark','speedup'):d.update(CPU_1R=None,CPU_4R=None,GPU_1R=None,speedup=None)
    write_json(R/(name+'.json'),d)
write_json(ROOT/'benchmarks/sv1_3j/status.json',load('benchmark'))
write_json(R/'stage_result.json',{'timestamp':now(),'status':'FAIL','reason':'SVMULTIPHYSICS_GPU_FAIL',
    'failure_class':'MPI_FORTRAN_PREDEFINED_DATATYPES_UNAVAILABLE','compatibility_winner':'CUDA_12_3_2',
    'PETSc_GPU_runtime':'PASS','svMultiPhysics_build':'PASS','svMultiPhysics_link':'PASS','official_fluid_smoke':'FAIL',
    'classification':'BLOCKED','production':'CPU_EARLY_STOP_PRODUCTION','production_changed':False,
    'GPU_full_production_started':False,'mesh_convergence_started':False,'matrix_exhausted':False})
print('CUDA12.3.2 winner; PETSc GPU PASS; stage FAIL: SVMULTIPHYSICS_GPU_FAIL (MPI datatype gap)')
