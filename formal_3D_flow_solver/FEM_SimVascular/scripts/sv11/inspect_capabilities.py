#!/usr/bin/env python3
import json,re,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from sv_validation.sv11 import REPORT,load
from sv_validation.provenance import sha256,write_json
source=ROOT/'external/svMultiPhysics'
commit=subprocess.check_output(['git','-C',str(source),'rev-parse','HEAD'],text=True).strip()
assert commit==load('flow_execution','sv1')['solver_commit']
assert not subprocess.check_output(['git','-C',str(source),'status','--porcelain'],text=True)
files=['Code/Source/solver/Parameters.cpp','Code/Source/solver/read_files.cpp','Code/Source/solver/petsc_impl.cpp','Code/Source/solver/ComMod.h','Code/Source/solver/consts.cpp','Code/Source/solver/CMakeLists.txt','Code/Source/solver/LinearAlgebra.cpp','Code/Source/solver/output.cpp','Code/Source/linear_solver/precond.cpp','Code/Source/linear_solver/ns_solver.cpp']
evidence={f:{'sha256':sha256(source/f),'url':f'https://github.com/SimVascular/svMultiPhysics/blob/{commit}/{f}'} for f in files}
impl=(source/files[2]).read_text();reader=(source/files[1]).read_text()
assert 'KSPSetFromOptions(psol[cEq].ksp)' in impl
assert '//            KSPGMRESSetRestart' in impl
assert 'lEq.FSILS.RI.relTol = linear_solver.tolerance.value()' in reader
assert 'lEq.ls.relTol = linear_solver.tolerance.value()' not in reader
cache=(ROOT/'external/build/native/svMultiPhysics-build/CMakeCache.txt').read_text()
result={'status':'SOURCE_VERIFIED','commit':commit,'source_unmodified':True,
 'FSILS_supported':True,'PETSc_supported':True,'Trilinos_supported':True,
 'baseline_build_PETSc':not bool(re.search(r'^SV_PETSC_DIR:STRING=\s*$',cache,re.M)),
 'baseline_build_Trilinos':'SV_USE_TRILINOS:BOOL=ON' in cache,
 'PETSc_enable_cmake':'SV_PETSC_DIR points to an include/lib prefix; USE_PETSC definition is then enabled',
 'xml':{'backend':{'path':'LS/Linear_algebra@type','supported':['fsils','petsc','trilinos']},
 'KSP_type':{'path':'LS@type','effective':True,'supported':['GMRES','CG','BICGS']},
 'PC_type':{'path':'LS/Linear_algebra/Preconditioner','effective':True,'supported':['petsc-jacobi','petsc-rcs']},
 'relative_tolerance':{'path':'LS/Tolerance','parsed':True,'effective_for_PETSc':False,'reason':'read_ls writes FSILS.RI.relTol; PetscImpl passes ls.relTol, whose fluid default is 0'},
 'absolute_tolerance':{'path':'LS/Absolute_tolerance','parsed':True,'effective_for_PETSc':False,'reason':'read_ls writes FSILS.RI.absTol; PetscImpl passes ls.absTol, default 1e-8'},
 'max_iterations':{'path':'LS/Max_iterations','parsed':True,'effective_for_PETSc':True,'reason':'read_ls assigns both ls.mItr and FSILS.RI.mItr'},
 'restart':{'path':'LS/Krylov_space_dimension','parsed':True,'effective_for_PETSc':False,'reason':'read_ls writes FSILS.RI.sD only; KSPGMRESSetRestart is also commented out'},
 'NS_inner_settings':{'effective_for_PETSc':False,'reason':'Only applied to FSILS NS subsolvers'}},
 'existing_runtime_interface':{'available':True,'entry':'PetscInitialize followed by KSPSetFromOptions, before KSPSetUp',
 'options_transport':'PETSC_OPTIONS environment variable; complete owned value frozen per run, no XML fiction',
 'KSP_prefix':'none for this single-equation fluid case; ns_ is added only for nEq > 1',
 'can_set':['ksp_type','pc_type','ksp_rtol','ksp_atol','ksp_max_it','ksp_gmres_restart','ksp_pc_side','ksp_norm_type','ksp_monitor_true_residual','ksp_converged_reason','ksp_view','pc_asm_overlap','sub_ksp_type','sub_pc_type','sub_pc_factor_levels'],
 'must_verify':'Compile-time linked PETSc, runtime KSP view/tolerances/reasons, and actual library option consumption'},
 'baseline_actual_settings':load('baseline_solver_diagnosis')['configuration'],
 'sources':evidence,
 'official_PETSc_references':['https://petsc.org/release/manualpages/Sys/PetscInitialize/','https://petsc.org/release/manualpages/KSP/KSPSetFromOptions/','https://petsc.org/release/manualpages/KSP/KSPConvergedDefault/'],
 'constraints':{'commit_change':False,'solver_logic_change':False,'new_backends_attempted':['petsc']}}
write_json(REPORT/'linear_solver_capabilities.json',result)
print('Pinned interface inspected: XML rtol/atol/restart are ineffective for PETSc; existing options interface required')
