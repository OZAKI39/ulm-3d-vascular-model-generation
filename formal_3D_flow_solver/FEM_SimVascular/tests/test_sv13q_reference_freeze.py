from sv13q_support import *
def test_stage_p_delivery_and_frozen_science():
 ref=read('reference_freeze');assert ref['Stage_P_files_verified']==len(json.loads((ROOT/'reports/sv1_3p/delivery_manifest.json').read_text())['files'])
 assert digest(ROOT/'reports/sv1_3p/delivery_manifest.json')==ref['Stage_P_delivery_sha256']
 assert digest(ROOT/policy()['production_policy_path'])==policy()['production_policy_sha256']
 assert ref['PETSc_build']['version']=='3.25.5'
 assert ref['PETSc_build']['CUDA_version']=='13.2'
 for d in cases():
  scientific_xml_gate(ROOT/'outputs/sv1_3p/REAL_VASCULAR_GPU_PC_WINNER/solver.xml',ROOT/'outputs/sv1_3q'/d['name']/'solver.xml')
  assert d['PETSc_library_sha256']==ref['GPU_baseline']['PETSc_library_sha256']
  for s in d['runtime_semantics']:semantics_gate(s,d['PETSC_OPTIONS'])
def test_no_new_packages_or_extra_solver_strategies():
 for d in cases():
  o=d['PETSC_OPTIONS'];assert '-pc_type asm' in o and '-sub_pc_factor_levels 2' in o and '-ksp_gmres_restart 100' in o
  assert not any(x in o.lower() for x in ('hypre','gamg','amgx','fieldsplit'))
