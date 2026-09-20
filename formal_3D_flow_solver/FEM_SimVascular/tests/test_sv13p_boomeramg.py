from sv13p_support import *
def test_one_gpu_supported_boomer_profile():
 d=candidate('P3');assert d['result'] and not d['scientific_parameters_changed']
 if d.get('config') is None:return
 s=d['config']['PETSC_OPTIONS'];assert '-pc_hypre_boomeramg_coarsen_type PMIS' in s
 assert '-pc_hypre_boomeramg_interp_type ext+i' in s and '-pc_hypre_boomeramg_relax_type_all l1scaled-Jacobi' in s
 assert '-pc_hypre_boomeramg_relax_type_coarse l1scaled-Jacobi' in s
 if d.get('window'):assert d['smoke']['status']=='PASS'
