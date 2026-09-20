from sv13p_support import *
def test_single_documented_gpu_ilu_profile():
 d=candidate('P2');assert d['result'] and d['scientific_parameters_changed'] is False
 if d.get('config') is None:return
 s=d['config']['PETSC_OPTIONS'];help=(ROOT/'logs/sv1_3p/remote/hypre_ilu_standalone_help.log').read_text()
 for x in re.findall(r'(?:^|\s)(-pc_\w+)',s):assert x in help
 assert '-pc_hypre_ilu_type Block-Jacobi-ILUk' in s and '-pc_hypre_ilu_level 0' in s
 assert '-pc_hypre_ilu_local_reordering false' in s
 assert '-pc_hypre_ilu_tri_solve false' in s
 if d.get('window'):assert d['smoke']['status']=='PASS'
