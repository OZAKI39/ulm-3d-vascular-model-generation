from sv13p_support import *
def test_existing_mixed_matrix_layout_is_unchanged():
 a=read('gamg_matrix_audit');assert a['status']=='AUDITED' and a['block_size']==4 and a['matrix_type']=='seqaijcusparse'
 assert not a['near_nullspace_presence'] and not a['DOF_reordered'] and not a['physical_near_nullspace_added']
 d=candidate('P4');assert d['config']['PETSC_OPTIONS'].endswith('-pc_type gamg')
 if d['smoke']['status']=='FAIL':assert d['result']=='GAMG_NOT_SUITABLE_CURRENT_MONOLITHIC_LAYOUT'
