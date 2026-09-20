from sv13p_support import *
def test_amgx_is_bounded_and_does_not_disrupt_stack():
 d=read('amgx_feasibility');assert d['PETSc_version']=='3.25.5' and d['CUDA_version']=='13.2'
 assert d['system_stack_modified'] is False and d['package_version']=='2.4.0'
 assert d['result'] and d['evidence']
 c=candidate('P5');assert c['result'] and not c['scientific_parameters_changed']
 if c.get('window'):assert c['smoke']['status']=='PASS'
