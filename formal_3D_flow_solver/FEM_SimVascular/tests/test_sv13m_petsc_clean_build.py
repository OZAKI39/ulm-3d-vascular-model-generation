from sv13m_support import *
from sv_validation.sv13m import *
def test_fresh_source_and_original_configuration():
 d=accepted('petsc_build');assert all(d['selftest_observed'].values())
 assert len(d['configuration_diff_from_L'])==2
 assert all(x['after'].startswith(('--prefix=','PETSC_ARCH=')) for x in d['configuration_diff_from_L'])
 for name in ('configure','make','install','selftest'):assert d[name]['exit_code']==0
 integrity=accepted('patch_integrity');assert integrity['fresh_source'] and not integrity['previous_objects_reused']
