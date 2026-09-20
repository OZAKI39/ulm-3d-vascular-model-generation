from sv13n_support import *
def test_actual_source_references_have_locations():
 d=accepted('petsc_api_inventory');assert d['count']>=100 and d['headers']
 for a in d['api']:assert a['locations'] and a['kind'] in ('function','type','enum','macro')
 assert {'VecGhostGetLocalForm','PetscBool','KSPSetTolerances'}<={a['name'] for a in d['api']}
def test_bool_layout_audited():
 d=read('petsc_bool_audit');assert d['status']=='NO_DEPENDENCY_FOUND' and len(d['findings'])==5
