from particle_3d.nearfield_regularization import contract,STATES

def test_scope():
 c=contract();assert c['near_field_scope']==['SPHERE_WALL','SPHERE_SPHERE']
 assert not c['particle7_started'] and c['no_cfd_executed'] and not c['full_production_hydrodynamics']
 assert {'FAR_FIELD','NEARFIELD_TRANSITION','FULL_LUBRICATION','CONTINUUM_HANDOFF_CONTACT','GEOMETRIC_HARD_CONTACT','NONSPHERICAL_NOT_SUPPORTED'}<=set(STATES)
