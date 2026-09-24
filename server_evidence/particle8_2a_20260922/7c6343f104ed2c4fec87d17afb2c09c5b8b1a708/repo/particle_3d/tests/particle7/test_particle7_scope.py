def test_scope(contract):
 assert contract['particle8_started'] is False and contract['no_cfd_executed'] is True
 assert not any(contract[k] for k in ['production_particle_timestep_frozen','production_neighbor_cutoff_frozen','production_neighbor_skin_frozen','nonspherical_lubrication_frozen','production_mb_pk_profile_frozen','tube_hct_enforced'])
 assert contract['real_rbc_passage']=='NOT_ESTABLISHED'
