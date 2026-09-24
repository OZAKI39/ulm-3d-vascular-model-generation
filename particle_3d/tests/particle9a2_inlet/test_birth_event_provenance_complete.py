from particle_3d.injection_population import C_MB

def test_complete_event_and_entering_rate(synthetic):
    e=synthetic.event(3)
    fields=['particle_id','birth_time_s','diameter_m','radius_m','diameter_source_draw','diameter_global_rejections','position_proposal_count','position_xyz','inlet_triangle_id','local_flux_weight_m_s','clearance_m','admission_status','source_distribution_contract_sha256','flow_sha256','geometry_sha256','seed','diameter_fixed_during_position_sampling']
    assert set(fields)<=e.keys()
    assert e['birth_time_s']==synthetic.clock.time_at(3)
    assert abs(e['birth_time_s']*C_MB*synthetic.sampler.Q_m3_s-3)<1e-14
    assert e['admission_status']=='ACCEPTED' and e['diameter_fixed_during_position_sampling']
