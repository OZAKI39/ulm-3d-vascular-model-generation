def test_actual_full_original_p5_and_lammps_replay(real_result):
    d=real_result
    assert d['status']=='PASS' and d['requested_steps']==256
    assert d['original_p5_state_max_error']==d['pair_gap_max_error_m']==d['wall_gap_max_error_m']==d['residual_max_error']==0
    assert d['nearfield_activation_equal'] and d['outlet_events_equal']
    assert all(e['exact_equal'] and e['neighbor_mismatch_count']==0 for e in d['errors'])
    init=d['initialization']['original_p4_initialization']
    assert d['dt_s']==init['validation_dt_s'] and d['horizon_s']==init['horizon_s']
    for row in d['bridge_world_states']:
        for p,old in zip(row['particles'],init['initial_state']['particles']):assert p['radius_m']==old['radius_m']
