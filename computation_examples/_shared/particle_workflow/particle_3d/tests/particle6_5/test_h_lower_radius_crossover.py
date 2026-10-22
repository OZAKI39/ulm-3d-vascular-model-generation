def test_crossover(policy):
 assert policy.lower_handoff_gap(1.999e-6)['which_component_dominates']=='MOLECULAR'
 assert policy.lower_handoff_gap(2.001e-6)['which_component_dominates']=='SUSPENSION'
 assert policy.lower_handoff_gap(2e-6)['h_lower_m']==2e-9
