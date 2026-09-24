from particle_3d.particle6_validation import resistance_case

def test_resistance_rhs_solution_contact_active_set(mu):
    d=resistance_case(mu)
    assert all(e==0 for e in d['max_errors'].values())
    assert d['sparsity_equal'] and d['active_constraints_equal']
    assert len(d['standalone']['active_constraints'])>0
