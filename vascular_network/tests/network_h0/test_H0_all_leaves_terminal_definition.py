import numpy as np

def test_exact_leaf_set_and_zero_gauge_pressure(domain,baseline):
    expected=np.flatnonzero((np.bincount(domain.edges.ravel())==1)&(domain.ids!=2410))
    np.testing.assert_array_equal(domain.terminals,expected)
    assert len(expected)==122
    np.testing.assert_array_equal(baseline['pressure'][expected],0.)
    assert domain.port_nodes[3] in expected and domain.ids[domain.port_nodes[3]]==4484
