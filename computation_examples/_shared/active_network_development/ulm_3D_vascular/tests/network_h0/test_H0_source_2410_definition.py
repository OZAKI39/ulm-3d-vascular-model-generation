import numpy as np
from network_1d0d.idealized_h0 import h0_boundary

def test_explicit_structural_source_needs_no_physiological_label(domain,baseline):
    source,terminals=h0_boundary(domain.ids,domain.edges)
    assert domain.ids[source]==2410
    assert source not in terminals
    assert baseline['unit'].pressure[source]==1.
    assert baseline['unit'].node_outflow[source]>0

def test_analysis_component_only_plus_exact_virtual_cuts(domain):
    assert len(domain.ids)==7422 and len(domain.edges)==7421
    assert np.count_nonzero(domain.ids>0)==7419
    assert set(domain.ids[domain.ids<0])=={-10001,-10002,-10003}
