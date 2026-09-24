from dataclasses import replace
import numpy as np
from network_1d0d.idealized_h0 import port_signs,solve_operating_point

def test_reversing_all_stored_edges_preserves_physical_port_flows(domain,baseline):
    edges=domain.edges[:,::-1].copy()
    signs=port_signs(domain.ids,domain.xyz_m,edges,domain.internal_roi,domain.ports)
    reversed_domain=replace(domain,edges=edges,signs=signs)
    result=solve_operating_point(reversed_domain)
    np.testing.assert_allclose(result['unit_roi_q'],baseline['unit_roi_q'],rtol=1e-10,atol=0)
    np.testing.assert_allclose(result['flow'],-baseline['flow'],rtol=1e-8,atol=1e-26)
    assert [r['coefficient'] for r in signs]==[-r['coefficient'] for r in domain.signs]
    assert all(r['to_roi_dot_outward']<-.99 for r in signs)
