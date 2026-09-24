from stage03_helpers import *

def test_actual_facet_integral_matches_target_at_fixed_direct_solver_tolerance():
    f=qc()['flux'];Q=config()['physics']['inlet_volume_flow_m3_s']
    assert f['Q_in_signed_m3_s']<0
    assert f['Q_in_m3_s']==-f['Q_in_signed_m3_s']
    error=abs(f['Q_in_m3_s']-Q)/Q
    assert error==f['relative_inlet_error'] and error<=1e-10
    assert 'integrated on physical tagged facets' in f['method']
