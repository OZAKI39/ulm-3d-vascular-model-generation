import math
from sv12_support import accepted
def test_accepted_outlet_split_from_integrals():
    q=accepted();assert q['accepted']
    assert math.isclose(sum(q['outlet_fractions'].values()),1.,abs_tol=1e-14)
    for role,flow in q['outlet_flows_m3_s'].items():assert q['outlet_fractions'][role]==flow/q['Q_out_total_m3_s']
