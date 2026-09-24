from stage03_helpers import *

def test_signed_outlet_flows_and_measured_fractions_are_reported():
    f=qc()['flux'];names={'outlet_01','outlet_02','outlet_03'}
    assert set(f['Q_outlets_signed_m3_s'])==set(f['flow_fractions'])==names
    for n in names:close(f['flow_fractions'][n],f['Q_outlets_signed_m3_s'][n]/f['Q_out_total_m3_s'])
    close(sum(f['flow_fractions'].values()),1)
    expected='MANUAL_PHYSICS_REVIEW' if any(v<0 for v in f['Q_outlets_signed_m3_s'].values()) else 'ALL_NET_OUTFLOWS_POSITIVE'
    assert f['outlet_sign_review']==expected
    assert config()['gates']['outlet_positive_flux_hard_gate'] is False
