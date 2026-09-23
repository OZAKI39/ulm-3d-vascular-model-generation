from stage03_helpers import *

def test_all_three_integrated_outlets_close_the_prescribed_inlet():
    f=qc()['flux'];out=sum(f['Q_outlets_signed_m3_s'].values())
    assert out==f['Q_out_total_m3_s']
    error=abs(out-f['Q_in_m3_s'])/f['Q_target_m3_s']
    assert error==f['relative_mass_closure'] and error<=1e-10
