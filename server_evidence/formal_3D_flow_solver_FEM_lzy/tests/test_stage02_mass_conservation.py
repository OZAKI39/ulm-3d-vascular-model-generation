import pytest
from stage02_helpers import read_case


@pytest.mark.parametrize("profile",["pipe_coarse","pipe_medium","pipe_fine"])
@pytest.mark.parametrize("mode",["natural","reference"])
def test_independent_inlet_outlet_integrals_close(profile,mode):
    q=read_case(profile+"_"+mode,"flux")
    assert q["Q_out_signed_m3_s"]>0
    assert abs(q["Q_out_signed_m3_s"]+q["Q_in_signed_m3_s"])/q["target_Q_m3_s"]<=1e-10
    assert abs(q["wall_flux_m3_s"])/q["target_Q_m3_s"]<1e-12
