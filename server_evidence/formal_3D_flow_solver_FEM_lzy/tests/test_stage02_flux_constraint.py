import pytest
from stage02_helpers import read_case


@pytest.mark.parametrize("profile",["pipe_coarse","pipe_medium","pipe_fine"])
@pytest.mark.parametrize("mode",["natural","reference"])
def test_integrated_inflow_matches_positive_target(profile,mode):
    q=read_case(profile+"_"+mode,"flux")
    assert q["Q_in_signed_m3_s"]<0 and q["actual_Q_in_m3_s"]>0
    assert abs(q["actual_Q_in_m3_s"]-q["target_Q_m3_s"])/q["target_Q_m3_s"]<=1e-10
    assert q["relative_inlet_constraint_error"]<=1e-10
