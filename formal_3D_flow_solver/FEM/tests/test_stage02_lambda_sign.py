import numpy as np
import pytest
from stage02_helpers import read_case


@pytest.mark.parametrize("profile",["pipe_coarse","pipe_medium","pipe_fine"])
@pytest.mark.parametrize("mode",["natural","reference"])
def test_multiplier_positive_for_positive_inflow(profile,mode):
    s=read_case(profile+"_"+mode,"solver")
    assert np.isfinite(s["lambda_pa"]) and s["lambda_pa"]>0
    assert s["real_global_dofs"]==sum(s["real_owned_dofs_per_rank"])==1
