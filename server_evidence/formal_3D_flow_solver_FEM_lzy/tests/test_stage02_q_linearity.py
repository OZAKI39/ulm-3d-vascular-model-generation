import pytest
from stage02_helpers import coefficient_scaling_error


@pytest.mark.parametrize("case,factor",[("q_half",.5),("q_double",2.)])
def test_all_velocity_and_pressure_coefficients_and_lambda_scale_with_Q(case,factor):
    assert max(coefficient_scaling_error(case,(factor,factor,factor)).values())<1e-11
