from stage02_helpers import coefficient_scaling_error


def test_mu_doubles_pressure_and_lambda_with_fixed_velocity():
    assert max(coefficient_scaling_error("mu_double",(1.,2.,2.)).values())<1e-11
