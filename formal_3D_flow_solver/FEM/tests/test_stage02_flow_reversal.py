from stage02_helpers import coefficient_scaling_error,read_case


def test_internal_signed_rhs_reverses_fields_without_absolute_lambda():
    assert max(coefficient_scaling_error("reverse",(-1.,-1.,-1.)).values())<1e-11
    assert read_case("reverse","solver")["lambda_pa"]<0
    assert read_case("reverse","flux")["Q_in_signed_m3_s"]>0
