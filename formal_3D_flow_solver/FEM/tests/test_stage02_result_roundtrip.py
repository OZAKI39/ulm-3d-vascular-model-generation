from stage02_helpers import read_report


def test_new_process_restores_primary_and_derived_results():
    result=read_report("result_roundtrip")
    assert result["status"]=="PASS" and result["case_count"]>=12
    for case in result["cases"].values():
        assert case["status"]=="PASS" and case["new_process"]
        assert case["relative_flux_difference"]<1e-11
        assert case["primary_coefficients_equal"] and case["derived_values_equal"]
