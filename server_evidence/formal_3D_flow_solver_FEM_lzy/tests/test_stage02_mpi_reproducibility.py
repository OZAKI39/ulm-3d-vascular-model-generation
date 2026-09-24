from stage02_helpers import read_report


def test_two_ranks_preserve_integrals_and_one_global_multiplier():
    result=read_report("mpi_reproducibility")
    assert result["status"]=="PASS" and result["ranks"]==[1,2]
    assert max(result["relative_differences"].values())<=result["tolerance"]
    assert result["real_global_dofs"]==1
