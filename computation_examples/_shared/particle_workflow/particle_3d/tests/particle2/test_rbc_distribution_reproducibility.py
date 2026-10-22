import json
import runpy
import numpy as np
import pytest
from particle_3d.rbc_distribution import sample_rbc_geometries


def test_two_independent_processes_match_float64_and_csv(p2_repo):
    code=runpy.run_path(str(p2_repo/"particle_3d/scripts/verify_rbc_reproducibility.py"))
    record=code["verify"](257,9917)
    assert record["status"]=="PASS"
    assert record["csv_byte_identical"] and record["float64_arrays_identical"]


def test_full_population_reproducibility_receipt_matches_saved_csv(p2_data,population):
    from particle_3d.rbc_distribution import digest
    receipt=json.loads((p2_data/"reproducibility.json").read_text())
    assert receipt["N"]==100000 and receipt["independent_process_count"]==2 and receipt["status"]=="PASS"
    for record in receipt["process_records"]:
        assert record["csv_sha256"]==digest(p2_data/"C57BL6_RBC_GEOMETRY_VALIDATION_100000.csv")
        assert record["sample_structured_array_sha256"]==population.metadata["sample_structured_array_sha256"]
    assert not receipt["cross_numpy_version_byte_identity_claimed"]


@pytest.mark.parametrize("n,seed",[(0,3),(-1,3),(True,3),(3,None),(3,-1),(3,True),(3,1.5)])
def test_rng_requires_explicit_valid_n_and_seed(n,seed):
    with pytest.raises(ValueError):sample_rbc_geometries(n,seed)


def test_sampler_does_not_touch_global_rng_and_owns_output():
    before=np.random.get_state()
    p=sample_rbc_geometries(5,23)
    after=np.random.get_state()
    np.testing.assert_array_equal(before[1],after[1])
    assert before[2:]==after[2:]
    with pytest.raises(ValueError):p.samples.setflags(write=True)
