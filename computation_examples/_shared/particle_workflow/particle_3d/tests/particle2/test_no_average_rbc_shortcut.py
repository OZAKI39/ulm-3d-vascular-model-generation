import numpy as np
from particle_3d.rbc_distribution import quantile_indices,stratified_indices


def test_biological_geometries_are_actual_population_members(population,geometries):
    ids=quantile_indices(population.samples)
    assert len(set(ids))==5
    for i,g in zip(ids,geometries):
        assert g.provenance.rbc_id==population.samples[i]["rbc_id"]
        assert g.provenance.role=="RBC_GEOMETRY_VALIDATION_ONLY"
        assert g.a_m==population.samples[i]["a_m"]
        assert g.c_m==population.samples[i]["c_m"]
    assert max(g.r for g in geometries)/min(g.r for g in geometries)>3
    wide=stratified_indices(population.samples,2026092064)
    assert len(set(wide))==64
    np.testing.assert_array_equal(wide,stratified_indices(population.samples,2026092064))
