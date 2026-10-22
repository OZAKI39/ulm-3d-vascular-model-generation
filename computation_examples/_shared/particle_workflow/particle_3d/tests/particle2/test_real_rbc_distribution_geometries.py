import numpy as np
from particle_3d.rbc_distribution import quantile_indices


def test_real_cases_use_five_actual_shapes_same_start_orientation_and_three_dt(real_cases,population):
    ids=quantile_indices(population.samples)
    assert len(real_cases)==15
    all_dts=set()
    for rows,s in real_cases:
        original=population.samples[ids[s["geometry_index"]]]
        assert s["rbc_id"]==original["rbc_id"]
        for a in "abc":assert s["geometry"][f"{a}_m"]==original[f"{a}_m"]
        p=np.array([rows[0][f"p_{a}"] for a in "xyz"])
        np.testing.assert_allclose(p,np.array([1,2,3])/np.sqrt(14),atol=16*np.finfo(float).eps,rtol=0)
        assert s["initial_orientation_role"]=="VALIDATION_INITIAL_ORIENTATION_ONLY"
        assert s["timestep_role"]=="VALIDATION_ONLY"
        assert not s["production_particle_timestep_frozen"]
        assert not s["production_orientation_distribution_frozen"]
        all_dts.add(s["validation_dt_s"])
    coarse,middle,fine=sorted(all_dts,reverse=True)
    assert coarse==2*middle==4*fine
