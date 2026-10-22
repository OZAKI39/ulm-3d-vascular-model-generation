import numpy as np
from particle_3d.rbc_distribution import classify_candidate


def test_thick_candidate_is_rejected_whole_without_clipping(contract,population):
    assert classify_candidate(4.,40.,contract)==("SHAPE",None)
    rejected=population.candidates[population.candidates["status"]=="SHAPE"]
    assert len(rejected)>0
    assert np.all(3*rejected["V_raw_fL"]/(np.pi*rejected["D_raw_um"]**2)>=rejected["D_raw_um"]/2)
    assert not np.isin(rejected["candidate_id"],population.samples["candidate_id"]).any()


def test_every_accepted_geometry_is_finite_oblate_and_guarded(population):
    s=population.samples
    assert all(np.isfinite(s[name]).all() for name in s.dtype.names)
    assert np.all((s["D_um"]>=4)&(s["D_um"]<=9.58))
    assert np.all((s["V_fL"]>=26.6324)&(s["V_fL"]<=69.1676))
    assert np.array_equal(s["a_m"],s["b_m"])
    assert np.all((s["c_m"]>0)&(s["c_m"]<s["a_m"]))
    assert np.all((s["r"]>0)&(s["r"]<1)&(s["jeffery_lambda"]>-1)&(s["jeffery_lambda"]<0))
