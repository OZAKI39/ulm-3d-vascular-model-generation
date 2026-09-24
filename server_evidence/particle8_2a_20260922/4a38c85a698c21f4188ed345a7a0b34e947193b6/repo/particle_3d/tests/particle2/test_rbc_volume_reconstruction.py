import numpy as np


def test_all_100000_ellipsoid_volumes_reconstruct_sampled_volume(population):
    s=population.samples
    np.testing.assert_allclose(4*np.pi*s["a_um"]*s["b_um"]*s["c_um"]/3,s["V_fL"],rtol=32*np.finfo(float).eps,atol=0)
    np.testing.assert_allclose(4*np.pi*s["a_m"]*s["b_m"]*s["c_m"]/3,s["volume_m3"],rtol=32*np.finfo(float).eps,atol=0)
