import numpy as np


def test_exact_si_conversions(population):
    s=population.samples
    for a in "abc":np.testing.assert_array_equal(s[f"{a}_m"],s[f"{a}_um"]*1e-6)
    np.testing.assert_array_equal(s["volume_m3"],s["V_fL"]*1e-18)
