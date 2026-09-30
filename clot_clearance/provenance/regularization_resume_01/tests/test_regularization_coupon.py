import json
from pathlib import Path
import numpy as np
from pd_clot.regularization.coupon import run_coupon, cyclic_increment


def test_coupon_measures_actual_force_work_and_complete_crack():
    c=json.loads((Path(__file__).parents[1]/'configs/streaming_regularized_demo_coarse.json').read_text())
    r=run_coupon(c,.000125,.012520992043877192)
    assert r['weak_plane_complete']
    assert r['relative_work_residual']<1e-7
    assert abs(r['Gc_measured_J_m2']/.01-1)<.01
    assert abs(r['final_stored_energy_after_unloading_J'])<1e-16
    assert np.isclose(r['dissipated_fracture_work_J']/r['created_crack_area_m2'],r['Gc_measured_J_m2'])


def test_energy_increment_extensive_and_exposure_additive():
    W=np.array([0.,1.,2.]);crit=np.array([1.,2.,3.])
    a=cyclic_increment(W,crit,1000,.001,1.5)
    np.testing.assert_allclose(a,2*cyclic_increment(W,crit,500,.001,1.5))
    np.testing.assert_allclose(a,cyclic_increment(W/8,crit/8,1000,.001,1.5))
