import numpy as np
from particle_3d.particle82a_visuals import entry_progress


def test_failed_entry_does_not_animate_unsearched_storage_suffix():
    path=np.array([[0,0,0,0],[1,0,0,10e-6]],float)
    result=dict(accepted=False,searched_distance_m=6e-6)
    center,distance,_=entry_progress(path,result,1.)
    assert distance==6e-6
    np.testing.assert_allclose(center,[0,0,6e-6],atol=1e-20)


def test_final_birth_distance_label_matches_saved_center():
    path=np.array([[0,0,0,0],[1,0,0,10e-6]],float)
    result=dict(accepted=True,s_birth_m=2e-6,birth_center_m=[0,0,2e-6])
    center,distance,_=entry_progress(path,result,.96)
    assert distance==2e-6
    assert np.array_equal(center,result['birth_center_m'])
