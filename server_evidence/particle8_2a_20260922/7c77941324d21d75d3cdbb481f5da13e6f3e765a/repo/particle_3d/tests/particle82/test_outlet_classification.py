import pytest
from particle_3d.particle82_integration import environment

def test_each_completed_natural_exit_crosses_original_cap(saved_scene):
    s=saved_scene;env=environment()
    for pid in s.completed:
        a=s.arrays[pid];e=s.entries[pid];hit=env.classifier.first_event(a[-2,1:4],a[-1,1:4])
        assert hit is not None and hit.role==e['exit_outlet']
        assert e['exit_time_s']==pytest.approx(e['birth_time_s']+a[-1,0],abs=1e-12)
        with pytest.raises(ValueError,match='extrapolation'):s.position(pid,a[-1,0]+1e-5)
